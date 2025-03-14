/*
 
COPIADO DEL LINK: http://ftp-archive.freebsd.org

http://www.google.es/url?sa=t&rct=j&q=&esrc=s&source=web&cd=2&ved=0CCoQFjAB&url=http%3A%2F%2Fftp-archive.freebsd.org%2Fpub%2FFreeBSD-Archive%2Fold-releases%2Fi386%2F1.0-RELEASE%2Fports%2Furt%2Flib%2Fcolorquant.c&ei=gmygU_nNFoay0QXouYD4DA&usg=AFQjCNFC72mmG4Ju4DDa6R55MRMei1qQJg&bvm=bv.68911936,d.bGE&cad=rja

 * This software is copyrighted as noted below.  It may be freely copied,
 * modified, and redistributed, provided that the copyright notice is 
 * preserved on all copies.
 * 
 * There is no warranty or other guarantee of fitness for this software,
 * it is provided solely "as is".  Bug reports or fixes may be sent
 * to the author, who may or may not act on them as he desires.
 *
 * You may not include this software in a program or other software product
 * without supplying the source, or without informing the end-user that the 
 * source is available for no extra charge.
 *
 * If you modify this software, you should include a notice giving the
 * name of the person performing the modification, the date of modification,
 * and the reason for such modification.
 */
/*
 * 
 *
 * Perform variance-based color quantization on a "full color" image.
 * Author:	Craig Kolb
 *		Department of Mathematics
 *		Yale University
 *		kolb@yale.edu
 * Date:	Tue Aug 22 1989
 * Copyright (C) 1989 Craig E. Kolb
 * $Id: colorquant.c,v 1.1.1.1 1993/10/05 08:51:50 ljo Exp $
 *
 * $Log: colorquant.c,v $
 * Revision 1.1.1.1  1993/10/05  08:51:50  ljo
 * Utah Raster Toolkit
 *
 * Revision 3.0.1.2  90/11/29  15:18:04  spencer
 * Remove a typo.
 * 
 * Revision 3.0.1.1  90/11/19  16:59:48  spencer
 * Use inv_cmap instead of find_colors.
 * Changes to process multiple files into one colormap (accum_hist arg).
 * Delete 'otherimages' argument -- unnecessary with faster inv_cmap code.
 * 
 * 
 * Revision 3.0  90/08/03  15:20:11  spencer
 * Establish version 3.0 base.
 * 
 * Revision 1.6  90/07/29  08:06:06  spencer
 * If HUGE isn't defined, make it HUGE_VAL (for ansi).
 * 
 * Revision 1.5  90/07/26  17:25:48  rgb
 * Added a parameter to colorquant for rgbmap construction.
 * 
 * Revision 1.4  90/07/13  14:53:31  spencer
 * Get rid of ARB_ARG sh*t.
 * Change a couple of vars to double so that HUGE won't cause problems.
 * 
 * Revision 1.3  90/06/28  21:42:56  spencer
 * Declare internal functions properly.
 * Delete unused global variable.
 * 
 * Revision 1.2  90/06/28  13:18:56  spencer
 * Make internal functions static.
 * Build entire RGB cube, not just those colors used by this image.  This
 * is so dithering will work.
 * 
 * Revision 1.1  90/06/18  20:45:17  spencer
 * Initial revision
 * 
 * Revision 1.3  89/12/03  18:27:16  craig
 * Removed bogus integer casts in distance calculation in makenearest().
 * 
 * Revision 1.2  89/12/03  18:13:12  craig
 * FindCutpoint now returns FALSE if the given box cannot be cut.  This
 * to avoid overflow problems in CutBox.
 * "whichbox" in GreatestVariance() is now initialized to 0.
 * 
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


#include <sys/time.h> //Añado pues da error de prototipo en gettimeofday

//para usar imágenes ppm
#include "../archivos_imagenes/header.h" //añadida ruta

#define ERR_NUM_DATOS  "ERROR: numero de datos leidos distinto del esperado"



/* Ansi uses HUGE_VAL. */
#ifndef HUGE
#define HUGE HUGE_VAL
#endif



#define slow_color_map 1

static void QuantHistogram();
static void BoxStats();
static void UpdateFrequencies();
static void ComputeRGBMap();
static void SetRGBmap();
#ifdef slow_color_map
static void find_colors();
static int  getneighbors();
static int  makenearest();
#else
extern void inv_cmap();
#endif
static int  CutBoxes();
static int  CutBox();
static int  GreatestVariance();
static int  FindCutpoint();



/* 

 * The colorquant() routine has been tested on an Iris 4D70 workstation,
 * a Sun 3/60 workstation, and (to some extent), a Macintosh.
 * 
 * Calls to bzero() may have to be replaced with the appropriate thing on
 * your system.  (bzero(ptr, len) writes 'len' 0-bytes starting at the location
 * pointed to by ptr.)
 * 
 * Although I've tried to avoid integer overflow problems where ever possible,
 * it's likely I've missed a spot where an 'int' should really be a 'long'.
 * (On the machine this was developed on, an int == long == 32 bits.)
 * 
 * Note that it's quite easy to optimize this code for a given value for
 * 'bits'.  In addition, if you assume bits is small and
 * that the total number of pixels is relatively small, there are several
 * places that integer arithmetic may be substituted for floating-point.
 * One such place is the loop in BoxStats -- mean and var need not necessary
 * be floats.
 * 
 * As things stand, the maximum number of colors to which an image may
 * be quantized is 256.  This limit may be overcome by changing rgbmap and
 * colormap from pointers to characters to pointers to something larger.
 */

// para escribir las dimensiones del fichero de salida
int FILAS_FIG, COLS_FIG;

/*
 * Maximum number of colormap entries.  To make larger than 2^8, the rgbmap
 * type will have to be changed from unsigned chars to something larger.
 */
#define MAXCOLORS		256
/*
 * Value corrersponding to full intensity in colormap.  The values placed
 * in the colormap are scaled to be between zero and this number.  Note
 * that anything larger than 255 is going to lead to problems, as the
 * colormap is declared as an unsigned char.
 */
#define FULLINTENSITY		255 
#define MAX(x,y)	((x) > (y) ? (x) : (y))

/*
 * Readability constants.
 */
#define REDI		0	
#define GREENI		1
#define BLUEI		2	
#define TRUE		1
#define FALSE		0

typedef struct {
	float	        weightedvar,	/* weighted variance */
				    mean[3];		/* centroid */
	unsigned long 	weight,			    /* # of pixels in box */
			        freq[3][MAXCOLORS];	/* Projected frequencies */
	int 		    low[3], high[3];	/* Box extent */
} Box;



static unsigned long	*Histogram,		/* image histogram */
			            NPixels,		/* # of pixels in an image*/
			            SumPixels;		/* total # of pixels */
static unsigned int	Bits,			/* # significant input bits */
			        ColormaxI;		/* # of colors, 2^Bits */
static Box		    *Boxes;			/* Array of color boxes. */

/*
colorquant(Ir,Ig,Ib,size,colormap,MAXCOLORS,bits,rgbmap,fast,accum_hist);
 * Perform variance-based color quantization on a 24-bit image.
 *
 * Input consists of:
 *	red, green, blue	Arrays of red, green and blue pixel
 *				intensities stored as unsigned characters.
 *				The color of the ith pixel is given
 *				by red[i] green[i] and blue[i].  0 indicates
 *				zero intensity, 255 full intensity.
 *	pixels			The length of the red, green and blue arrays
 *				in bytes, stored as an unsigned long int.
 *	colormap		Points to the colormap.  The colormap
 *				consists of red, green and blue arrays.
 *				The red/green/blue values of the ith
 *				colormap entry are given respectively by
 *				colormap[0][i], colormap[1][i] and
 *				colormap[2][i].  Each entry is an unsigned char.
 *	colors			The number of colormap entries, stored
 *				as an integer.	
 *	bits			The number of significant bits in each entry
 *				of the red, green and blue arrays. An integer.
 *	rgbmap		An array of unsigned chars of size (2^bits)^3.
 *				This array is used to map from pixels to
 *				colormap entries.  The 'prequantized' red,
 *				green and blue components of a pixel
 *				are used as an index into rgbmap to retrieve
 *				the index which should be used into the colormap
 *				to represent the pixel.  In short:
 *				index = rgbmap[(((r << bits) | g) << bits) | b];
 * 	fast		If non-zero, the rgbmap will be constructed
 *				quickly.  If zero, the rgbmap will be built
 *				much slower, but more accurately.  In most
 *				cases, fast should be non-zero, as the error
 *				introduced by the approximation is usually
 *				small.  'Fast' is stored as an integer.
 *	accum_hist		If non-zero the histogram will accumulate and 
 *				reflect pixels from multiple images.
 *				when 1, the histogram will be initialized and
 *				summed, but not thrown away OR processed. when 
 *				2 the image RGB will be added to it.  When 3 
 *				Boxes are cut and a colormap and rgbmap
 *				are be returned, Histogram is freed too.
 *				When zero, all code is executed as per normal.
 *
 * colorquant returns the number of colors to which the image was
 * quantized.
 */
#define INIT_HIST 1
#define USE_HIST 2
#define PROCESS_HIST 3

unsigned char  *Ir, *Ig, *Ib;
int  *Ir_i, *Ig_i, *Ib_i; 

/* --------------------------------------------------------------------
 * Calls to bzero() may have to be replaced with the appropriate thing on
 * your system.  (bzero(ptr, len) writes 'len' 0-bytes starting at the location
 * pointed to by ptr.)
 * tipo = 1: he pasado un vector dinamico de colors elementos
 *      = 0: he pasado un elemento de una matriz
-------------------------------------------------------------------- */
void MI_bzero(Box *ptr, int len,  int tipo, int colors)
{

	Box *ptr_aux=ptr;	
	int i, j;
	
	if(tipo == 1)
	{
	  for(i=0; i< colors; i++)
	  {
		  //ptr_aux->weightedvar = 0;
		  //ptr_aux->mean[0] = ptr_aux->mean[1] = ptr_aux->mean[2] = 0;
		  //ptr_aux->weight = 0;
		 // ptr_aux->low[0]  =  ptr_aux->low[1]  =  ptr_aux->low[2] = 0;
		  //ptr_aux->high[0] =  ptr_aux->high[1] =  ptr_aux->high[2] = 0;
		  
	      for(j=0; j<MAXCOLORS; j++)
	      {
	        ptr_aux->freq[0][j] = ptr_aux->freq[1][j] = ptr_aux->freq[2][j] = 0;
	      }
		  
		  ptr_aux++;
	  }
	}else
	{
		//  ptr->weightedvar = 0;
//		  ptr->mean[0] = ptr->mean[1] = ptr->mean[2] = 0;
	//	  ptr->weight = 0;
	//	  ptr->low[0]  =  ptr->low[1]  =  ptr->low[2] = 0;
		//  ptr->high[0] =  ptr->high[1] =  ptr->high[2] = 0;
	      for(j=0; j<MAXCOLORS; j++)
	      {
		    ptr->freq[0][j] = ptr->freq[1][j] = ptr->freq[2][j] = 0;
	      }

	}
}




/* --------------------------------------------------------------
* Lee de un fichero PPM los pixels de la imagen original y los
* copia en 3 vectores (para cada una de las 3 componentes de color RGB).
*
* PARAMETROS:
*    nombre_fich: nombre del fichero que contiene los datos
* -------------------------------------------------------------- */
void leer_pixels_imagen_ppm(char nombre_fich[])
{
Image *org_img;  // imagen PPM de entrada
uchar *org_data; // valores de los pixels que se extraen de org_min
int num_elems,   // número de elementos de org_data (3*nº de pixels de la imagen)
  index,         // contador de pixels de la imagen
  ik;            // contador de elementos de org_data
int n_puntos; // total de pixels de la imagen
int n_leidos = 0;  // numero de pixels realmente leidos del fichero



   org_img = read_img (nombre_fich);


   // copio en las variables globales las dimensiones de la imagen    
   FILAS_FIG = org_img->num_rows;  // - altura de la imagen           
   COLS_FIG = org_img->num_cols;   // - anchura de la imagen
         
   
  n_puntos = FILAS_FIG * COLS_FIG;
   
  // reservar espacio para los vectores de colores
   Ir = (unsigned char *) malloc( sizeof(unsigned char) * n_puntos); // para rojo
   Ig = (unsigned char *) malloc( sizeof(unsigned char) * n_puntos); // para verde
   Ib = (unsigned char *) malloc( sizeof(unsigned char) * n_puntos); // para azul

   if( (Ir == NULL) || (Ig == NULL) || (Ib == NULL) )
   {
	printf("\n ERROR al reservar espacio para los colores. FIN.");
	exit(1);
   }   

   
          
   // elementos del vector org_img
   num_elems = COLS_FIG *FILAS_FIG * 3;  
        
 
   // se pasan los datos de la imagen a forma de vector
   /* los datos se guardan en un vector, colocando consecutivos
   los valores R G y B de cada pixel */
   org_data = get_data_ptr ( org_img ); 
       

   /* Se recorren los num_lemens elmeentos de org_data, tomando cada
    grupo de 3 valores consecutivos como las componente RGB de un pixel.
    Dichas componentes se almacenan en los vectores globales rojo, verde
    y azul que usa el programa */		    
   index = 0; // primer pixel de la imagen
   for (ik = 0; ik < num_elems; ik += 3 )
   {
       Ir[index] = (unsigned char) org_data[ik];
       Ig[index] = (unsigned char) org_data[ik + 1];
       Ib[index] = (unsigned char) org_data[ik + 2];
          
      index++; // siguiente pixel de la imagen
   }


   // se libera la memoria asociada a la imagen PPM
   free_img ( org_img ); 
   
   /* si no se pudieron leer tantos valores como se esperaba, se
   avisa al usuario y el programa acaba */
   if (index != n_puntos)  
   {
      printf ("\n%s (index=%d de %d)", ERR_NUM_DATOS, index, n_puntos) ; 
      exit(-1);	
   }
}


/* ---------------------------------------------------------------	  
* pruebo a volcar imagen en formato PPM
*
* PARÁMETROS:
--------------------------------------------------------------- */	  
void volcar_imagen_ppm(char nombre_fich[], unsigned char *rgbmap, 
                       unsigned char colormap[3][MAXCOLORS], int bits)
{	  
Image *out_img;
uchar *out_data;
int ik;
int index; // para copiar datos del formato PPM al formato que usa Wu
int num_elems = FILAS_FIG * COLS_FIG * 3;	
int i;   
   
   out_img = alloc_img ( FILAS_FIG, COLS_FIG, MAX_RGB );
   out_data = get_data_ptr ( out_img );
   

   i = 0;
   for ( ik = 0; ik < num_elems; ik += 3 )
   {      
      //de variance-based:
      index = rgbmap[(((Ir[i] << bits) | Ig[i]) << bits) | Ib[i]];
      
      out_data[ik]     = (int) colormap[0][index];
      out_data[ik + 1] = (int) colormap[1][index];
      out_data[ik + 2] = (int) colormap[2][index];
            
      i++;
   }
      //printf("\n%s\n",nombre_fich);
   write_img ( out_img, nombre_fich);
 
   free_img(out_img);
 }







int colorquant(red,green,blue, pixels,colormap,colors,bits,rgbmap,fast,accum_hist)
unsigned char *red, *green, *blue; //int *red, *green, *blue; //unsigned char *red, *green, *blue;
long int pixels ; //unsigned long pixels;
unsigned char colormap[3][MAXCOLORS]; //int colormap[3][MAXCOLORS];  //unsigned char *colormap[3];
int colors, bits;
unsigned char *rgbmap;
int fast, accum_hist;
{
    int	i,			/* Counter */
    OutColors,			/* # of entries computed */
    Colormax;			/* quantized full-intensity */ 
    float	Cfactor;	/* Conversion factor */
int j; // contador
int aux1, aux2, aux3;
char char_aux;	
	unsigned char uc;



    if (accum_hist < 0 || accum_hist > PROCESS_HIST)
	   fprintf(stderr, "colorquant: bad value for accum_hist\n");

    ColormaxI = 1 << bits;	/* 2 ^ Bits */   // MAXCOLORS
    Colormax = ColormaxI - 1;
    Bits = bits;
    NPixels = pixels;
    Cfactor = (float)FULLINTENSITY / Colormax;

	// muestro algunos valores
	//printf("\nBits=%d, NPixels=%ld, Cfactor=%lf", bits, NPixels, (double) Cfactor);
	//printf("\nRango de colores: (%d, %d)", ColormaxI, Colormax);
    
	
    if (! accum_hist || accum_hist == INIT_HIST) 
	{
	//	printf("\nTam. del histograma: %d",ColormaxI*ColormaxI*ColormaxI );
	//	printf("\nTam. de Boxes: %d",colors );
	    Histogram = (unsigned long *)calloc(ColormaxI*ColormaxI*ColormaxI, sizeof(long)); //, 
	    Boxes = (Box *)malloc(colors * sizeof(Box));
		
		//YO: compruebo que la memoria dinámica se haya asignado correctamente
	  if( (Histogram==NULL) || (Boxes==NULL) )
	  {
			printf("\nERROR al reservar espacio para Boxes o Histogram. FIN");
			exit(1);
	  }
		
	  /*
	   * Zero-out the projected frequency arrays of the largest box.
	   */
		//printf("\nColormaxI=%d", ColormaxI);
		//printf("\ncolormaxI*sizeof(unsigned long) es: %ld", ColormaxI*sizeof(unsigned long));
	 /* bzero(Boxes->freq[0], ColormaxI * sizeof(unsigned long));
	  bzero(Boxes->freq[1], ColormaxI * sizeof(unsigned long));
	  bzero(Boxes->freq[2], ColormaxI * sizeof(unsigned long));*/
		
		MI_bzero(Boxes,ColormaxI * sizeof(unsigned long), 1, colors);
		//printf("\ncolors es %d", colors); fflush(stdout);
		/*for(j=0; j<colors;j++)
		{
			Boxes[j].freq[0] = Boxes[j].freq[1] = Boxes[j].freq[2] = 0;
		}*/
		//Boxes.freq[0]= Boxes.freq[1]=Boxes.freq[2]=0; //YO
	  SumPixels = 0;
		//printf("\nPaso 1-B");fflush(stdout);
    }

	
//printf("\nPaso 2");fflush(stdout);
    SumPixels += NPixels;

    if ( accum_hist != PROCESS_HIST ) 
	{
	  QuantHistogram(red, green, blue, &Boxes[0]);
//printf("\nPaso 2.1");fflush(stdout);
	}
	
//printf("\nPaso 3");fflush(stdout);
    if ( !accum_hist || accum_hist == PROCESS_HIST) 
	{
		
	  OutColors = CutBoxes(Boxes, colors);
//printf("\nPaso 4. OutColors es %d; Cfactor=%lf, colors es %d",  OutColors, Cfactor, colors);fflush(stdout);
		
	  /*
	  * We now know the set of representative colors.  We now
	  * must fill in the colormap and convert the representatives
	  * from their 'prequantized' range to 0-FULLINTENSITY.
	  */
	  for (i = 0; i < OutColors; i++) 
	  {
		//printf("\n[i=%d], %lf,", i+1, Boxes[i].mean[0]); fflush(stdout);
		/*printf("\n[i=%d], boxes[i]=%lf", i+1, (Boxes[i].mean[REDI])); fflush(stdout);
		printf("\n[i=%d], boxes[i]=%lf", i+1, (Boxes[i].mean[GREENI])); fflush(stdout);
		printf("\n[i=%d], boxes[i]=%lf", i+1, (Boxes[i].mean[BLUEI])); fflush(stdout);
*/
/*		  		printf("\n[i=%d], boxes[i]=%lf", i+1, (Boxes[i].mean[REDI] * Cfactor + 0.5)); fflush(stdout);
		printf("\n[i=%d], boxes[i]=%c", i+1, (unsigned char)(Boxes[i].mean[GREENI] * Cfactor + 0.5)); fflush(stdout);
		printf("\n[i=%d], boxes[i]=%c", i+1,  (unsigned char)(Boxes[i].mean[BLUEI] * Cfactor + 0.5)); fflush(stdout);
*/
		  //aux = (Boxes[i].mean[REDI] * Cfactor + 0.5)
		aux1 = (int) (Boxes[i].mean[REDI] * Cfactor + 0.5);
		aux2 = (int) (Boxes[i].mean[GREENI] * Cfactor + 0.5);
		aux3 = (int) (Boxes[i].mean[BLUEI] * Cfactor + 0.5);
		

		  colormap[0][i] = (unsigned char) aux1;
		  colormap[1][i] = (unsigned char) aux2;
		  colormap[2][i] = (unsigned char) aux3;
		  
		 // printf("\nPaso 4-%d", i+1);fflush(stdout);
	  }

	//	printf("\nPaso 4 terminado");fflush(stdout);
	  ComputeRGBMap(Boxes, OutColors, rgbmap, bits, colormap, fast);
//printf("\nPaso 4 -b terminado");fflush(stdout);
	  free((char *)Histogram);
	  free((char *)Boxes);
	  return OutColors;	/* Return # of colormap entries */
    }
    return 0;
}

 

/* --------------------------------------------------------------------
 * Compute the histogram of the image as well as the projected frequency
 * arrays for the first world-encompassing box.
 * -------------------------------------------------------------------- */
static void QuantHistogram(r, g, b, box)
register unsigned char *r, *g, *b;
Box *box;
{
	unsigned long *rf, *gf, *bf, i;
	
	rf = box->freq[0];
	gf = box->freq[1];
	bf = box->freq[2];
   // printf("\nEn QuantHistogram "); fflush(stdout); //%d %d %d)", (int)rf, (int)gf, (int)bf);
	/*
	 * We compute both the histogram and the proj. frequencies of
	 * the first box at the same time to save a pass through the
	 * entire image. 
	 */
	for (i = 0; i < NPixels; i++) 
	{
		rf[*r]++;
		gf[*g]++;
		bf[*b]++;
		Histogram[((((*r++)<<Bits)|(*g++))<<Bits)|(*b++)]++;
	}
	 //  printf("\nEn QuantHistogram "); fflush(stdout);
}


/* --------------------------------------------------------------------
 * Interatively cut the boxes.
 -------------------------------------------------------------------- */
static int CutBoxes(boxes, colors) 
Box	*boxes;
int	colors;
{
	int curbox;

	boxes[0].low[REDI] = boxes[0].low[GREENI] = boxes[0].low[BLUEI] = 0;
	boxes[0].high[REDI] = boxes[0].high[GREENI] = boxes[0].high[BLUEI] = ColormaxI;
	boxes[0].weight = SumPixels;

	BoxStats(&boxes[0]);

	for (curbox = 1; curbox < colors; curbox++) 
	{
		if (CutBox(&boxes[GreatestVariance(boxes, curbox)],  &boxes[curbox]) == FALSE)
				break;
	}

	return curbox;
}


/* --------------------------------------------------------------------
 * Return the number of the box in 'boxes' with the greatest variance.
 * Restrict the search to those boxes with indices between 0 and n-1.
 -------------------------------------------------------------------- */
static int GreatestVariance(boxes, n)
Box *boxes;
int n;
{
	register int i, whichbox = 0;
	float max;

	max = -1;
	for (i = 0; i < n; i++) 
	{
		if (boxes[i].weightedvar > max) 
		{
			max = boxes[i].weightedvar;
			whichbox = i;
		}
	}
	return whichbox;
}

/* --------------------------------------------------------------------
 * Compute mean and weighted variance of the given box.
 --------------------------------------------------------------------*/
static void BoxStats(box)
register Box *box;
{
	register int i, color;
	unsigned long *freq;
	float mean, var;

	if(box->weight == 0) 
	{
		box->weightedvar = 0;
		return;
	}

	box->weightedvar = 0.;
	for (color = 0; color < 3; color++) 
	{
		var = mean = 0;
		i = box->low[color];
		freq = &box->freq[color][i];
		
		for (; i < box->high[color]; i++, freq++) 
		{
			mean += i * (*freq);
			var += i*i* (*freq);
		}
		box->mean[color] = mean / (float)box->weight;
		//printf("\nMedia de color %d es %lf", color, box->mean[color]);

		
		box->weightedvar += var - box->mean[color]*box->mean[color]*
					(float)box->weight;
	}
	box->weightedvar /= SumPixels;
}

/* --------------------------------------------------------------------
 * Cut the given box.  Returns TRUE if the box could be cut, FALSE otherwise.
 -------------------------------------------------------------------- */
static int CutBox(box, newbox)
Box *box, *newbox;
{
	int i;
	double totalvar[3];
	Box newboxes[3][2];

	if (box->weightedvar == 0. || box->weight == 0)
		/*
		 * Can't cut this box.
		 */
		return FALSE;

	/*
	 * Find 'optimal' cutpoint along each of the red, green and blue
	 * axes.  Sum the variances of the two boxes which would result
	 * by making each cut and store the resultant boxes for 
	 * (possible) later use.
	 */
	for (i = 0; i < 3; i++) 
	{
		if (FindCutpoint(box, i, &newboxes[i][0], &newboxes[i][1]))
			totalvar[i] = newboxes[i][0].weightedvar +
				newboxes[i][1].weightedvar;
		else
			totalvar[i] = HUGE;
	}

	/*
	 * Find which of the three cuts minimized the total variance
	 * and make that the 'real' cut.
	 */
	if (totalvar[REDI] <= totalvar[GREENI] &&
	    totalvar[REDI] <= totalvar[BLUEI]) 
	{
		*box = newboxes[REDI][0];
		*newbox = newboxes[REDI][1];
	} else if (totalvar[GREENI] <= totalvar[REDI] &&
		 totalvar[GREENI] <= totalvar[BLUEI]) 
	{
		*box = newboxes[GREENI][0];
		*newbox = newboxes[GREENI][1];
	} else 
	{
		*box = newboxes[BLUEI][0];
		*newbox = newboxes[BLUEI][1];
	}

	return TRUE;
}

/* --------------------------------------------------------------------
 * Compute the 'optimal' cutpoint for the given box along the axis
 * indcated by 'color'.  Store the boxes which result from the cut
 * in newbox1 and newbox2.
 -------------------------------------------------------------------- */ 
static int FindCutpoint(box, color, newbox1, newbox2)
Box *box, *newbox1, *newbox2;
int color;
{
	float u, v, max;
	int i, maxindex, minindex, cutpoint;
	unsigned long optweight, curweight;

	if (box->low[color] + 1 == box->high[color])
		return FALSE;	/* Cannot be cut. */
	minindex = (int)((box->low[color] + box->mean[color]) * 0.5);
	maxindex = (int)((box->mean[color] + box->high[color]) * 0.5);

	cutpoint = minindex;
	optweight = box->weight;

	curweight = 0;
	for (i = box->low[color] ; i < minindex ; i++)
		curweight += box->freq[color][i];
	u = 0.;
	max = -1;
	for (i = minindex; i <= maxindex ; i++) 
	{
		curweight += box->freq[color][i];
		if (curweight == box->weight)
			break;
		u += (float)(i * box->freq[color][i]) /
					(float)box->weight;
		v = ((float)curweight / (float)(box->weight-curweight)) *
				(box->mean[color]-u)*(box->mean[color]-u);
		if (v > max) 
		{
			max = v;
			cutpoint = i;
			optweight = curweight;
		}
	}
	cutpoint++;
	*newbox1 = *newbox2 = *box;
	newbox1->weight = optweight;
	newbox2->weight -= optweight;
	newbox1->high[color] = cutpoint;
	newbox2->low[color] = cutpoint;
	UpdateFrequencies(newbox1, newbox2);
	BoxStats(newbox1);
	BoxStats(newbox2);

	return TRUE;	/* Found cutpoint. */
}

/* --------------------------------------------------------------------
 * Update projected frequency arrays for two boxes which used to be
 * a single box.
 -------------------------------------------------------------------- */
static void UpdateFrequencies(box1, box2)
register Box *box1, *box2;
{
	register unsigned long myfreq, *h;
	register int b, g, r;
	int roff;
int j; // contador
	
/*	bzero(box1->freq[0], ColormaxI * sizeof(unsigned long));
	bzero(box1->freq[1], ColormaxI * sizeof(unsigned long));
	bzero(box1->freq[2], ColormaxI * sizeof(unsigned long)); */
//box1->freq[0]=box1->freq[1]=box1->freq[2]=0; //YO
	MI_bzero(box1, ColormaxI * sizeof(unsigned long), 0, 0);
	
/*		for(j=0; j<colors;j++)
		{
			box1[j]->freq[0] = box1[j]->freq[1] = box1[j]->freq[2] = 0;
		}
	*/
	for (r = box1->low[0]; r < box1->high[0]; r++) 
	{
		roff = r << Bits;
		for (g = box1->low[1];g < box1->high[1]; g++) 
		{
			b = box1->low[2];
			h = Histogram + (((roff | g) << Bits) | b);
			for ( ; b < box1->high[2]; b++) 
			{
				if ((myfreq = *h++) == 0)
					continue;
				
				box1->freq[0][r] += myfreq;
				box1->freq[1][g] += myfreq;
				box1->freq[2][b] += myfreq;
				
				box2->freq[0][r] -= myfreq;
				box2->freq[1][g] -= myfreq;
				box2->freq[2][b] -= myfreq;
			}
		}
	}
}

/* --------------------------------------------------------------------
 * Compute RGB to colormap index map.
 -------------------------------------------------------------------- */
static void ComputeRGBMap(boxes, colors, rgbmap, bits, colormap, fast)
Box *boxes;
int colors;
unsigned char *rgbmap, *colormap[3];
int bits, fast;
{
	register int i;

	if (fast) 
	{
		//printf("\nRAPIDO...");
		/*
		 * The centroid of each box serves as the representative
		 * for each color in the box.
		 */
		for (i = 0; i < colors; i++)
			SetRGBmap(i, &boxes[i], rgbmap, bits);
	} 
	else
	{
		/*
		 * Find the 'nearest' representative for each
		 * pixel.
		 */
#ifdef slow_color_map
		find_colors(boxes, colors, rgbmap, bits, colormap, 0);
	//	printf("\nNORAPIDO...1");
#else
		inv_cmap(colors, colormap, bits, Histogram, rgbmap);
		//printf("\nNORAPIDO...2");
#endif
		//printf("\nNORAPIDO...");
	}
}

/* --------------------------------------------------------------------
 * Make the centroid of "boxnum" serve as the representative for
 * each color in the box.
 -------------------------------------------------------------------- */
static void SetRGBmap(boxnum, box, rgbmap, bits)
int boxnum;
Box *box;
unsigned char *rgbmap;
int bits;
{
	register int r, g, b;
	
	for (r = box->low[REDI]; r < box->high[REDI]; r++) 
	{
		for (g = box->low[GREENI]; g < box->high[GREENI]; g++) 
		{
			for (b = box->low[BLUEI]; b < box->high[BLUEI]; b++) 
			{
				rgbmap[(((r<<bits)|g)<<bits)|b]=(char)boxnum;
			}
		}
	}
}

#ifdef slow_color_map
/* --------------------------------------------------------------------
 * Form colormap and NearestColor arrays.
 --------------------------------------------------------------------*/
static void find_colors(boxes, colors, rgbmap, bits, colormap, otherimages)
int colors;
Box *boxes;
unsigned char *rgbmap;
int bits;
unsigned char *colormap[3];
int otherimages;
{
	register int i;
	int num, *neighbors;

	neighbors = (int *)malloc(colors * sizeof(int));

	/*
	 * Form map of representative (nearest) colors.
	 */
	for (i = 0; i < colors; i++) 
	{
		/*
		 * Create list of candidate neighbors and
		 * find closest representative for each
		 * color in the box.
		 */
		num = getneighbors(boxes, i, neighbors, colors, colormap);
		makenearest(boxes, i, num, neighbors, rgbmap, bits, colormap, 
			    otherimages);
	}
	free((char *)neighbors);
}

/* --------------------------------------------------------------------
 * In order to minimize our search for 'best representative', we form the
 * 'neighbors' array.  This array contains the number of the boxes whose
 * centroids *might* be used as a representative for some color in the
 * current box.  We need only consider those boxes whose centroids are closer
 * to one or more of the current box's corners than is the centroid of the
 * current box. 'Closeness' is measured by Euclidean distance.
 -------------------------------------------------------------------- */
static int getneighbors(boxes, num, neighbors, colors, colormap)
Box *boxes;
int num, colors, *neighbors;
int *colormap[3]; //unsigned char *colormap[3];
{
	register int i, j;
	Box *bp;
	float dist, LowR, LowG, LowB, HighR, HighG, HighB, ldiff, hdiff;

	bp = &boxes[num];

	ldiff = bp->low[REDI] - bp->mean[REDI];
	ldiff *= ldiff;
	hdiff = bp->high[REDI] - bp->mean[REDI];
	hdiff *= hdiff;
	dist = MAX(ldiff, hdiff);

	ldiff = bp->low[GREENI] - bp->mean[GREENI];
	ldiff *= ldiff;
	hdiff = bp->high[GREENI] - bp->mean[GREENI];
	hdiff *= hdiff;
	dist += MAX(ldiff, hdiff);

	ldiff = bp->low[BLUEI] - bp->mean[BLUEI];
	ldiff *= ldiff;
	hdiff = bp->high[BLUEI] - bp->mean[BLUEI];
	hdiff *= hdiff;
	dist += MAX(ldiff, hdiff);

#ifdef IRIS
	dist = fsqrt(dist);
#else
	dist = (float)sqrt((double)dist);
#endif

	/*
	 * Loop over all colors in the colormap, the ith entry of which
	 * corresponds to the ith box.
	 *
	 * If the centroid of a box is as close to any corner of the
	 * current box as is the centroid of the current box, add that
	 * box to the list of "neighbors" of the current box.
	 */
	HighR = (float)bp->high[REDI] + dist;
	HighG = (float)bp->high[GREENI] + dist;
	HighB = (float)bp->high[BLUEI] + dist;
	LowR = (float)bp->low[REDI] - dist;
	LowG = (float)bp->low[GREENI] - dist;
	LowB = (float)bp->low[BLUEI] - dist;
	for (i = j = 0, bp = boxes; i < colors; i++, bp++) 
	{
		if (LowR <= bp->mean[REDI] && HighR >= bp->mean[REDI] &&
		    LowG <= bp->mean[GREENI] && HighG >= bp->mean[GREENI] &&
		    LowB <= bp->mean[BLUEI] && HighB >= bp->mean[BLUEI])
			neighbors[j++] = i;
	}

	return j;	/* Return the number of neighbors found. */
}


/* --------------------------------------------------------------------
 * Assign representative colors to every pixel in a given box through
 * the construction of the NearestColor array.  For each color in the
 * given box, we look at the list of neighbors passed to find the
 * one whose centroid is closest to the current color.
 -------------------------------------------------------------------- */
static int makenearest(boxes, boxnum, nneighbors, neighbors, rgbmap, bits, colormap, 
	    otherimages)
Box *boxes;
int boxnum;
int nneighbors, *neighbors, bits;
unsigned char *rgbmap, *colormap[3];
int otherimages;
{
	register int n, b, g, r;
	double rdist, gdist, bdist, dist, mindist;
	int which, *np;
	Box *box;

	box = &boxes[boxnum];

	for (r = box->low[REDI]; r < box->high[REDI]; r++) 
	{
		for (g = box->low[GREENI]; g < box->high[GREENI]; g++) 
		{
			for (b = box->low[BLUEI]; b < box->high[BLUEI]; b++) 
			{
/*
 * The following "if" is to avoid doing extra work when the RGBmap is
 * not going to be used for other images.
 */
			        if ((!otherimages) && 
				    (Histogram[(((r<<bits)|g)<<bits)|b] == 0))
					continue;

				mindist = HUGE;
				/*
				 * Find the colormap entry which is
				 * closest to the current color.
				 */
				np = neighbors;
				for (n = 0; n < nneighbors; n++, np++) 
				{
					rdist = r-boxes[*np].mean[REDI];
					gdist = g-boxes[*np].mean[GREENI];
					bdist = b-boxes[*np].mean[BLUEI];
					dist = rdist*rdist + gdist*gdist + bdist*bdist;
					if (dist < mindist) 
					{
						mindist = dist;
						which = *np; 
					}
				}
				/*
				 * The colormap entry closest to this
				 * color is used as a representative.
				 */
				rgbmap[(((r<<bits)|g)<<bits)|b] = which;
			}
		}
	}
}
#endif /* slow_color_map */


/* ------------------------------------------------
 * Dado el nombre de un fichero de datos, que puede incluir un path y una
 * extensión, se queda sólo con el nombre (omite path y .extension)
 * (se hace para utilizar dicha cadena como parte de un nombre de fichero de salida)
 *
 *	 OJO: valido para linux, pues considero directorios separados con '/'
 * PARÁMETROS:
 *   cad: nombre de fichero, que puede contener ruta y extensión
 *   fich_E: nombre extraido de la cadena anterior (sin path previo ni .txt final)
 ------------------------------------------------ */
void extraer_nombre_fich(char cad[], char fich_E[])
{
int j, i;              // contadores
char fich_E_reves[30]; // solo el nombre del fichero de entrada, sin extensión

	
	j=0; // inicio a cero el contador
	
	/* recorro la cadena cad desde el final, y voy copiando en fich_E_reves, 
	 hasta que encuentro la barra que separa del nombre de un directorio */
	for(i=strlen(cad)-1; i>=0; i--) 
	{
			if (cad[i]!= '/') 
			{
				fich_E_reves[j]=cad[i]; 
				j++;
			}else
				break;
	}
	
	fich_E_reves[j]='\0'; // añado el caracter de fin de cadena

	
	/* recorro la cadena fich_E_reves desde el final, y voy copiando en fich_E,
	 hasta que encuentro el punto que separa la extensión */
	i=0;
	for(j=strlen(fich_E_reves)-1; j>=0; j--)  
	{
		if(fich_E_reves[j] !='.')
		{
		   fich_E[i] = fich_E_reves[j];
		   i++;
		}else
			break;
	}
	
 fich_E[i]='\0'; // añado el carácter de fin de cadena
	
#ifdef DEBUG	
	    printf("\nExtracccion del nombre del fichero que contiene la imagen: ");
	    printf("\nderecho: |%s|, reves: |%s|", fich_E, fich_E_reves);
#endif	    
}




void main(int argc, char *argv[])
{

/* añadido para indicar el fichero de datos (pixels) */
char nombre_fich[50]="";
long int size;
unsigned char colormap[3][MAXCOLORS];
//int colormap[3][MAXCOLORS];

unsigned char	*rgbmap; //			An array of unsigned chars of size (2^bits)^3.
int bits=8; //bits significativos de cada array R , G, B
//unsigned char rgbmap[256]; //256*256*256]; //pow(pow(2,8),3)];
int fast =0; //0 o 1
int accum_hist=0; //0; //sólo funciona con cero
int i, j;
	int rojo, verde, azul;
	int indice;
FILE *fs;
	int ret;
int size_int;


		// para calcular el tiempo empleado en los cálculos
// --- para calcular tiempo en segundos
time_t T_ini_a, T_fin_a; // instantes inicial y final de la aplicación del algoritmo 
// --- para calcular tiempo en milisegundos
struct timeval ti_a, tf_a;   // para la aplicación del algoritmo 

unsigned long long tiem_a, tiem_m;  // para almacenar la diferencia entre los dos tiempos previos
double tiem_double;
		
char nombre_fin_sin_path[100]; /* nombre de la imagen original, quitando el posible
	path, que usaré para construir el nombre del fichero de salida */



int colores; // nuevo, para leer desde el terminal el tamaño de la paleta cuantizada

	
	/* ---- input R,G,B components into Ir, Ig, Ib;
	set size to width*height ---- */
       // Código para leer el nombre del fich. de datos y cargar los datos
	/* si se paso al menos un argumento, se toma como nombre del fichero de
	datos sobre los articulos. Es obligatorio */
	if(argc >=3)
	{
		strcpy(nombre_fich, argv[1]);
		colores= atoi(argv[2]);
	}
	else
	{
		printf("\nDebe indicar el nombre del fichero que contiene los valores RGB de los pixels y el tam. de la paleta cuantizada.\n");
		exit(1);
	}
	
	
	/* quito al nombre del fichero la ruta y las multiples extensiones,
	para usarlo como parte del nombre del fichero de salida */
	extraer_nombre_fich(nombre_fich, nombre_fin_sin_path);

	(void) time(&T_ini_a); 
        gettimeofday(&ti_a, NULL);   // para calcularlo en usegundos



    leer_pixels_imagen_ppm(nombre_fich); 


   (void) time(&T_ini_a); 
   gettimeofday(&ti_a, NULL);   // para calcularlo en usegundos

	
	
NPixels = FILAS_FIG * COLS_FIG; 
size = NPixels;  



   // reservo memoria para el rgbmap
   rgbmap = (unsigned char *) malloc( sizeof(unsigned char)*(pow(pow(2, bits), 3))); // para rojo
  
   if( rgbmap == NULL) 
   {
	printf("\n ERROR al reservar espacio para rgbmap. FIN.");
	exit(1);
   }

	fast=1;
     ret = colorquant(Ir,Ig,Ib,size,colormap, colores,
                   bits,rgbmap,fast,accum_hist);
	                              


   char nombre_fich_s[150];
   

    sprintf(nombre_fich_s, "VB_%d_%s.ppm", colores, nombre_fin_sin_path);   
	
   (void) time(&T_fin_a);    
   gettimeofday(&tf_a, NULL);   // para calcularlo en usegundos


   

 	
   tiem_a = (tf_a.tv_sec - ti_a.tv_sec)*1000 +(tf_a.tv_usec - ti_a.tv_usec)/1000;

   printf("%s %s %llu\n", argv[2], argv[1], tiem_a);
   
 	
   //VUELCO LA IMAGEN A DISCO - formato PPM
   volcar_imagen_ppm(nombre_fich_s, rgbmap, colormap, bits);
	
	

	//Instrucciones para liberar memoria dinámica
	free(Ig); free(Ib); free(Ir); free(rgbmap);

}
