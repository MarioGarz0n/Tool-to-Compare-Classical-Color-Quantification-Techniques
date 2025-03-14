/*
 * COMPILAR:
 *   gcc TGA221_YO_ppm.c MEDIAN.c  -lm ppm_io.c mdmalloc.c -o mc_ppm
 *
 * EJECUTAR:
 *   ./mc_ppm <imagen.ppm> <numero colores paleta cuantizada>
 *
 *
** File:        tga221.c             Copyright (c) Truda Software
** Author:      Anton Kruger         215 Marengo Rd, #2, 
** Date:        March 1992           Oxford, IA 52322-9383
** Revision:    1.0                  March 1992
** 
** Description: A driver for "median". Reads in Targa type 2
**              (true-color) image files, applies the median-cut
**              algorithm to find 256 good colors for the image,
**              and writes out Targa type 1 (color-mapped) file.
**
** Compilers:  MSC 5.1, 6.0.
**
** Note:       Compile in large memory model.


*/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <math.h>    
#include <sys/time.h> 


#include "../archivos_imagenes/header.h"  //para usar imágenes ppm

#define ERR_NUM_DATOS  "ERROR: numero de datos leidos distinto del esperado"




#define IOBUFFSIZE 8192          /* size of i/o buffers        */

int MAXCOLORS =32; // LE DOY UN VALOR POR DEFECTO, POR SI ACASO


#define HSIZE     1240000 //262144 //AMPLIO YO // 32768         /* size of image histogram    */

typedef unsigned char byte;      /* range 0-255                */
typedef unsigned short word;     /* range 0-65,535             */
typedef unsigned long dword;     /* range 0-4,294,967,295      */

/* Macros for converting between (r,g,b)-colors and 15-bit     */
/* colors follow.                                              */
#define RGB(r,g,b) (word)(((b)&~7)<<7)|(((g)&~7)<<2)|((r)>>3)
#define RED(x)     (byte)(((x)&31)<<3)
#define GREEN(x)   (byte)((((x)>>5)&255)<< 3)
#define BLUE(x)    (byte)((((x)>>10)&255)<< 3)

word Hist[HSIZE];                /* 15-bit image histogram     */

struct ImageHeader {   /* Offset          Meaning              */ 
   byte  IdSize;       /*  0   id field size at end of header  */
   byte  ColMapType;   /*  1   0=>no color map, 1=>color map   */
   byte  TypeCode;     /*  2   1=>type 1, 2=>type 2 images,... */
   word  ColMapOrigin; /*  3   index of first color map entry  */
   word  ColMapLength; /*  5   number of entries in color map  */
   byte  ColMapBits;   /*  7   number of bits/color map entry  */
   word  XOrigin;      /*  8   lower left corner of image      */
   word  YOrigin;      /*  10  lower left corner of image      */
   word  Width;        /*  12  width of image in pixels        */
   word  Height;       /*  14  height of image in pixels       */
   byte  BitsPerPixel; /*  16  see note 1: below               */
   byte  Descriptor;   /*  17  see note 2: below               */
   byte  IdField[256]; /*  18  image identification field      */
   /*                                                          */
   /* Note 1: The "BitsPerPixel" field depends of the type of  */
   /*         Targa file. For Targa type 1 images, it is the   */
   /*         number of bits per pixel index.  For a true-     */
   /*         color Targa image this is 24, since pixels need  */
   /*         1 byte each for r, g, and b.                     */
   /* Note 2: Individual bits have meaning here. For example,  */
   /*         toggling bit 5 flips the image upside down.      */
};

unsigned char  *Ir, *Ig, *Ib;

void ReadTargaHeader(FILE *fp1, struct ImageHeader *Image);
void WriteTargaHeader(FILE *fp1, struct ImageHeader Image);
void getrgb(byte * r,byte * g,byte * b, FILE *fp1);
word MedianCut(word Hist[],byte ColMap[][3], int maxcubes);
void inform(char *mess);
void main(int argc, char *argv[]);


int FILAS_FIG, COLS_FIG; /* número de filas y columnas de la imagen original,
                          para usarlas al volcar la imagen cuantizada */




/* --------------------------------------------------------------
* Lee de un fichero PPM los pixels de la imagen original y los
* copia en 3 vectores (para cada una de las 3 componentes de color RGB).
*
* PARAMETROS:
*    nombre_fich: nombre del fichero que contiene los datos
* -------------------------------------------------------------- */
void leer_pixels_imagen_ppm(char nombre_fich[], int *n_f)
{
Image *org_img;  // imagen PPM de entrada
uchar *org_data; // valores de los pixels que se extraen de org_min
int num_elems,   // número de elementos de org_data (3*nº de pixels de la imagen)
  index,         // contador de pixels de la imagen
  ik;            // contador de elementos de org_data
int n_puntos; // total de pixels de la imagen
int n_leidos = 0;  // numero de pixels realmente leidos del fichero



   // se lee la imagen original (un fichero PPM)
   org_img = read_img (nombre_fich);


   // copio en las variables globales las dimensiones de la imagen    
   FILAS_FIG = org_img->num_rows;  // - altura de la imagen           
   COLS_FIG = org_img->num_cols;   // - anchura de la imagen

        
   // de Variance-based >>>     
   n_puntos = FILAS_FIG * COLS_FIG;
   
   
	/* el parámetro esperaba el número total de puntos, así que copio
	el producto de filas y columnas */
	*n_f = FILAS_FIG * COLS_FIG;
	   
  // reservar espacio para los vectores de colores
   Ir = (unsigned char *) malloc( sizeof(unsigned char) * (*n_f) ); // para rojo
   Ig = (unsigned char *) malloc( sizeof(unsigned char) * (*n_f) ); // para verde
   Ib = (unsigned char *) malloc( sizeof(unsigned char) * (*n_f) ); // para azul

   if( (Ir == NULL) || (Ig == NULL) || (Ib == NULL) )
   {
	printf("\n ERROR al reservar espacio para los colores. FIN.");
	exit(1);
   }   
   // de Variance-based <<<
   
          
   // elementos del vector org_img
   num_elems = COLS_FIG *FILAS_FIG * 3;  
        
 

   org_data = get_data_ptr ( org_img ); 
       
    
   index = 0; 
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




/* ------------------------------------------------
 * Dado el nombre de un fichero de datos, que puede incluir un path y una
 * extensión .txt, se queda sólo con el nombre (omite path y .txt)
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

	
	j=0; 
	

	for(i=strlen(cad)-1; i>=0; i--) 
	{
			if (cad[i]!= '/') 
			{
				fich_E_reves[j]=cad[i]; 
				j++;
			}else
				break;
	}
	
	fich_E_reves[j]='\0'; 

	

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
	
 fich_E[i]='\0';     
}




/* ---------------------------------------------------------------	  
* Volcar imagen en formato PPM
*
* PARÁMETROS:
--------------------------------------------------------------- */	  
void volcar_imagen_ppm(char nombre_fich[], byte ColMap[256][3]  )
{	  
Image *out_img;
uchar *out_data;
int ik;
int num_elems = FILAS_FIG * COLS_FIG * 3;	
int j;  


byte r_yo, g_yo, b_yo;
byte index; 
byte               r,g,b;
word             color;

   
   out_img = alloc_img ( FILAS_FIG, COLS_FIG, MAX_RGB );
   out_data = get_data_ptr ( out_img );
   

   j = 0;
   for ( ik = 0; ik < num_elems; ik += 3 )
   {      
      g=Ig[j]; b=Ib[j]; r=Ir[j]; 
      color = RGB(r,g,b);
      index = (byte)Hist[color];    


      r_yo = ColMap[index][0]; 
      g_yo = ColMap[index][1]; 
      b_yo = ColMap[index][2]; 

      out_data[ik]     = (int) r_yo;
      out_data[ik + 1] = (int) g_yo;
      out_data[ik + 2] = (int) b_yo;
            
      j++;
   }
      
            
   write_img ( out_img, nombre_fich);
 
   free_img(out_img);
 }




void main(int argc, char *argv[])
{
   byte               r,g,b;
   byte               index,ColMap[256][3];
   word               i,n,color;
   struct ImageHeader Image;	


int tam=0;    // numero de pixels comparados
byte r_yo, g_yo, b_yo;
char nombre_fich[50]; 
int size;
FILE *fd_salida;
	int j; 

// --- para calcular tiempo en milisegundos
struct timeval ti_a, tf_a;   // para la aplicación del algoritmo 

unsigned long long tiem_a, tiem_m;  // para almacenar la diferencia entre los dos tiempos previos
	
char nombre_fs[150], // nombre del fichero de salida
     nombre_fich_sin_path[100]; // nombre del fich de entrada sin path ni parte final separada por .



   MAXCOLORS = 32; //un valor por defecto    // tamaño de la paleta

   
   if (argc >=3)
   {
      strcpy(nombre_fich, argv[1]);  //fichero PPM de entrada
      MAXCOLORS = atoi(argv[2]);     // tamaño de la paleta
   }
   else
   {
	printf("\nEJECUCION: ./a.out fich_datos.ppm <tamaño-paleta>");
	exit(1);
   }

	/* tomo la parte del fichero de entrada que excluye el path y los trozos
	finales separados por puntos */
	extraer_nombre_fich(nombre_fich, nombre_fich_sin_path);

	//defino el nombre del fichero de salida
	sprintf(nombre_fs, "MC_%d_%s.ppm", MAXCOLORS, nombre_fich_sin_path);

	
// <---- fin
   gettimeofday(&ti_a, NULL);   // para calcularlo en usegundos

   
   //para leer imagen desde PPM
   leer_pixels_imagen_ppm(nombre_fich, &size);
   
   
   //printf("\n%s COLORES=%d", nombre_fich, MAXCOLORS);   

  
   /*
   ** Clear the histogram, read the image header and then build
   ** the image histogram. Next call "MedianCut" to quantize.
   */

 gettimeofday(&ti_a, NULL);   // para calcularlo en usegundos NO CUENTO TIEMPO DE LECTURA 
	
   i = (word)HSIZE;


   while (i--)
   {
      Hist[i] = (word)0;
   }
	


   for(j=0; j<size; j++)
   {
      g=(byte) Ig[j]; b=(byte)Ib[j]; r=(byte)Ir[j]; //      color = RGB(r,g,b);
      color = RGB(r,g,b);
      Hist[color] = Hist[color] + 1;
   }


// <-------------

   n = MedianCut(Hist,ColMap,(int)MAXCOLORS);


	
   /*
   ** Rewind the input file, and skip over the header. We are
   ** going write out a Targa type 1 image, so make the proper
   ** changes to the image header. Then write it out.
   */


   /* Now write out the color map, 1 byte/color */
/*   printf("\nMAPA...:");
   for (j=0;j<n;j++) 
   {
      b = ColMap[j][2]; //fputc(b,fp2); 
      g = ColMap[j][1]; //fputc(g,fp2);
      r = ColMap[j][0]; //fputc(r,fp2);

      printf("\n[%d] (%d, %d, %d)", j+1, r, g, b);
   }
*/
	
   /*
   ** Finally, we can remap the input image. For each pixel in
   ** the input image, convert pixel's color to a 15-bit color.
   ** The pixel's index in the color map is looked up in "Hist"
   ** which now functions as an inverse color map.
   */



	gettimeofday(&tf_a, NULL);



   /* ###  VOLCADO A DISCO DE LA IMAGEN CUANTIZADA ### */
   volcar_imagen_ppm(nombre_fs, ColMap);
   



   // calculo el error medio y muestro resultados
   tiem_a = (tf_a.tv_sec - ti_a.tv_sec)*1000 +(tf_a.tv_usec - ti_a.tv_usec)/1000;

	
	// muestro el numero de colores de la paleta y el tiempo de ejecución
	printf("%s %s %llu\n", argv[2], argv[1], tiem_a);
	
	//Instrucciones para liberar memoria dinámica
	free(Ig); free(Ib); free(Ir);
}









void getrgb(byte * r,byte * g,byte * b, FILE *fp1)
{
   *b = (byte)getc(fp1); 
   *g = (byte)getc(fp1);
   *r = (byte)getc(fp1);
   return;
}


void WriteTargaHeader(FILE *fp1, struct ImageHeader Image)
{
   /*  
   ** Writes the Targa header information to a file. Note the
   ** two macros below work reliable ONLY for 2-byte variables.
   ** "LSB" returns the lower byte, "MSB" returns the upper byte.
   */
#define LSB(x)  ((x) & 255)
#define MSB(x)  (((x)>>8)&255)
   
   byte        header[18];
   size_t      nbytes;

   header[0]  = (byte)Image.IdSize;
   header[1]  = (byte)Image.ColMapType;
   header[2]  = (byte)Image.TypeCode;
   header[3]  = (byte)LSB(Image.ColMapOrigin);
   header[4]  = (byte)MSB(Image.ColMapOrigin);
   header[5]  = (byte)LSB(Image.ColMapLength);
   header[6]  = (byte)MSB(Image.ColMapLength);
   header[7]  = (byte)Image.ColMapBits;
   header[8]  = (byte)LSB(Image.XOrigin);
   header[9]  = (byte)MSB(Image.XOrigin);
   header[10] = (byte)LSB(Image.YOrigin);
   header[11] = (byte)MSB(Image.YOrigin);
   header[12] = (byte)LSB(Image.Width);
   header[13] = (byte)MSB(Image.Width);
   header[14] = (byte)LSB(Image.Height);
   header[15] = (byte)MSB(Image.Height);
   header[16] = (byte)Image.BitsPerPixel;
   header[17] = (byte)Image.Descriptor;
   nbytes = 18;
   fwrite((void *)header,sizeof(byte),nbytes,fp1);

   /* Write out image identification field, if any */

   if ((nbytes = Image.IdSize) != 0)
      fwrite((void *)Image.IdField,sizeof(byte),nbytes,fp1);
}


void ReadTargaHeader(FILE *fp1, struct ImageHeader *Image)
{
   /*
   ** Reads a Targa file header, and fills "Image" structure.
   */

   byte        header[18];
   size_t      nbytes;

   nbytes = 18;
   if (fread((void *)header,sizeof(byte),nbytes,fp1) != nbytes){
      inform("ERROR: Unexpected EOF\n");
      exit(-1);
   }
   Image->IdSize       = header[0];
   Image->ColMapType   = header[1];
   Image->TypeCode     = header[2];
   Image->ColMapOrigin = header[3] + 256*header[4];
   Image->ColMapLength = header[5] + 256*header[6];
   Image->ColMapBits   = header[7];
   Image->XOrigin      = header[8]  + 256*header[9];
   Image->YOrigin      = header[10] + 256*header[11];
   Image->Width        = header[12] + 256*header[13];
   Image->Height       = header[14] + 256*header[15]; 
   Image->BitsPerPixel = header[16];
   Image->Descriptor   = header[17];

   /* Read in image identification field, if any */

   if ((nbytes = Image->IdSize) != 0)
      fread((void*)Image->IdField,sizeof(byte),nbytes,fp1);
}


void inform(char *mess)
{
   fprintf(stderr,"%s",mess);
}

