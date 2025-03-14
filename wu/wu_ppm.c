/* **********************************************************************
 *                           wu_ppm.c

Método de Wu aplicado tomando como datos de entrada imágenes en formato PPM.

/*Having received many constructive comments and bug reports about my previous
C implementation of my color quantizer (Graphics Gems vol. II, p. 126-133),
I am posting the following second version of my program (hopefully 100%
healthy) as a reply to all those who are interested in the problem.
*/



/**********************************************************************
	    C Implementation of Wu's Color Quantizer (v. 2)
	    (see Graphics Gems vol. II, pp. 126-133)

Author:	Xiaolin Wu
	Dept. of Computer Science
	Univ. of Western Ontario
	London, Ontario N6A 5B7
	wu@csd.uwo.ca

Algorithm: Greedy orthogonal bipartition of RGB space for variance
	   minimization aided by inclusion-exclusion tricks.
	   For speed no nearest neighbor search is done. Slightly
	   better performance can be expected by more sophisticated
	   but more expensive versions.

The author thanks Tom Lane at Tom_Lane@G.GP.CS.CMU.EDU for much of
additional documentation and a cure to a previous bug.

Free to distribute, comments and suggestions are appreciated.
**********************************************************************/	

#include<stdio.h>

#include <stdlib.h> 
#include <string.h>
#include <math.h>  
#include <time.h>
#include <sys/time.h> 


#define MAXCOLOR	256

#define RED     2
#define GREEN	1
#define BLUE	0




#define N_PUNTOS 512*768 //CELEBI MAYO 2022

/*Defino 2 ctes para las filas y columnas de las fig. cuantizadas
que vuelco a disco */
int FILAS_DISCO, COLS_DISCO;

int    N;                 // número de pixels reales de una imagen




//AÑADO PARA PROCESAR IMÁGENES PPM
#include "../archivos_imagenes/header.h" //hay que poner la ruta para llegar


#define verboso 0 //1: mostrar mensajes; 0: no mostrarlos 

struct box {
    int r0;			 /* min value, exclusive */
    int r1;			 /* max value, inclusive */
    int g0;  
    int g1;  
    int b0;  
    int b1;
    int vol;
};


/* --------------------------------------------------------------
 * Histogram is in elements 1..HISTSIZE along each axis,
 * element 0 is for base or marginal value
 * NB: these must start out 0!
 -------------------------------------------------------------- */
float		m2[33][33][33];
long int	wt[33][33][33], 
            mr[33][33][33],	mg[33][33][33],	mb[33][33][33];
unsigned char   *Ir, *Ig, *Ib; 

int	   size;             /* image size. Número de pixels de la imagen */
int		K;                /* color look-up table size*/
unsigned short int *Qadd;




void Mark( struct box *cube, int label, unsigned char *tag);




/* ----------------------------------------------------------------
  At conclusion of the histogram step, we can interpret
 *   wt[r][g][b] = sum over voxel of P(c)
 *   mr[r][g][b] = sum over voxel of r*P(c)  ,  similarly for mg, mb
 *   m2[r][g][b] = sum over voxel of c^2*P(c)
 * Actually each of these should be divided by 'size' to give the usual
 * interpretation of P() as ranging from 0 to 1, but we needn't do that here.
 * ---------------------------------------------------------------- */
void Hist3d(vwt, vmr, vmg, vmb, m2) 
/* build 3-D color histogram of counts, r/g/b, c^2 */
long int *vwt, *vmr, *vmg, *vmb;
float	*m2;
{
register int ind, 
	          r, g, b; 
int	     inr, ing, inb, 
	        table[256]; 
register long int i;

	
	for(i=0; i<256; ++i)
		table[i]=i*i;

	Qadd = (unsigned short int *)malloc(sizeof(short int)*size);
	if (Qadd==NULL)
	{
		printf("Not enough space\n");
		exit(1);
		} 

	if(verboso == 1)
	   printf("\nsize vale ...%d\n", size); fflush(stdout);

	for(i=0; i<size; ++i)
	{
	    r = Ir[i]; g = Ig[i]; b = Ib[i];

	    inr=(r>>3)+1; 
	    ing=(g>>3)+1; 
	    inb=(b>>3)+1; 

 	    Qadd[i]= ind = (inr<<10)+(inr<<6)+inr+(ing<<5)+ing+inb;
		
	    ++vwt[ind];

	    vmr[ind] += r;
	    vmg[ind] += g;
	    vmb[ind] += b;
    	m2[ind] += (float)(table[r]+table[g]+table[b]);
	}
}




/* ----------------------------------------------------------------
 * We now convert histogram into moments so that we can rapidly calculate
 * the sums of the above quantities over any desired box.
 ---------------------------------------------------------------- */
void M3d(vwt, vmr, vmg, vmb, m2) /* compute cumulative moments. */
long int *vwt, *vmr, *vmg, *vmb;
float	*m2;
{
register unsigned short int ind1, ind2;
register unsigned char i, r, g, b;
long int line, line_r, line_g, line_b,
	 area[33], 
	area_r[33], 
	area_g[33], 
	area_b[33];
float    line2, area2[33];

	
   for(r=1; r<=32; ++r)
	{
	   for(i=0; i<=32; ++i) 
	      area2[i] = area[i] = area_r[i] = area_g[i] = area_b[i] = 0;
		
	   for(g=1; g<=32; ++g)
	   {
	       line2 = line = line_r = line_g = line_b = 0;
		   
	       for(b=1; b<=32; ++b)
		   {
		      ind1 = (r<<10) + (r<<6) + r + (g<<5) + g + b; /* [r][g][b] */
			   
		      line   += vwt[ind1];
		      line_r += vmr[ind1]; 
		      line_g += vmg[ind1]; 
		      line_b += vmb[ind1];
			   
		      line2 += m2[ind1];
			   
		      area[b]   += line;
		      area_r[b] += line_r;
		      area_g[b] += line_g;
		      area_b[b] += line_b;
			   
		      area2[b] += line2;
		      ind2 = ind1 - 1089; /* [r-1][g][b] */
				
		      vwt[ind1] = vwt[ind2] + area[b];
		      vmr[ind1] = vmr[ind2] + area_r[b];
		      vmg[ind1] = vmg[ind2] + area_g[b];
		      vmb[ind1] = vmb[ind2] + area_b[b];
				
		      m2[ind1] = m2[ind2] + area2[b];
	       }
	   }
    }
}


long int Vol(cube, mmt) 
/* Compute sum over a box of any given statistic */
struct box *cube;
long int mmt[33][33][33];
{
    return( mmt[cube->r1][cube->g1][cube->b1] - mmt[cube->r1][cube->g1][cube->b0]
	       -mmt[cube->r1][cube->g0][cube->b1] + mmt[cube->r1][cube->g0][cube->b0]
	       -mmt[cube->r0][cube->g1][cube->b1] + mmt[cube->r0][cube->g1][cube->b0]
	       +mmt[cube->r0][cube->g0][cube->b1] - mmt[cube->r0][cube->g0][cube->b0] );
}


/* ----------------------------------------------------------------
 The next two routines allow a slightly more efficient calculation
 * of Vol() for a proposed subbox of a given box.  The sum of Top()
 * and Bottom() is the Vol() of a subbox split in the given direction
 * and with the specified new upper bound.
 ---------------------------------------------------------------- */

long int Bottom(cube, dir, mmt)
/* Compute part of Vol(cube, mmt) that doesn't depend on r1, g1, or b1 */
/* (depending on dir) */
struct box *cube;
unsigned char dir;
long int mmt[33][33][33];
{
   switch(dir)
	{
	   case RED:   return( - mmt[cube->r0][cube->g1][cube->b1] 
	                       + mmt[cube->r0][cube->g1][cube->b0]
		                    + mmt[cube->r0][cube->g0][cube->b1] 
	                       - mmt[cube->r0][cube->g0][cube->b0] );
	      break;
			
	   case GREEN:   return( - mmt[cube->r1][cube->g0][cube->b1]
		                      + mmt[cube->r1][cube->g0][cube->b0]
		                      + mmt[cube->r0][cube->g0][cube->b1]
		                      - mmt[cube->r0][cube->g0][cube->b0] );
	      break;
			
	   case BLUE:     return( - mmt[cube->r1][cube->g1][cube->b0]
		                       + mmt[cube->r1][cube->g0][cube->b0]
		                       + mmt[cube->r0][cube->g1][cube->b0]
		                       - mmt[cube->r0][cube->g0][cube->b0] );
	      break;
    }
}


long int Top(cube, dir, pos, mmt)
/* Compute remainder of Vol(cube, mmt), substituting pos for */
/* r1, g1, or b1 (depending on dir) */
struct box *cube;
unsigned char dir;
int   pos;
long int mmt[33][33][33];
{
    switch(dir)
	{
	 case RED:
	    return(  mmt[pos][cube->g1][cube->b1] 
		        - mmt[pos][cube->g1][cube->b0]
		        - mmt[pos][cube->g0][cube->b1]
		        + mmt[pos][cube->g0][cube->b0] );
	    break;
	 case GREEN:
	    return( mmt[cube->r1][pos][cube->b1] 
		   -mmt[cube->r1][pos][cube->b0]
		   -mmt[cube->r0][pos][cube->b1]
		   +mmt[cube->r0][pos][cube->b0] );
	    break;
	 case BLUE:
	    return( mmt[cube->r1][cube->g1][pos]
		   -mmt[cube->r1][cube->g0][pos]
		   -mmt[cube->r0][cube->g1][pos]
		   +mmt[cube->r0][cube->g0][pos] );
	    break;
    }
}


float Var(cube)
/* Compute the weighted variance of a box */
/* NB: as with the raw statistics, this is really the variance * size */
struct box *cube;
{
float dr, dg, db, xx;

    dr = Vol(cube, mr); 
    dg = Vol(cube, mg); 
    db = Vol(cube, mb);
	
    xx =  m2[cube->r1][cube->g1][cube->b1] 
	 -m2[cube->r1][cube->g1][cube->b0]
	 -m2[cube->r1][cube->g0][cube->b1]
	 +m2[cube->r1][cube->g0][cube->b0]
	 -m2[cube->r0][cube->g1][cube->b1]
	 +m2[cube->r0][cube->g1][cube->b0]
	 +m2[cube->r0][cube->g0][cube->b1]
	 -m2[cube->r0][cube->g0][cube->b0];

    return( xx - (dr*dr+dg*dg+db*db)/(float)Vol(cube,wt) );    
}

/* We want to minimize the sum of the variances of two subboxes.
 * The sum(c^2) terms can be ignored since their sum over both subboxes
 * is the same (the sum for the whole box) no matter where we split.
 * The remaining terms have a minus sign in the variance formula,
 * so we drop the minus sign and MAXIMIZE the sum of the two terms.
 */


float Maximize(cube, dir, first, last, cut, whole_r, whole_g, whole_b, whole_w)
struct box *cube;
unsigned char dir;
int first, last, *cut;
long int whole_r, whole_g, whole_b, whole_w;
{
register long int half_r, half_g, half_b, half_w;
long int          base_r, base_g, base_b, base_w;
register int i;
register float temp, max;

	
    base_r = Bottom(cube, dir, mr);
    base_g = Bottom(cube, dir, mg);
    base_b = Bottom(cube, dir, mb);
    base_w = Bottom(cube, dir, wt);
	
    max = 0.0;
    *cut = -1;
	
    for(i=first; i<last; ++i)
	{
	   half_r = base_r + Top(cube, dir, i, mr);
	   half_g = base_g + Top(cube, dir, i, mg);
	   half_b = base_b + Top(cube, dir, i, mb);
	   half_w = base_w + Top(cube, dir, i, wt);
		
      /* now half_x is sum over lower half of box, if split at i */
      if (half_w == 0)       /* subbox could be empty of pixels! */
		{
			continue;             /* never split into an empty box */
	   } else
          temp = ( (float)half_r*half_r + (float)half_g*half_g +
                   (float)half_b*half_b)/half_w;

	   half_r = whole_r - half_r;
	   half_g = whole_g - half_g;
	   half_b = whole_b - half_b;
	   half_w = whole_w - half_w;
		
      if (half_w == 0)       /* subbox could be empty of pixels! */
		{
		   continue;             /* never split into an empty box */
	   } else
          temp += ((float)half_r*half_r + (float)half_g*half_g +
                   (float)half_b*half_b)/half_w;

      if (temp > max) 
		{
			max=temp; 
			*cut=i;
		}
    }
	
    return(max);
}


int Cut(set1, set2)
struct box *set1, *set2;
{
unsigned char dir;
int      cutr, cutg, cutb;
float    maxr, maxg, maxb;
long int whole_r, whole_g, whole_b, whole_w;


	 
    whole_r = Vol(set1, mr);
    whole_g = Vol(set1, mg);
    whole_b = Vol(set1, mb);
    whole_w = Vol(set1, wt);

	
    maxr = Maximize(set1, RED, set1->r0+1, set1->r1, &cutr,
		                whole_r, whole_g, whole_b, whole_w);
    maxg = Maximize(set1, GREEN, set1->g0+1, set1->g1, &cutg,
		            whole_r, whole_g, whole_b, whole_w);
    maxb = Maximize(set1, BLUE, set1->b0+1, set1->b1, &cutb,
		            whole_r, whole_g, whole_b, whole_w);


   if( (maxr>=maxg) && (maxr>=maxb) ) 
	{
	   dir = RED;
		
	   if (cutr < 0) 
			return 0; /* can't split the box */
    }
    else

       if( (maxg>=maxr) && (maxg>=maxb) ) 
	      dir = GREEN;

       else
	      dir = BLUE; 


    set2->r1 = set1->r1;
    set2->g1 = set1->g1;
    set2->b1 = set1->b1;


   switch (dir)
	{
	  case RED:
	     set2->r0 = set1->r1 = cutr;
	     set2->g0 = set1->g0;
	     set2->b0 = set1->b0;
	     break;
			
	  case GREEN:
	     set2->g0 = set1->g1 = cutg;
	     set2->r0 = set1->r0;
	     set2->b0 = set1->b0;
	     break;
			
	  case BLUE:
	     set2->b0 = set1->b1 = cutb;
	     set2->r0 = set1->r0;
	     set2->g0 = set1->g0;
	     break;
    }
	
    set1->vol= (set1->r1-set1->r0) * (set1->g1-set1->g0) * (set1->b1-set1->b0);
    set2->vol= (set2->r1-set2->r0) * (set2->g1-set2->g0) * (set2->b1-set2->b0);
	
    return 1;
}


void Mark(cube, label, tag)
struct box *cube;
int label;
unsigned char *tag;
{
register int r, g, b;

    for(r=cube->r0+1; r<=cube->r1; ++r)
       for(g=cube->g0+1; g<=cube->g1; ++g)
	      for(b=cube->b0+1; b<=cube->b1; ++b)
	         tag[(r<<10) + (r<<6) + r + (g<<5) + g + b] = label;
}






/* ************************************************ 
 * Lee los parámetros pasados en la linea de comandos. 
 * ************************************************ */
void leer_parametros_ejecucion_CELEBI(int argc, char *argv[], char nombre_fich[])
{
	// formato de la orden adecuado para pasar todos los parámetros al programa
       //printf("\n  FORMATO LLAMADA: ./<ejecutable> <fich_pixels> <colores>\n");
	
	/* --- primer argumento (obligatorio): nombre del fichero de datos.
	 Si no se indica, el programa acaba */
	if(argc >= 2)
		strcpy(nombre_fich, argv[1]);
	else
	{
		printf("Debe indicar el nombre del fichero que contiene los valores RGB de los pixels.");
		exit(1);
	}
	
	/* --- segundo argumento (obligatorio): número de colores de la paleta para
	 el algoritmo de Wu. Si no se indica o se da un valor mayor que 256, el
	 programa acaba */
	if(argc >= 3) 
	{
		K= atoi(argv[2]);
		
		if(K<1 || K>256)
		{
			printf("Tamaño paleta para alg. Wu: debe estar en el intervalo [1, 256].");
			exit(1);
		}
	}
	else
	{
		printf("Debe indicar el numero de colores para Wu.");
		exit(1);
	}
}


/* --------------------------------------------------
 * Tomo el nombre del fichero de datos y le quito tanto el path como la extensión
 * -------------------------------------------------- */
void extraer_nombre_fich(char nombre_fich[], char nombre_solo[])
{
int i, j, k;
int seguir =1;

	i=strlen(nombre_fich);



	while(seguir == 1)
	{
		if( nombre_fich[i]=='/') 
			seguir=0;
                else
                 {
                 if(i ==0)
                 { i--;
                   seguir =0; 
                 }else
                   i--;

		
		if (i<0) 
		{
		   break;
		}
		}
	}


	k=0;
	for(j=i+1; j<strlen(nombre_fich); j++, k++)  //for(j=i+1; j<strlen(nombre_fich); j++, k++)
	{
		nombre_solo[k] = nombre_fich[j];
	    if (nombre_fich[j] == '\0')
	    { j++; k++;
	      break;
	      }

        }
        nombre_solo[k] ='\0';
        
return;
	//printf("\n|%s|-> %d elementos", nombre_solo, strlen(nombre_solo));

	//recorro la cadena hacia atras para borrar el .txt
	//for(j=strlen(nombre_solo)-1; j>=0; j--)
	seguir=1;

	j=strlen(nombre_solo)-1;
	while( (seguir == 1) || (j <= 0) )
	{
		if(nombre_solo[j] == '.')
		{
			nombre_solo[j]='\0';
			seguir = 0;
		}
		j--;
	}
}


// >>>>>>>>>>>>>>>>>>>>>>>>>>> MAIN >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
void main(int argc, char *argv[])
{
struct box	        cube[MAXCOLOR];
unsigned char	    *tag;
register long int	i, weight;
register int	    k;
int		          next;
float		          vv[MAXCOLOR], temp;

// colores de la paleta cuantizada
unsigned char	lut_r[MAXCOLOR], lut_g[MAXCOLOR], lut_b[MAXCOLOR];

	
char nombre_fich[50]="";  /* nombre del fichero de datos (pixels) */

char nombre_fich_salida[100];
char nombre_solo[80]; //nombre del fichero de entrada, sin ruta



// para calcular tiempo de ejecución en milisegundos
struct timeval ti_a, tf_a;          // instante inicial/final de ejecución del algoritmo 
unsigned long long tiem_a, tiem_m;  // diferencia entre los dos tiempos previos


struct timeval ti_a_extra_WU, tf_a_extra_WU; 
   unsigned long long tiem_a_extra_WU;




// NUEVAS VARIABLES PARA LEER Y EXCRIBIR IMÁGENES ppm
Image *org_img;
int num_rows, num_cols;
int num_elems, num_pixels ;
uchar *org_data;
int index; // para copiar datos del formato PPM al formato que usa Wu
	  Image *out_img;
	   uchar *out_data;
int ik;	   
	   
	   
	
   /* #############################################
      ### LEE INFORMACIÓN DE ENTRADA 
    * ############################################# */ 
	leer_parametros_ejecucion_CELEBI(argc, argv, nombre_fich);


 /* extraigo el nombre del fichero de datos, para usarlo luego al generar 
   el nombre del fichero que almacenará la imagen cuantizada */
   extraer_nombre_fich(nombre_fich, nombre_solo);


	/* +++++++++++++++++++++++++++++++++++++++++++++
	AÑADO EL CODIGO QUE LEE UN FICHERO PPM
	+++++++++++++++++++++++++++++++++++++++++++++ */
        // se lee la imagen original (un fichero PPM)
        org_img = read_img (nombre_fich);


        FILAS_DISCO = org_img->num_rows;
        COLS_DISCO = org_img->num_cols;
        
       num_pixels = FILAS_DISCO*COLS_DISCO;
        

        num_elems = num_pixels * 3;

        // se pasan los datos de la imagen a forma de vector
       /* los datos se guardan en un vector, colocando consecutivos
       los valores R G y B de cada pixel */
       org_data = get_data_ptr ( org_img ); 
       
       

	
	/* -2- lee los valores RGB de los puntos de la imagen y los almacena en
	 los vectores Ir, Ig, Ib, que son globales y se les reserva espacio
	 dinámicamente en la función de lectura */
    
size =num_pixels;




   Ir = (unsigned char *) malloc( sizeof(unsigned char)*num_pixels); // para rojo
   Ig = (unsigned char *) malloc( sizeof(unsigned char)*num_pixels); // para verde
   Ib = (unsigned char *) malloc( sizeof(unsigned char)*num_pixels); // para azul
   
   if( (Ir == NULL) || (Ig == NULL) || (Ib == NULL) )
   {
	  printf("ERROR al reservar espacio para los colores. FIN."); 
	  exit(1);
	  }  
 

    /* sustituyo la llamada a la función por el bucle que copia los datos
    de la imagen de los vectores globales que usa el método de Wu (que
    son unsigned char) */		    
    index = 0;
 for (ik = 0; ik < num_elems; ik += 3 )
  {

   Ir[index] = (unsigned char)org_data[ik];
   Ig[index] = (unsigned char)org_data[ik + 1];
   Ib[index] = (unsigned char)org_data[ik + 2];
   
        
   index++;
  }


free_img ( org_img );   



   /* #############################################
      ### -A- OPERACIONES DEL MÉTODO DE WU 
      ############################################# */  
	
   gettimeofday(&ti_a, NULL);   // tiempo de comienzo de las operaciones
   
	Hist3d(wt, mr, mg, mb, m2); 
	
	M3d(wt, mr, mg, mb, m2);
	//printf("Moments done\n");

	cube[0].r0 = cube[0].g0 = cube[0].b0 = 0;
	cube[0].r1 = cube[0].g1 = cube[0].b1 = 32;
	next = 0;

	// repetir... (K es el número de colores de la paleta cuantizada)
   for(i=1; i<K; ++i)
	{
       if (Cut(&cube[next], &cube[i])) 
	   {
              /* volume test ensures we won't try to cut one-cell box */
              vv[next] = (cube[next].vol>1) ? Var(&cube[next]) : 0.0;
              vv[i] = (cube[i].vol>1) ? Var(&cube[i]) : 0.0;
	   }
	   else
	   {
              vv[next] = 0.0;   /* don't try to split this box again */
              i--;              /* didn't create box i */
	   }

		
      next = 0;     
      temp = vv[0]; 

      for(k=1; k<=i; ++k)
         if (vv[k] > temp) 
		   {
            temp = vv[k];  
				next = k;      
		   }


	   if (temp <= 0.0)
	   {
              K = i+1;
              //fprintf(stderr, "Only got %d boxes\n", K);
              break;
	   }
	}
	
    // printf("Partition done\n");

	/* the space for array m2 can be freed now */

	tag = (unsigned char *)malloc(33*33*33);
	if (tag==NULL)
	{
		printf("Not enough space\n"); 
		exit(1);
	}

	for(k=0; k<K; ++k)
	{
	    Mark(&cube[k], k, tag);
	    weight = Vol(&cube[k], wt);
		
	    if (weight) 
		{
			lut_r[k] = Vol(&cube[k], mr) / weight;
			lut_g[k] = Vol(&cube[k], mg) / weight;
			lut_b[k] = Vol(&cube[k], mb) / weight;
	    }
	    else
		{
	      fprintf(stderr, "bogus box %d\n", k);
	      lut_r[k] = lut_g[k] = lut_b[k] = 0;		
	    }
	}

	for(i=0; i<size; ++i)
	{
		Qadd[i] = tag[Qadd[i]];
	}



   gettimeofday(&tf_a, NULL);   // instante final de ejecución de WU






	  
	  
	  
	  // +++++++++++++++++++++++++++++++
	  // +++++++++++++++++++++++++++++++
	  // >>>>>>>>>>>>> VUELCO IMAGEN PPM generada por el método de Wu 
   
	  out_img = alloc_img ( FILAS_DISCO, COLS_DISCO, MAX_RGB );
	  out_data = get_data_ptr ( out_img );
	 

	  index = 0;
          for ( ik = 0; ik < num_elems; ik += 3 )
         {
            out_data[ik] = lut_r[Qadd[index]];
            out_data[ik + 1] = lut_g[Qadd[index]];
            out_data[ik + 2] = lut_b[Qadd[index]];            
            index++;
        }

 

        sprintf(nombre_fich_salida, "WU_%s_%s", argv[2], nombre_solo);
        write_img ( out_img, nombre_fich_salida );
 
        free_img(out_img);
 
	  
	  
	  
   // se calcula el tiempo de volcado de la imagen
   tiem_a_extra_WU = (tf_a_extra_WU.tv_sec - ti_a_extra_WU.tv_sec)*1000 +
                     (tf_a_extra_WU.tv_usec - ti_a_extra_WU.tv_usec)/1000;


 
   // se calcula y muestra el tiempo de ejecución
   tiem_a = (tf_a.tv_sec - ti_a.tv_sec)*1000 +(tf_a.tv_usec - ti_a.tv_usec)/1000;

	printf("%s %s %llu\n", argv[2], argv[1], tiem_a);

	//Instrucciones para liberar memoria dinámica
	free(Ig); free(Ib); free(Ir);
}
