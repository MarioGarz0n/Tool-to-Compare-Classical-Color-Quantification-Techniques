/* *********************************************************************
 *         Método Neuquant
 *
 * COMPILAR: 
 *    gcc NQ_PPM.c  neuquant32.c mdmalloc.c ppm_io.c -lm
 *
 * EJECUTAR:
 *   ./a.out <imagen original.ppm> <numColores paleta cuantizada>
 *
 ********************************************************************* */

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <sys/time.h>


#include "neuquant32.h"


//AÑADO PARA PROCESAR IMÁGENES PPM
#include "../archivos_imagenes/header.h" //hay que poner la ruta para llegar



// redefino el tipo que había en el código Neuquant, para desvincularlo del formato PNG
typedef struct _rwpng_color_struct {  
   unsigned char red;
   unsigned char green;
   unsigned char blue;
} rwpng_color;

struct timeval ti_global, tf_global, ti_a, tf_a;       
unsigned long long tiem_a;  

int width, height;


// pixels de la imagen original
unsigned char *row_pointers_MOD = NULL;




/* ---------------------------------------
* Lee la imagen de disco, del fichero indicado como parámetro (un fichero PPM)
*
* ---------------------------------------*/
void leer_datos_PPM(char *nombre_fich)
{
int num_elems;
int i;
Image *org_img;
uchar *org_data;



   org_img = read_img (nombre_fich);

   height = org_img->num_rows; // height es global
   width = org_img->num_cols;  // width es global
        

   // se pasan los datos de la imagen a forma de vector
   /* los datos se guardan en un vector, colocando consecutivos
   los valores R G y B de cada pixel */
   org_data = get_data_ptr ( org_img ); 
       
   num_elems = height* width*3;
   //printf("\n (%d x %d x 3= %d valores)", height, width, num_elems);
    
              
   row_pointers_MOD = (unsigned char*) malloc(sizeof(unsigned char) * num_elems );


   /* se copian las componentes RGB de todos los pixels en el vector
   con el formato que requiere el algoritmo NEUQUANT */
   for(i = 0; i < num_elems; i++)
      row_pointers_MOD[i] = (unsigned char)org_data[i];
   

   free_img ( org_img );   
}




static void remap_simple_PARA_PPM2(unsigned int cols, unsigned int rows, unsigned char map[MAXNETSIZE][3], 
                        unsigned int* remap,  unsigned char *row_pointers, rwpng_color palette[256],
                        char nombre_fich_salida[100])
{
    int color; //numero de un color
    
    unsigned int i,row;

int ik;
Image *out_img;
uchar *out_data;
    

   // para volcar en PPM  	         
   out_img = alloc_img ( rows, cols, MAX_RGB); 
   out_data = get_data_ptr ( out_img );
	 
   ik = 0;
          
    /* Do each image row */
    for(i=0; i < cols*rows*3; i+=3)
    {
       // siguiente pixel de la lista
       unsigned char *punto = &(row_pointers[i]); // (se extiende en 3 posiciones del vector row_pointers)

       color = remap[ inxsearch_MOD( punto[2], punto[1], punto[0] )  ];            

       out_data[ik]     = palette[color].red;  
       out_data[ik + 1] = palette[color].green; 
       out_data[ik + 2] = palette[color].blue;  
          
       ik += 3; 
    }


   write_img ( out_img, nombre_fich_salida );
   free_img(out_img);   
}





/* 
****************************
Program Skeleton
   ----------------
   	** [select samplefac in range 1..30]
   	**pic = (unsigned char*) malloc(4*width*height);
   	**[read image from input file into pic]
	**initnet(pic,4*width*height,samplefac,colors);
	**learn();
	unbiasnet();
	[write output image header, using writecolourmap(f),
	possibly editing the loops in that function]
	**inxbuild();
	[write output image using inxsearch(a,b,g,r)]	

*/
/*
*PARÁMETROS:
*   sample_factor, 
*   rows, cols: filas y columnas de la imagen original 
*   newcolors: número de colores de la paleta cuantizada
*   quantization_gamma
*/
		
void aplicar_neuquant(int sample_factor, int rows, int cols , int newcolors, int verbose, double quantization_gamma,
                       char nombre_img_cuan[100])
{

int bot_idx, top_idx; /* for remapping of indices */
int x;

unsigned char map[MAXNETSIZE][3];  //[4];
unsigned int remap[MAXNETSIZE];

rwpng_color palette[256];	        /* write */  
int i, j;

  
  /* Start neuquant */
  initnet_MOD(//(unsigned char*)row_pointers, //(unsigned char*)rwpng_info.rgba_data, 
            row_pointers_MOD, rows*cols*3, //    rows*cols*4, 
            newcolors, quantization_gamma);            
  //printf("\ninicializado..."); fflush(stdout);
  
              
  learn_MOD(sample_factor,verbose);  
  
  inxbuild_MOD();    
  
  getcolormap_MOD((unsigned char*)map);   
  


   /* Remap indexes so all tRNS chunks are together */  
   // mi versión del if anterior ( desvinculado de la exietencia del canal alfa)
   for (top_idx = newcolors-1, bot_idx = x = 0;  x < newcolors;  ++x)        
         remap[x] = bot_idx++;
    

           
    /* Remap and make palette entries */     
    for (x = 0; x < newcolors; ++x)
    {
       palette[remap[x]].red  = map[x][0];  
       palette[remap[x]].green = map[x][1];
       palette[remap[x]].blue = map[x][2];
    }
  
   remap_simple_PARA_PPM2(cols,rows,map,remap,row_pointers_MOD, palette, nombre_img_cuan);
    
     
   if(row_pointers_MOD != NULL)    
      free(row_pointers_MOD); 
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



int main(int argc, char *argv[]) {

  int sample_factor = 1;   /* will be set depending on image size */
  double force_gamma = 1;  //0.1;
  int n_colours;     /* number of colours to quantize to. Default 256 */
  int verbose = 1;
  char nombre_fich_salida[100];
  char nombre_fich_sin_path[100]; // nombre del fich de entrada sin path ni parte final separada por .
            
                  
  if(argc != 3) 
    abort();
   
  n_colours = atoi(argv[2]); //se guarda el número de colores de la paleta cuantizada
  
  leer_datos_PPM (argv[1]);
  
  /* tomo la parte del fichero de entrada que excluye el path y los trozos
	finales separados por puntos */
  extraer_nombre_fich(argv[1], nombre_fich_sin_path);
  
  gettimeofday(&ti_a, NULL); 
  
  sprintf(nombre_fich_salida, "NQ_%d_%s.ppm", n_colours, nombre_fich_sin_path);
  aplicar_neuquant(sample_factor, height, width, n_colours, verbose, force_gamma, nombre_fich_salida);
  
  gettimeofday(&tf_a, NULL); 
  tiem_a = (tf_a.tv_sec - ti_a.tv_sec)*1000 +(tf_a.tv_usec - ti_a.tv_usec)/1000;
  printf("%s %s %llu\n", argv[2], argv[1], tiem_a);
  
  return 0;
}
