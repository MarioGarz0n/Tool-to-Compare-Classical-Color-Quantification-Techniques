/* **********************************************************
*                        kmeans_imagenes.c 
*
* **********************************************************
*
* COMPILACION:
*	 gcc kmeans_imagenes.c -lm -o kmeans 
*
* EJECUCION:  
*	 ./kmeans <imagen_original.ppm> <numero_colores> 
* 
********************************************************** */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <sys/time.h>
#include <malloc.h>
#include <math.h>
#include <float.h> 

//PARA PROCESAR IMÁGENES PPM
#include "../archivos_imagenes/header.h" //hay que poner la ruta para llegar

// ------------- CONSTANTES -----------------------
// mensajes de error => tras presentarlos el programa acaba
#define ERR_PEDIR_MEM  "ERROR: no hay memoria disposible"
#define ERR_FOPEN      "ERROR al abrir el fichero"
#define ERR_L_FICH     "ERROR al leer del fichero"
#define ERR_NUM_DATOS  "ERROR: numero de datos leidos distinto del esperado"


#define MAX_PASADAS 10 //25 //20 //3 //10 // número máximo de iteraciones del algoritmo
                       
#define MAX_TEST 10    // numero de pruebas independientes

#define ERROR 0.01 //0.1 // 1e-4
// ------------- VARIABLES -----------------------
int N;      // numero real de documentos (= numero de hormigas)
int R=3;    // numero real de terminos considerados por documento
float **d=NULL;  // valores RGB de los pixels de la imagen (1 pixel por fila)
float **cuant=NULL; //valores RGB de los pixels de la imagen cuantizada (1 pixel por fila)

char nombre_fich_salida[200];    /* fichero que contiene los pixels de la imagen cuantizada
                                 (resultado del algoritmo) */
int n_clusters_quiero= 32; //256; //32 64 128 256;

// ------------- PROTOTIPOS -----------------------
void aviso_y_fin(char *cad_ini, char *cad_fin);
float **RESERVAR_MAT_REAL(int n_filas, int n_col);
void generar_imagen(float **imagen,int num_filas, int num_cols);




/* ############################################
 *     operaciones sobre ficheros
 * ############################################ */
/* -------------------------------------------------------------
* Abre un fichero en el modo indicado. En caso de error el programa termina.
-------------------------------------------------------------- */
FILE * Abrir_fichero(char *nombre_fich, char *modo)
{
FILE *fp;   

	
   fp = fopen(nombre_fich, modo);

   if (fp == NULL)
   {
      printf("\n(%s %s, modo: %s\n)\n", ERR_FOPEN, nombre_fich, modo);
      exit(-2);
   }

   return fp;
}


/* --------------------------------------------------------------
* Lee de un fichero una matriz de reales.
* La primera linea del fichero contiene dos enteros que indican las filas y las
* columnas, respectivamente. El resto de las lineas contienen los datos de las
* sucesivas filas de la matriz
* -------------------------------------------------------------- */
int LEER_M_REAL(char *nombre_fich, float ***mat, int *n_f, int *n_c)
{
FILE *fd;          
int i, j,          
    n_leidos = 0;  


   fd = Abrir_fichero(nombre_fich, "rt");

	if( fscanf(fd,"%d %d", n_f, n_c) != 2)
		aviso_y_fin("Tamaño del array", ERR_L_FICH);

	*mat = RESERVAR_MAT_REAL(*n_f, *n_c);

   for(i=0; i<(*n_f); i++)
      for(j=0; j<(*n_c); j++)

         if(!feof(fd))
         {
            fscanf(fd,"%f", &(*mat)[i][j]); 
            n_leidos++;
         }

   fclose(fd);  


   if (n_leidos != (*n_f)*(*n_c))
   {
      printf ("\n%s (%d de %d)", ERR_NUM_DATOS, n_leidos, (*n_f)*(*n_c));
      exit(-1);	
   }

   return n_leidos;  
}



/* ############################################
 *     operaciones sobre matrices
 * ############################################ */

/* -----------------------------------------------------------
* Crea una matriz dinamica de reales con n_filas filas y n_col columnas.
* Si se produce un error, el programa acaba.
------------------------------------------------------------- */
float **RESERVAR_MAT_REAL(int n_filas, int n_col)
{
float **p;  
int i;       

	

	p= (float **) calloc(n_filas,sizeof(float *));

	if (p == NULL)
		aviso_y_fin("Matriz de reales (punteros a filas)", ERR_PEDIR_MEM);


	for(i=0; i<n_filas; i++)
	{
	     p[i] = (float *) calloc (n_col, sizeof(float));

     	     if (p[i] == NULL)
		   aviso_y_fin("Matriz de reales (filas)", ERR_PEDIR_MEM);
	}


	return p;
}



/* -----------------------------------------------------------
* Libera la memoria ocupada por una matriz dinamica de reales
----------------------------------------------------------- */
void LIBERAR_MAT_REAL(float ***p, int n_filas)
{
int i; 


	for(i=0; i<n_filas; i++)
		free((*p)[i]); 


	free(*p);
	*p = NULL; 
}



/* ############################################
 *     operaciones miscelaneas
 * ############################################ */

/* -------------------------------------------------------------
 * Muestra un aviso al usuario y fuerza el final del programa
 ------------------------------------------------------------- */
void aviso_y_fin(char *cad_ini, char *cad_fin)
{
	printf("\n\n %s -- %s\n", cad_ini, cad_fin);
	exit(1);
}





/* ############################################
 *     operaciones relacionadas con k-medias
 * ############################################ */

/* ---------------------------------------------------------------
 * Copia de k_means en la que paso yo los valores de los centroides.
 *
 PARAMETROS:
 *   data: array que contiene los items a clasificar
 *   n_items: numero de items a clasificar
 *   m: número de claves consideradas para definir los items a clasificar
 *   k: numero de clusters que se van a generar en la partición
 *   t: error permitido
 *   centroids: 
 *   centroidsXXX: centroides iniciales que paso para la ejecución de k-means
 *    es_alea: 1 => es aleatorio (tomar los centroides pasados como parámetros
               0 => nolo es (tomas los de dentro de esta función
 ---------------------------------------------------------------- */
int *k_means_MOD_X(float **data, int n_items, int m, int k, float t, 
   float **centroids,  float **centroidsXXX, int es_alea)
{
int h, i, j;                      
float old_error, error = DBL_MAX; /* sum of squared euclidean distance */

struct timeval ti_global, tf_global, ti_a, tf_a;       
unsigned long long tiem_a, tiem_a0;  


/* output cluster label for each data point */
int *labels = (int*)calloc(n_items, sizeof(int));
int *counts = (int*)calloc(k, sizeof(int)); 
		
float **c =   centroids ? centroids : (float**)calloc(k, sizeof(float*));  
float **c1 =  (float**)calloc(k, sizeof(float*)); 

int pasada =0; 
double r1, r2, r3;


int contador =0;	
int t_total=0;


	
   gettimeofday(&ti_global, NULL);	
   gettimeofday(&ti_a, NULL);   



 
   /****
   ** initialization. 
   Se definen los centroides iniciales */
   for (h = i = 0; i < k; h += n_items / k, i++) 
   {
      c1[i] = (float*)calloc(m, sizeof(float));

      if (!centroids) 
         c[i] = (float*)calloc(m, sizeof(float));

      for (j = m; j-- > 0; c[i][j] = data[h][j]);
   }


   /* Redefino los centroides de los bucles previos, copiando los que he
   pasado como parámetro. Dejo esos bucles porque reservan memoria  */
   if(es_alea == 1)
   {
      for(i=0; i<k; i++) 
	for(j=0; j<m; j++)
		c[i][j]=centroidsXXX[i][j];
   }
	


   gettimeofday(&tf_a, NULL);   
   tiem_a0 = (tf_a.tv_sec - ti_a.tv_sec)*1000 +(tf_a.tv_usec - ti_a.tv_usec)/1000;
   t_total += tiem_a0;


	
   /****
   ** bucle principal */
   do {
	gettimeofday(&ti_a, NULL); 
		
	pasada++;
		
        if (pasada>MAX_PASADAS)
		   break;



      /* save error from last step */
      old_error = error, error = 0;

      /* clear old counts and temp centroids */
      for (i = 0; i < k; counts[i++] = 0) 
      {
         for (j = 0; j < m; c1[i][j++] = 0);
      }


      for (h = 0; h < n_items; h++) 
      {
 			
         float min_distance = DBL_MAX; 

         for (i = 0; i < k; i++) 
         {
            float distance = 0;
            for (j = m; j-- > 0; distance += pow(data[h][j] - c[i][j], 2));

            if (distance < min_distance) 
            {
               labels[h] = i;
               min_distance = distance;
            }
         }


         /* update size and temp centroid of the destination cluster */
         for (j = m; j-- > 0; c1[labels[h]][j] += data[h][j]);
         counts[labels[h]]++;
         
         /* update standard error */
         error += min_distance;
      }



      for (i = 0; i < k; i++)   
      { /* update all centroids */
         for (j = 0; j < m; j++)  
         {
            c[i][j] = counts[i] ? c1[i][j] / counts[i] : c1[i][j];
         }
      }


      gettimeofday(&tf_a, NULL);   
      tiem_a0 = (tf_a.tv_sec - ti_a.tv_sec)*1000 +(tf_a.tv_usec - ti_a.tv_usec)/1000;
      t_total+=tiem_a0; 



		
     
		
   } while (fabs(error - old_error) > t);



    centroids = c;


    gettimeofday(&tf_global, NULL);
    gettimeofday(&tf_a, NULL);   

    tiem_a0 = (tf_global.tv_sec - ti_global.tv_sec)*1000 +(tf_global.tv_usec - ti_global.tv_usec)/1000;

    gettimeofday(&ti_a, NULL);   



   // ---VUELCO A DISCO ---


   for(i=0; i<n_items; i++) 
   {
	   // sustituyo por un bucle que escriba de 3 en 3, por si es más rapido
	   r1 = d[i][0] - (int)centroids[labels[i]][0] ;
	   r2 = d[i][1] - (int)centroids[labels[i]][1] ;
	   r3 = d[i][2] - (int)centroids[labels[i]][2] ;
	   
	  
	  //se guarda en cuant los valores de la nueva imagen cuantizada
	  cuant[i][0] = (int)centroids[labels[i]][0] ;
	  cuant[i][1] = (int)centroids[labels[i]][1] ;
	  cuant[i][2] = (int)centroids[labels[i]][2] ;
   }


	
	
   //printf("Tiempo KMeans (milisegundos): %llu\n", tiem_a0); 


   gettimeofday(&tf_a, NULL);   
   tiem_a = (tf_a.tv_sec - ti_a.tv_sec)*1000 +(tf_a.tv_usec - ti_a.tv_usec)/1000;

	
	
	
   /****
   ** housekeeping */
   for (i = 0; i < k; i++) 
   {
      if (!centroids) 
         free(c[i]);

		free(c1[i]);
   }

   if (!centroids) 
      free(c);

   free(c1);
   free(counts);

   return labels;
}



/* -------------------------------------------------------------------
 * llama al algoritmo k-means pasándole como centroides iniciales varios
 * puntos elegidos al azar
 --------------------------------------------------------------------- */
void k_medias_aleatorio_X(int n_items, int n_clusters, int n_claves, int X,int num_filas, int num_cols)
{
float **centroids;      /* centroides para pasar a k-medias */	
int *candidato =NULL;    /* indices de los items que se tomarán como centroides */
  int i, k;
int alea;



    // --1-- DEFINIR CENTROIDES INICIALES PARA K-MEANS
    centroids = RESERVAR_MAT_REAL(n_clusters, n_claves);

   for(i=0; i<n_clusters; i++)
   {
  	alea = rand() % n_items;
        for (k = n_claves; k-- > 0; centroids[i][k] = d[alea][k]);
   }

	
	
   // --2-- LLAMAR A K-MEANS
   int *c = k_means_MOD_X(d, n_items, n_claves, n_clusters, ERROR, NULL, centroids, 1); 
   
   // --3-- LIBERAR MEMORIA DINÁMICA
   free(c);  
	
   LIBERAR_MAT_REAL(&centroids, n_clusters); 
	
   if(candidato != NULL)
		free(candidato);
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


void main(int argc, char *argv[])
{
char nombre_fich[200];           /* fichero que contiene los pixels de la imagen original */
time_t t;                        // para iniciar el generador de números aleatorios

// para calcular el tiempo empleado en los cálculos
// --- para calcular tiempo en milisegundos
struct timeval ti_a, tf_a;   // para la aplicación del algoritmo de hormigas
unsigned long long tiem_a;   // para almacenar la diferencia entre los dos tiempos previos


struct tm *instante_actual;
time_t momento;
char cad[42];    // para el nombre del fichero de salida
int pruebas;     // contador de ejecuciones independientes realizadas

// NUEVAS VARIABLES PARA LEER Y EXCRIBIR IMÁGENES ppm
Image *org_img;	
int num_rows, num_cols;
int num_elems, num_pixels;
int ik,index;
uchar *org_data;
int FILAS_DISCO, COLS_DISCO;
Image *out_img;
uchar *out_data;
char nombre_fich_sin_path[100]; // nombre del fich de entrada sin path ni parte final separada por .

   srand((unsigned) time(&t));


   time(&momento);


   instante_actual = localtime(&momento);

   strftime(cad, 42,"%d%b%y_%H%M%S.txt", instante_actual);
	
   // -------------------------------------------------------
   // -- SE DETERMINAN LOS VALORES A UTILIZAR PARA LOS CÁLCULOS --
   // -------------------------------------------------------	
	if(argc >=2){
		strcpy(nombre_fich, argv[1]);
		n_clusters_quiero= atoi(argv[2]);
	}
	else
	{
		printf("\nDebe indicar el nombre del fichero que contieneos pixels");
		exit(1);
	}
	
	/* tomo la parte del fichero de entrada que excluye el path y los trozos
	finales separados por puntos */
	extraer_nombre_fich(nombre_fich, nombre_fich_sin_path);
	
   /* +++++++++++++++++++++++++++++++++++++++++++++
	AÑADO EL CODIGO QUE LEE UN FICHERO PPM
	+++++++++++++++++++++++++++++++++++++++++++++ */
        // se lee la imagen original (un fichero PPM)
        org_img = read_img (nombre_fich);

        
        FILAS_DISCO = org_img->num_rows; //numero de filas
        COLS_DISCO = org_img->num_cols; //numero de columnas
        
        num_pixels=FILAS_DISCO*COLS_DISCO;
        num_elems=num_pixels*3; //cada pixel tiene 3 elementos (rojo, verde y azul)
        
        // se pasan los datos de la imagen a forma de vector
       // los datos se guardan en un vector, colocando consecutivos
       //los valores R, G y B de cada pixel 
       org_data = get_data_ptr ( org_img );
       
       //variables globales que usa el kmeans
       N=num_pixels; //numero total de pixeles que tiene la imagen
       R=3; //hay 3 componentes RGB en cada pixel
       
       d = RESERVAR_MAT_REAL(N, R); //se reserva memoria para guardar los datos de la imagen
       cuant = RESERVAR_MAT_REAL(N, R); //se reserva memoria para guardar los datos de la imagen cuantizada
       
       index=0;
       for(ik=0; ik<num_elems;ik+=3)
       {
       	  d[index][0]=org_data[ik];
       	  d[index][1]=org_data[ik+1];
       	  d[index][2]=org_data[ik+2];
       	  index++;
       }
  	
  	free_img ( org_img );
   
	

   // -------------------------------------------------------
   // -- APLICO K-MEDIAS --
   // -------------------------------------------------------	
	gettimeofday(&ti_a, NULL); 
	
	for(pruebas=0; pruebas <MAX_TEST; pruebas++)
        {
	   k_medias_aleatorio_X( N, n_clusters_quiero, R, 5, FILAS_DISCO, COLS_DISCO);
       }
	
	gettimeofday(&tf_a, NULL); 
	tiem_a = (tf_a.tv_sec - ti_a.tv_sec)*1000 +(tf_a.tv_usec - ti_a.tv_usec)/1000;
	printf("%s %s %llu\n", argv[2], argv[1], tiem_a);
	
   // -------------------------------------------------------
   // -- VOLCADO DE LA IMAGEN --
   // -------------------------------------------------------
   out_img = alloc_img ( FILAS_DISCO, COLS_DISCO, MAX_RGB ); //MAX_RGB: maximos colores que se pueden utilizar
   out_data = get_data_ptr ( out_img );
	 
   index = 0;//modificar una vez que se tenga todo lo anterior
   for ( ik = 0; ik < num_elems; ik += 3 )
   {

            //modificar una vez que se tenga todo lo anterior
            out_data[ik] = cuant[index][0];
            out_data[ik + 1] = cuant[index][1];
            out_data[ik + 2] = cuant[index][2];            
            index++;
   }

 
   sprintf(nombre_fich_salida, "KM_%s_%s.ppm", argv[2], nombre_fich_sin_path);
   write_img ( out_img, nombre_fich_salida );
 
   free_img(out_img);
   
   // -------------------------------------------------------
   // -- LIBERACIÓN DE LA MEMORIA DINAMICA UTILIZADA --
   // -------------------------------------------------------
   LIBERAR_MAT_REAL(&d, N);     

   //printf("\n\nValores para calculos: MAX. COLORES: %d epsilon=%.2lf\n", n_clusters_quiero, ERROR); 
   //printf("\nMAX_PASADAS: %d\n", MAX_PASADAS);
   
}
