/**********************************************************/
/***                     "octree.c"                     ***/
/***          Octree Color Quantizer Program            ***/
/***             produces 256 color images              ***/
/***           based on code from the book,             ***/
/***           "Practical Ray Tracing in C"             ***/
/***          written in Borland/Turbo C 	            ***/
/***                         by                         ***/
/***     Wolfgang Stuerzlinger and Craig A. Lindley     ***/
/***		simpified by F.S.Hill,Jr. 4/30/94		    ***/
/***                                                    ***/
/**********************************************************/


/*
This code implements the algorithm described in the article
"A Simple Method for Color Quantization: Octree Quantization" by
Michael Gervautz, Werner Purgathofer from the book "Graphic Gems"
edited by Andrew Glassner. This code was originally written
and placed in the public domain by Wolfgang Stuerzlinger. The
code was subsequently modified for use in the PC environment by
the author.

Then Sandy clobbered it for ECE661 in April, 1994. It is now
reasonably system and ANSI compiler independent.

The main() program can be used to see how the routines fit into
your own program, and to try out the tools here.

Note: this is a fairly memory hungry utility: the octree uses
a lot of storage. On a PC type machine you may have to use the
compact, large, or huge memory model.

If HLdebug is defined, no image files are opened. Instead a 'fake'
image of 10 rows and 50 columns of pseudorandom pixel triples is 
generated and processed.

If HLverbose is defined, various internal values are printed
on the screen during the process, so you can see what is happening.

In normal use (HLdebug not defined), the user would give the filename
of the original 24 bit/pixel image to be processed, and the name of the
new file to be created.

It is assumed that from this file the program can capture:
1. The number of rows, ImageHeight, of the filed image.
2. The number of columns, ImageWidth, of the filed image.
3. The sequence of RGB bytes for each pixel.

After processing, the new output file should contain:

0). header information, e.g. a BMP header: up to you to handle this;
1). a 256 entry LUT of RGB BYTE triples;
2). a sequence of BYTE indices into the LUT, one for each pixel.

*/

#include <time.h>   
#include <sys/time.h> 
#include "../archivos_imagenes/header.h"  //Para imágenes ppm

/* comment out one or both of the following as you wish */
//#define HLdebug
#define HLverbose

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <malloc.h> 
#include <math.h>

#define TRUE 1
#define FALSE 0
#define MAX(a,b) (((a) > (b)) ? (a) : (b))
#define TESTBIT(a,i) (((a)&(1<<(i)))>>(i))
#define MAXDEPTH 7           /* max depth of octree - 1 */

#define verboso 0 //1=> MOSTRAR MENS; 0=>NO


int MAXCOLOR=256; //256 como valor inicial
typedef unsigned char BYTE;

typedef struct {  		/* LUT palette is an array of these */
	BYTE Red;
	BYTE Green;
	BYTE Blue;
} ColorRegister;

struct ColorSum {
  unsigned long R;
  unsigned long G;
  unsigned long B;
};

typedef struct Node * OCTREE;

struct Node {                /* an octree node data structure */
  unsigned Leaf;     		 /* 1 if this node is a leaf */
  unsigned Level;       	 /* the level this node is at */
  unsigned ColorIndex;  	 /* its LUT index -- assigned in pass 1 */
  unsigned Children;    	 /* number of children this node has */
  unsigned long ColorCount;  /* how many times this color has occurred */
  struct ColorSum RGBSum;    /* sum of color values accumulated */
  OCTREE NextReduceable;     /* link to next reduceable node at this level */
  OCTREE Next[8];       	 /* the 8 kin of this node */
};


/* Global variables used in the octree method */
static unsigned ImageWidth;
static unsigned ImageHeight;
static unsigned Size;
static unsigned ReduceLevel;
static unsigned LeafLevel;
static ColorRegister Palette[256];   /* the LUT palette - global */
static OCTREE Tree;                       /* the octree */
static OCTREE ReduceList[MAXDEPTH + 1];
static unsigned NumNodes = 0;             /* num.nodes in octree */
static unsigned MaxNodes = 0;             /* max. num. ever in octree */
static unsigned seed;


int FILAS_IMG, COLS_IMG; //para volcar el tamaño correcto 


/*<<<<<<<<<<<<<<<< prototypes >>>>>>>>>>>>>>>>>>>>>>>*/
unsigned Quantize(OCTREE Tree, ColorRegister TheColor);
void FillPalette(OCTREE Tree, unsigned *Index);
void NewandInit(OCTREE *Tree, unsigned Depth);
void GetReduceable(OCTREE *Node);
void MakeReduceable(unsigned Level, OCTREE Node);
void KillTree(OCTREE *Node);
void ReduceTree(void);
void InsertTree(OCTREE *Tree, ColorRegister RGB, unsigned Depth);
void GenOctree(char *FileName, OCTREE *Tree);
void ProcessPass1(char *FileName);
//void
int ProcessPass2(char *FileName, char *OutFileName);

/*<<<<<<<<<<<<<<<<<<<<<<< Quantize >>>>>>>>>>>>>>>>>>>>>>>>*/
/* This function returns the color index of the color that most closely
matches the color passed to it. This function can only be called after
the Octree is built.*/

unsigned Quantize(OCTREE Tree, ColorRegister TheColor)
{

  if (Tree->Leaf)
     return(Tree->ColorIndex);
  else
     return(Quantize(Tree->Next[
               TESTBIT(TheColor.Red,  MAXDEPTH - Tree->Level) * 4 +
               TESTBIT(TheColor.Green,MAXDEPTH - Tree->Level) * 2 +
               TESTBIT(TheColor.Blue, MAXDEPTH - Tree->Level)],TheColor));
}

/*<<<<<<<<<<<<<<<<<<<<<<< FillPalette >>>>>>>>>>>>>>>>>>>>>>*/
/* This function traverses the octree assigning a unique color index
value to each leaf node. It also calculate the RGB color components
to be used in the leaf nodes by averaging. Finally, it copies the
RGB color value for the leaf node into the Palette array.*/

void FillPalette(OCTREE Tree, unsigned *Index)
{

  unsigned Ind;

  if (Tree != NULL)  
  {
     if (Tree->Leaf || Tree->Level == LeafLevel)  
	 {
        Palette[*Index].Red   = (BYTE)(Tree->RGBSum.R / Tree->ColorCount);
        Palette[*Index].Green = (BYTE)(Tree->RGBSum.G / Tree->ColorCount);
        Palette[*Index].Blue  = (BYTE)(Tree->RGBSum.B / Tree->ColorCount);
#ifdef HLverbose
		 /* esto muestra la paleta,creo
		printf("\n Pal[%3u]=(%u,%u,%u)",*Index,Palette[*Index].Red,
		 Palette[*Index].Green, Palette[*Index].Blue);
		*/
#endif
        Tree->ColorIndex = *Index;
        Tree->Leaf = 1;
        *Index = *Index + 1;
     }
     else   
	 {
        for (Ind=0; Ind < 8; Ind++)
		   FillPalette(Tree->Next[Ind],Index);
     }
  }
}

/*<<<<<<<<<<<<<<<<<<<<<< NewandInit >>>>>>>>>>>>>>>>>>>>>>>*/
/* This function is called whenever a new node needs to be allocated. */

void NewandInit(OCTREE *Tree, unsigned Depth)
{

	*Tree = (struct Node *) calloc(1,sizeof(struct Node));
	if (*Tree == NULL)
	{
		 printf("Error: out of memory"); exit(-1);
	}


  NumNodes++;
  MaxNodes = MAX(MaxNodes,NumNodes);

  (*Tree)->Level = Depth;
  (*Tree)->Leaf = (Depth >= LeafLevel);
  if ((*Tree)->Leaf)
     Size++;
}

/* <<<<<<<<<<<<<<<<<<<<<< GetReduceable >>>>>>>>>>>>>>>>>>*/
/* Get a pointer to the next node to be reduced. */
void GetReduceable(OCTREE *Node)
{

  unsigned NewReduceLevel;

  NewReduceLevel = ReduceLevel;
  while (ReduceList[NewReduceLevel] == NULL) /* find a level with a node */
     NewReduceLevel--;                       /* to reduce */
  *Node = ReduceList[NewReduceLevel];        /* return ptr to that node */
  ReduceList[NewReduceLevel] =               /* pop from list of reducibles */
         ReduceList[NewReduceLevel]->NextReduceable;
}



/* <<<<<<<<<<<<<<<<<<<<<< KillTree >>>>>>>>>>>>>>>>>>>>>*/
/* This function recursively frees memory associated with the node
it is pointing at and all nodes under it (its siblings).*/

void KillTree(OCTREE *Tree)   {

  register unsigned Index;

  if (*Tree == NULL)
    return;
  for (Index=0; Index < 8; Index++)
    KillTree(&((*Tree)->Next[Index]));

  NumNodes--;
  free(*Tree);
  *Tree = NULL;
}

/*<<<<<<<<<<<<<<<<<<<< ReduceTree >>>>>>>>>>>>>>>>>>>>>>*/
/* This function performs the reduction of the octree. */

void ReduceTree(void)
{

  OCTREE Node;
  register unsigned Index, Depth;

  GetReduceable(&Node);                /* get ptr to deepest node needing reduction */
  Node->Leaf = 1;                      /* change it from intermediate node to leaf */
  Size = Size - Node->Children + 1;    /* reduce number of leaves by the number of children. */
  Depth = Node->Level;

  for (Index=0; Index < 8; Index++)    /* free memory occupied by children */
    KillTree(&(Node->Next[Index]));

  if (Depth < ReduceLevel)  {          /* if depth of reduced node is less */
	 ReduceLevel = Depth;			  /* than current reduce level, reduce */
#ifdef HLverbose
	  if(verboso ==1)
	    printf("\n ReduceLevel -> %u",ReduceLevel);
	// getchar(); //getch();
#endif
	 LeafLevel = ReduceLevel + 1;      /* reduce level and leaf level to */
  }                                    /* deepest intermediate (non leaf) node. */
}

/* <<<<<<<<<<<<<<<<<< MakeReduceable >>>>>>>>>>>>>>>>>>>>>*/
/* This function links a node that has been determined to be
reducible into a linked list of reducible nodes at the same
depth level in the octree. */

void MakeReduceable(unsigned Level, OCTREE Node)
{

  Node->NextReduceable = ReduceList[Level];
  ReduceList[Level] = Node;
}

/* <<<<<<<<<<<<<<<<<<< InsertTree >>>>>>>>>>>>>>>>>>>>>>*/
/* This function builds the octree from the pixel color values passed
to it. */

void InsertTree(OCTREE *Tree, ColorRegister RGB, unsigned Depth)
{

  unsigned Branch;

  if (*Tree == NULL)
     NewandInit(Tree,Depth);
  (*Tree)->ColorCount++;
  (*Tree)->RGBSum.R += (unsigned long) RGB.Red;
  (*Tree)->RGBSum.G += (unsigned long) RGB.Green;
  (*Tree)->RGBSum.B += (unsigned long) RGB.Blue;
  if (((*Tree)->Leaf == FALSE) && (Depth < LeafLevel))  {
     Branch = TESTBIT(RGB.Red,  MAXDEPTH - Depth) * 4 +
              TESTBIT(RGB.Green,MAXDEPTH - Depth) * 2 +
              TESTBIT(RGB.Blue, MAXDEPTH - Depth);
     if ((*Tree)->Next[Branch] == NULL)  {
        (*Tree)->Children++;
        if ((*Tree)->Children == 2)    /* any node with 2 or more children */
           MakeReduceable(Depth,*Tree);/* is candidate for reduction */
     }
     InsertTree(&((*Tree)->Next[Branch]),RGB,Depth + 1);
  }
}




/* <<<<<<<<<<<<<<<< GenOctree >>>>>>>>>>>>>>>>>>>>>>*/
/* This function makes a pass over the image data and calls
InsertTree which builds the octree from each pixel. */
void GenOctree(char *FileName, OCTREE *Tree)
{
   FILE *ImageFile;              
   register unsigned Row, Col;
   unsigned RowNum;
   int ReturnCode;
   ColorRegister RGB;


int ro, ve, az; //VALORES RGB
int n_f, n_c; //NUMERO DE FILAS Y COLUMNAS

// >>>>>>>>>>>>>>>>
Image *org_img;  // PARA LEER IMAGENES PPM. imagen PPM de entrada
uchar *org_data; //  PARA LEER IMAGENES PPM.valores de los pixels que se extraen de org_min
int num_elems,   //  PARA LEER IMAGENES PPM.número de elementos de org_data (3*nº de pixels de la imagen)
  index,         //  PARA LEER IMAGENES PPM.contador de pixels de la imagen
  ik =0;         //  PARA LEER IMAGENES PPM.contador de elementos de org_data

// <<<<<<<<<<<<<<<<


   Size = 0;                           /* initialize variables */
   ReduceLevel = MAXDEPTH;
   LeafLevel = ReduceLevel + 1;

#ifndef HLdebug
   /* Attempt to open the file */
  /*ImageFile = fopen(FileName,"rb");
   if(!ImageFile)
   {
      printf("\n can't open %s", FileName); 
      exit(0);
   }
   */


   //  PARA LEER IMAGENES PPM.
   // se lee la imagen original (un fichero PPM)
   org_img = read_img (FileName);   
   
   // se pasan los datos de la imagen a forma de vector
   /* los datos se guardan en un vector, colocando consecutivos
   los valores R G y B de cada pixel */
   org_data = get_data_ptr ( org_img );    
   
     /* aviso: las siguientes instrucciones son válidas
          cuando se lee una imagen de 512 x 512 pixels */
	ImageWidth = 512;   
	ImageHeight = 512; //para imagenes Kodak:  768; 

	
	
#else
   puts("\n debug mode: generate image randomly..");
   ImageWidth = 50;
   ImageHeight = 10;
   puts("\n indica una semilla : ");
   scanf("%d",&seed);  	/* set global seed for random pixels. */
   srand(seed);
#endif

   /*	Now scan complete image building the octree.  */
   /* leo la primera fila del fichero de imagen, que contiene el numero de pixel y
   un 3, que indica el número de colores que hay por pixel */
   //fscanf(ImageFile,"%d %d", &n_f, &n_c);


   //  PARA LEER IMAGENES PPM.
   // copio en las variables globales las dimensiones de la imagen    
   FILAS_IMG = org_img->num_rows;  // - altura de la imagen           
   COLS_IMG = org_img->num_cols;   // - anchura de la imagen
   

   n_f = FILAS_IMG ;
   n_c = COLS_IMG ;
   
	ImageHeight = n_f;
	ImageWidth = n_c;   

   
   	
	
   for (Row=0; Row < ImageHeight; Row++)
	for (Col=0; Col < ImageWidth; Col++)
   	{

#ifndef HLdebug
            // lo que he leido como enteros, lo copio como caracter sin signo
	    /*RGB.Blue =  (unsigned char)az; 
	    RGB.Green = (unsigned char)ve; 
            RGB.Red = (unsigned char) ro;  */

            RGB.Red = (unsigned char)org_data[ik];
	    ik++; 
	    RGB.Green = (unsigned char)org_data[ik];
	    ik++;	    
    	    RGB.Blue = (unsigned char)org_data[ik];   //PARA LEER IMAGENES PPM.  
	    ik++;
       
#else
	   RGB.Red = 	(BYTE)(rand() % 256); /* pseudorandom pixels */
	   RGB.Green =  (BYTE)(rand() % 256);
	   RGB.Blue = 	(BYTE)(rand() % 256);
#endif
#ifdef HLverbose
/*	   printf("\n %u,%u,%u,%u -> (%u,%u,%u)", Row, Col,Size,
				NumNodes, RGB.Red, RGB.Green, RGB.Blue);
*/
#endif

	   InsertTree(Tree,RGB,0);       /* insert color into octree */
       if (Size > MAXCOLOR - 1)      /* max number of colors */
         ReduceTree();               /* reduce tree if more than 256 colors */
   }

	if(verboso==1)
	{
       printf("\nPass 1 complete.\n");
	   printf("\nMaxNodes = %d\n\n",MaxNodes);
	}

#ifndef HLdebug
   //fclose(ImageFile);
   
      // se libera la memoria asociada a la imagen PPM
   free_img ( org_img );    //  PARA LEER IMAGENES PPM.
#endif
}





/*<<<<<<<<<<<<<<<<<< ProcessPass1 >>>>>>>>>>>>>>>>>>>*/
/* Pass one: build the octree and load the color palette.*/
void ProcessPass1(char *FileName)
{/* just pass this file name on to GenOctree() */
unsigned Index;
int      ReturnCode;

	if(verboso == 1)
	{
      printf("\nPass1\n");
      printf("\nfichero que contiene la imagen: |%s|", FileName);
	}

   /* Initialize all of the Palette memory to zeros. */
   memset(Palette,'\0',sizeof(Palette));

   Tree = NULL;
   GenOctree(FileName,&Tree);     /* read through the file */

   Index = 1;                     /* entry 0 is left black */
   /* scan whole octree, and for each leaf compute average color,
	  assign index to leaf, and load avg intopoalette */

#ifdef HLdebug
	puts("\n ready to load LUT, <CR> to go on");
#endif
   FillPalette(Tree,&Index);
   /* on return Index is # of colors */

}


/*<<<<<<<<<<<<<<<<<<< ProcessPass2 >>>>>>>>>>>>>>>>>>>>>>>*/
/* This function maps the original image pixel data from the
image file into the new color map. */

/* OJO: AQUI EL VOLCADO A PPM LO VOY A INTEGRAR EN ESTA FUNCIÓN
   volcar_imagen_ppm(OutFileName, centros_R, centros_G, centros_B,
        cluster_cada_item);

*/

//void 
int ProcessPass2(char *FileName, char *OutFileName)
{
   FILE *ImageFile, *OutFile;
   register unsigned Row, Col, PixVal;
   unsigned RowNum;
   ColorRegister RGB;

int ro, ve, az,           // valores RGB de un pixels
    ro_sa, ve_sa, az_sa;  // valores RGB de otro pixel
int n_f, n_c;   // número de filas y columnas de una imagen


// --- para calcular tiempo en milisegundos
struct timeval ti_a, tf_a;   // para la aplicación del algoritmo 

unsigned long long tiem_a, tiem_m;  // para almacenar la diferencia entre los dos tiempos previos
int T_entero =0; 


// >>>>>>>>>>>>>>>>  PARA VOLCAR IMAGEN PPM
Image *out_img;   // PARA VOLCAR IMAGEN PPM
uchar *out_data;  // PARA VOLCAR IMAGEN PPM
int ik =0;           // PARA VOLCAR IMAGEN PPM
int num_elems = FILAS_IMG * COLS_IMG * 3; // PARA VOLCAR IMAGEN PPM	
int aux;   // PARA VOLCAR IMAGEN PPM
// <<<<<<<<<<<<<<<<<

// >>>>>>>>>>>>>>>>
Image *org_img;  // PARA LEER IMAGENES PPM. imagen PPM de entrada
uchar *org_data; //  PARA LEER IMAGENES PPM.valores de los pixels que se extraen de org_min
int        //  PARA LEER IMAGENES PPM.contador de pixels de la imagen
  ik_in=0;            //  PARA LEER IMAGENES PPM.contador de elementos de org_data

// <<<<<<<<<<<<<<<<


    // PARA VOLCAR IMAGEN PPM
   out_img = alloc_img ( FILAS_IMG, COLS_IMG, MAX_RGB );
   out_data = get_data_ptr ( out_img );


	/* Rescan the raw image file so that it can now be displayed.*/

	if(verboso ==1)   printf("\nPass 2\n");

#ifndef HLdebug

   
   //PARA LEER IMAGENES PPM
   org_img = read_img (FileName);
      
   // se pasan los datos de la imagen a forma de vector
   /* los datos se guardan en un vector, colocando consecutivos
   los valores R G y B de cada pixel */
   org_data = get_data_ptr ( org_img );        

	/* now write header to output file, and Palette[] contents */


#else
	puts("\n debug still: generate same image randomly ");
	srand(seed);  /* same seed as before */
#endif


   for (Row=0; Row < ImageHeight; Row++)
      for (Col=0; Col < ImageWidth; Col++)
      {
         /* Read the raw image file (or fake pixel values, into three
		variables, Red, Green & Blue a pixel at a time. */
         #ifndef HLdebug
	    
            ro = RGB.Red = (unsigned char)org_data[ik_in];
	    ik_in++;            
  	    ve = RGB.Green = (unsigned char)org_data[ik_in];
	    ik_in++;          
	    az =RGB.Blue = (unsigned char)org_data[ik_in];   //PARA LEER IMAGENES PPM.  
	    ik_in++;
	    
         #else
	  RGB.Red = (BYTE)(rand()%256);
	  RGB.Green = (BYTE)(rand()%256);
	  RGB.Blue = (BYTE)(rand()%256);
         #endif
        
         PixVal = Quantize(Tree,RGB);/* select best color index in LUT */

         /* write PixVal to outfile */

         #ifndef HLdebug

        // PARA VOLCAR IMAGEN PPM
        out_data[ik]     =  (unsigned char)(Palette[PixVal].Red);
        ik++;
        out_data[ik] =  (unsigned char)(Palette[PixVal].Green);
         ik++;
        out_data[ik] =  (unsigned char) (Palette[PixVal].Blue);
         ik++;
  
	    
         #endif
         
         #ifdef HLverbose
	   //	 printf("\n %u,%u -> index = %u",Row,Col,PixVal);
         #endif

         // se acumula un término para calcular los errores MSE y MAE
         ro_sa = (int) Palette[PixVal].Red;
         ve_sa = (int) Palette[PixVal].Green;
         az_sa = (int) Palette[PixVal].Blue;
      }












#ifndef HLdebug
   free_img ( org_img );    //  PARA IMAGENES PPM.
   
   // PARA VOLCAR IMAGEN PPM
   write_img ( out_img, OutFileName);
   free_img(out_img);   
#endif



	return T_entero;
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



  
       
        
        
/* <<<<<<<<< Main Octree Quantizer Program >>>>>>>>>>>>> */
void main(int argc, char *argv[])
{
char    ImageFileName[80]; //Amplio, porque al incluir un path me corta la cadena
char    OutputFileName[80];


// para calcular el tiempo empleado en los cálculos
// --- para calcular tiempo en segundos
time_t T_ini_a, T_fin_a; // instantes inicial y final de la aplicación del algoritmo 
// --- para calcular tiempo en milisegundos
struct timeval ti_a, tf_a;   // para la aplicación del algoritmo 

unsigned long long tiem_a, tiem_m;  // para almacenar la diferencia entre los dos tiempos previos
double tiem_double;

int T_proceso2;
int K; //numero de colores de la paleta cuantizada



char nombre_fs[100], // nombre del fichero de salida
     nombre_fich_sin_path[100]; // nombre del fich de entrada sin path ni parte final separada por .
char nombre_fich_salida_ppm[100];

	
      (void) time(&T_ini_a); 
  
        
        /* al ejecutar desde el terminal se debe indicar el nombre del fichero de entrada.
	Si no se proporciona, el programa concluye */
	if(argc >=3) 
		strcpy(ImageFileName, argv[1]);
	else
	{
		printf("\nDebe indicar el nombre del fichero que contiene la imagen original y el numero de colores de la nueva paleta");
		exit(1);
	}

	 
	/* tomo la parte del fichero de entrada que excluye el path y los trozos
	finales separados por puntos */
	extraer_nombre_fich(ImageFileName, nombre_fich_sin_path);
	
	//se guarda el numero de colores de la paleta cuantizada
	MAXCOLOR= atoi(argv[2]);

	//defino el nombre del fichero de salida (será una imagen PPM)
	sprintf(OutputFileName, "OC_%d_%s.ppm", MAXCOLOR, nombre_fich_sin_path);
	

   
   if(verboso==1)
	puts("\n******* Begin octree quantization program ********");


        gettimeofday(&ti_a, NULL);   // para calcular tiempo en usegundos

	
	/* Pass one preprocesses the image data to build the Octree
	for quantization. */
	ProcessPass1(ImageFileName);

	
	/* Pass two scans the octree and assigns indices to each octree
	   color. It then scans the image again and maps the pixel
	   codes from the image to the closest color found in the Octree. */
	T_proceso2= ProcessPass2(ImageFileName,OutputFileName);

        (void) time(&T_fin_a);    
        gettimeofday(&tf_a, NULL);   // para calcularlo en usegundos
	
        tiem_a = (tf_a.tv_sec - ti_a.tv_sec)*1000 +(tf_a.tv_usec - ti_a.tv_usec)/1000;
   
	printf("%s %s %llu\n", argv[2], argv[1], tiem_a);

	

	//printf("\n OJO: definir nº colores (cte. MAXCOLOR) y compilar\n");


#ifdef HLverbose
/*	  printf("\n final statistics:");
	  printf("\ntotal pixels = %d (%dx%d)", ImageWidth *ImageHeight, ImageWidth, ImageHeight);
	  printf(" max. nodes in octree at any time = %u", MaxNodes);
	  printf("  num. nodes in octree now = %u", NumNodes);
	  printf("  ReduceLevel = %u", ReduceLevel);
	  printf("  LeafLevel = %u", LeafLevel);
	  */


#endif

	/* Free the memory occupied by the Octree */
	KillTree(&Tree);
	//puts("\n done! <CR> to exit ");
}
