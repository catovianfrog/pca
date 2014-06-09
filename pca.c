/******************************************************************************  
 *  pca.c  
 *  Performs Principal Component Analysis (PCA) on a set of data 
 */
    #define	VERSION "2.0"
/*  Author:	Bruno S. Charrière
 *  Date:	9.06.2014
 *****************************************************************************/
#include <stdio.h>
#define  _GNU_SOURCE
#include <stdlib.h>
#include <math.h>
#include "dbg.h"	// error handling functions
#include "libal.h"	// linear algebra & statistics library

#define	PROGNAME	"pca - Principal component analysis"
#define DEFAULT_OUT	"pca_out.dat"	// default name of output data file.
#define DEFAULT_IN	"pca_input.dat"	// default name of input  data file.
#define LINE_LENGTH	10000		// max length of data file lines
#define NAME_LENGTH	255		// length of strings & labels
#define MAX_ROWS	40000		// maximum number of data lines (arbitrary);
#define MAX_COLS	10000		// maximum number of data columns (parameters)

typedef char	    t_fname  [NAME_LENGTH];
typedef char	    t_header [NAME_LENGTH];
typedef	char	    *t_tag;                     // label type, variable length str
typedef t_tag	    *t_labels;			// dynamic array of labels
typedef struct 	{   t_matrix	data;		// numerical datasets
		    t_matrix	data_cr;	// centered reduced dataset
		    t_header	title;		// project title
		    t_labels	headers;	// parameters names
		    t_labels	abbrevs;	// abbreviated parameters names
		    t_labels	obs_id; 	// observations labels (lines)
		    int		n_obs;		// number of lines (observations)
		    int		ncols;		// number of columns (parameters) 
		} t_dataset;
typedef struct 	{   int		order;		// order of the dataset (ncols)
		    t_matrix    means;		// means of columns
		    t_matrix    sdevs;		// standard deviations of columns
		    t_matrix    covar;		// covariance matrix
		    t_matrix    correl;		// correlation matrix
		} t_datastats;
typedef struct	{   /**** ACP data structure ****/
		    t_dataset   *data;		// pointer to dataset (dont forget to assign it;
		     				// we are not in C++ !) 
    		    t_datastats *stats;		// pointer to data base statistics
		    t_matrix    e_values;	// eigenvalues
		    t_matrix    inertia;	// inertia of eigenvalues
		    t_matrix    e_vectors;	// eigenvectors
		    t_matrix    princ_comp;	// principal components (coordinates
		    //of observations in eigenvectors rotated coordinate system)
		} t_acp_data;
const char	*g_version =VERSION;
const char	*g_progname=PROGNAME;
//---------------------------function prototypes--------------------------------
int	read_dataset(const t_fname fname, t_dataset *dataset);
int	str2headers(char *linestr, t_labels headers);
int	str2hvector(char *linestr, t_matrix *v);

/**********************************************************************
 *		    M A I N
 **********************************************************************/
int main( int argc, char *argv[]) {        // args not used so far

    FILE	*outputf;
    char	in_fname[NAME_LENGTH] = DEFAULT_IN;
    char	outfname[NAME_LENGTH] = DEFAULT_OUT;
    t_dataset	*dataset;	// structure containing input data
    t_datastats	*datastats;     // structure containing dataset base statistics
    t_acp_data	*pca;           // structure containing dataset PCA results
                 
    printf("%s, v%s\n",g_progname, g_version);
    dataset   =	malloc(sizeof(*dataset));
    datastats =	malloc(sizeof(*datastats));
    pca	      =	malloc(sizeof(*pca));
    // didn't put memory checks here. Shouldn't be needed.
    pca->data=dataset;
    pca->stats=datastats;

    read_dataset(in_fname,dataset);

    
    free(pca);
    free(datastats);
    free(dataset);
    return 0;
}
/**********************************************************************
* read_dataset: read data file and return a t_dataset structure
**********************************************************************/
int	read_dataset(const t_fname fname, t_dataset *dataset) {

    FILE	*fichier_data; 
    char	s[LINE_LENGTH];
    int		nrows, ncols;
    char	*ptr;
    t_matrix	*V;		// hVector

    fichier_data=fopen(fname, "r");
    if (fichier_data == NULL) {
	fprintf(stderr, "Erreur: impossible d'ouvrir le fichier %s.\n", fname);
	exit(8);
    }
    printf("Reading data from file: %s\n", fname);
    // read first line (headings, separated by spaces or tabs)
    // first string is the name of data set (no spaces)
    // following strings are headings (names of parameters)
    // number of headings defines number of characters (columns)
    // Underscores in headings are turned into spaces
    // TODO if a line starts with ABBREV, the strings are abbreviated headings

    ncols=0;
    while (fgets(s, sizeof(s), fichier_data) != NULL) {
	if(s[0]!='#') break;
    }
    if(s == NULL) {
	fprintf(stderr, "Error: no headers line in data file %s.\n",fname);
	exit(8);
    }                    
    s[strlen(s)-1]='\0';  // remove NL at the end
    // read first line = name + headers
    ncols=str2headers(s,dataset->headers)-1;  // first colunm are labels
    
    V=matrix_new(1,ncols);
    nrows=0;
    while (fgets(s, sizeof(s), fichier_data) != NULL) {
	if(nrows>=MAX_ROWS) {
	    fprintf(stderr,"Error: too many lines in data file (max %d)\n",MAX_ROWS);
	    exit(8);
	}
	if(s[0]=='#') continue;
        if(sscanf(s,"%s",dataset->obs_id[nrows])!=1) {
	    fprintf(stderr, "Error: observation label n°%d cannot be read.\n",nrows+1);
	    exit(8);
	}

    // ptr shift: eliminates the first label from s
	ptr=s+strlen(dataset->obs_id[nrows]);  

	if (str2hvector(ptr, V)==0) {
	    fprintf(stderr, "Error: reading line %d of input file %s.\n",nrows, fname);
		 exit(8);
	    }
	    matrix_set_block(V,&dataset->data,nrows,0);
	nrows++;
    }
    // A longer title should be planned. For the moment long title= short title (Headers[0])
    strcpy(dataset->title,dataset->headers[0]);
    dataset->n_obs=nrows;
    dataset->ncols=ncols;   
    matrix_free(V);
    fclose(fichier_data);
    return 0;


 }
/**********************************************************************
 * str2headers: converts str in a vector of headers 
 *	    headers are separated by spaces or tabs
 *	    underscores are turned into spaces in headers
 *	    returns the number of headers read
 * TODO	    manage quotes and escape sequences in text data
 **********************************************************************/
int str2headers(char *linestr, t_labels headers) {
    t_header	 s;     // string used to build the header substr
    int		 j;     // indice de position du caractère  du header
    int		 k;	// indice de position dans la chaine
    int		cols;	// number of headers/columns in data set
    
    k=0;                        
    cols=0;
    strcat(linestr," ");  // add a space at the end to terminate the headers counting
    printf("=%s=\n",linestr);
    while(k<strlen(linestr)){
        if(cols>=MAX_COLS) {
	    fprintf(stderr,"Error: to many parameters in data file (max comlumns: %d).\n",MAX_COLS);
	    exit(8);
	}
        j=0;
	s[0]='\0';
        while (linestr[k]==' ' || linestr[k]=='\t') k++; /* skip leading spaces */
        while (linestr[k]!=' ' && linestr[k]!='\t' && linestr[k]!='\0') /* read next header */
        {
    		 s[j]=linestr[k]; 
		 if(s[j]=='_') {s[j]=' ';}; //transform underscores into spaces
    		 k++;
    		 j++;
        }
        s[j]='\0';
	strcpy(headers[cols],s);
	cols++;
    }
    return cols-1;
}                                                    
/**********************************************************************
 * str2hvector: converts the string into a vector of N=v->nrows reals
 *		argument v already exists, and is just filled in
 **********************************************************************/
int str2hvector(char *linestr, t_matrix *v) {
    char s[LINE_LENGTH];// string used to build the number
    int		 j;     // indice de position du caractère numérique dans le nombre      
    int		 k;	// indice de position dans la chaine
    int		 l;     // indice de position dans le vecteur horizontal
    
//printf("___________________\n%s\n",linestr);

    s[0]='\0'; 
    k=0;                        
    l=0;
    while(k<LINE_LENGTH && l<v->nrows){
        j=0;
        if (linestr[k]=='\0') break;
        while (linestr[k]==' ' || linestr[k]=='\t') k++; /* skip leading spaces */
        while (linestr[k]!=' ' && linestr[k]!='\t' && linestr[k]!='\0') /* read next float value */
        {
    		 s[j]=linestr[k];
    		 k++;
    		 j++;
        }
        s[j]='\0';
        if (sscanf(s, "%lf", &v->data[l]) != 1) {
// ATTENTION, if non digit characters follow digits, their are ignored.
		     printf("Error: %s is not a floating point number.\n", s);
    		 return(1);
        }
//printf("s: %s\tx:%4.2f\n",s,v[l]);
        l++;
    };
if(l==v->nrows) return l;
return 0;  /* the line doesn't have enough float values */
}
         
