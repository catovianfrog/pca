/******************************************************************************  
 *  pca.c  
 *  Performs Principal Component Analysis (PCA) on a set of data 
 */
    #define	VERSION "2.0"
/*  Author:	Bruno S. Charrière
 *  Date:	15.06.2014
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
#define LINE_LENGTH	100000		// max length of data file lines
#define NAME_LENGTH	255		// length of descriptive strings
#define MAX_ROWS	40000		// maximum number of data lines (arbitrary);
#define MAX_COLS	10000		// maximum number of data columns (parameters)

typedef char	    t_fname [NAME_LENGTH];
typedef char	    t_title [NAME_LENGTH];
typedef struct 	{   t_matrix	*data;		// numerical datasets
		    t_matrix	*data_cr;	// centered reduced dataset
		    t_title	title;		// project title
		    t_headers	*headers;	// parameters names
		    t_headers   *abbrevs;        // abbreviated parameters names
		    t_headers	*obs_id; 	// observations labels (lines)
		    int		n_obs;		// number of lines (observations)
		    int		ncols;	        // number of columns (parameters) 
		} t_dataset;
typedef struct 	{   int		order;		// order of the dataset (ncols)
		    t_matrix    *means;		// means of columns
		    t_matrix    *sdevs;		// standard deviations of columns
		    t_matrix    *covar;		// covariance matrix
		    t_matrix    *correl;		// correlation matrix
		} t_datastats;
typedef struct	{   /**** ACP data structure ****/
		    t_dataset   *data;		// pointer to dataset (dont forget to assign it;
		     				// we are not in C++ !) 
    		    t_datastats *stats;		// pointer to data base statistics
		    t_matrix    *e_values;	// eigenvalues
		    t_matrix    *inertia;	// inertia of eigenvalues
		    t_matrix    *e_vectors;	// eigenvectors
		    t_matrix    *princ_comp;	// principal components (coordinates
		    //of observations in eigenvectors rotated coordinate system)
		} t_acp_data;
const char	*g_version =VERSION;
const char	*g_progname=PROGNAME;
//---------------------------function prototypes--------------------------------
int	read_dataset(const t_fname fname, t_dataset *dataset);
t_headers *str2headers(char *linestr);
int	str2hvector(char *linestr, t_matrix *v);
void	dataset_free(t_dataset *d);
void	datastats_free(t_datastats *d);
void	print_dataset(FILE *fptr, t_dataset *dataset);
void	compute_datastats(t_dataset *D,t_datastats *S);
void	print_datastats(FILE *fptr, t_dataset *D, t_datastats *S);

/**********************************************************************
 *		    M A I N
 **********************************************************************/
int main( int argc, char *argv[]) {        // args not used so far

//    FILE	*outputf;
    char	in_fname[NAME_LENGTH] = DEFAULT_IN;
//    char	outfname[NAME_LENGTH] = DEFAULT_OUT;
    t_dataset	*dataset;	// structure containing input data
    t_datastats	*datastats;     // structure containing dataset base statistics
    t_acp_data	*pca;           // structure containing dataset PCA results
    int		iterations;	// number of iterations to computec eigenvalues
                 
    printf("%s, v%s\n",g_progname, g_version);
    // didn't put memory checks here. Shouldn't be needed.
    dataset   =	malloc(sizeof(*dataset));
    datastats =	malloc(sizeof(*datastats));
    pca	      =	malloc(sizeof(*pca));
    pca->data =dataset;
    pca->stats=datastats;
    read_dataset(in_fname,dataset);
    print_dataset(stdout,dataset);
    compute_datastats(dataset, datastats);
    print_datastats(stdout,dataset,datastats);
    pca->e_values=matrix_new(1,dataset->ncols);
    iterations=matrix_eigenvalues(datastats->correl, pca->e_values,PRECISION);
    printf("%d iterations.\n",iterations);
    if(iterations==0) goto error;


error:
    matrix_free(pca->e_values);
    free(pca);
    dataset_free(dataset);
    datastats_free(datastats);
    return 0;
}
/**********************************************************************
 * dataset_free: free dynamic memory of argument *dataset
 **********************************************************************/
 void	dataset_free(t_dataset *d) {
     free(d->obs_id->tags);
     free(d->obs_id);
     matrix_free_headers(d->headers);
     matrix_free(d->data);
     matrix_free(d->data_cr);
     free(d);
 }
/**********************************************************************
 * datastats_free: free dynamic memory of argument *datastats
 **********************************************************************/
 void	datastats_free(t_datastats *d) {
     matrix_free(d->means);
     matrix_free(d->sdevs);
     matrix_free(d->covar);
     matrix_free(d->correl);
     free(d);
 }
/**********************************************************************
* print_dataset(FILE *fptr, t_dataset *dataset)
**********************************************************************/
void print_dataset(FILE *fptr, t_dataset *dataset) {

    int		i,j;
    int		ncols, n_obs;
    t_tag	*title;
    t_title	*long_title;

    //not efficient but better legibility
    ncols=dataset->ncols;
    n_obs=dataset->n_obs;
    title=&dataset->headers->tags[0];
    long_title=&dataset->title;

    fprintf(fptr,"Data set: %s (%s)\n",*title, *long_title);
    fprintf(fptr,"Number of parameters:\t%d\n",ncols);
    fprintf(fptr,"Number of observations:\t%d\n",n_obs);
    fprintf(fptr,"\nVariables:  \t");
    for(i=1; i<=ncols; i++) {
	fprintf(fptr,"%s\t",dataset->headers->tags[i]);
    }
    fprintf(fptr,"\n");
    for(i=0;i<n_obs;i++) {
	fprintf(fptr,"#%-3d %-8s\t",i+1,dataset->obs_id->tags[i]);
	for(j=0;j<ncols;j++) {
	    fprintf(fptr,"%6.3f\t",dataset->data->data[i*ncols+j]);
	}
	fprintf(fptr,"\n");
    }
}                        

/**********************************************************************
* print_datastats(FILE *fptr, t_dataset *D, t_datastats *S)
**********************************************************************/
void print_datastats(FILE *fptr, t_dataset *D, t_datastats *S) {
    int		i,j;
    int		ncols;

    ncols=D->ncols;
    fprintf(fptr,"\n#-----------\t");
    for(i=1; i<=ncols; i++) {
	fprintf(fptr,"%s\t",D->headers->tags[i]);
    }
    fprintf(fptr,"\n");
    fprintf(fptr,"# Mean:     \t");
    for(j=0;j<ncols;j++) {
	fprintf(fptr,"%6.3f\t",S->means->data[j]);
    }
    fprintf(fptr,"\n");
    fprintf(fptr,"# Std Dev:  \t");
    for(j=0;j<ncols;j++) {
	fprintf(fptr,"%6.3f\t",S->sdevs->data[j]);
    }
    fprintf(fptr,"\n");

/**************************** OPTIONAL. Prints centered reduced dataset
    // print centered reduced data matrix
    fprintf(fptr,"\n# Centered reduced data");
    fprintf(fptr,"\n------------\t");
    for(i=1; i<=ncols; i++) {
	fprintf(fptr,"%s\t",D->headers->tags[i]);
    }
    fprintf(fptr,"\n");

    int	 n_obs;
    n_obs=D->n_obs;
    for(i=0;i<n_obs;i++) {
	fprintf(fptr,"#%-3d %s\t",i+1,D->obs_id->tags[i]);
	for(j=0;j<ncols;j++) {
	    fprintf(fptr,"% 4.2f\t",D->data_cr->data[i*ncols+j]);
	}
	fprintf(fptr,"\n");
    }
// ------------- end optional printing ---------------------*/
    // print correlations matrix
    fprintf(fptr,"\n#Correlations\t");
    for(i=1; i<=ncols; i++) {
	fprintf(fptr,"%s\t",D->headers->tags[i]);
    }
    fprintf(fptr,"\n");
    for(i=0;i<ncols;i++) {
	fprintf(fptr,"# %-10s\t",D->headers->tags[i+1]);
	for(j=0;j<ncols;j++) {
	    fprintf(fptr,"%6.3f\t",S->correl->data[i*ncols+j]);
	}
	fprintf(fptr,"\n");
    }
}
 
 /**********************************************************************
 * compute_datastats(t_dataset *D,t_datastats *S)
 **********************************************************************/
void compute_datastats(t_dataset *D, t_datastats *S) {


    t_matrix	*V;
    t_matrix	*T; // transposed intermediate
    int		i,j;
    double	norm;	//to divide matrix by number of data lines

    S->means  =matrix_new(1,D->ncols);
    S->sdevs  =matrix_new(1,D->ncols);
    D->data_cr=matrix_new(D->n_obs, D->ncols);
    norm=1/(double)D->n_obs;
    //compute means and standard deviation of each vector
    for(j=0;j<D->ncols;j++){
	V=matrix_get_vector(D->data,j);
	S->means->data[j]=mean(V);
	S->sdevs->data[j]=sdev(V);
	matrix_free(V);
	for(i=0;i<D->n_obs;i++) {
	    // put centered reduced data in D->data_cr
	    D->data_cr->data[i*D->ncols+j] = (D->data->data[i*D->ncols+j] - S->means->data[j])/S->sdevs->data[j];
	}
    }
    T=matrix_transpose(D->data_cr);
    S->covar=matrix_prod(T, D->data_cr);
    matrix_free(T);
    S->correl=matrix_scale(S->covar,norm);
}
/**********************************************************************
* read_dataset: read data file and return a t_dataset structure 
*   The first string on first line is the name of data set (no spaces)
*   the rest of the first line iare headings, separated by spaces or tabs
*   The number of headings defines the number of characters (columns)
*   Underscores in headings are turned into spaces
*   TODO if a line starts with ABBREV, the strings are abbreviated headings
*	 manage escape characters such as \' (apostrophe)
**********************************************************************/
int	read_dataset(const t_fname fname, t_dataset *dataset) {
    FILE	*fichier_data; 
    char	s[LINE_LENGTH];
    t_tag	tag;
    int		nrows;
    char	*ptr;
    t_matrix	*V;		// hVector

    fichier_data=fopen(fname, "r");
    if (fichier_data == NULL) {
	fprintf(stderr, "Erreur: impossible d'ouvrir le fichier %s.\n", fname);
	exit(8);
    }
    printf("Reading data from file: %s\n", fname);
    // read first line (headings, separated by spaces or tabs)
    while (fgets(s, sizeof(s), fichier_data) != NULL) {
	if(s[0]!='#') break;
    }
    if(s == NULL) {
	fprintf(stderr, "Error: no headers line in data file %s.\n",fname);
	exit(8);
    }                    
    s[strlen(s)-1]='\0';  // remove NL at the end
    // read first line = name + headers
    dataset->headers=str2headers(s);  // first colunm are labels
    V=matrix_new(1,dataset->headers->ntags-1);
    dataset->data=matrix_new(0,dataset->headers->ntags-1);
    dataset->ncols=dataset->headers->ntags-1;
    nrows=0;
    dataset->obs_id=calloc(1,sizeof(t_headers));
    dataset->obs_id->tags=NULL;
    while (fgets(s, sizeof(s), fichier_data) != NULL) {
	if(nrows>=MAX_ROWS) {
	    fprintf(stderr,"Error: too many lines in data file (max %d)\n",MAX_ROWS);
	    exit(8);
	}
	if(s[0]=='#') continue;
        if(sscanf(s,"%s",tag)!=1) {
	    fprintf(stderr, "Error: observation label n°%d cannot be read.\n",nrows+1);
	    exit(8);
	}
	dataset->obs_id->tags=realloc(dataset->obs_id->tags,sizeof(t_tag)*(nrows+1));
	strcpy(dataset->obs_id->tags[nrows],tag);
	// ptr shift: eliminates the first label from s before reading numerical values
	ptr=s+strlen(tag);
	s[strlen(s)-1]='\0';  // remove NL at the end
	if (str2hvector(ptr, V)==0) {
	    fprintf(stderr, "Error: reading line %d of input file %s.\n",nrows, fname);
		 exit(8);
	    }
	    matrix_add_rows(dataset->data,1);
	    matrix_set_block(V,dataset->data,nrows,0);
	nrows++;
    }
    // A longer title should be planned. 
    // For the moment long title=short title (Headers[0])
    strcpy(dataset->title,dataset->headers->tags[0]);
    dataset->n_obs=nrows;
    dataset->obs_id->ntags=nrows;
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
t_headers *str2headers(char *linestr) {
    t_headers	 *headers;
    t_tag	 s;     // string used to build the header substr
    int		 j;     // indice de position du caractère  du header
    int		 k;	// indice de position dans la chaine
    int		cols;	// number of headers/columns in data set
    
    headers=calloc(1,sizeof(t_headers));
    headers->tags=NULL;
    // error handling here
    k=0;                        
    cols=0;
    strcat(linestr," ");  // add a space at the end to terminate the headers counting
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
	cols++;
	headers->tags=realloc(headers->tags,sizeof(t_tag)*cols);
	strcpy(headers->tags[cols-1],s);
    }
    headers->ntags=cols-1;
    return headers;
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
    
    s[0]='\0'; 
    k=0;                        
    l=0;
    while(k<LINE_LENGTH && l<v->ncols){
        j=0;
        if (linestr[k]=='\0') break;
        while (linestr[k]==' ' || linestr[k]=='\t') k++; /* skip leading spaces */
        while (linestr[k]!=' ' && linestr[k]!='\t' && linestr[k]!='\0') /* read next float*/
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
	//printf("s: %s\tx:%4.2f\n",s,v->data[l]);       //for debug
        l++;
    };
if(l==v->ncols) return l;
return 0;  /* the line doesn't have enough float values */
}
         
