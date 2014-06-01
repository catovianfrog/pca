/******************************************************************************  
 *  pca.c  
 *  Performs Principal Component Analysis (PCA) on a set of data */
 #define	VERSION "1.2" 
/*  Author:	Bruno S. Charrière
 *  Date:	10.05.2014
 ******************************************************************************
 * VERSION LOG
 * v0.9  start from initial linear algebra stub (qr.c)
 * v0.91 not there yet. Trying to read data file!
 * v0.92 read and store data set. Still debugging
 * v0.93 replaced MAX_COLS by MAX_LINES and MAX_COLS for matrices
 * v0.94 added size controls for dataset
 * v0.95 all basic stats calculated, including correlation matrix
 * v0.96 eigenvalues and eigenvectors calculated and displayed
 * v0.97 principal components (rotated axes) calculated and displayed
 * v0.98 outputs coordinates of observations in principal plane and correlation circle
 * v1.0	 working version. Limited to 150 lines of data (static data, in the stack)
 * v1.1	 all dynamic memory allocations in Main. Works up to 245 lines      
 * v1.2  changed from syslinQR to syslinGauss (muchg less memory usage)
 *	 changed from stack storage to dynamic allocation for matrices in functions
 *	 works up to 350 lines.
 *	 Corrected a bug in str2headers (added a space at the end of linestr)
 *
 *****************************************************************************/
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#define DEFAULT_NAME   "matrice.txt"// default name of input data file.
#define MAX_ITERATION  1e5	    // max number of iterations for eigenvalues
#define PRECISION      1e-6	    // stop iterations when residue < PRECISION
#define LINE_LENGTH    1000	    // max length of data file lines
#define NAME_LENGTH    255	    // length of strings & labels
#define MAX_LINES      350	    // max number of observations (individus)
#define MAX_COLS       MAX_LINES    // max number of parameters
/****************************************************************************
 *			W A R N I N G
 *                      ¨¨¨¨¨¨¨¨¨¨¨¨¨
 *  MAX_COLS must be the same as MAX_LINES, because the transposed matrix
 *  used in calculating data covariance have as many columns as there are
 *  lines in data set.
 *  Ideally, matrices should be allocated dynamically based on their
 *  dimensions, but the programme has to be completely rewritten.
 *  TODO: An alternative is to declare "double t_hmatrix[MAC_COLS][MAX_LINES]"
 *  for containers that receive transposed data matrices (not many). But then
 *  specific matrix functions (product, etc) will be required with the reversed
 *  two-dimensional array type.
 *
 *  v1.2 with static allocations for vectors runs with 350 rowsx350 cols maximum.
 *  Above this, valgrind detects invalid memory reads and writes.
 *
 ***************************************************************************/

typedef double	    t_matrix [MAX_LINES][MAX_COLS];
typedef double	    t_hmatrix[MAX_COLS][MAX_LINES]; // unused: storage of transposed t_matrix
typedef char	    t_fname  [NAME_LENGTH];
typedef char	    t_header [NAME_LENGTH];
typedef double	    t_vector [MAX_LINES];
typedef t_header    t_labels [MAX_LINES];
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
		    t_vector    means;		// means of columns
		    t_vector    sdevs;		// standard deviations of columns
		    t_matrix    covar;		// covariance matrix
		    t_matrix    correl;		// correlation matrix
		} t_datastats;
typedef struct	{   /**** ACP data structure ****/
		    t_dataset   *data;		// pointer to dataset (dont forget to assign it;
		     				// we are not in C++ !) 
    		    t_datastats *stats;		// pointer to data base statistics
		    t_vector    e_values;	// eigenvalues
		    t_vector    inertia;	// inertia of eigenvalues
		    t_matrix    e_vectors;	// eigenvectors
		    t_matrix    princ_comp;	// principal components (coordinates
		    //of observations in eigenvectors rotated coordinate system)
		} t_acp_data;

char	*version=VERSION;
/*------------------------ functions ------------------------------------------*/
int	sign(double x); 
int	str2hvector(char *s, t_vector v, int order);
int	str2headers(char *linestr, t_labels headers);
int	normalize(t_vector V, int n);
double	tracemat(t_matrix M, int n);                                
double	precision=PRECISION;
double	norm(t_vector V, int n);
int	qr_decomp(t_matrix M, t_matrix Q, t_matrix R, int n);
void	Id(t_matrix I, int n);
void	transpose(t_matrix M, t_matrix T,int rows, int cols);
void	copymat(t_matrix M, t_matrix C, int rows, int cols);
void	copyvect(t_vector U, t_vector V, int n);
void	matrix_line_from_vector(t_matrix M, t_vector V, int lineNr, int order);
void	vector_from_matrix(t_matrix M,t_vector V,int nrows,int coln);
void	print_matrix(t_matrix M, int order, char *msg);
void	prodmat(t_matrix A, t_matrix B, t_matrix P, int order);
void	print_matrixMN(t_matrix M, int n, int p, char *msg);
void	prodmat_scalar(t_matrix M, double a, t_matrix aM, int n, int p); 
void	prodVtV(t_vector V, t_matrix P, int order);
void	prodmatMN(t_matrix A, t_matrix B, t_matrix P, int m, int n, int p);
int	syslinQR(t_matrix A, t_vector x, t_vector b, int order); // solve Ax=b using QR decomposition
int	syslinGauss(t_matrix A, t_vector x, t_vector b, int order); // same with gauss pivot method
int	eigenvalues(t_matrix A, t_vector EV, int order, double eps);
void	eigenvector(t_matrix A, t_vector EV, double r, int order);
void	eigenvectors(t_matrix A, t_matrix EV, t_vector e_values, int n);
void    ev_inertia(t_vector e_values,t_vector inertia, int order);
double	mean(t_vector V, int n);
double	sumsquares(t_vector V, int n);
double	variance(t_vector V, int n);	    //variance of V as an entire population (/n)
double	variance_sample(t_vector V, int n); //variance of V, as a sample (/n-1)
double  sdev(t_vector V, int n);	    //std deviation of V, entire population(/n)
double  sdev_sample(t_vector V, int n);	    //std deviation of V as a sample (/n-1)
int	read_dataset(t_fname fname, t_dataset *dataset);
void	compute_datastats(t_dataset *D,t_datastats *S);
void	print_dataset(FILE *fptr, t_dataset *dataset);
void	print_datastats(FILE *fptr, t_dataset *D, t_datastats *S);
void	print_princ_comp(FILE *fptr, t_acp_data acp);
void	print_cercle_corr(FILE *fptr, t_acp_data acp);
void	print_eigenvalues(FILE *fptr, t_acp_data acp);
void    print_eigenvectors(FILE *fptr, t_acp_data acp);
void	print_acp_results(FILE* fptr, t_acp_data acp);

//------------------------------------------------------------------------------
//					M A I N 
//  TODO:
//	    - static data on stack in subroutines. Needs dynamic allocation
//	    - arguments reading routine
//	    - ask for data file name if not on command line
//	    - set some options
//	    - help routine
//  NOTES:
//	    - Maximum size of objects are hardcoded. Takes up a lot of memory
//	    - Maximum size of objects should be set as compile-time parameters
//------------------------------------------------------------------------------
int main( int argc, char *argv[]) {        // args not used so far

    FILE	*outputf;
    char	 data_fname[NAME_LENGTH] = DEFAULT_NAME;
    char	*outfname="dataset.txt";
    t_dataset	*dataset;
    t_datastats	*datastats;
    t_acp_data	*acp;

    // Initialization of memory is not necessary at this stage (checked)
    //dataset	=calloc(1,sizeof(*dataset));
    //datastats	=calloc(1,sizeof(*datastats));
    //acp	=calloc(1,sizeof(*acp));
    dataset	=malloc(sizeof(*dataset));
    datastats	=malloc(sizeof(*datastats));
    acp		=malloc(sizeof(*acp));

    acp->data=dataset;			// CRITICAL or core-dump!
    acp->stats=datastats;		// need to find a smarter way !!
    printf("%s, v%s\n", argv[0],version);   // should be on request only (arg -V)

    read_dataset(data_fname, dataset);
    print_dataset(stdout,dataset);
    compute_datastats(dataset, datastats);
    print_datastats(stdout,dataset,datastats);
    eigenvalues(datastats->correl, acp->e_values,dataset->ncols,precision);
    ev_inertia(acp->e_values, acp->inertia, dataset->ncols);
    eigenvectors(datastats->correl, acp->e_vectors, acp->e_values, dataset->ncols);
    prodmatMN(dataset->data_cr,acp->e_vectors,acp->princ_comp,dataset->n_obs,
		dataset->ncols,dataset->ncols);
    print_acp_results(stdout,*acp);
//---------- write results in output file --------------------------------------
    outputf=fopen(outfname, "w");
    print_dataset(outputf, dataset);
    print_datastats(outputf,dataset,datastats);
    print_acp_results(outputf, *acp);
    fclose(outputf);
    // Cleanup
    free(acp);
    free(datastats);
    free(dataset);
    return 0;
}

/**********************************************************************
 * print_acp_results(FILE* fptr, t_acp_data acp);
 *	prints the results of the analysis: eigenvalues, inertia,
 *	eigenvectors (optional),principal components (i.e. data points
 *	coordinates in eigenvectors coordinates system), and variables 
 *	coordinates for the correlation circles
 **********************************************************************/
 void	print_acp_results(FILE* fptr, t_acp_data acp) {

    print_eigenvalues(fptr, acp);
    print_eigenvectors(fptr, acp);
    print_princ_comp(fptr, acp);
    print_cercle_corr(fptr, acp);
    fprintf(fptr,"\n\n");
 }

/**********************************************************************
 * print_eigenvalues(FILE *fptr, t_acp_data acp);
 *	prints eigenvalues, their inertia and cumulated inertia in file *fptr
 **********************************************************************/
 void print_eigenvalues(FILE *fptr, t_acp_data acp) {
     int	j,n;
     double	ci;	// cumulated iniertia

     n=acp.data->ncols;
     fprintf(fptr,"\n#---- PCA ---\t");
     for(j=0;j<n;j++) {
	 fprintf(fptr,"PC%-3d\t",j+1);
     }
     fprintf(fptr,"\n# Eigenvalues:\t");
     for(j=0;j<n;j++) { 
	 fprintf(fptr,"%6.3f\t",acp.e_values[j]);
     }
     fprintf(fptr,"\n# Inertia:\t"); 
     for(j=0;j<n;j++) {
	 fprintf(fptr,"%3d%c\t",(int)(acp.inertia[j]*100),'%');
     }
     fprintf(fptr,"\n# Cum. inertia:\t");
     ci=0;
     for(j=0;j<n;j++) {
 	 ci+=acp.inertia[j];
 	 fprintf(fptr,"%3d%c\t",(int)(ci*100),'%');
     }
 }     

/**********************************************************************
 * print_eigenvectors(FILE *fptr, t_acp_data acp);
 *	prints the eigenvectors coordinates
 *	(directions of the principal components axes)
 *	This print out is optional, not used in further analysis
 **********************************************************************/
void    print_eigenvectors(FILE *fptr, t_acp_data acp){

    int i,j,n;
    n=acp.data->ncols;
    fprintf(fptr,"\n\n#---- Eigenvectors ----\t");
    for(i=0;i<n;i++) {
	fprintf(fptr,"\nPC%d\t\t",i+1);
	for(j=0;j<n;j++) {
	    fprintf(fptr,"%6.3f\t",acp.e_vectors[i][j]);
	}
    }
}

/**********************************************************************
* read_dataset: read data file and return a t_dataset structure
*
**********************************************************************/
int	read_dataset(t_fname fname, t_dataset *dataset) {

    FILE	*fichier_data; 
    char	s[LINE_LENGTH];
    int		nrows, ncols;
    char	*ptr;
    t_vector	*V;			//

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
    ncols=str2headers(s,dataset->headers)-1;  // first colunm are labels
    V=malloc(sizeof(*V));
    nrows=0;
    while (fgets(s, sizeof(s), fichier_data) != NULL) {
	if(nrows>=MAX_LINES) {
	    fprintf(stderr,"Error: too many lines in data file (max %d)\n",MAX_LINES);
	    exit(8);
	}
	if(s[0]=='#') continue;
        if(sscanf(s,"%s",dataset->obs_id[nrows])!=1) {
	    fprintf(stderr, "Error: observation label n°%d cannot be read.\n",nrows+1);
	    exit(8);
	}

    // ptr shift: eliminates the first label from s
	ptr=s+strlen(dataset->obs_id[nrows]);  

	if (str2hvector(ptr, *V, ncols) !=0) {
	    fprintf(stderr, "Erreur: lecture de la line %d de la matrice.\n",nrows);
		 exit(8);
	    }
	    matrix_line_from_vector(dataset->data, *V, nrows, ncols);
	nrows++;
    }
    // A longer title should be planned. For the moment long title= short title (Headers[0])
    strcpy(dataset->title,dataset->headers[0]);
    dataset->n_obs=nrows;
    dataset->ncols=ncols;   
    free(V);
    fclose(fichier_data);
    return 0;
 }

/**********************************************************************
 * compute_datastats(t_dataset *D,t_datastats *S)
 **********************************************************************/
void compute_datastats(t_dataset *D, t_datastats *S) {


    t_vector	*V;
    t_matrix	*T; // transposed intermediate
    int		i,j;
    double	norm;	//to divide matrix by number of data lines

    V=malloc(sizeof(*V));
    norm=1/(double)D->n_obs;
    //compute means dans standard deviation of each vector
    for(j=0;j<D->ncols;j++){
	vector_from_matrix(D->data,*V,D->n_obs,j);  
	S->means[j]=mean(*V,D->n_obs);
	S->sdevs[j]=sdev(*V,D->n_obs);
	for(i=0;i<D->n_obs;i++) {
	    // put centered reduced data in D->data_cr
	    D->data_cr[i][j] = (D->data[i][j] - S->means[j])/S->sdevs[j];
	}
    }
    T=malloc(sizeof(*T));
    if(T==NULL) {
	fprintf(stderr,"Error: not enough memory in compute_datastats. Abort.\n");
	exit(16);
    }
    transpose(D->data_cr,*T,D->n_obs,D->ncols);
    prodmatMN(*T,D->data_cr, S->covar,D->ncols,D->n_obs,D->ncols);
    free(T);
    prodmat_scalar(S->covar,norm,S->correl,D->ncols,D->ncols);
    free(V);
}
/**********************************************************************
* print_dataset(FILE *fptr, t_dataset *dataset)
**********************************************************************/
void print_dataset(FILE *fptr, t_dataset *dataset) {

    int		i,j;
    int		ncols, n_obs;
    t_header	*title, *long_title;

    //not efficient but better legibility
    ncols=dataset->ncols;
    n_obs=dataset->n_obs;
    title=&dataset->headers[0];
    long_title=&dataset->title;

    fprintf(fptr,"Data set: %s (%s)\n",*title, *long_title);
    fprintf(fptr,"Number of parameters:\t%d\n",ncols);
    fprintf(fptr,"Number of observations:\t%d\n",n_obs);
    fprintf(fptr,"\nVariables:  \t");
    for(i=1; i<=ncols; i++) {
	fprintf(fptr,"%s\t",dataset->headers[i]);
    }
    fprintf(fptr,"\n");
    for(i=0;i<n_obs;i++) {
	fprintf(fptr,"#%-3d %-8s\t",i+1,dataset->obs_id[i]);
	for(j=0;j<ncols;j++) {
	    fprintf(fptr,"%6.3f\t",dataset->data[i][j]);
	}
	fprintf(fptr,"\n");
    }
}
/**********************************************************************
* print_princ_comp(FILE *fptr, t_acp_data acp)
*	prints the coordinates of the observations in the eigenvectors 
*	coordinate system
**********************************************************************/
void print_princ_comp(FILE *fptr, t_acp_data acp) {

    int		i,j;
    int		ncols, n_obs;

    ncols=acp.data->ncols;
    n_obs=acp.data->n_obs;
                                                  
    fprintf(fptr,"\n\n#--- Principal Components -------");
    fprintf(fptr,"\n    axes     \t");
    for(j=0;j<ncols;j++) {
	 fprintf(fptr,"PC%-3d\t",j+1);
    }
    for(i=0;i<n_obs;i++) {
       fprintf(fptr,"\n%-8s\t",acp.data->obs_id[i]);
       for(j=0;j<ncols;j++) {
           fprintf(fptr,"%6.3f\t",acp.princ_comp[i][j]);
       }
    }
    fprintf(fptr,"\n");
} 
/**********************************************************************
 * print_cercle_corr(FILE *fptr, t_acp_data acp);
 **********************************************************************/
void	print_cercle_corr(FILE *fptr, t_acp_data acp) {

    int	    i,j;
    int	    n;

    n=acp.data->ncols;
    fprintf(fptr,"\n#--- Correlations Circle -------");
    for(i=0;i<n;i++) {
	fprintf(fptr,"\n%-8s\t",acp.data->headers[i+1]);
	for(j=0;j<n;j++) {
	    fprintf(fptr,"%6.3f\t",acp.e_vectors[i][j]*sqrt(acp.e_values[j]));
	}
    }
    fprintf(fptr,"\n");
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
	fprintf(fptr,"%s\t",D->headers[i]);
    }
    fprintf(fptr,"\n");
    fprintf(fptr,"# Mean:     \t");
    for(j=0;j<ncols;j++) {
	fprintf(fptr,"%6.3f\t",S->means[j]);
    }
    fprintf(fptr,"\n");
    fprintf(fptr,"# Std Dev:  \t");
    for(j=0;j<ncols;j++) {
	fprintf(fptr,"%6.3f\t",S->sdevs[j]);
    }
    fprintf(fptr,"\n");

/**************************** OPTIONAL. Prints centered reduced dataset
    // print centered reduced data matrix
    fprintf(fptr,"\n# Centered reduced data");
    fprintf(fptr,"\n------------\t");
    for(i=1; i<=ncols; i++) {
	fprintf(fptr,"%s\t",D->headers[i]);
    }
    fprintf(fptr,"\n");

    int	 n_obs;
    n_obs=D->n_obs;
    for(i=0;i<n_obs;i++) {
	fprintf(fptr,"#%-3d %s\t",i+1,D->obs_id[i]);
	for(j=0;j<ncols;j++) {
	    fprintf(fptr,"% 4.2f\t",D->data_cr[i][j]);
	}
	fprintf(fptr,"\n");
    }
// ------------- end optional printing ---------------------*/
    // print correlations matrix
    fprintf(fptr,"\n#Correlations\t");
    for(i=1; i<=ncols; i++) {
	fprintf(fptr,"%s\t",D->headers[i]);
    }
    fprintf(fptr,"\n");
    for(i=0;i<ncols;i++) {
	fprintf(fptr,"# %-10s\t",D->headers[i+1]);
	for(j=0;j<ncols;j++) {
	    fprintf(fptr,"%6.3f\t",S->correl[i][j]);
	}
	fprintf(fptr,"\n");
    }
}
 
 /**********************************************************************
 * eigenvalues:	calculates eigenvalues using QR algorithm
 *		return 0 if error
 *	    or
 *		returns number of iterations k
 **********************************************************************/
int eigenvalues(t_matrix A, t_vector EV, int order, double eps) {

    int		i,j,k,err_code;
    t_matrix	*Q,*R,*Ak;
    double	trace;
    double	residue;                
    
    Q=malloc(sizeof(t_matrix));
    R=malloc(sizeof(t_matrix));
    Ak=malloc(sizeof(t_matrix));

    err_code=qr_decomp(A,*Q,*R,order);
    if (err_code!= 0) {
	fprintf(stderr, "Error.\tSingular matrix in 'eigenvalues'.\n");
	return 0;}
    // Calculate lim Ak+1=RkQk
    // trace=norm
    trace=tracemat(*R,order);
    residue=trace;
    // residue= norm of bottom left triangle below the diagonal
    // stop when residue/trace < eps  (precision epsilon)
    for(k=0;fabs(residue/trace)>precision;k++){
	prodmat(*R,*Q,*Ak,order);
	qr_decomp(*Ak,*Q,*R,order);
	residue=0;
	for(j=0;j<order;j++) {
	    for(i=j+1;i<order;i++){
		residue+=fabs((*Ak)[i][j]);}
	}
    }
    for(i=0;i<order;i++) {
	EV[i]=(*Ak)[i][i];
    }
    free(Ak);
    free(R);
    free(Q);
    return k; // returns number of iterations
}
/**********************************************************************
 * ev_inertia: calculate the inertia of each eigenvalue
 **********************************************************************/ 
void    ev_inertia(t_vector e_values,t_vector inertia, int order) {
    int	i;
    for(i=0;i<order;i++) {
	inertia[i]=e_values[i]/order;
    }
}

/**********************************************************************
 * eigenvector: calculate eigenvector for eigenvalue r
 **********************************************************************/
void eigenvector(t_matrix A, t_vector EV, double r, int order) {

    t_matrix	*AA;	// AA=A-rId
    t_matrix	*M;     // (n-1) order matrix to resolve  M.X=B
    t_vector	*b;     // second member (order-1)
    t_vector	*X;	// t(EV)={1, {X}}
    int		i,j, err_code;

    AA=malloc(sizeof(t_matrix));
    M=malloc(sizeof(t_matrix));
    b=malloc(sizeof(t_vector));
    X=malloc(sizeof(t_vector));

    if(AA==NULL || M==NULL || b==NULL || X==NULL) {
	fprintf(stderr,"Error: not enough memory in eigenvector. Abort.\n");
	exit(16);
    }
    copymat(A,*AA,order,order);
    for(i=0;i<order;i++) { (*AA)[i][i] -=  r;};
    for(i=0;i<order-1;i++) {
	(*b)[i]=-(*AA)[i][0];
	for(j=0;j<order-1;j++) {
	    (*M)[i][j]=(*AA)[i][j+1];
	}
    }
//-----Code below to use QR decomposition for linear system resolution---
//  err_code=syslinQR(M,X,b,order-1);
//--------------Code below to use Gauss pivot system resolution-------
    err_code=syslinGauss(*M,*X,*b,order-1);
    if(err_code !=0) {
	fprintf(stderr, "Error:\tthe covariance matrix is singular, and has no Eigenvalues.\n ");
	fprintf(stderr, "\tAt least one variable is a linear combination of others.\n");
	fprintf(stderr, "\tEliminate correlated variables before running a PCA.\n");
	exit(8);
    }
    EV[0]=1;
    for(i=0;i<order-1;i++) {EV[i+1]=(*X)[i];};
    normalize(EV,order);
    free(X);
    free(b);
    free(M);
    free(AA);
    return;
}

/**********************************************************************
 * eigenvectors: calculate eigenvectors of A and stores them in matrix EV
 **********************************************************************/
void	eigenvectors(t_matrix A, t_matrix EV, t_vector e_values, int n) {

    int		i,j;
    t_vector	*V;

    V=malloc(sizeof(*V));
    for(j=0;j<n;j++){
	eigenvector(A,*V,e_values[j],n);
	for(i=0;i<n;i++){
	    EV[i][j]=(*V)[i];
	}
    }
    free(V);
}

/**********************************************************************
 * syslinQR:	solve linear system Ax=b using QR decomposition
 *		return non-zero int if error (matrix not regular)
 **********************************************************************/
int syslinQR(t_matrix A, t_vector x, t_vector b, int order) {

    t_matrix	*Q,*R,*tQ,*y,*B;
    double	s;
    int		i,j, err_code;

    Q=malloc(sizeof(t_matrix));
    R=malloc(sizeof(t_matrix));
    tQ=malloc(sizeof(t_matrix));
    y=malloc(sizeof(t_matrix));
    B=malloc(sizeof(t_matrix));

    if(Q==NULL || R==NULL || Q==NULL || y==NULL || B==NULL) {
	fprintf(stderr,"Error: not enough memory in syslinQR. Abort.\n");
	exit(16);
    } 
    err_code=qr_decomp(A,*Q,*R,order);
    if (err_code!= 0) {
	fprintf(stderr, "Error: La matrice n'est pas régulière.\n");
	return 1;}
    transpose(*Q,*tQ,order,order);     // calculate tQ
    // Change vector b into a matrix for prodmatMN
    for(i=0;i<order;i++){ (*B)[i][0]=b[i];  }
    prodmatMN(*tQ,*B,*y,order,order,1); // calculate y=tQ.b
    x[order-1]=(*y)[order-1][0]/(*R)[order-1][order-1];   // Xn = Yn / Rnn
    for(i=order-2;i>=0;i--) {
	s=0;
	for(j=i+1;j<order;j++) {
	    s+=(*R)[i][j]*x[j];
	}
	x[i]=((*y)[i][0]-s)/(*R)[i][i];
    }
    free(B);
    free(y);
    free(tQ);
    free(R);
    free(Q);
    return 0;
}

/**********************************************************************
 * syslinGauss:	solve linear system Ax=b using Gauss-Jordan pivot method
 *		return non-zero int if error (matrix not regular)
 **********************************************************************/
int syslinGauss(t_matrix A, t_vector x, t_vector b, int order) {

    t_matrix	*M;      // working matrix (A remains untouched)
    double	amax;	// largest pivot in column below line k
    double	alpha;	// normalisation factor
    double	tmp;
    int		i,j,k,l,km;

    M=malloc(sizeof(*M));
    if(M==NULL) {
	fprintf(stderr,"Error: not enough memory in syslinGauss. Abort.\n");
	exit(16);
    }
    copymat(A,*M,order,order);
    copyvect(b,x,order);
    for(k=0;k<order-1;k++) {
	// look for line l0>=k with largest pivotamax
	amax=0;
	km=0;
	for(l=k;l<order;l++) {
	    if(fabs((*M)[l][k])>amax) {
		amax=fabs((*M)[l][k]);
		km=l;
	    }
	}
	if (amax < precision*precision) {  //pivot close to zero: irregular matrix
	    fprintf(stderr, "Error:\tThe pivot is too small: %e - ", amax);
	    fprintf(stderr,         "\tThe linear system has no solution.\n");
	    return 2;
	}
	if (amax < precision) {  //pivot small: WARNING
	    fprintf(stderr, "\n**** WARNING: The pivot is very small: %-6.0e\n", amax);
	    fprintf(stderr,   "****          The matrix is close to singular.\n");
	    fprintf(stderr,   "****          The linear system may have no solution.\n");
	    /*return 1;*/
	}
	// switch lines k and km (in M and x)
	if( km!=k) {
	    for(j=k;j<order;j++) {                                                          
    		tmp=(*M)[km][j];
    		(*M)[km][j]=(*M)[k][j];
    		(*M)[k][j]=tmp;
    	    }
    	    tmp=x[km];x[km]=x[k];x[k]=tmp;
	}
	for(i=k+1;i<order;i++) {
	    alpha=(*M)[i][k]/(*M)[k][k];
	    // normalise line i to zero coeffs in column k except pivot
	    for(j=k+1;j<order;j++) {
		(*M)[i][j]=(*M)[i][j]-alpha*(*M)[k][j];
		}
	    (*M)[i][k]=0;
	    x[i]=x[i]-alpha*x[k];
	}
    }
    // M is now an upper triangular matrix 
    // calculate vector solution X starting from the bottom
    x[order-1]=x[order-1]/(*M)[order-1][order-1];
    for(k=order-2;k>=0;k--) {
    	for(j=k+1;j<order;j++){
	    x[k]-=(*M)[k][j]*x[j];
	}
	x[k]=x[k]/(*M)[k][k];
    }                                 
    free(M);
    return 0;
}

/**********************************************************************
 * print_matrix
 **********************************************************************/
void print_matrix(t_matrix M, int order, char *msg) {
    int i;
    int j;
    char spacer[]={"_____________________________"};
    printf("%s %8s     %s\n",spacer, msg, spacer);
    for(i=0;i<order;i++) {
    		 for(j=0; j<order;j++) {
    		     printf("%6.3f\t", M[i][j]);}
    		 printf("\n");}
}
/**********************************************************************
 * print_matrixMN (print non square matrix)
 **********************************************************************/
void print_matrixMN(t_matrix M, int n, int p, char *msg) {
    int i;
    int j;
    char spacer[]={"_____________________________"};
    printf("%s %8s     %s\n",spacer, msg, spacer);
    for(i=0;i<n;i++) {
    		 for(j=0; j<p;j++) {
    		     printf("%6.3f\t", M[i][j]);}
    		 printf("\n");}
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
 * str2hvector: converts the string into a vector of N=order floating point numbers
 **********************************************************************/
int str2hvector(char *linestr, t_vector v, int order) {
    char s[LINE_LENGTH];// string used to build the number
    int		 j;     // indice de position du caractère numérique dans le nombre      
    int		 k;	// indice de position dans la chaine
    int		 l;     // indice de position dans le vecteur horizontal
    
//printf("___________________\n%s\n",linestr);

    s[0]='\0';
    k=0;                        
    l=0;
    while(k<LINE_LENGTH){
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
        if (sscanf(s, "%lf", &v[l]) != 1) {
// ATTENTION, if non digit characters follow digits, their are ignored.
		     printf("Error: %s is not a floating point number.\n", s);
    		 return(1);
        }
//printf("s: %s\tx:%4.2f\n",s,v[l]);
        l++;
};
return 0;  /* the line doesn't have enough float values */
}
 
/**********************************************************************
 *  matrix_line_from_vector: copie le vecteur V dans la ligne 'lineNr'
 *		 	     de la matrice d'ordre 'order'
 **********************************************************************/
void matrix_line_from_vector(t_matrix M, t_vector V, int lineNr, int order)
{
    int j;
    for (j=0; j<order; j++) {
		 M[lineNr][j]=V[j];
    }
    return;
}
/**********************************************************************
 * vector_from_matrix(Matrix, Vector, NbRows, ColNr)
 *	    copy column coln of M into vector V
 **********************************************************************/
void	vector_from_matrix(t_matrix M,t_vector V,int nrows,int coln){

    int	    i;
    for(i=0;i<nrows;i++)  V[i]=M[i][coln];
    return;
}

/**********************************************************************
 * prodmat: produit de deux matrices carrées d'ordre 'order'
 * *******************************************************************/
void prodmat(t_matrix A, t_matrix B, t_matrix P, int order)
{
    int	    i, j, k;
    double  aij;
    for (i=0; i<order; i++) {
		 for (j=0; j<order; j++) {
		     aij=0;
		     for (k=0;k<order;k++) {
		 		 aij=aij+A[i][k]*B[k][j];
		     }
		     P[i][j]=aij;
		 }
    }
    return;
}
/**********************************************************************
 * prodmat_scalar: product of a scalar a with a matrix M  into aM
 **********************************************************************/
void prodmat_scalar(t_matrix M, double a, t_matrix aM, int n, int p) {
    int i,j;
    for(i=0;i<n;i++){
	for(j=0;j<p;j++) {
	    aM[i][j]=a*M[i][j];
        }
    }
}
/**********************************************************************
 * prodmatMN: product of two matrices A(m,n) x B(n,p) = P(m,p)
 * *******************************************************************/

void prodmatMN(t_matrix A, t_matrix B, t_matrix P, int m, int n, int p)
{
    int	i, j, k;
    double  aij;

    for (i=0; i<m; i++) {
		 for (j=0; j<p; j++) {
		     aij=0;
		     for (k=0;k<n;k++) {
		 		 aij=aij+A[i][k]*B[k][j];
		     }
		     P[i][j]=aij;
		 }
    }
    return;
}

/**********************************************************************
 * prodVtV: product of a vector with its transpose P = V x t(V)
 * *******************************************************************/
void   prodVtV(t_vector V, t_matrix P, int order) {
    int i;
    t_matrix *U,*tV;

    U=calloc(1,sizeof(t_matrix));
    tV=calloc(1,sizeof(t_matrix));

    if (U==NULL || tV == NULL) {
	fprintf(stderr,"Error: not enough memory in prodVtV\n");
	exit(16);
    }
    for (i=0;i<order;i++) {
	(*U)[i][0]=V[i];
	(*tV)[0][i]=V[i];
    }
    prodmatMN(*U,*tV,P,order,1,order);

    free(tV);
    free(U);
    return;
}


/**********************************************************************
 * Transpose: transpose matrix M and returns into matrix T
 * *******************************************************************/
void transpose(t_matrix M, t_matrix T,int rows, int cols) {
    int i,j;
    for(i=0;i<cols;i++) {
		 for(j=0; j<rows; j++) {
		     T[i][j]=M[j][i];
		 }
    }
}

/**********************************************************************
 * copy matrix M into C
 * *******************************************************************/
void copymat(t_matrix M, t_matrix C, int rows, int cols) {
    int i,j;

    for(i=0;i<rows;i++) {
		 for(j=0;j<cols;j++) {
		     C[i][j]=M[i][j];
		 }
    }
}

/**********************************************************************
 * copyvect copy vector U into V
 * *******************************************************************/
void copyvect(t_vector U, t_vector V, int n){
    int	    i;
    for (i=0;i<n;i++) { V[i]=U[i];}
}

 /**********************************************************************
 * Id: returns Identity matrix of order n
 * *******************************************************************/
void Id(t_matrix I, int n) {
    int i,j;
    for(i=0;i<n;i++) {
		 for(j=0;j<n;j++) {
		     if(i==j) {I[i][j]=1;}
		     else {I[i][j]=0;}
		 }
    }
}

/**********************************************************************
 * sign: returns the sign(integer +1 or -1) of parameter
 * 	 requires <math.h> and compilation with -lm
 * *******************************************************************/
int sign(double x) {
	if (x==0) { return 1;}
	else {return nearbyint(fabs(x)/x);}
}

/**********************************************************************
 * norm: returns the norm of a vector, zero if null vector
 * 	 requires <math.h> and compilation with -lm
 * *******************************************************************/
double norm(t_vector V, int n) {
	int i;
	double s=0;
	for(i=0;i<n;i++) {s=s+V[i]*V[i];}
	return sqrt(s);
}

/**********************************************************************
 * normalize:	vector is normalized to modulus 1.
 *		returns 0 if OK
 *		returns 2 if modulus is zero (null vector)
 *		If null vector, error management must be done by caller
 * *******************************************************************/
int normalize(t_vector V, int n) {
	int i;
	double s;
	s=norm(V,n);
	if (s==0) { return 2;} 
	else {
		for(i=0;i<n;i++) {
			V[i]=V[i]/s;
		}
	return 0;
	}
}

/**********************************************************************  
 * tracemat: returns the trace of matrix M
 *********************************************************************/ 
double tracemat(t_matrix M, int n) {
    int i;
    double t=0;
    for(i=0;i<n;i++) { t+=M[i][i];}
    return t;
}

/**********************************************************************
 * qr_decomp: decompose matrix A into Q & R such as A=QR
 *	using Householder elimination method
 *	Q is Hermitian tQ x Q = I
 *	R is a triangular superior matrix
 *	R eigenvalues are the eigenvalues of A
 * *******************************************************************/
int  qr_decomp(t_matrix A, t_matrix Q, t_matrix R, int n) {
    int i,j,k;
    double   alpha;		// signed norm of vector x1
    t_vector *vector_e1;	// first vector of orthonormal basis
    t_vector *vector_x1;		// first column of Householder sub-matrix
    t_vector *vector_u;		// u=x1-alpha.e1, normalized
    t_matrix *UtU;		// Hk=Id - 2.U.t(U)
    t_matrix *H;			// Householder intermediate matrices
    t_matrix *HH;		// temporary matrix to calculate and copy H
    t_matrix *I;			// Identity matrix

    vector_e1=calloc(1,sizeof(t_vector));
    vector_x1=calloc(1,sizeof(t_vector));
    vector_u=malloc(sizeof(t_vector));
    UtU	= malloc(sizeof(t_matrix));
    H	= malloc(sizeof(t_matrix));
    HH	= malloc(sizeof(t_matrix));
    I	= malloc(sizeof(t_matrix));

    if(vector_e1==NULL ||vector_x1==NULL || vector_u==NULL || UtU==NULL
	    || H==NULL || HH==NULL || I==NULL){
	fprintf(stderr,"Error: not enough memory in qr_decomp. Abort.\n");
	exit(16);
    }
    (*vector_e1)[0]=1;
    copymat(A,R,n,n);
    Id(Q,n);	// 1st iteration Q starts as Identity
    for(k=0;k<n-1;k++) {
	//N.B.	indices for Hk submatrix range from k to n-1
	//	the order of vectors of this matrix is (n-k)
	for(i=0;i<n-k;i++) {
	    (*vector_x1)[i]=R[(k+i)][k];
	}
	alpha=norm(*vector_x1, n-k) * sign((*vector_x1)[0]);
	for(i=0; i<n-k;i++) {
	    (*vector_u)[i]=(*vector_x1)[i] - alpha * (*vector_e1)[i];
	}
	if(normalize(*vector_u,n-k)!=0) {
	    fprintf(stderr,"Error: Null vector during QR decomposition. The matrix must be singular.\n");
	    return 2;
	}
	// vector_u x Transpose(vector_u) gives a quare matrix 
	prodVtV(*vector_u,*UtU,n-k);
	// Fill the bottom right of H with calculated Hk submatrix
	Id(*H,n);
	Id(*I,n-k);
	for(i=0;i<n-k;i++) {
	    for(j=0;j<n-k;j++) {
		(*H)[i+k][j+k]=(*I)[i][+j] - 2*(*UtU)[i][j];
	    }
	}
	prodmat(*H,R,*HH,n);
        copymat(*HH,R,n,n);
        prodmat(Q,*H,*HH,n);
	copymat(*HH,Q,n,n);
    }
    free(vector_e1);
    free(vector_x1);
    free(vector_u);
    free(I);
    free(HH);
    free(H);
    free(UtU);
 return 0;
}

/**********************************************************************
 *     BASIC STATISTIC FUNCTIONS
 *
 **********************************************************************/
//------------------------------ mean(t_vector V, int n);
double  mean(t_vector V, int n) {
    int	    i;
    double  s=0;
    for(i=0;i<n;i++) s+=V[i];
    return s/n;
}
//------------------------------ sumsquares(t_vector V, int n);
double	sumsquares(t_vector V, int n) {
    int i;
    double s=0;
    for(i=0;i<n;i++) s+=V[i]*V[i];
    return s;
}
//--------------------------------variance of V (entire population)
double	variance(t_vector V, int n){
    double  m;
    m=mean(V,n);
    return sumsquares(V,n)/n-m*m;
}
//-------------------------- std deviation of V (entire population)
double  sdev(t_vector V, int n){
    return sqrt(variance(V,n));
}
//-------------------------------variance of V, as a sample (/n-1)
double	variance_sample(t_vector V, int n){
    return variance(V,n)*n/(n-1);
}
//-------------------------------std deviation of V  as a sample
double  sdev_sample(t_vector V, int n){
    return sqrt(variance_sample(V,n));
}




/************************** END OF PROGRAM *******************************/
