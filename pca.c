/******************************************************************************  
 *  pca.c  
 *  Performs Principal Component Analysis (PCA) on a set of data */
 #define	VERSION "2.0"
/*  Author:	Bruno S. Charrière
 *  Date:	9.06.2014
 *****************************************************************************/
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "dbg.h"	// error handling functions
#include "libal.h"	// linear algebra & statistics library


#define DEFAULT_NAME   "matrice.txt"// default name of input data file.
#define LINE_LENGTH    1000	    // max length of data file lines
#define NAME_LENGTH    255	    // length of strings & labels

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

char	*g_version=VERSION;

/**********************************************************************
 *		    M A I N
 **********************************************************************/
int main( int argc, char *argv[]) {        // args not used so far

    FILE	*outputf;
    char	 data_fname[NAME_LENGTH] = DEFAULT_NAME;
    char	*outfname="dataset.txt";
    t_dataset	*dataset;
    t_datastats	*datastats;
    t_acp_data	*acp;


}
