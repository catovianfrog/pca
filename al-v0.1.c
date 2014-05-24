/******************************************************************************  
 *  qr.c  
 *  version:	1.07
 *  Author:	Bruno S. Charrière
 *  Date:	07.04.2014
 ****************************************************************************
 * CHANGE LOG
 * v0.1  One dimentional matrices
 * v0.7  QR decomposition finally works!
 * v0.8  Iterative qr decomposition gives eigenvalues
 * v1.01 Creation of t_matrix & t_vector types
 * v1.02 rewrite from one dimension to two-dimension matrices
 * v1.03 Linear system resolution using QC decomposition (gauss pivot yet to do)
 * v1.04 Calculates also eigenvectors
 * v1.05 Clean-up main - move activity blocks into functions
 * v1.06 Linear system resolution using Gauss pivot 
 *
 *****************************************************************************/
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#define	VERSION	       "1.07"
#define DEFAULT_NAME   "matrice.txt"
#define LINE_LENGTH    1000
#define NAME_LENGTH    255
#define MAX_ORDER      20
#define MAX_ITERATION  1e5
#define PRECISION      1e-6

typedef double t_matrix[MAX_ORDER][MAX_ORDER];
typedef double t_vector[MAX_ORDER];

char	*version=VERSION;
int	sign(double x); 
int	str2hvector(char *s, t_vector v, int order);
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
void	print_matrix(t_matrix M, int order, char *msg);
void	prodmat(t_matrix A, t_matrix B, t_matrix P, int order);
void	prodVtV(t_vector V, t_matrix P, int order);
void	prodmatMN(t_matrix A, t_matrix B, t_matrix P, int m, int n, int p);
int	syslinQR(t_matrix A, t_vector x, t_matrix b, int order); // solve Ax=b using QR decomposition
int	syslinGauss(t_matrix A, t_vector x, t_vector b, int order);
void	eigenvector(t_matrix A, t_vector EV, double r, int order);

//------------------------------------------------------------------------------
//					M A I N 
//------------------------------------------------------------------------------
//
int main( int argc, char *argv[]) {
    FILE	*fichier_data;
    char	nom_fichier_data[NAME_LENGTH] = DEFAULT_NAME;
    int		order;
    char	line_read[LINE_LENGTH];
    int		i,j,k;
    int		err_code;
    int		iterations;
    t_vector	V,e_values;
    t_matrix	matrix, EVmatrix;


    printf("version:\t %s\n",version);
    fichier_data=fopen(nom_fichier_data, "r");
    if (fichier_data == NULL) {
	fprintf(stderr, "Erreur: impossible d'ouvrir le fichier %s.\n", nom_fichier_data);
	return(1);
    }
    printf("Data read from %s\n", nom_fichier_data);
// read matrix order, and then read matrix line by line
    fgets(line_read, sizeof(line_read), fichier_data);
    sscanf(line_read,"%d", &order);
    if (order == 0) {
	fprintf(stderr, "Erreur: la matrice est vide");
	return(1);
    }
    // read the matrix line by line
    for(i=0; i<order; i++){
	fgets(line_read, sizeof(line_read), fichier_data);
	if(line_read[0] == '#') { i--; continue;}  //ignore comment lines, with # in first column
	if (str2hvector(line_read, V, order) !=0) {
	    fprintf(stderr, "Erreur: lecture de la line %d de la matrice.\n",i);
		 exit(1);
	    }
	    matrix_line_from_vector(matrix, V, i, order);
	}

    iterations=eigenvalues(matrix,e_values,order,precision);
    if(iterations == 0) {return 1;}
    // quits returning 1 if error (e.g. matrix is singular)
    // returns residue in 'precision'
    // returns 0 if error, or number of iterations


    printf(  "\nIterations: %7d\nResidue:\t%12.8f\n",iterations,precision); //precision=residue/trace
    char spacer[]={"_____________________________"};
    printf("%s %8s     %s\n",spacer, "Eigenvalues ", spacer);
    for(i=0;i<order;i++) {
	printf("%8.3f\t", e_values[i]);
    }
    printf("\n");
    
    for(k=0;k<order;k++) {
	eigenvector(matrix,V,e_values[k],order);
	for(i=0; i<order;i++) {EVmatrix[i][k]=V[i];}
    }
    print_matrix(EVmatrix,order,"Eigenvectors");
    printf("\n");


// TRY GAUSS PIVOT linear system resolution
    t_matrix u={{1},{2},{3},{4}};
    t_matrix prod;
    t_vector b;
    prodmatMN(matrix,u,prod,order,order,1);

    for(i=0;i<order;i++){b[i]=prod[i][0];}
    printf("\n____ Equation second member________\n");
    for(i=0;i<order;i++){printf("b(%d)=%f\t",i,b[i]);}

    err_code=syslinGauss(matrix,V,b,order);
    if(err_code!=0){return 1;}
    printf("\n\n____ Equation solution__(Gauss)__\n");
    for(i=0;i<order;i++){printf("x(%d)=%f\t",i,V[i]);}
    printf("\n\n");

    err_code=syslinQR(matrix,V,prod,order);
    if(err_code!=0){return 1;}
    printf("\n\n____ Equation solution__(QR)_____\n");
    for(i=0;i<order;i++){printf("x(%d)=%17.9f\t",i,V[i]);}
    printf("\n\n");



    fclose(fichier_data);
    return 0;
}

 /**********************************************************************
 * eigenvalues:	calculates eigenvalues using QR algorithm
 *		return 0 if error
 *	    or
 *		returns number of iterations k
 **********************************************************************/
int eigenvalues(t_matrix A, t_vector EV, int order, double eps) {

    int	i,j,k,err_code;
    t_matrix	Q,R,Ak;
    double	trace;
    double	residue;                
    
    print_matrix(A, order,"A");
    err_code=qr_decomp(A,Q,R,order);
    if (err_code!= 0) {
	fprintf(stderr, "La matrice n'est pas régulière.\n");
	return 1;}
    print_matrix(Q, order,"Q");
    print_matrix(R, order,"R");
    // Calculate lim Ak+1=RkQk
    // trace=norm
    trace=tracemat(R,order);
    residue=trace;
    // residue= norm of bottom left triangle below the diagonal
    // stop when residue/trace < eps  (precision epsilon)
    for(k=0;fabs(residue/trace)>precision;k++){
	prodmat(R,Q,Ak,order);
	qr_decomp(Ak,Q,R,order);
	residue=0;
	for(j=0;j<order;j++) {
	    for(i=j+1;i<order;i++){
		residue+=fabs(Ak[i][j]);}
	}
    }
    for(i=0;i<order;i++) {
	EV[i]=Ak[i][i];
    }
    return k; // returns number of iterations
}
 


/**********************************************************************
 * eigenvector: calculate eigenvector for eigenvalue r
 **********************************************************************/
void eigenvector(t_matrix A, t_vector EV, double r, int order) {

    t_matrix	AA;	// AA=A-rId
    t_matrix	M;      // (n-1) order matrix to resolve  M.X=B
    t_matrix	B;	// second member (order-1)
    t_vector	X;	// t(EV)={1, {X}}
    int		i,j, err_code;

    copymat(A,AA,order,order);
    for(i=0;i<order;i++) { AA[i][i] -=  r;};
    for(i=0;i<order-1;i++) {
	B[i][0]=-AA[i][0];
	for(j=0;j<order-1;j++) {
	    M[i][j]=AA[i][j+1];
	}
    }
    err_code=syslinQR(M,X,B,order-1);
    EV[0]=1;
    for(i=0;i<order-1;i++) {EV[i+1]=X[i];};
    normalize(EV,order);
    return;
}



/**********************************************************************
 * syslinQR:	solve linear system Ax=b using QR decomposition
 *		return non-zero int if error (matrix not regular)
 **********************************************************************/
int syslinQR(t_matrix A, t_vector x, t_matrix b, int order) {

    t_matrix	Q,R,tQ,y;
    double	s;
    int		i,j, err_code;

    err_code=qr_decomp(A,Q,R,order);
    if (err_code!= 0) {
	fprintf(stderr, "La matrice n'est pas régulière.\n");
	return 1;}
    transpose(Q,tQ,order,order);     // calculate tQ
    prodmatMN(tQ,b,y,order,order,1); // calculate y=tQ.b
    x[order-1]=y[order-1][0]/R[order-1][order-1];   // Xn = Yn / Rnn
    for(i=order-2;i>=0;i--) {
	s=0;
	for(j=i+1;j<order;j++) {
	    s+=R[i][j]*x[j];
	}
	x[i]=(y[i][0]-s)/R[i][i];
    }
    return 0;
}

/**********************************************************************
 * syslinGauss:	solve linear system Ax=b using Gauss-Jordan pivot method
 *		return non-zero int if error (matrix not regular)
 **********************************************************************/
int syslinGauss(t_matrix A, t_vector x, t_vector b, int order) {

    t_matrix	M;      // working matrix (A remains untouched)
    double	amax;	// largest pivot in column below line k
    double	alpha;	// normalisation factor
    double	tmp;
    int		i,j,k,l,km, err_code=0;

    copymat(A,M,order,order);
    copyvect(b,x,order);
    for(k=0;k<order-1;k++) {
	// look for line l0>=k with largest pivotamax
	amax=0;
	km=0;
	for(l=k;l<order;l++) {
	    if(fabs(M[l][k])>amax) {
		amax=fabs(M[l][k]);
		km=l;
	    }
	}
	if (amax < precision) {  //pivot close to zero: irregular matrix
	    fprintf(stderr, "Error: irregular matrix - no solution.\n");
	    return 1;}
	// switch lines k and km (in M and x)
	if( km!=k) {
	    for(j=k;j<order;j++) {
    		tmp=M[km][j];
    		M[km][j]=M[k][j];
    		M[k][j]=tmp;
    	    }
    	    tmp=x[km];x[km]=x[k];x[k]=tmp;
	}
	for(i=k+1;i<order;i++) {
	    alpha=M[i][k]/M[k][k];
	    // normalise line i to zero coeffs in column k except pivot
	    for(j=k+1;j<order;j++) {
		M[i][j]=M[i][j]-alpha*M[k][j];
		}
	    M[i][k]=0;
	    x[i]=x[i]-alpha*x[k];
	}
    }
    // M is now an upper triangular matrix 
    // calculate vector solution X starting from the bottom
    x[order-1]=x[order-1]/M[order-1][order-1];
    for(k=order-2;k>=0;k--) {
    	for(j=k+1;j<order;j++){
	    x[k]-=M[k][j]*x[j];
	}
	x[k]=x[k]/M[k][k];
    }
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
    		     printf("%8.3f\t", M[i][j]);}
    		 printf("\n");}
}



/**********************************************************************
 * str2hvector: convertit la chaîne en une série de N=order nombres double précision
 **********************************************************************/
int str2hvector(char *linestr, t_vector v, int order) {
    char s[LINE_LENGTH];// string used to build the number
    int		 j;     // indice de position du caractère numérique dans le nombre      
    int		 k;	// indice de position dans la chaine
    int		 l;     // indice de position dan sle vecteur horizontal
    
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
        l++;
		 if (l == order) return 0; /* read the right number of float values */
};
return 1;  /* the line doesn't have enough float values */
}
 
/**********************************************************************
 *  matrix_line_from_vector: copie le vecteur V dans la ligne 'lineNr'
 *		 		 		      de la matrice d'ordre 'order'
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
 * prodmatMN: produit de deux matrices A(m,n) x B(n,p) = P(m,p)
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
    int i,j;
    double aij;
    t_matrix U,tV;

    for (i=0;i<order;i++) {
	U[i][0]=V[i];
	tV[0][i]=V[i];
    }
    prodmatMN(U,tV,P,order,1,order);
    return;
}


/**********************************************************************
 * Transpose: transpose matrix M and returns into matrix T
 * *******************************************************************/
void transpose(t_matrix M, t_matrix T,int rows, int cols) {
    int i,j;
    for(i=0;i<rows;i++) {
		 for(j=0; j<cols; j++) {
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
 * normalize:	returns vector of modulus 1 (no change if null vector)
 * *******************************************************************/
 
int normalize(t_vector V, int n) {
	int i;
	double s;
	s=norm(V,n);
	// need error management if s=0
	if (s==0) { return 1;} 
	else {
		for(i=0;i<n;i++) {
			V[i]=V[i]/s;
		}
	return 0;
	}
}

/******************************************************************************  
 * tracemat: returns the trace of matrix M
 *****************************************************************************/ 
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
    int i,j,k,sgn,err;
    double   alpha;		// signed norm of vector x1
    t_vector vector_e1={1,0};	// first vector of orthonormal basis
    t_vector vector_x1;		// first column of Householder sub-matrix
    t_vector vector_u;		// u=x1-alpha.e1, normalized
    t_matrix UtU;		// Hk=Id - 2.U.t(U)
    t_matrix H;			// Householder intermediate matrices
    t_matrix HH;		// temporary matrix to calculate and copy H
    t_matrix I;			// Identity matrix

    copymat(A,R,n,n);
    Id(Q,n);	// 1st iteration Q starts as Identity
    for(k=0;k<n-1;k++) {
	//N.B.	indices for Hk submatrix range from k to n-1
	//	the order of vectors of this matrix is (n-k)
	for(i=0;i<n-k;i++) {
	    vector_x1[i]=R[(k+i)][k];
	}
	alpha=norm(vector_x1, n-k) * sign(vector_x1[0]);
	for(i=0; i<n-k;i++) {
	    vector_u[i]=vector_x1[i] - alpha * vector_e1[i];
	}
	normalize(vector_u,n-k);
	// vector_u x Transpose(vector_u) gives a quare matrix 
	prodVtV(vector_u,UtU,n-k);
	// Fill the bottom right of H with calculated Hk submatrix
	Id(H,n);
	Id(I,n-k);
	for(i=0;i<n-k;i++) {
	    for(j=0;j<n-k;j++) {
		H[i+k][j+k]=I[i][+j] - 2*UtU[i][j];
	    }
	}
	prodmat(H,R,HH,n);
        copymat(HH,R,n,n);
        prodmat(Q,H,HH,n);
	copymat(HH,Q,n,n);
    }
 return 0;
}

	
	
