#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h> 

#define MAX_ITERATION  1e5	    // max number of iterations for eigenvalues
#define PRECISION      1e-6	    // stop iterations when residue < PRECISION
#define LINE_LENGTH    1000	    // max length of data file lines
#define NAME_LENGTH    255	    // length of strings & labels
#define MAX_LINES      350	    // max number of observations (individus)
#define MAX_COLS       MAX_LINES    // max number of parameters
//*********************************************************************
typedef	struct  {   int	    nrows;
		    int	    ncols;
		    double  *data; }	t_matrix;
/**********************************************************************
 * TODO
 *	    - reset relevant function parameters as 'const'
 **********************************************************************/
const double	g_precision=PRECISION;

//---------- Functions prototypes ----------------------------------------------
t_matrix*   matrix_new(int nrows, int ncols);
void	    matrix_free(t_matrix *m);
t_matrix*   matrix_new_vector(int n, void *p);
inline double matrix_get(t_matrix *M, int i, int j); // body must be in header
inline int  matrix_set(t_matrix *M, const int i, const int j, const double r);
void	    matrix_assign(t_matrix *m, void *p); 
void	    matrix_printf(FILE *f, t_matrix *M, char* format, char* msg);
void	    matrix_print(FILE *f, t_matrix *m, char* msg);
t_matrix*   matrix_prod(const t_matrix *A, const t_matrix *B);
t_matrix*   matrix_add(const t_matrix *A, const t_matrix *B, const double scale);
t_matrix*   matrix_copy(const t_matrix *M);
int	    matrix_move(t_matrix *Dest, const t_matrix *Source);
t_matrix*   matrix_scale(const t_matrix *M, const double scalar);
t_matrix*   matrix_Id(const int nrows, const int ncols);
t_matrix*   matrix_transpose(const t_matrix *M);
double	    matrix_norm(const t_matrix *M);
t_matrix*   matrix_normalize(t_matrix *M);
double	    matrix_trace(const t_matrix* M);
double	    matrix_lower_residue(const t_matrix *M);
t_matrix*   matrix_get_vector(const t_matrix *M, const int n);
int	    matrix_set_vector(t_matrix *M,const int n, const t_matrix *V);
int	    matrix_qr_decomp(t_matrix *A, t_matrix *Q, t_matrix *R);
int	    syslinQR(t_matrix *A, t_matrix *X, t_matrix *B) ;
int	    syslinGauss(t_matrix const *A, t_matrix *X, t_matrix const *B) ;
double	    sign(double x);
//------------------------------------------------------------------------------

/**********************************************************************
 *  matrix_new: returns a new matrix, or NULL if not enough memory
 *		all new matrices and vectors are initialized with 0
 **********************************************************************/
t_matrix*    matrix_new(int nrows, int ncols) {
    t_matrix	*matrix;
    matrix=calloc(1,sizeof(t_matrix));
    matrix->data=calloc(nrows*ncols, sizeof(double));
    if(matrix->data==NULL) return NULL;
    matrix->nrows=nrows;
    matrix->ncols=ncols;
    return matrix;
}
/**********************************************************************
 * matrix_free	    (works also for vectors)
 **********************************************************************/
void	matrix_free(t_matrix *m){   
	if(m==NULL) return;
        free(m->data);
	free(m);
}
/**********************************************************************
 * matrix_new_vector: creates a new vector(n). returns a t_matrix[*,1]
 *	       fills it from *p if p not NULL
 **********************************************************************/
t_matrix*   matrix_new_vector(int n, void *p) {
        t_matrix    *v;
        v=matrix_new(n,1);
	if(p !=NULL) memcpy(v->data, p, sizeof(*v->data));
	return v;
}
/**********************************************************************
 * matrix_set(M,i,j,r)  ==>  M[i,j]=r                    
 *  in case you can't remember the normal syntax M->data[i*M->ncols+j]=r; 
 **********************************************************************/
inline int matrix_set(t_matrix *M, const int i, const int j, const double r) {
    if(i>=M->nrows || j>=M->ncols) return 1;
    M->data[i*M->ncols+j]=r;
    return 0;
}
/**********************************************************************
 * matrix_get: aij = matrix_get(M,i,j);
 **********************************************************************/
inline double matrix_get(t_matrix *M, int i, int j) {
    if(i>=M->nrows || j>=M->ncols) return nan("");
    return M->data[i*M->ncols+j];
}
/**********************************************************************
 * matrix_assign       SHOULD BE IN THE HEADER (inline)
 **********************************************************************/
void	matrix_assign(t_matrix *m, void *p) {
    memcpy(m->data,p,m->nrows*m->ncols*sizeof(double));
}
/**********************************************************************
 *  matrix_printf(FILE *f, t_matrix *M, char* format, char* msg);
 *
 *   Prints out the matrix m with message 'msg' as a header line.
 *   The string 'format' is used if not empty and valid formating 
 *   string for a doubleprecision real. If 'format' is empty or invalid,
 *   default numbers printout format is "% 8.2g\t"
 **********************************************************************/
void matrix_printf(FILE *f, t_matrix *M, char* format, char* msg) {
    int	    i,j;
    char    teststr[50];
    //char    form[20]="% 8.2g";
    char    form[20]="% 9.3f";
    char    *spacer=" __________________________ ";

    //To check validity of the format, use snprintf to compute the result
    //in teststr, and then test the code returned by snprintf is positive
    if(strlen(format)>0 && snprintf(teststr,40,format,(double)1.23-5)>0) {
	strncpy(form, format, sizeof(form)-1);
    }
    strcat(form, "\t");
    fprintf(f,"%s%12.12s%s\n",spacer,msg,spacer);
    for(i=0;i<M->nrows;i++) {
	for(j=0;j<M->ncols;j++) {
	    printf(form,matrix_get(M,i,j));
	}
	printf("\n");
    }
}
/**********************************************************************
 *  matrix_print
 *   prints out the matrix with message 'msg' as header line
 *   Numbers printout format is "% 8.2g\t"
 **********************************************************************/
void matrix_print(FILE *f, t_matrix *m, char* msg) {
    matrix_printf(f,m,"",msg);
}
/**********************************************************************
 * matrix_prod
 * returns NULL if nrows(B) != ncols(A), or return matrix product A x B
 **********************************************************************/
t_matrix* matrix_prod(const t_matrix *A, const t_matrix *B) {
    int	    i,j,k,m,n,p;
    double  aij;
    t_matrix	*P;

    m=A->nrows;
    n=A->ncols;
    p=B->ncols;
    if(B->nrows !=n) return NULL;
    P=matrix_new(m,p);
    if(P==NULL) return P; 
    for(i=0;i<m;i++) {
	for(j=0;j<p;j++) {
	    aij=0;
	    for(k=0;k<n;k++) {
		aij=aij + A->data[i*n+k] * B->data[k*p+j]; // A[i,k].B[k,j]
	    }
	    P->data[i*p+j]=aij;
	}
    }
    return P;
}
/**********************************************************************
 * matrix_add
 * returns matrix S=A+scale.B
 **********************************************************************/
t_matrix*   matrix_add(const t_matrix *A, const t_matrix *B, const double scale) {
    t_matrix	*S;
    int		i,j;

    if(A->nrows!=B->nrows || A->ncols!=B->ncols) return NULL;
    S=matrix_new(A->nrows, A->ncols);
    if(S!=NULL){
	for(i=0;i<A->nrows;i++) {
	    for(j=0;j<A->ncols;j++) {
		S->data[i*A->ncols+j] = A->data[i*A->ncols+j] + scale * B->data[i*A->ncols+j];
	    }
	}
    }
    return S;
}
/**********************************************************************
 * matrix_copy
 **********************************************************************/
t_matrix*   matrix_copy(const t_matrix *M) {

    t_matrix	*C;
    C=matrix_new(M->nrows,M->ncols);
    if(C!=NULL) memcpy(C->data,M->data,M->nrows*M->ncols*sizeof(double));
    return C;
}
/**********************************************************************
 * matrix_move(*Dest, *Source): copy the content of *Source into *Dest,
 * Source and Dest already exist. Memory is neither allocated nor freed
 * returns 0 if error, or byte copied
 **********************************************************************/
int    matrix_move(t_matrix *Dest, const t_matrix *Source) {
       int  nbytes=0;
       if(Source->nrows == Dest->nrows && Source->ncols==Dest->ncols) {
	   nbytes=Source->nrows*Source->ncols*sizeof(double);
	   memcpy(Dest->data, Source->data,nbytes); 
	   // no need to copy ncols and nrows, they are the same
       }
       return nbytes;
}
/**********************************************************************
 * matrix_scale  (multiply matrix M by scalar:)
 **********************************************************************/    
t_matrix*   matrix_scale(const t_matrix *M, const double scalar) {
    int		i,j;
    t_matrix	*sM;

    sM=matrix_copy(M);
    if(sM!=NULL) {
	for(i=0;i<sM->nrows;i++) {
	    for(j=0;j<sM->ncols;j++){
		sM->data[i*sM->ncols+j]=scalar*sM->data[i*sM->ncols+j];
	    }
	}
    }
    return sM;
}
/**********************************************************************
 * Matrix_Id(nrows, ncols) returns pseudo-identity matrix
 **********************************************************************/
t_matrix*   matrix_Id(const int nrows, const int ncols) {
    int		i;
    t_matrix	*I;

    I=matrix_new(nrows,ncols);
    if(I!= NULL) {
	for(i=0;i<nrows && i<ncols;i++) {
	    I->data[i*ncols+i]=1;
	}
    }
    return I;
}
/**********************************************************************
 * matrix_transpose
 **********************************************************************/
t_matrix* matrix_transpose(const t_matrix *M){
    int	    i,j;
    t_matrix	*T;

    T=matrix_new(M->ncols,M->nrows);
    if(T==NULL) return T;
    for(i=0;i<M->nrows;i++) {
	for (j=0;j<M->ncols;j++) {
	    T->data[j*T->ncols+i]=M->data[i*M->ncols+j];
	}
    }
    return T;
}
/**********************************************************************
 * matrix_norm: returns the norm of the matrix/vector
 **********************************************************************/
double	matrix_norm(const t_matrix *M) {
    int	    i,j;
    double  s;
    s=0;
    for(i=0;i<M->nrows;i++){
	for(j=0;j<M->ncols;j++) {
	    s+=M->data[i*M->ncols+j]*M->data[i*M->ncols+j];
	}
    }
    return sqrt(s);
}
/**********************************************************************
 * matrix_normalize: normalize the argument matrix or vector
 **********************************************************************/
t_matrix    *matrix_normalize(t_matrix *M) {
    double	norm;
    norm=matrix_norm(M);
    if(norm == 0) return NULL;
    return matrix_scale(M,1/norm);
}
/**********************************************************************
 * matrix_trace
 **********************************************************************/
double	matrix_trace(const t_matrix* M) {
    int	    i;
    double  t;
    t=0;
    for(i=0;i<M->nrows && i<M->ncols;i++) {
	t+=M->data[i*M->ncols+i];
    }
    return t;
}
/**********************************************************************
 * matrix_lower_residue
 *     returns the sum of absolute values of elements under the diagonal
 **********************************************************************/
double	matrix_lower_residue(const t_matrix *M) {
    int	    i,j;
    double  r;
    r=0;
    for(i=0;i<M->nrows;i++) {
	for(j=0;j<i;j++) {
	    r+=fabs(M->data[i*M->ncols+j]);
	}
    }
    return r;
}
/**********************************************************************
 * matrix_get_vector:	returns a vector from a matrix column, col index
 *			is like an array: zero to M->ncols-1
 **********************************************************************/
t_matrix*    matrix_get_vector(const t_matrix *M, const int n) {
    t_matrix	*V;
    int		i;
    if( n>=M->ncols) return NULL;
    V=matrix_new_vector(M->nrows,NULL);
    for(i=0;i<M->nrows;i++)  {
	V->data[i]=M->data[i*(M->ncols)+n];
    }
    return V;
}
/**********************************************************************
 * matrix_set_vector: sets matrix column n (0.. ncols-1) from vector V
 **********************************************************************/
int	matrix_set_vector(t_matrix *M,const int n, const t_matrix *V) {
    int		i;
    if( n>=M->ncols) return 1;
    for(i=0;i<M->nrows;i++) {
	M->data[i*M->ncols+n]=V->data[i];
    }
    return 0;
}
/**********************************************************************
 * matrix_qr_decomp: decompose matrix A into Q & R such as A=QR using
 *		     householder elimination method
 *		     Q is hermitian ( tQ . Q = I )
 *		     R is a triangular superior matrix
 *		     R eigenvalues are the eigenvalues of A
 **********************************************************************/
int matrix_qr_decomp(t_matrix *A, t_matrix *Q, t_matrix *R) {

    int		n;	    // order of A
    int		i,j,k;      // matrix indices
    double	alpha;	    // signed norm of vector x1
    t_matrix	*vector_e1; // first vector of orthogonal basis
    t_matrix	*vector_x1; // first column of Householder sub-matrix
    t_matrix	*vector_u;  // u=x1 - alpha.e1, normalized
    t_matrix	*tU;	    // transpose of vector_u (intermediate)
    t_matrix	*UtU;	    // Hk=Id - 2.U.t(U)  householder submatrix
    t_matrix	*H;	    // Householder intermediate matrices
    t_matrix	*I;	    // Identity matrix (order n-k)
    t_matrix	*matrix_ptr;// temporary matrix pointer for normalisation

    n=A->ncols;
    if(n!=A->nrows) {       // error - matrix is not square)
	fprintf(stderr,"Error in QR Decompsition: matrix is not square ");
	fprintf(stderr,"(%d rows, %d cols.\n Aborting.\n",A->nrows,n);
	exit(1);
    }
    vector_e1=matrix_new_vector(n,NULL);
    if(vector_e1==NULL) {
	fprintf(stderr,"Error: not enough memory in qr_decomp. Abort.\n");
	exit(16);
    }
    vector_e1->data[0]=1;		// define unit vector                                            !
    matrix_move(R,A);			// R stats as A
    matrix_move(Q, I=matrix_Id(n,n));	// 1st iteration, Q starts as Identity
    matrix_free(I);
    if(R==NULL || Q==NULL) {
	fprintf(stderr,"Error: not enough memory in qr_decomp. Abort.\n");
	exit(16);
    }
    for(k=0;k<n-1;k++) {
 	// N.B. indices for Hk submatrix range from k to n-1
	// the order of vectors in this matrix is (n-k)
	vector_x1=matrix_new_vector(n-k,NULL);
	vector_u =matrix_new_vector(n-k,NULL);
	if(vector_x1==NULL || vector_u==NULL) {
	    fprintf(stderr,"Error: not enough memory in qr_decomp. Abort.\n");
	    exit(16);
	}
	for(i=0;i<n-k;i++) {
	    vector_x1->data[i]=R->data[(k+i)*R->ncols+k];
	}
	alpha=matrix_norm(vector_x1)*sign(vector_x1->data[0]);
	for(i=0;i<n-k;i++) {
	    vector_u->data[i]=vector_x1->data[i]- alpha*vector_e1->data[i];
	    // N.B. vector_e1 has a higher dimension, but it doesn't matter here
	}
	if(matrix_norm(vector_u)==0) {
	    fprintf(stderr,"Error: null vector during QR decomposition. The matrix is singular.\n");
	    exit(2); // dynamic memory should be freed here!
	}
	matrix_ptr=matrix_normalize(vector_u);
	matrix_free(vector_u);
	vector_u=matrix_ptr;
	// vector_U x transpose(vector_u) gives a square matrix UtU of order (n-k)
	tU=matrix_transpose(vector_u);
	UtU=matrix_prod(vector_u,tU);
	matrix_free(tU);
	// fill in the bottom right of H with calculated Hk submatrix
	H=matrix_Id(n,n);
        I=matrix_Id(n-k,n-k);
	for(i=0;i<n-k;i++) {
	    for(j=0;j<n-k;j++) {
		// H[i+k,j+k]=I[i,j] -2 * UtU[i,j]
		H->data[(i+k)*n+j+k]=I->data[i*(n-k)+j]-2*UtU->data[i*(n-k)+j];
	    }
	}
	matrix_free(UtU);
	matrix_free(I);
	matrix_free(vector_x1);
	matrix_free(vector_u);
	// return R(k) and Q(k) with R=H.R  and Q=Q.H
	matrix_ptr=matrix_prod(H,R);
	matrix_move(R,matrix_ptr);
	matrix_free(matrix_ptr);
	matrix_ptr=matrix_prod(Q,H);
	matrix_move(Q,matrix_ptr);
	matrix_free(matrix_ptr);
	matrix_free(H);
    }
    matrix_free(vector_e1);
    return 0;
}
/**********************************************************************
 * syslinQR:	Solve linear system AX=B using QR decomposition algorithm
 *		Return non-zero int if error (matrix not regular)
 *		TODO correct error mgmt (return values)
 **********************************************************************/
int syslinQR(t_matrix *A, t_matrix *X, t_matrix *B) {
    t_matrix	*Q,*R,*tQ,*Y;
    double	s;
    int		n,i,j,err_code;

    n=A->ncols;	    //order of the system
    Q=matrix_new(n,n);
    R=matrix_new(n,n);
    if(Q==NULL || R==NULL) {
	fprintf(stderr,"Error: not enough memory in SyslinQR. Abort\n");
	exit(16);
    }
    err_code=matrix_qr_decomp(A,Q,R);
    if(err_code!=0) {
	fprintf(stderr,"Error: The matrix is not regular.\n");
	fprintf(stderr,"       the linear system has no solution. ABort.\n");
	exit(1);
    }
    Y=matrix_prod(tQ=matrix_transpose(Q),B); // calculate Y=tQ.B
    matrix_free(tQ);
    X->data[n-1]=Y->data[n-1]/R->data[n*n-1];	    // Xn = Yn /Rnn   
    //calculate backwards Xn-1 .. X0 by detriangularization
    for(i=n-2;i>=0;i--){
	s=0;
	for(j=i+1;j<n;j++) {
	    s+=R->data[i*n+j]*X->data[j];	    // s=SIGMA( Rij / Xj )
	}
	X->data[i]=(Y->data[i]-s)/R->data[i*n+i];	    // Xi = [ Yi-SIGMA( Rij / Xj ) ] / Rii
    }
    matrix_free(Q);
    matrix_free(R);
    matrix_free(Y);
    return 0;
}
/**********************************************************************
 * syslinGauss:	Solve linear system AX=B using Gauss pivot method
 *		Return non-zero int if error (matrix not regular)
 *		TODO correct error mgmt / return values!
 **********************************************************************/
int syslinGauss(t_matrix const *A, t_matrix *X, t_matrix const *B) {
   
    int		n;		    // order of A
    int		i,j,k,l,km;	    // matrix indices
    double	tmp,amax,alpha;	    // maximum coeff = pivot
    t_matrix	*M;		    //working matrix (A remains untouched)
 
    M=matrix_copy(A);       // M new matrix, starts equal to A
    n=M->ncols;		    // order of the system
    matrix_move(X,B);	    // X alreday exists, starts as B in Gauss transformation
    for(k=0;k<n-1;k++){
	// look for line l0>=k with largest pivotmax
	amax=0;
	km=0;
	for(l=k;l<n;l++){
	    if((tmp=fabs(M->data[l*n+k]))>amax) {
		amax=tmp;
		km=l;
	    }
	}
	// pivot amax should not be too small
	if(amax < g_precision*g_precision) {  //pivot close to zero: irregular matrix 
	    fprintf(stderr,"Error: Pivot too small(%-6.0e).  The matrix is not regular.\n",amax);
	    fprintf(stderr,"       the linear system has no solution.\n");
	    return 2;
	}
	if (amax < g_precision) {  //pivot small: WARNING
	    fprintf(stderr, "\n**** WARNING: The pivot is very small: %-6.0e\n", amax);
	    fprintf(stderr,   "****          The matrix is close to singular.\n");
	    fprintf(stderr,   "****          The linear system may have no solution.\n");
	}
    	// switch lines k and km (in M and X)
	if( km!=k) {
	    for(j=k;j<n;j++) {                                                          
    		tmp=M->data[km*n+j];
    		M->data[km*n+j]= M->data[k*n+j];
    		M->data[k*n+j]=tmp;
    	    }
    	    tmp=X->data[km];X->data[km]=X->data[k];X->data[k]=tmp;
	}
	for(i=k+1;i<n;i++) {
	    alpha=M->data[i*n+k]/M->data[k*(n+1)]; // alpha = Mik / Mkk
	    // normalise line i to zero coeffs in column k except pivot
	    for(j=k+1;j<n;j++) {
		M->data[i*n+j]=M->data[i*n+j]-alpha*M->data[k*n+j];
	    }
	    M->data[i*n+k]=0;
	    X->data[i]=X->data[i]-alpha*X->data[k];
	}
    }
    // M is now an upper triangular matrix 
    // calculate vector solution X starting from the bottom
    X->data[n-1]=X->data[n-1]/M->data[n*n-1];
    for(k=n-2;k>=0;k--) {
    	for(j=k+1;j<n;j++){
	    X->data[k]-=M->data[k*n+j] * X->data[j];
	}
	X->data[k]=X->data[k]/M->data[k*(n+1)];
    }                                 
    matrix_free(M);
    return 0;           
}                   
/**********************************************************************
 * sign:  returns signe of x (-1, +1 or 0) as a double
 **********************************************************************/
double	    sign(double x) {
    if(x==0) return (double) 0;
    if(x<0)  return (double)-1;
	     return (double) 1;
}

/**********************************************************************
 * 			    
 * 			    MAIN
 *
 **********************************************************************/
int main(int argc, char *argv[]) {
    t_matrix	*A,*X,*B;
    double	m[25]={-1,3,4,2,2,
			2,1,1,1,0,
			-4,4,0,-5,
			5,-3,2,-3,
			-3,0,-1,1,
			-3,-3,0.01};
    double	b[5]={4,1,-2,3,2};

    int n=sqrt(sizeof(m)/sizeof(double));
    A=matrix_new(n,n);
    X=matrix_new(n,1);
    B=matrix_new(n,1);
    
    matrix_assign(A,m);
    matrix_assign(B,b);
    syslinQR(A,X,B);
    syslinGauss(A,X,B);

    matrix_free(A);
    matrix_free(B);
    matrix_free(X);
return 0;
}

