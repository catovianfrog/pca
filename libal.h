/**********************************************************************
 *	al.h 
 *	Headers for al.c, Linear algrebra library
 *	(C) Bruno S. Charrière  2014 under GPL
 **********************************************************************/
#ifndef __libal_h__
#define __libal_h__

typedef	struct  {   int	    nrows;
		    int	    ncols;
		    double  *data; }	t_matrix;
//------------------------------------------------------------------------------

t_matrix*   matrix_new(int nrows, int ncols);
void	    matrix_free(t_matrix *m);
t_matrix*   matrix_new_vector(int n, void *p);
void*	    matrix_add_rows(t_matrix *m, int n);
void	    matrix_printf(FILE *f, const t_matrix *M, char* format, char* msg);
void	    matrix_print(FILE *f, const t_matrix *m, char* msg);
t_matrix*   matrix_prod(const t_matrix *A, const t_matrix *B);
t_matrix*   matrix_add(const t_matrix *A, const t_matrix *B, const double scale);
t_matrix*   matrix_copy(const t_matrix *M);
int	    matrix_move(t_matrix *Dest, const t_matrix *Source);
t_matrix*   matrix_scale(const t_matrix *M, const double scalar);
t_matrix*   matrix_Id(const int nrows, const int ncols);
t_matrix*   matrix_transpose(const t_matrix *M);
double	    matrix_norm(const t_matrix *M);
t_matrix*   matrix_normalize(const t_matrix *M);
double	    matrix_trace(const t_matrix* M);
double	    matrix_lower_residue(const t_matrix *M);
t_matrix*   matrix_get_vector(const t_matrix *M, const int n);
int	    matrix_set_vector(t_matrix *M,const int n, const t_matrix *V);
t_matrix*   matrix_get_block(const t_matrix *M,int n,int m,int i0,int j0);
int         matrix_set_block(const t_matrix *S,t_matrix *D,int i0,int j0);
int	    matrix_qr_decomp(t_matrix *A, t_matrix *Q, t_matrix *R);
int	    matrix_eigenvalues(const t_matrix *A, t_matrix *EV, double precision);
t_matrix*   matrix_ev_inertia(const t_matrix *eigenvalues);
t_matrix*   matrix_eigenvector(const t_matrix *M, const double r);
t_matrix*   matrix_eigenvectors(const t_matrix *M, const t_matrix *e_values);
int	    syslinQR(t_matrix *A, t_matrix *X, t_matrix *B);
int	    syslinGauss(t_matrix const *A, t_matrix *X, t_matrix const *B);
//--------- basic statistics functions  ----------------------------------------
double	    sign(double x);
double	    mean(t_matrix *V);
double	    sumsquares(t_matrix *V);
double	    variance(t_matrix *V);
double	    sdev(t_matrix *V);
double	    variance_sample(t_matrix *V);
double	    sdev_sample(t_matrix *V);

/**********************************************************************
 *	    Inline functions  
 **********************************************************************
inline 
double	    matrix_get(const t_matrix *M, int i, int j) {
	    if(i>=M->nrows || j>=M->ncols) return nan(""); return M->data[i*M->ncols+j];}
inline 
int	    matrix_set(t_matrix *M,const int i, const int j, const double r){
	    if(i>=M->nrows || j>=M->ncols) return 1; M->data[i*M->ncols+j]=r; return 0;}
inline 
void	    matrix_assign(t_matrix *m, void *p){memcpy(m->data,p,m->nrows*m->ncols*sizeof(double));}
*********************************************************************/
#endif

