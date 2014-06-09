pca
===

Principal component analysis software, and linear algebra library

	    V E R S I O N   H I S T O R Y
	    =============================
v0.9 

    Start from initial linear algebra stub (qr.c)
v1.0

    Fully working PCA version.
    All the linear algebra (matrix) functions are integrated
    All the data are stored statically, in the stack.
    So each matriux or vector takes up MAX_COLxMAX_COL size in stack memory
    Data sets are limited to 150 lines
v1.2

    All data sets (in main and in functions) are allocated dynamically
    and stored in the heap. Also using syslinGauss instead of syslinQR uses much 
    less memory (clearly visible with valgrind). However, size of matrices and data 
    set is still fixed (inevitable because matrices are two-dimensional arrays).
    This versions runs v1.2 with 350 rows x 350 cols maximum.
    Above this, valgrind detects invalid memory reads and writes. 
    Also corrected a bug in str2headers (added a space at the end of linestr) 

 
