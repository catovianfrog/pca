pca
===

Principal component analysis software, and linear algebra library

#### CHANGE LOG

v0.9 

Start from initial linear algebra stub (qr.c)


v1.0

Fully working PCA version.
    All the linear algebra (matrix) functions are integrated
    All the data are stored statically, in the stack.
    the mlatrices are stored as two dimensional arrays, therefore
    requiring fixed size allocation.
    So each matrix or vector takes up MAX_COLxMAX_COL size in stack memory
    Data sets are limited to 150 linesa
    running the programm typically requires 2 Gb of memory!

v1.2

All data sets (in main and in functions) are allocated dynamically
    and stored in the heap. Also using syslinGauss instead of syslinQR uses much 
    less memory (clearly visible with valgrind). However, size of matrices and data 
    set is still fixed (inevitable because matrices are two-dimensional arrays).
    This versions runs v1.2 with 350 rows x 350 cols maximum.
    Above this, valgrind detects invalid memory reads and writes. 
    Also corrected a bug in str2headers (added a space at the end of linestr)

v2.0

__Complete overhaul__
    Programme split between a linear algebra library (libal.c) and a main pca.c
    principal componet analysis programme.
    *libal.c* groups  basic statistics functions (mean, standard deviation, 
    sumsquares,..), matrix structure definition (including dimensions n,m and 
    pointer to one-dimensional array, to be allocated dynamically), and all 
    the matrix operations including transpose, scaling, porduct, inversion,
    linear system resolution, QR transformation, eigenvalues, eigevctors, etc..

**pca.c** is the main PCA resolution programme. It reads data from  a 
    default input programme "pca_input.dat", calculates the pca and displays it
    screen and text file.
    
v2.1 -- 
    Uses **libstring** library function 'tokenize()' instead of internal 'str2headers()'
    to read headers. v2.1.2 is the first fully working version (without memory
    leaks) that uses _libstring.c_.

