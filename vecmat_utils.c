/*----------------------------------------------------------------------------
 *
 * vecmat_utils.c
 *
 * Author:        Anthony Brown
 * Last updated:  20.03.2006
 * First created: 10.05.2005
 *
 * Define some utility functions for allocating memory for vectors and
 * matrices. The matrices can be allocated column-major by calling the
 * appropriate routines. The code is based on the Numerical Recipes in C
 * nrutil.c code.
 *
 *----------------------------------------------------------------------------
 */

#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>

/*==========================================================================*/

/* allocate a char vector with nelements elements v[0..nelements-1] */
char *cvector(long nelements)
{
  char *v;

  v = (char*)malloc((size_t) (nelements*sizeof(char)));
  if (!v) {
    printf("Allocation failure in cvector()!\n");
    exit(1);
  }

  return v;
}

/*==========================================================================*/

/* allocate a float vector with nelements elements v[0..nelements-1] */
float *vector(long nelements)
{
  float *v;

  v = (float*)malloc((size_t) (nelements*sizeof(float)));
  if (!v) {
    printf("Allocation failure in vector()!\n");
    exit(1);
  }

  return v;
}

/*==========================================================================*/

/* allocate an integer vector with nelements elements v[0..nelements-1] */
int *ivector(long nelements)
{
  int *v;

  v = (int*)malloc((size_t) (nelements*sizeof(int)));
  if (!v) {
    printf("Allocation failure in ivector()!\n");
    exit(1);
  }

  return v;
}

/*==========================================================================*/

/* allocate a long vector with nelements elements v[0..nelements-1] */
long *lvector(long nelements)
{
  long *v;

  v = (long*)malloc((size_t) (nelements*sizeof(long)));
  if (!v) {
    printf("Allocation failure in lvector()!\n");
    exit(1);
  }

  return v;
}

/*==========================================================================*/

/* allocate a double vector with nelements elements v[0..nelements-1] */
double *dvector(long nelements)
{
  double *v;

  v = (double*)malloc((size_t) (nelements*sizeof(double)));
  if (!v) {
    printf("Allocation failure in dvector()!\n");
    exit(1);
  }

  return v;
}

/*==========================================================================*/

/* allocate a char matrix in COLUMN-MAJOR format with ncolumns*nrows
   elements c[0..ncolumns-1][0..nrows-1] */
char **cmatrix_cm(long ncolumns, long nrows)
{
  long i;
  char **m;

  /* allocate pointers to columns */
  m = (char**) malloc((size_t) (ncolumns*sizeof(char*)));
  if (!m) {
    printf("Allocation failure in cmatrix_cm()!\n");
    exit(1);
  }

  /* allocate columns and set pointers to them */
  m[0] = (char*) malloc((size_t) (ncolumns*nrows*sizeof(char)));
  if (!m[0])  {
    printf("Allocation failure in cmatrix_cm()!\n");
    exit(1);
  }

  for (i=1;i<ncolumns;i++) m[i]=m[i-1]+nrows;

  /* return pointer to array of pointers to columns */
  return m;
}

/*==========================================================================*/

/* allocate a float matrix in COLUMN-MAJOR format with ncolumns*nrows
   elements m[0..ncolumns-1][0..nrows-1] */
float **matrix_cm(long ncolumns, long nrows)
{
  long i;
  float **m;

  /* allocate pointers to columns */
  m = (float**) malloc((size_t) (ncolumns*sizeof(float*)));
  if (!m) {
    printf("Allocation failure in matrix_cm()!\n");
    exit(1);
  }

  /* allocate columns and set pointers to them */
  m[0] = (float*) malloc((size_t) (ncolumns*nrows*sizeof(float)));
  if (!m[0])  {
    printf("Allocation failure in matrix_cm()!\n");
    exit(1);
  }

  for (i=1;i<ncolumns;i++) m[i]=m[i-1]+nrows;

  /* return pointer to array of pointers to columns */
  return m;
}

/*==========================================================================*/

/* allocate a double matrix in COLUMN-MAJOR format with ncolumns*nrows
   elements m[0..ncolumns-1][0..nrows-1] */
double **dmatrix_cm(long ncolumns, long nrows)
{
  long i;
  double **m;

  /* allocate pointers to columns */
  m = (double**) malloc((size_t) (ncolumns*sizeof(double*)));
  if (!m) {
    printf("Allocation failure in dmatrix_cm()!\n");
    exit(1);
  }

  /* allocate columns and set pointers to them */
  m[0] = (double*) malloc((size_t) (ncolumns*nrows*sizeof(double)));
  if (!m[0])  {
    printf("Allocation failure in dmatrix_cm()!\n");
    exit(1);
  }

  for (i=1;i<ncolumns;i++) m[i]=m[i-1]+nrows;

  /* return pointer to array of pointers to columns */
  return m;
}

/* allocate a double matrix in ROW-MAJOR format with nrows*ncolumns
   elements m[0..nrows-1][0..ncolumns-1] */
double **dmatrix(long nrows, long ncolumns)
{
  long i;
  double **m;

  /* allocate pointers to rows */
  m = (double**) malloc((size_t) (nrows*sizeof(double*)));
  if (!m) {
    printf("Allocation failure in dmatrix()!\n");
    exit(1);
  }

  /* allocate rows and set pointers to them */
  m[0] = (double*) malloc((size_t) (nrows*ncolumns*sizeof(double)));
  if (!m[0])  {
    printf("Allocation failure in dmatrix()!\n");
    exit(1);
  }

  for (i=1;i<nrows;i++) m[i]=m[i-1]+ncolumns;

  /* return pointer to array of pointers to rows */
  return m;
}

/*==========================================================================*/

/* Allocate a float 3D-matrix as a stack of 2D column-major format matrices. The
   first index refers to the 3rd dimension or 'z-axis' of the cube. The cube
   has ndepth*ncolumns*nrows elements
   c[0..ndepth-1][0..ncolumns-1][0..nrows-1]. */

float ***matstack_cm(long ndepth, long ncolumns, long nrows)
{
  long i,j;
  float ***c;

  /* allocate pointers to pointers to columns */
  c = (float***) malloc((size_t) (ndepth*sizeof(float**)));
  if (!c) {
    printf("Allocation failure in matstack_cm()!\n");
    exit(1);
  }

  /* allocate pointers to columns and set pointers to them */
  c[0] = (float**) malloc((size_t) (ndepth*ncolumns*sizeof(float*)));
  if (!c[0]) {
    printf("Allocation failure in matstack_cm()!\n");
    exit(1);
  }

  /* allocate columns and set pointers to them */
  c[0][0] = (float *) malloc((size_t) (ndepth*ncolumns*nrows*sizeof(float)));
  if (!c[0][0]) {
    printf("Allocation failure in matstack_cm()!\n");
    exit(1);
  }

  for (j=1;j<ncolumns;j++) c[0][j]=c[0][j-1]+nrows;
  for (i=1;i<ndepth;i++) {
    c[i]=c[i-1]+ncolumns;
    c[i][0]=c[i-1][0]+ncolumns*nrows;
    for (j=1;j<ncolumns;j++) {
      c[i][j]=c[i][j-1]+nrows;
    }
  }
  
  /* return pointer to array of pointers to columns */
  return c;
}

/*==========================================================================*/

/* Allocate a double 3D-matrix as a stack of 2D column-major format matrices. The
   first index refers to the 3rd dimension or 'z-axis' of the cube. The cube
   has ndepth*ncolumns*nrows elements
   c[0..ndepth-1][0..ncolumns-1][0..nrows-1]. */

double ***dmatstack_cm(long ndepth, long ncolumns, long nrows)
{
  long i,j;
  double ***c;

  /* allocate pointers to pointers to columns */
  c = (double***) malloc((size_t) (ndepth*sizeof(double**)));
  if (!c) {
    printf("Allocation failure in dmatstack_cm()!\n");
    exit(1);
  }

  /* allocate pointers to columns and set pointers to them */
  c[0] = (double**) malloc((size_t) (ndepth*ncolumns*sizeof(double*)));
  if (!c[0]) {
    printf("Allocation failure in dmatstack_cm()!\n");
    exit(1);
  }

  /* allocate columns and set pointers to them */
  c[0][0] = (double *) malloc((size_t) (ndepth*ncolumns*nrows*sizeof(double)));
  if (!c[0][0]) {
    printf("Allocation failure in dmatstack_cm()!\n");
    exit(1);
  }

  for (j=1;j<ncolumns;j++) c[0][j]=c[0][j-1]+nrows;
  for (i=1;i<ndepth;i++) {
    c[i]=c[i-1]+ncolumns;
    c[i][0]=c[i-1][0]+ncolumns*nrows;
    for (j=1;j<ncolumns;j++) {
      c[i][j]=c[i][j-1]+nrows;
    }
  }
  
  /* return pointer to array of pointers to columns */
  return c;
}

/*==========================================================================*/

/* resize a double vector to nelements; v[0..nelements-1] */
double *dvector_resize(double *vector, long nelements)
{
  double *v;

  v = (double*)realloc(vector, (size_t) (nelements*sizeof(double)));
  if (!v) {
    printf("Re-allocation failure in dvector_resize()!\n");
    exit(1);
  }

  return v;
}

/*==========================================================================*/

/* resize a char vector to nelements; v[0..nelements-1] */
char *cvector_resize(char *vector, long nelements)
{
  char *v;

  v = (char*)realloc(vector, (size_t) (nelements*sizeof(char)));
  if (!v) {
    printf("Re-allocation failure in cvector_resize()!\n");
    exit(1);
  }

  return v;
}

/*==========================================================================*/

/* free a char vector allocated with cvector() */
void free_cvector(char *v)
{
  free((char*) v);
}

/*==========================================================================*/

/* free a float vector allocated with vector() */
void free_vector(float *v)
{
  free((char*) v);
}

/*==========================================================================*/

/* free an integer vector allocated with ivector() */
void free_ivector(int *v)
{
  free((char*) v);
}

/*==========================================================================*/

/* free a long vector allocated with lvector() */
void free_lvector(long *v)
{
  free((char*) v);
}

/*==========================================================================*/

/* free a double vector allocated with dvector() */
void free_dvector(double *v)
{
  free((char*) v);
}

/*==========================================================================*/

/* free a char matrix allocated with cmatrix_cm() */
void free_cmatrix_cm(char **m)
{
  free((char*) m[0]);
  free((char*) m);
}

/*==========================================================================*/

/* free a float matrix allocated with matrix_cm() */
void free_matrix_cm(float **m)
{
  free((char*) m[0]);
  free((char*) m);
}

/*==========================================================================*/

/* free a double matrix allocated with dmatrix_cm() */
void free_dmatrix_cm(double **m)
{
  free((char*) m[0]);
  free((char*) m);
}

/* free a double matrix allocated with dmatrix() */
void free_dmatrix(double **m)
{
  free((char*) m[0]);
  free((char*) m);
}

/*==========================================================================*/

/* free a float 3D matrix allocated with matstack_cm() */
void free_matstack_cm(float ***c)
{
  free((char*) c[0][0]);
  free((char*) c[0]);
  free((char*) c);
}

/*==========================================================================*/

/* free a double 3D matrix allocated with dmatstack_cm() */
void free_dmatstack_cm(double ***c)
{
  free((char*) c[0][0]);
  free((char*) c[0]);
  free((char*) c);
}

/*==========================================================================*/

/* Float compare of *p1 and *p2. For use with the C library function qsort. */
int fcmp(const void *p1, const void *p2)
{
  float v1,v2;

  v1 = *(float*)p1;
  v2 = *(float*)p2;
  if (v1 < v2)
    return -1;
  else if (v1 == v2)
    return 0;
  else
    return 1;
}

/*==========================================================================*/

/* Double compare of *p1 and *p2. For use with the C library function qsort. */
int dcmp(const void *p1, const void *p2)
{
  double v1,v2;

  v1 = *(double*)p1;
  v2 = *(double*)p2;
  if (v1 < v2)
    return -1;
  else if (v1 == v2)
    return 0;
  else
    return 1;
}

/*==========================================================================*/

/* Sort an array of floats using the C-library qsort function. */
void sort(float *array, long n)
{
  qsort(array, n, sizeof(array[0]), fcmp);
}

/*==========================================================================*/

/* Sort an array of doubles using the C-library qsort function. */
void dsort(double *array, long n)
{
  qsort(array, n, sizeof(array[0]), dcmp);
}

/*==========================================================================*/

/* Interchange floats v[i] and v[j] */
void swap(float *v, long i, long j)
{
  float temp;

  temp=v[i];
  v[i]=v[j];
  v[j]=temp;
}

/*==========================================================================*/

/* Interchange doubles v[i] and v[j] */
void dswap(double *v, long i, long j)
{
  double temp;

  temp=v[i];
  v[i]=v[j];
  v[j]=temp;
}

/*==========================================================================*/

/* Sort an array of floats using a hand-coded quicksort function (from The
   Practice of Programming book by Kernighan & Pike). */
void quicksort(float *array, long n)
{
  long i, last;

  if (n<=1) return;

  swap(array, 0, rand() % n);
  last=0;
  for (i=1;i<n;i++) {
    if (array[i] < array[0]) swap(array, ++last, i);
  }
  swap(array, 0, last);
  quicksort(array,last);
  quicksort(array+last+1, n-last-1);
}

/*==========================================================================*/

/* Sort an array of doubles using a hand-coded quicksort function (from The
   Practice of Programming book by Kernighan & Pike). */
void dquicksort(double *array, long n)
{
  long i, last;

  if (n<=1) return;

  dswap(array, 0, rand() % n);
  last=0;
  for (i=1;i<n;i++) {
    if (array[i] < array[0]) dswap(array, ++last, i);
  }
  dswap(array, 0, last);
  dquicksort(array,last);
  dquicksort(array+last+1, n-last-1);
}

/*==========================================================================*/

/* Calculate the sum of all elements of the float vector v[0..n-1]. Sort the
   elements first and then carry out the sum in order to ensure numerical
   accuracy. */
float vtotal(float *v, long n)
{
  long i;
  float *vcopy;
  float total=0.0;

  vcopy=vector(n);
  for (i=0;i<n;i++) {
    vcopy[i]=v[i];
  }

  quicksort(vcopy, n);

  for (i=0;i<n;i++) {
    total+=vcopy[i];
  }

  return total;
}

/*==========================================================================*/

/* Calculate the sum of all elements of the double vector v[0..n-1]. Sort the
   elements first and then carry out the sum in order to ensure numerical
   accuracy. */
double dvtotal(double *v, long n)
{
  long i;
  double *vcopy;
  double total=0.0;

  vcopy=dvector(n);
  for (i=0;i<n;i++) {
    vcopy[i]=v[i];
  }

  dquicksort(vcopy, n);

  for (i=0;i<n;i++) {
    total+=vcopy[i];
  }

  return total;
}
