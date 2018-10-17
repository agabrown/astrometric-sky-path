/*----------------------------------------------------------------------------
 *
 * vecmat_utils.h
 *
 * Header file for vecmat_utils.c
 *
 * Author:        Anthony Brown
 * Last updated:  20.03.2006
 * First created: 11.05.2005
 *
 *----------------------------------------------------------------------------
 */

char *cvector(long nelements);
float *vector(long nelements);
int *ivector(long nelements);
long *lvector(long nelements);
double *dvector(long nelements);
char **cmatrix_cm(long ncolumns, long nrows);
float **matrix_cm(long ncolumns, long nrows);
double **dmatrix_cm(long ncolumns, long nrows);
double **dmatrix(long nrows, long ncolumns);
float ***matstack_cm(long ndepth, long ncolumns, long nrows);
double ***dmatstack_cm(long ndepth, long ncolumns, long nrows);

char *cvector_resize(char *vector, long nelements);
double *dvector_resize(double *vector, long nelements);

void free_cvector(char *v);
void free_vector(float *v);
void free_ivector(int *v);
void free_lvector(long *v);
void free_dvector(double *v);
void free_cmatrix_cm(char **m);
void free_matrix_cm(float **m);
void free_dmatrix_cm(double **m);
void free_dmatrix(double **m);
void free_matstack_cm(float ***c);
void free_dmatstack_cm(double ***c);

int fcmp(const void *p1, const void *p2);
int dcmp(const void *p1, const void *p2);

void sort(float *array, long n);
void dsort(double *array, long n);
void swap(float *v, long i, long j);
void dswap(double *v, long i, long j);
void quicksort(float *array, long n);
void dquicksort(double *array, long n);
float vtotal(float *v, long n);
double dvtotal(double *v, long n);
