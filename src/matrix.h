/*************************************************************************
    > File Name: .matrix.h
    > Author: Ge Tan
    > Mail: gtan@me.com 
    > Created Time: Wed Jul 30 22:57:00 2014
 ************************************************************************/

#include<stdio.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#define Eps 1e-12
#include "Rdefines.h"
#include "Rinternals.h"

/*********************************************************************
 * Print a GSL matrix in a proper shape
 ********************************************************************/
void printGSLMatrix(const gsl_matrix *m);
void printGSLMatrixComplex(const gsl_matrix_complex *m);

/********************************************************************
 * Make diagonal matrix from vector for Real and Complex data
 * *****************************************************************/
void asDiagonalComplex(const gsl_vector_complex *X, gsl_matrix_complex *mat);

/********************************************************************
 * Create identity complex matrix
 * *****************************************************************/
void create_identity_matrix(gsl_matrix_complex *I);

/*********************************************************************
 * Matrix add/sub, taken from http://mimo-coherent-with-c.googlecode.com/svn-history/r60/trunk/moperations.h
 * Not verified yet.
 * ******************************************************************/
void matrix_add(gsl_matrix *c, const gsl_matrix *a, const gsl_matrix *b);
void matrix_sub(gsl_matrix *c, const gsl_matrix *a, const gsl_matrix *b);
void matrix_add_constant(gsl_matrix* c, gsl_matrix *a, const double x);
void matrix_complex_add(gsl_matrix_complex *c, const gsl_matrix_complex *a,
  const gsl_matrix_complex *b);
void matrix_complex_sub(gsl_matrix_complex *c, const gsl_matrix_complex *a,
  const gsl_matrix_complex *b);
void matrix_complex_add_constant(gsl_matrix_complex *c, gsl_matrix_complex *a,
    gsl_complex x);

/********************************************************************
 * Matrix inverse for Real and Complex
 * *****************************************************************/
void gsl_matrix_inverse(const gsl_matrix *m, gsl_matrix *minvert);
void gsl_matrix_complex_inverse(const gsl_matrix_complex *m,
    gsl_matrix_complex *minvert);

/*********************************************************************
 * Matrix conjug for Real and Complex
 * Not verified yet.
 * ******************************************************************/
void gsl_matrix_complex_conjug(gsl_matrix_complex *c,
    gsl_matrix_complex *a);

/********************************************************************
 * Generate mutation probability matrix for PAM distance from PAM1 matrix
 * *****************************************************************/
void PAMn(const gsl_matrix *m, const double distance, gsl_matrix *mPAM);

