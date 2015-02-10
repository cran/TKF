/*************************************************************************
    > File Name: matrix.c
    > Author: Ge Tan
    > Mail: gtan@me.com 
    > Created Time: Wed Jul 30 22:55:24 2014
 ************************************************************************/

#include<stdio.h>
#include "matrix.h"
#include "assert.h"

/*********************************************************************
 * Print a GSL matrix in a proper shape
 ********************************************************************/
void printGSLMatrix(const gsl_matrix *m){
  //int status, n = 0;
  for(size_t i = 0; i < m->size1; i++){
    for(size_t j = 0; j < m->size2; j++){
      //if((status = printf("%g ", gsl_matrix_get(m, i, j))) < 0)
      //  return -1;
      // n += status;
      Rprintf("%g ", gsl_matrix_get(m, i, j));
    }
    //if((status = printf("\n")) < 0)
    Rprintf("\n");
    //  return -1;
    //n += status;
  }
  //return n;
}

void printGSLMatrixComplex(const gsl_matrix_complex *m){
  //int status, n = 0;
  for(size_t i = 0; i < m->size1; i++){
    for(size_t j = 0; j < m->size2; j++){
      //if((status = printf("%g ", GSL_REAL(gsl_matrix_complex_get(m, i, j)))) < 0)
      //    return -1;
      //n += status;
      Rprintf("%g ", GSL_REAL(gsl_matrix_complex_get(m, i, j)));
    }
    //if((status = printf("\n")) < 0)
    //  return -1;
    //n += status;
    Rprintf("\n");
  }
  //return n;
}

/********************************************************************
 * Make diagonal matrix from vector for Real and Complex data
 * *****************************************************************/
void asDiagonalComplex(const gsl_vector_complex *X, gsl_matrix_complex *mat){
  int i, j;
  int nrow = mat->size1;
  for(i = 0; i < nrow; i++){
    for(j = 0; j < nrow; j++){
      if(i == j){
        gsl_matrix_complex_set(mat, i, j,
            gsl_vector_complex_get(X, i));
      }
    }
  }
}

/********************************************************************
 * Create identity complex matrix
 * *****************************************************************/
void create_identity_matrix(gsl_matrix_complex *I){
  int i,j;
  int nrow = I->size1;
  for(i = 0; i < nrow; i++){
    for(j = 0;j < nrow; j++){
      if(i == j){
        gsl_matrix_complex_set(I, i, j, gsl_complex_rect(1,0));
      }else{
        gsl_matrix_complex_set(I, i, j, gsl_complex_rect(0,0));
      }
    }
  }
}

/*********************************************************************
 * Matrix add/sub, taken from http://mimo-coherent-with-c.googlecode.com/svn-history/r60/trunk/moperations.h
 * Not verified yet.
 * ******************************************************************/
//matrix + matrix
void matrix_add(gsl_matrix *c, const gsl_matrix *a, const gsl_matrix *b){
  gsl_matrix_memcpy(c, a);
  gsl_matrix_add(c, b);
}

//matrix - matrix
void matrix_sub(gsl_matrix *c, const gsl_matrix *a, const gsl_matrix *b){
  gsl_matrix_memcpy(c, a);
  gsl_matrix_sub(c, b);
}

//matrix + constant
void matrix_add_constant(gsl_matrix* c, gsl_matrix *a, const double x) {
  gsl_matrix_memcpy(c, a);
  gsl_matrix_add_constant(c, x);
}

//complex matrix + complex matrix
void matrix_complex_add(gsl_matrix_complex *c, const gsl_matrix_complex *a, 
  const gsl_matrix_complex *b){
  gsl_matrix_complex_memcpy(c, a);
  gsl_matrix_complex_add(c, b);
}

//complex matrix - complex matrix
void matrix_complex_sub(gsl_matrix_complex *c, const gsl_matrix_complex *a, 
  const gsl_matrix_complex *b){
  gsl_matrix_complex_memcpy(c, a);
  gsl_matrix_complex_sub(c, b);
}

//complex matrix + constant
void matrix_complex_add_constant(gsl_matrix_complex *c, gsl_matrix_complex *a, 
    gsl_complex x){
  gsl_matrix_complex_memcpy(c, a);
  gsl_matrix_complex_add_constant(c, x);
}


/********************************************************************
 * Matrix inverse for Real and Complex
 * *****************************************************************/
void gsl_matrix_inverse(const gsl_matrix *m, gsl_matrix *minvert){
  int s;
  gsl_matrix *b = gsl_matrix_alloc(m->size1, m->size2);
  gsl_permutation *p = gsl_permutation_alloc(m->size1);
  gsl_matrix_memcpy(b, m);
  gsl_linalg_LU_decomp(b, p, &s);
  gsl_linalg_LU_invert(b, p, minvert);
  gsl_matrix_free(b);
  gsl_permutation_free(p);
}

void gsl_matrix_complex_inverse(const gsl_matrix_complex *m, 
    gsl_matrix_complex *minvert){
  // validated.
  int s;
  gsl_matrix_complex *b = gsl_matrix_complex_alloc(m->size1, m->size2);
  gsl_permutation *p = gsl_permutation_alloc(m->size1);
  gsl_matrix_complex_memcpy(b, m);
  gsl_linalg_complex_LU_decomp(b, p, &s);
  gsl_linalg_complex_LU_invert(b, p, minvert);
  gsl_matrix_complex_free(b);
  gsl_permutation_free(p);
}

/*********************************************************************
 * Matrix conjug for Real and Complex
 * Not verified yet.
 * ******************************************************************/
void gsl_matrix_complex_conjug(gsl_matrix_complex *c, 
    gsl_matrix_complex *a){
  gsl_matrix_complex_memcpy(c, a);
  int i,j;
  for(i = 0; i< c->size1; i++){
    for(j = 0; j < c->size2; j++){
      gsl_complex *p = gsl_matrix_complex_ptr(c, i, j);
      *p = gsl_complex_conjugate(*p);
    }
  }
}

/********************************************************************
 * Generate mutation probability matrix for PAM distance from PAM1 matrix
 * *****************************************************************/

void PAMn(const gsl_matrix *m, const double distance, gsl_matrix *mPAM){
  int nrow;
  int i, j;
  nrow = m->size1;
  gsl_matrix *m2 = gsl_matrix_alloc(nrow, nrow);
  gsl_matrix_memcpy(m2, m);
  gsl_vector_complex *eval = gsl_vector_complex_alloc(nrow);
  gsl_matrix_complex *evec = gsl_matrix_complex_alloc(nrow, nrow);
  gsl_eigen_nonsymmv_workspace * w = gsl_eigen_nonsymmv_alloc(nrow);
  gsl_eigen_nonsymmv(m2, eval, evec, w);
  
  // free the eigen nonsymmv workspace
  gsl_eigen_nonsymmv_free(w);

  gsl_eigen_nonsymmv_sort(eval, evec,
      GSL_EIGEN_SORT_ABS_DESC);

  // build the diagonal matrix with pow(eval, distance)
  gsl_matrix_complex *diagMatrix = gsl_matrix_complex_calloc(nrow, nrow);
  for(i = 0; i < nrow; i++){
    for(j = 0; j < nrow; j++){
      if(i == j){
        gsl_matrix_complex_set(diagMatrix, i, j, 
            gsl_complex_pow_real(gsl_vector_complex_get(eval, i), distance));
        //gsl_complex_log(gsl_vector_complex_get(eval, i)));
      }
    }
  }
  // invert the eigen vector matrix
  gsl_matrix_complex *evecinv = gsl_matrix_complex_calloc(nrow, nrow);
  gsl_matrix_complex_inverse(evec, evecinv);

  // matrix multiplication
  gsl_matrix_complex *mPAMComplex = gsl_matrix_complex_alloc(nrow, nrow);
  gsl_matrix_complex *mPAMComplex2 = gsl_matrix_complex_alloc(nrow, nrow);
  gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, GSL_COMPLEX_ONE, 
      evec, diagMatrix, GSL_COMPLEX_ZERO, mPAMComplex);
  gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, GSL_COMPLEX_ONE,
      mPAMComplex, evecinv, GSL_COMPLEX_ZERO, mPAMComplex2);

  for(i = 0; i < nrow; i++){
    for(j = 0; j < nrow; j++){
      // The matrix should be almost real!
      assert(GSL_IMAG(gsl_matrix_complex_get(mPAMComplex2, i, j)) < Eps);
      gsl_matrix_set(mPAM, i, j, 
          GSL_REAL(gsl_matrix_complex_get(mPAMComplex2, i, j))   );
    }
  }

  // free the allocated eigen values and vectors
  gsl_vector_complex_free(eval);
  gsl_matrix_complex_free(evec);
  gsl_matrix_complex_free(diagMatrix);
  gsl_matrix_complex_free(mPAMComplex);
  gsl_matrix_complex_free(mPAMComplex2);
  gsl_matrix_complex_free(evecinv);
  gsl_matrix_free(m2);
}

