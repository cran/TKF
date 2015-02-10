/*************************************************************************
    > File Name: TKF91.h
    > Author: Ge Tan
    > Mail: gtan@me.com 
    > Created Time: Mon Aug  4 19:05:09 2014
 ************************************************************************/

#include<stdio.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_math.h>
#include<gsl/gsl_min.h>
#include<gsl/gsl_multimin.h>
#include "Rdefines.h"
#include "Rinternals.h"
#include <R_ext/Rdynload.h>
#include "matrix.h"
#include "MathFunctions.h"

/********************************************************************
 * TKF91 stuff
 * *****************************************************************/

struct TKF91LikelihoodFunction1D_params
{
      double len, mu;
      gsl_matrix *substModel;
      gsl_vector *eqFrequencies;
      int *seq1Int, *seq2Int;
      int SA, SB;
};

struct TKF91LikelihoodFunction2D_params
{
  double len;
  gsl_matrix *substModel;
  gsl_vector *eqFrequencies;
  int *seq1Int, *seq2Int;
  int SA, SB;
};

/**** The actual function to do the computation ****/
double TKF91LikelihoodFunction(int *seq1Int, int *seq2Int, double len,
    double mu, double distance, gsl_matrix *substModel,
    gsl_vector *eqFrequencies, int SA, int SB);

/**** TKF91 1D and 2D main function for gsl ****/
double TKF91LikelihoodFunction1D(double distance, void *params);
double TKF91LikelihoodFunction2D(const gsl_vector *v,  void *params);

/**** purely calculate the likelihodd given the distance and mu ****/
SEXP TKF91LikelihoodFunctionWrapper(SEXP seq1IntR, SEXP seq2IntR, SEXP distanceR, SEXP muR, SEXP expectedLength, SEXP probMatR, SEXP eqFrequenciesR);

/**** df function for BFGS2 implementation ****/
void TKF91LikelihoodFunction2D_df(const gsl_vector *v, void *params,
    gsl_vector *df);
void TKF91LikelihoodFunction2D_fdf(const gsl_vector *v, void *params,
    double *f, gsl_vector *df);

/**** TKF91 1D with Brent implementation ****/
SEXP TKF91LikelihoodFunction1DMain(SEXP seq1IntR, SEXP seq2IntR, SEXP muR,
    SEXP expectedLength, SEXP probMatR, SEXP eqFrequenciesR);

/**** TKF91 with BFGS2 implementation ****/
SEXP TKF91LikelihoodFunction2DMain(SEXP seq1IntR, SEXP seq2IntR,
    SEXP expectedLength, SEXP probMatR, SEXP eqFrequenciesR);

/**** TKF91 with NM simplex implementation ****/
SEXP TKF91LikelihoodFunction2DMainNM(SEXP seq1IntR, SEXP seq2IntR,
    SEXP expectedLength, SEXP probMatR, SEXP eqFrequenciesR);



/********************************************************************
 * TKF92 stuff
 * *****************************************************************/


