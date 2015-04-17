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
#include <math.h> // pow

/********************************************************************
 * TKF92HG stuff
 * *****************************************************************/

struct TKF92HGLikelihoodFunction1D_params
{
      double len, mu, r, ps, kf;
      gsl_matrix *substModel;
      gsl_vector *eqFrequencies;
      int *seq1Int, *seq2Int;
      int SA, SB;
};

struct TKF92HGLikelihoodFunction5D_params
{
  double len;
  gsl_matrix *substModel;
  gsl_vector *eqFrequencies;
  int *seq1Int, *seq2Int;
  int SA, SB;
};

/**** The actual function to do the computation ****/
double TKF92HGLikelihoodFunction(int *seq1Int, int *seq2Int, double len,
    double mu, double r, double ps, double kf,
    double distance, gsl_matrix *substModelS, gsl_matrix *substModelF,
    gsl_vector *eqFrequencies, int SA, int SB);

/**** TKF92HG 1D and 5D main function for gsl ****/
double TKF92HGLikelihoodFunction1D(double distance, void *params);
double TKF92HGLikelihoodFunction5D(const gsl_vector *v,  void *params);

/**** purely calculate the likelihodd given the distance, mu and r ****/
SEXP TKF92HGLikelihoodFunctionWrapper(SEXP seq1IntR, SEXP seq2IntR, SEXP distanceR, SEXP muR, SEXP rR, SEXP psR, SEXP kfR, SEXP expectedLength, SEXP probMatR, SEXP eqFrequenciesR);


/**** TKF92HG 1D with Brent implementation ****/
SEXP TKF92HGLikelihoodFunction1DMain(SEXP seq1IntR, SEXP seq2IntR, 
    SEXP muR, SEXP rR, SEXP psR, SEXP kfR,
    SEXP expectedLength, SEXP probMatR, SEXP eqFrequenciesR);

/**** TKF92HG with NM simplex implementation ****/
SEXP TKF92HGLikelihoodFunction5DMainNM(SEXP seq1IntR, SEXP seq2IntR,
    SEXP expectedLength, SEXP probMatR, SEXP eqFrequenciesR);


