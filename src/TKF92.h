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
 * TKF92 stuff
 * *****************************************************************/

struct TKF92LikelihoodFunction1D_params
{
      double len, mu, r;
      gsl_matrix *substModel;
      gsl_vector *eqFrequencies;
      int *seq1Int, *seq2Int;
      int SA, SB;
};

struct TKF92LikelihoodFunction3D_params
{
  double len;
  gsl_matrix *substModel;
  gsl_vector *eqFrequencies;
  int *seq1Int, *seq2Int;
  int SA, SB;
};

/**** The actual function to do the computation ****/
double TKF92LikelihoodFunction(int *seq1Int, int *seq2Int, double len,
    double mu, double r, double distance, gsl_matrix *substModel,
    gsl_vector *eqFrequencies, int SA, int SB);

/**** TKF92 1D and 3D main function for gsl ****/
double TKF92LikelihoodFunction1D(double distance, void *params);
double TKF92LikelihoodFunction3D(const gsl_vector *v,  void *params);

/**** purely calculate the likelihodd given the distance, mu and r ****/
SEXP TKF92LikelihoodFunctionWrapper(SEXP seq1IntR, SEXP seq2IntR, SEXP distanceR, SEXP muR, SEXP rR, SEXP expectedLength, SEXP probMatR, SEXP eqFrequenciesR);


/**** TKF92 1D with Brent implementation ****/
SEXP TKF92LikelihoodFunction1DMain(SEXP seq1IntR, SEXP seq2IntR, 
    SEXP muR, SEXP rR,
    SEXP expectedLength, SEXP probMatR, SEXP eqFrequenciesR);

/**** TKF92 with NM simplex implementation ****/
SEXP TKF92LikelihoodFunction3DMainNM(SEXP seq1IntR, SEXP seq2IntR,
    SEXP expectedLength, SEXP probMatR, SEXP eqFrequenciesR);


