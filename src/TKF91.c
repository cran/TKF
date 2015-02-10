#include "TKF91.h"

double TKF91LikelihoodFunction(int *seq1Int, int *seq2Int, double len,
    double mu, double distance, gsl_matrix *substModel, 
    gsl_vector *eqFrequencies, int SA, int SB){
  double lambda = len / (len + 1.0) * mu;
  double alpha = -mu * distance;
  double lmt = exp((lambda-mu)*distance);
  double lbeta = gsl_log1p(-exp((lambda-mu)*distance)) - (log(mu) + gsl_log1p(-lambda/mu * exp((lambda-mu)*distance)));
  double beta = exp(lbeta); // beta is  not a very small number.
  double P1t = exp(- mu * distance) * (1.0 - lambda * beta);
  double lP12t = gsl_log1p(-lambda * beta);
  double lP01t = log(mu) + lbeta;
  double P11t = (-gsl_expm1(-mu*distance) - mu * beta) * (1.0 - lambda * beta);
  // initialize the entries tables, only log-likelihood is stored in thie table
  gsl_matrix *L0 = gsl_matrix_alloc(SA+1, SB+1);
  gsl_matrix *L1 = gsl_matrix_alloc(SA+1, SB+1);
  gsl_matrix *L2 = gsl_matrix_alloc(SA+1, SB+1);
  // initialize the boundary conditions
  gsl_matrix_set(L0, 0, 0, -INFINITY);
  gsl_matrix_set(L2, 0, 0, -INFINITY);
  gsl_matrix_set(L1, 0, 0, lP12t + gsl_log1p(-lambda/mu));


  int i, j;
  double temp;
  temp = 0;

 
  for(i = 1; i <= SA; i++){
    temp = temp + log(gsl_vector_get(eqFrequencies, seq1Int[i-1])) + lP01t;
    gsl_matrix_set(L0, i, 0, gsl_log1p(-lambda/mu) + i * (log(lambda) - log(mu)) + lP12t + temp);
    gsl_matrix_set(L1, i, 0, -INFINITY);
    gsl_matrix_set(L2, i, 0, -INFINITY);
  }
  temp = 0;
  for(j = 1; j <= SB; j++){
    temp = temp + log(gsl_vector_get(eqFrequencies, seq2Int[j-1]));
    gsl_matrix_set(L2, 0, j, gsl_log1p(-lambda/mu) + lP12t + j * (log(lambda) + lbeta) + temp);
    gsl_matrix_set(L1, 0, j, -INFINITY);
    gsl_matrix_set(L0, 0, j, -INFINITY);
  }


  //recursive iteration
  for(i = 1; i <= SA; i++){
    for(j = 1; j <= SB; j++){
      temp = GSL_MAX(GSL_MAX(gsl_matrix_get(L0, i-1, j),
            gsl_matrix_get(L1, i-1, j)),
          gsl_matrix_get(L2, i-1, j));
      gsl_matrix_set(L0, i, j, log(lambda) - log(mu) + log(gsl_vector_get(eqFrequencies, seq1Int[i-1])) + lP01t + temp + log(exp(gsl_matrix_get(L0, i-1, j) - temp) + exp(gsl_matrix_get(L1, i-1, j) - temp) + exp(gsl_matrix_get(L2, i-1, j) - temp)));
      temp = GSL_MAX(GSL_MAX(gsl_matrix_get(L0, i-1, j-1), gsl_matrix_get(L1, i-1, j-1)), gsl_matrix_get(L2, i-1, j-1));
      gsl_matrix_set(L1, i, j, log(lambda) - log(mu) + log(gsl_vector_get(eqFrequencies, seq1Int[i-1])) + log(gsl_matrix_get(substModel, seq1Int[i-1], seq2Int[j-1]) * P1t + gsl_vector_get(eqFrequencies, seq2Int[j-1]) * P11t) + temp + log(exp(gsl_matrix_get(L0, i-1, j-1) - temp) + exp(gsl_matrix_get(L1, i-1, j-1) - temp) + exp(gsl_matrix_get(L2, i-1, j-1) - temp)));

      if(isfinite(gsl_matrix_get(L1, i, j-1)) || isfinite(gsl_matrix_get(L2, i, j-1))){
        temp = GSL_MAX(gsl_matrix_get(L1, i, j-1), gsl_matrix_get(L2, i,j-1));
        gsl_matrix_set(L2, i, j, log(gsl_vector_get(eqFrequencies, seq2Int[j-1])) + log(lambda) + lbeta + temp + log(exp(gsl_matrix_get(L1, i, j-1) - temp) + exp(gsl_matrix_get(L2, i, j-1) - temp)));
      }else{
        gsl_matrix_set(L2, i, j, -INFINITY);
      }
    }
  }

  temp = GSL_MAX(GSL_MAX(gsl_matrix_get(L0, SA, SB), gsl_matrix_get(L1, SA, SB)), gsl_matrix_get(L2, SA, SB));
  double likelihood;
  likelihood = -(temp + log(exp(gsl_matrix_get(L0, SA, SB) - temp) + exp(gsl_matrix_get(L1, SA, SB) - temp) + exp(gsl_matrix_get(L2, SA, SB) - temp)));


  // free the allocated matrix
  gsl_matrix_free(L0);
  gsl_matrix_free(L1);
  gsl_matrix_free(L2);
  return likelihood;
}


double TKF91LikelihoodFunction1D(double distance, void *params){
  struct TKF91LikelihoodFunction1D_params *p = (struct TKF91LikelihoodFunction1D_params *) params;
  double len = p->len;
  double mu = p->mu;
  gsl_matrix *substModel = gsl_matrix_alloc(p->substModel->size1, p->substModel->size2);
  PAMn(p->substModel, distance, substModel);
  gsl_vector *eqFrequencies = p->eqFrequencies;
  int *seq1Int = p->seq1Int;
  int *seq2Int = p->seq2Int;
  int SA = p->SA;
  int SB = p->SB;
  double likelihood;
  likelihood = TKF91LikelihoodFunction(seq1Int, seq2Int, len, mu, distance, 
      substModel, eqFrequencies, SA, SB);
 
  // free the allocated matrix
  gsl_matrix_free(substModel);
  return likelihood;
}

SEXP TKF91LikelihoodFunctionWrapper(SEXP seq1IntR, SEXP seq2IntR, SEXP distanceR, SEXP muR, SEXP expectedLength, SEXP probMatR, SEXP eqFrequenciesR){
  double len = REAL(expectedLength)[0];
  double mu = REAL(muR)[0];
  double distance = REAL(distanceR)[0];
  int ncol, nrow;
  ncol = INTEGER(GET_DIM(probMatR))[1];
  nrow = INTEGER(GET_DIM(probMatR))[0];
  int i, j;
  // probMat
  gsl_matrix *probMat = gsl_matrix_alloc(nrow, ncol);
  for(i = 0; i < nrow; i++)
    for(j = 0; j < ncol; j++)
      gsl_matrix_set(probMat, i, j, REAL(probMatR)[i+j*ncol]);
  gsl_matrix *substModel = gsl_matrix_alloc(probMat->size1, probMat->size2);
  PAMn(probMat, distance, substModel);

  // eqFrequenciesR
  gsl_vector *eqFrequencies = gsl_vector_alloc(GET_LENGTH(eqFrequenciesR));
  for(i = 0; i < GET_LENGTH(eqFrequenciesR); i++){
    gsl_vector_set(eqFrequencies, i, REAL(eqFrequenciesR)[i]);
  }

  // seqInt preparation
  int *seq1Int, *seq2Int;
  seq1Int = (int *) R_alloc(GET_LENGTH(seq1IntR), sizeof(int));
  seq2Int = (int *) R_alloc(GET_LENGTH(seq2IntR), sizeof(int));
  for(i = 0; i < GET_LENGTH(seq1IntR); i++){
    seq1Int[i] = INTEGER(seq1IntR)[i];
  }
  for(i = 0; i < GET_LENGTH(seq2IntR); i++){
    seq2Int[i] = INTEGER(seq2IntR)[i];
  }
  int SA = GET_LENGTH(seq1IntR);
  int SB = GET_LENGTH(seq2IntR);
  
  double likelihood;
  likelihood = TKF91LikelihoodFunction(seq1Int, seq2Int, len, mu, distance,
      substModel, eqFrequencies, SA, SB);

  SEXP ans, ansNames;
  PROTECT(ans = NEW_NUMERIC(3)); // a vector of distance, mu and the negative log likelihood
  PROTECT(ansNames = NEW_CHARACTER(3));
  REAL(ans)[0] = REAL(distanceR)[0];
  REAL(ans)[1] = REAL(muR)[0];
  REAL(ans)[2] = likelihood;
  SET_STRING_ELT(ansNames, 0, mkChar("PAM"));
  SET_STRING_ELT(ansNames, 1, mkChar("Mu"));
  SET_STRING_ELT(ansNames, 2, mkChar("negLogLikelihood"));
  SET_NAMES(ans, ansNames);

  // free the allocated matrix
  gsl_matrix_free(substModel);
  gsl_matrix_free(probMat);
  gsl_vector_free(eqFrequencies);
  
  UNPROTECT(2);
  return ans;
}

double TKF91LikelihoodFunction2D(const gsl_vector *v,  void *params){
  if(gsl_vector_ispos(v) != 1){
    return GSL_POSINF;
  }
  double distance, mu;
  struct TKF91LikelihoodFunction2D_params *p = (struct TKF91LikelihoodFunction2D_params *) params;
  double len = p->len;
  distance = gsl_vector_get(v, 0);
  mu = gsl_vector_get(v, 1);
  gsl_matrix *substModel = gsl_matrix_alloc(p->substModel->size1, p->substModel->size2);
  PAMn(p->substModel, distance, substModel);
  gsl_vector *eqFrequencies = p->eqFrequencies;
  int *seq1Int = p->seq1Int;
  int *seq2Int = p->seq2Int;
  int SA = p->SA;
  int SB = p->SB;
  double likelihood;
  likelihood = TKF91LikelihoodFunction(seq1Int, seq2Int, len, mu, distance,
      substModel, eqFrequencies, SA, SB);
  // free the allocated matrix
  gsl_matrix_free(substModel);
  //Rprintf("The mu, distance and 2D likelihood is \n");
  //Rprintf("%f\t%f\t%f\n", mu, distance, likelihood);
  return likelihood;
}

void TKF91LikelihoodFunction2D_df(const gsl_vector *v, void *params,
    gsl_vector *df){
  gsl_vector *mEps = gsl_vector_alloc(2);
  gsl_vector_set(mEps, 0, 0.1);  // from Tools/aligner.cpp, tkf2d.setEps(Vector2d(0.1,0.01));
  gsl_vector_set(mEps, 1, 0.01); // from Tools/aligner.cpp
  //double mEps = 1e-4;
  gsl_vector *tempV = gsl_vector_alloc(2);
  gsl_vector_memcpy(tempV, v);

  double temp_dLikelihood;
  // df/d_distance
  gsl_vector_set(tempV, 0, gsl_vector_get(tempV, 0) + 
      gsl_vector_get(mEps, 0)/2);
  temp_dLikelihood = TKF91LikelihoodFunction2D(tempV, params);
  gsl_vector_set(tempV, 0, gsl_vector_get(tempV, 0) - gsl_vector_get(mEps, 0));
  temp_dLikelihood -= TKF91LikelihoodFunction2D(tempV, params);
  gsl_vector_set(tempV, 0, gsl_vector_get(tempV, 0) + 
    gsl_vector_get(mEps, 0)/2);
  temp_dLikelihood /= gsl_vector_get(mEps, 0);
  gsl_vector_set(df, 0, temp_dLikelihood);

  // df/d_mu
  gsl_vector_set(tempV, 1, gsl_vector_get(tempV, 1) + 
    gsl_vector_get(mEps, 1)/2);
  temp_dLikelihood = TKF91LikelihoodFunction2D(tempV, params);
  gsl_vector_set(tempV, 1, gsl_vector_get(tempV, 1) - gsl_vector_get(mEps, 1));
  temp_dLikelihood -= TKF91LikelihoodFunction2D(tempV, params);
  gsl_vector_set(tempV, 1, gsl_vector_get(tempV, 1) + 
    gsl_vector_get(mEps, 1)/2);
  temp_dLikelihood /= gsl_vector_get(mEps, 1);
  gsl_vector_set(df, 1, temp_dLikelihood);

  gsl_vector_free(tempV);
  gsl_vector_free(mEps);
}

void TKF91LikelihoodFunction2D_fdf(const gsl_vector *v, void *params,
    double *f, gsl_vector *df){
  *f = TKF91LikelihoodFunction2D(v, params);
  TKF91LikelihoodFunction2D_df(v, params, df);
}


SEXP TKF91LikelihoodFunction1DMain(SEXP seq1IntR, SEXP seq2IntR, SEXP muR,
    SEXP expectedLength, SEXP probMatR, SEXP eqFrequenciesR){
  int ncol, nrow;
  ncol = INTEGER(GET_DIM(probMatR))[1];
  nrow = INTEGER(GET_DIM(probMatR))[0];
  int i, j; 
  
  // probMat
  gsl_matrix *probMat = gsl_matrix_alloc(nrow, ncol);
  for(i = 0; i < nrow; i++) 
    for(j = 0; j < ncol; j++)
      gsl_matrix_set(probMat, i, j, REAL(probMatR)[i+j*ncol]);
  
  // eqFrequenciesR
  gsl_vector *eqFrequencies = gsl_vector_alloc(GET_LENGTH(eqFrequenciesR));
  for(i = 0; i < GET_LENGTH(eqFrequenciesR); i++){
    gsl_vector_set(eqFrequencies, i, REAL(eqFrequenciesR)[i]);
  }

  // seqInt preparation
  int *seq1Int, *seq2Int;
  seq1Int = (int *) R_alloc(GET_LENGTH(seq1IntR), sizeof(int));
  seq2Int = (int *) R_alloc(GET_LENGTH(seq2IntR), sizeof(int));
  for(i = 0; i < GET_LENGTH(seq1IntR); i++){
    seq1Int[i] = INTEGER(seq1IntR)[i];
  }
  for(i = 0; i < GET_LENGTH(seq2IntR); i++){
    seq2Int[i] = INTEGER(seq2IntR)[i];
  }

  // GSL minimizer 
  int status;
  int iter = 0, max_iter = 100;
  const gsl_min_fminimizer_type *T;
  gsl_min_fminimizer *s;
  gsl_function F;
  struct TKF91LikelihoodFunction1D_params params;
  params.len = REAL(expectedLength)[0];
  params.mu = REAL(muR)[0];
  params.substModel = probMat;
  params.eqFrequencies = eqFrequencies;
  params.seq1Int = seq1Int;
  params.seq2Int = seq2Int;
  params.SA = GET_LENGTH(seq1IntR);
  params.SB = GET_LENGTH(seq2IntR);
  F.function = &TKF91LikelihoodFunction1D;
  F.params = &params;
  double x_lo = 0.0494497, x_hi = 2000; 
  double x = 100;
  double mEps = 1e-5;

  T = gsl_min_fminimizer_brent;
  s = gsl_min_fminimizer_alloc (T);
  gsl_min_fminimizer_set (s, &F, x, x_lo, x_hi);
  
  Rprintf("using %s method\n", 
      gsl_min_fminimizer_name (s));
  Rprintf("%5s [%9s, %9s] %9s %10s %9s\n",
      "iter", "lower", "upper", "min", 
      "err", "err(est)");
  
  do
    {
      iter++;
      status = gsl_min_fminimizer_iterate (s);
      x = gsl_min_fminimizer_x_minimum (s);
      x_lo = gsl_min_fminimizer_x_lower (s);
      x_hi = gsl_min_fminimizer_x_upper (s);
      status = gsl_min_test_interval(x_lo, x_hi,
                                     mEps*mEps, mEps);
      if (status == GSL_SUCCESS)
        Rprintf ("Converged:\n");
      Rprintf ("%5d [%.7f, %.7f] "
              "%.7f %.7f\n",
              iter, x_lo,  x_hi,
              x, x_hi - x_lo);
    }
  while (status == GSL_CONTINUE && iter < max_iter);
 
  SEXP ans, ansNames;
  PROTECT(ans = NEW_NUMERIC(3)); // a vector of distance, mu and the negative log likelihood
  PROTECT(ansNames = NEW_CHARACTER(3));
  REAL(ans)[0] = x;
  REAL(ans)[1] = REAL(muR)[0];
  REAL(ans)[2] = gsl_min_fminimizer_f_minimum(s); 
  SET_STRING_ELT(ansNames, 0, mkChar("PAM"));
  SET_STRING_ELT(ansNames, 1, mkChar("Mu"));
  SET_STRING_ELT(ansNames, 2, mkChar("negLogLikelihood"));
  SET_NAMES(ans, ansNames);
  
  // free the allocation
  gsl_matrix_free(probMat);
  gsl_min_fminimizer_free (s);
  gsl_vector_free(eqFrequencies);

  UNPROTECT(2);
  return ans;
}

SEXP TKF91LikelihoodFunction2DMain(SEXP seq1IntR, SEXP seq2IntR,
    SEXP expectedLength, SEXP probMatR, SEXP eqFrequenciesR){
  int ncol, nrow;
  ncol = INTEGER(GET_DIM(probMatR))[1];
  nrow = INTEGER(GET_DIM(probMatR))[0];
  int i, j; 

  // probMat
  gsl_matrix *probMat = gsl_matrix_alloc(nrow, ncol);
  for(i = 0; i < nrow; i++)
    for(j = 0; j < ncol; j++)
      gsl_matrix_set(probMat, i, j, REAL(probMatR)[i+j*ncol]);

  // eqFrequenciesR
  gsl_vector *eqFrequencies = gsl_vector_alloc(GET_LENGTH(eqFrequenciesR));
  for(i = 0; i < GET_LENGTH(eqFrequenciesR); i++){
    gsl_vector_set(eqFrequencies, i, REAL(eqFrequenciesR)[i]);
  }

  // seqInt preparation
  int *seq1Int, *seq2Int;
  seq1Int = (int *) R_alloc(GET_LENGTH(seq1IntR), sizeof(int));
  seq2Int = (int *) R_alloc(GET_LENGTH(seq2IntR), sizeof(int));
  for(i = 0; i < GET_LENGTH(seq1IntR); i++){
    seq1Int[i] = INTEGER(seq1IntR)[i];
  }
  for(i = 0; i < GET_LENGTH(seq2IntR); i++){
    seq2Int[i] = INTEGER(seq2IntR)[i];
  }
  
  // GSL minimizer 
  int status;
  int iter = 0, max_iter = 200; // from Tools/aligner.cpp, optim2D.setMaxIterations(200);
  const gsl_multimin_fdfminimizer_type *T = 0;
  gsl_multimin_fdfminimizer *s;
  gsl_multimin_function_fdf F;
  struct TKF91LikelihoodFunction2D_params params;
  params.len = REAL(expectedLength)[0];
  params.substModel = probMat;
  params.eqFrequencies = eqFrequencies;
  params.seq1Int = seq1Int;
  params.seq2Int = seq2Int;
  params.SA = GET_LENGTH(seq1IntR);
  params.SB = GET_LENGTH(seq2IntR);
  
  // starting points
  gsl_vector *x;
  x = gsl_vector_alloc(2);
  gsl_vector_set(x, 0, 100); // distance, from Tools/aligner.cpp,  Vector2d(-3.0, 2.0); in Exp scale.
  gsl_vector_set(x, 1, exp(-3)); // mu, same

  double mAccuracy = 1e-3; // from Tools/aligner.cpp, optim2D.setAccuracy(1e-3);
  double mInitStepSize = 0.01; // from Tools/aligner.cpp, optim2D.setInitialStepSize(0.85);
  // Set initial step sizes 
  //ss = gsl_vector_alloc (2);
  //gsl_vector_set_all(ss, mInitStepSize);

  // Initialize method and iterate
  F.n = 2;
  F.f = &TKF91LikelihoodFunction2D;
  F.df = &TKF91LikelihoodFunction2D_df;
  F.fdf = &TKF91LikelihoodFunction2D_fdf;
  F.params = &params;
  T = gsl_multimin_fdfminimizer_vector_bfgs2;
  double accuracy = 0.1;
  s = gsl_multimin_fdfminimizer_alloc(T, 2);
  gsl_multimin_fdfminimizer_set(s, &F, x, mInitStepSize, accuracy);

  Rprintf("using %s method\n",
      gsl_multimin_fdfminimizer_name(s));
  do
  {
    iter++;
    status = gsl_multimin_fdfminimizer_iterate(s);
    if(status){
      Rprintf("The status code: %d\n", status);
      break;
    }
    status = gsl_multimin_test_gradient(s->gradient, mAccuracy);
    if(status == GSL_SUCCESS){
      Rprintf("converged to minimu at distance  mu  neg-log-likelihood\n");
    }
    Rprintf("%5d %.5f %.5f %10.5f\n",
        iter,
        gsl_vector_get(s->x, 0),
        gsl_vector_get(s->x, 1),
        s->f
        );
  }
  while(status == GSL_CONTINUE && iter < max_iter);
  gsl_vector_free(x);
  //gsl_vector_free(ss);
  gsl_multimin_fdfminimizer_free(s);

  return R_NilValue;
}

SEXP TKF91LikelihoodFunction2DMainNM(SEXP seq1IntR, SEXP seq2IntR,
    SEXP expectedLength, SEXP probMatR, SEXP eqFrequenciesR){
  int ncol, nrow;
  ncol = INTEGER(GET_DIM(probMatR))[1];
  nrow = INTEGER(GET_DIM(probMatR))[0];
  int i, j; 

  // probMat
  gsl_matrix *probMat = gsl_matrix_alloc(nrow, ncol);
  for(i = 0; i < nrow; i++)
    for(j = 0; j < ncol; j++)
      gsl_matrix_set(probMat, i, j, REAL(probMatR)[i+j*ncol]);

  // eqFrequenciesR
  gsl_vector *eqFrequencies = gsl_vector_alloc(GET_LENGTH(eqFrequenciesR));
  for(i = 0; i < GET_LENGTH(eqFrequenciesR); i++){
    gsl_vector_set(eqFrequencies, i, REAL(eqFrequenciesR)[i]);
  }

  // seqInt preparation
  int *seq1Int, *seq2Int;
  seq1Int = (int *) R_alloc(GET_LENGTH(seq1IntR), sizeof(int));
  seq2Int = (int *) R_alloc(GET_LENGTH(seq2IntR), sizeof(int));
  for(i = 0; i < GET_LENGTH(seq1IntR); i++){
    seq1Int[i] = INTEGER(seq1IntR)[i];
  }
  for(i = 0; i < GET_LENGTH(seq2IntR); i++){
    seq2Int[i] = INTEGER(seq2IntR)[i];
  }
  
  // GSL minimizer 
  int status;
  double size;
  int iter = 0, max_iter = 1000; 
  const gsl_multimin_fminimizer_type *T = 0;
  gsl_multimin_fminimizer *s;
  gsl_multimin_function F;
  struct TKF91LikelihoodFunction2D_params params;
  params.len = REAL(expectedLength)[0];
  params.substModel = probMat;
  params.eqFrequencies = eqFrequencies;
  params.seq1Int = seq1Int;
  params.seq2Int = seq2Int;
  params.SA = GET_LENGTH(seq1IntR);
  params.SB = GET_LENGTH(seq2IntR);
  
  // starting points
  gsl_vector *x;
  x = gsl_vector_alloc(2);
  gsl_vector_set(x, 0, 100); // distance, from Tools/aligner.cpp,  Vector2d(-3.0, 2.0); in Exp scale.
  gsl_vector_set(x, 1, exp(-3)); // mu, same

  // Set initial step sizes 
  gsl_vector *ss;
  ss = gsl_vector_alloc(2);
  gsl_vector_set(ss, 0, 1);
  gsl_vector_set(ss, 1, 0.01);

  // Initialize method and iterate
  F.n = 2;
  F.f = &TKF91LikelihoodFunction2D;
  F.params = &params;
  T = gsl_multimin_fminimizer_nmsimplex2;
  double accuracy = 0.1;
  s = gsl_multimin_fminimizer_alloc(T, 2);
  gsl_multimin_fminimizer_set(s, &F, x, ss);

  Rprintf("using %s method\n",
      gsl_multimin_fminimizer_name(s));
  do
  {
    iter++;
    status = gsl_multimin_fminimizer_iterate(s);
    if(status){
      Rprintf("The status code: %d\n", status);
      break;
    }
    size = gsl_multimin_fminimizer_size(s);
    status = gsl_multimin_test_size (size, 1e-5);
    if(status == GSL_SUCCESS){
      Rprintf("converged to minimu at \n");
    }
    Rprintf("%5d %10.3e %10.3e f() = %7.3f size = %.5f\n", 
        iter,
        gsl_vector_get (s->x, 0), 
        gsl_vector_get (s->x, 1), 
        s->fval, size);
  }
  while(status == GSL_CONTINUE && iter < max_iter);

  SEXP ans, ansNames;
  PROTECT(ans = NEW_NUMERIC(3)); // a vector of distance, mu and the negative log likelihood
  PROTECT(ansNames = NEW_CHARACTER(3));
  REAL(ans)[0] = gsl_vector_get (s->x, 0);
  REAL(ans)[1] = gsl_vector_get (s->x, 1);
  REAL(ans)[2] = s->fval;
  SET_STRING_ELT(ansNames, 0, mkChar("PAM"));
  SET_STRING_ELT(ansNames, 1, mkChar("Mu"));
  SET_STRING_ELT(ansNames, 2, mkChar("negLogLikelihood"));
  SET_NAMES(ans, ansNames);

  // free everything
  gsl_vector_free(x);
  gsl_vector_free(ss);
  gsl_multimin_fminimizer_free(s);

  UNPROTECT(2);
  return ans;
}


