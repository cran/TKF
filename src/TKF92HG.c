#include "TKF92HG.h"

double TKF92HGLikelihoodFunction(int *seq1Int, int *seq2Int, double len,
    double mu, double r, double ps, double kf,
    double distance, gsl_matrix *substModelS, gsl_matrix *substModelF, 
    gsl_vector *eqFrequencies, int SA, int SB){
  //double lambda = len / (len + 1.0) * mu; TKF91
  // The way to calculate the lambda, from Page 9, formular 10.
  double lambda = mu * (len * (1.0 - r) - 2.0 * r + sqrt(len * len * (1.0 - r) * (1.0 - r) + 4.0 * r)) / (2.0 * (len + 1.0) * (1.0 - r));
  //double alpha = -mu * distance;
  //double lmt = exp((lambda-mu)*distance);
  double lbeta = gsl_log1p(-exp((lambda-mu)*distance)) - (log(mu) + gsl_log1p(-lambda/mu * exp((lambda-mu)*distance)));
  double beta = exp(lbeta); // beta is  not a very small number.
  double P1t = exp(- mu * distance) * (1.0 - lambda * beta);
  double lP12t = gsl_log1p(-lambda * beta);
  double lP01t = log(mu) + lbeta;
  double P11t = (-gsl_expm1(-mu*distance) - mu * beta) * (1.0 - lambda * beta);
  double lP22t = gsl_log1p(-lambda * beta) + log(lambda) + lbeta;
  double pf = 1.0 - ps;

  // initialize the entries tables, only log-likelihood is stored in thie table
  gsl_matrix *L1 = gsl_matrix_alloc(SA+1, SB+1);
  gsl_matrix *L1S = gsl_matrix_alloc(SA+1, SB+1); // slow fragment
  gsl_matrix *L1F = gsl_matrix_alloc(SA+1, SB+1); // fast fragment
  gsl_matrix *L2 = gsl_matrix_alloc(SA+1, SB+1);
  gsl_matrix *L3 = gsl_matrix_alloc(SA+1, SB+1);
  gsl_matrix *L4 = gsl_matrix_alloc(SA+1, SB+1);
  gsl_matrix *L5 = gsl_matrix_alloc(SA+1, SB+1);
  gsl_matrix *L6 = gsl_matrix_alloc(SA+1, SB+1);

  double temp;
  // initialize the boundary conditions
  gsl_matrix_set(L1S, 0, 0, gsl_log1p(-lambda/mu) + lP12t);
  gsl_matrix_set(L1F, 0, 0, gsl_log1p(-lambda/mu) + lP12t);
  temp = GSL_MAX(gsl_matrix_get(L1S, 0, 0), gsl_matrix_get(L1F, 0, 0));
  gsl_matrix_set(L1, 0, 0, 
      temp + log(exp(gsl_matrix_get(L1S, 0, 0) - temp) + 
                 exp(gsl_matrix_get(L1F, 0, 0) - temp)));
  gsl_matrix_set(L2, 0, 0, -INFINITY);
  gsl_matrix_set(L3, 0, 0, -INFINITY);
  gsl_matrix_set(L4, 0, 0, -INFINITY);
  gsl_matrix_set(L5, 0, 0, -INFINITY);
  gsl_matrix_set(L6, 0, 0, -INFINITY);
  
  gsl_matrix_set(L2, 1, 0, log((1-lambda/mu)*lambda/mu*(1-r)) + lP12t + log(gsl_vector_get(eqFrequencies, seq1Int[0])) + lP01t);
  gsl_matrix_set(L1S, 1, 0, -INFINITY);
  gsl_matrix_set(L1F, 1, 0, -INFINITY);
  gsl_matrix_set(L1, 1, 0, -INFINITY);
  gsl_matrix_set(L3, 1, 0, -INFINITY);
  gsl_matrix_set(L4, 1, 0, -INFINITY);
  gsl_matrix_set(L5, 1, 0, -INFINITY);
  gsl_matrix_set(L6, 1, 0, -INFINITY);
  int i, j;
  temp = 0;
  for(i = 2; i <= SA; i++){
    temp = temp + log(gsl_vector_get(eqFrequencies, seq1Int[i-1])) + log(r) + gsl_log1p(exp(lP01t) * (1-r) / r * lambda / mu);
    gsl_matrix_set(L2, i, 0, log((1-lambda/mu)*lambda/mu*(1-r)) + lP12t + log(gsl_vector_get(eqFrequencies, seq1Int[0])) + lP01t + temp);
    gsl_matrix_set(L1S, i, 0, -INFINITY);
    gsl_matrix_set(L1F, i, 0, -INFINITY);
    gsl_matrix_set(L1, i, 0, -INFINITY);
    gsl_matrix_set(L3, i, 0, -INFINITY);
    gsl_matrix_set(L4, i, 0, -INFINITY);
    gsl_matrix_set(L5, i, 0, -INFINITY);
    gsl_matrix_set(L6, i, 0, -INFINITY);
  }

  gsl_matrix_set(L3, 0, 1, gsl_log1p(-lambda/mu) + lP22t + gsl_log1p(-r) + log(gsl_vector_get(eqFrequencies, seq2Int[0])));
  gsl_matrix_set(L1S, 0, 1, -INFINITY);
  gsl_matrix_set(L1F, 0, 1, -INFINITY);
  gsl_matrix_set(L1, 0, 1, -INFINITY);
  gsl_matrix_set(L2, 0, 1, -INFINITY);
  gsl_matrix_set(L4, 0, 1, -INFINITY);
  gsl_matrix_set(L5, 0, 1, -INFINITY);
  gsl_matrix_set(L6, 0, 1, -INFINITY);
  temp = 0;
  for(j = 2; j <= SB; j++){
    temp = temp + log(gsl_vector_get(eqFrequencies, seq2Int[j-1])) + log(r) + gsl_log1p(lambda/r*beta*(1-r));
    gsl_matrix_set(L3, 0, j, gsl_log1p(-lambda/mu) + lP22t + gsl_log1p(-r) + log(gsl_vector_get(eqFrequencies, seq2Int[0])) + temp);
    gsl_matrix_set(L1S, 0, j, -INFINITY);
    gsl_matrix_set(L1F, 0, j, -INFINITY);
    gsl_matrix_set(L1, 0, j, -INFINITY);
    gsl_matrix_set(L2, 0, j, -INFINITY);
    gsl_matrix_set(L4, 0, j, -INFINITY);
    gsl_matrix_set(L5, 0, j, -INFINITY);
    gsl_matrix_set(L6, 0, j, -INFINITY);
  }


  //recursive iteration
  double Km, Kn;
  double A, B;
  for(i = 1; i <= SA; i++){
    for(j = 1; j <= SB; j++){
      if(i == 1){
        Km = 0.0;
      }else{
        Km = 1.0;
      }
      if(j == 1){
        Kn = 0.0;
      }else{
        Kn = 1.0;
      }
      // we need some tricks to calculate the log of sums of exponentials
      A = r * Km * Kn;
      B = P1t * (1 - r) * lambda / mu;
      temp = GSL_MAX(GSL_MAX(GSL_MAX(GSL_MAX(GSL_MAX(GSL_MAX(
                  gsl_matrix_get(L1, i-1, j-1),
                  gsl_matrix_get(L2, i-1, j-1)),
                  gsl_matrix_get(L3, i-1, j-1)),
                  gsl_matrix_get(L4, i-1, j-1)),
                  gsl_matrix_get(L5, i-1, j-1)),
                  gsl_matrix_get(L6, i-1, j-1)),
                  gsl_matrix_get(L1S, i-1, j-1)
                  );
      gsl_matrix_set(L1S, i, j, 
          log(gsl_vector_get(eqFrequencies, seq1Int[i-1])) + 
          log(gsl_matrix_get(substModelS, seq1Int[i-1], seq2Int[j-1])) +
          temp +
          log(A*exp(gsl_matrix_get(L1S, i-1, j-1) - temp) +
              ps*B*(exp(gsl_matrix_get(L1, i-1, j-1) - temp) +
                    exp(gsl_matrix_get(L2, i-1, j-1) - temp) +
                    exp(gsl_matrix_get(L3, i-1, j-1) - temp) +
                    exp(gsl_matrix_get(L4, i-1, j-1) - temp) +
                    exp(gsl_matrix_get(L5, i-1, j-1) - temp) +
                    exp(gsl_matrix_get(L6, i-1, j-1) - temp))
            )
          );

      temp = GSL_MAX(GSL_MAX(GSL_MAX(GSL_MAX(GSL_MAX(GSL_MAX(
                    gsl_matrix_get(L1, i-1, j-1),
                    gsl_matrix_get(L2, i-1, j-1)),
                    gsl_matrix_get(L3, i-1, j-1)),
                    gsl_matrix_get(L4, i-1, j-1)),
                    gsl_matrix_get(L5, i-1, j-1)),
                    gsl_matrix_get(L6, i-1, j-1)),
                    gsl_matrix_get(L1F, i-1, j-1)
                    );
      gsl_matrix_set(L1F, i, j,
          log(gsl_vector_get(eqFrequencies, seq1Int[i-1])) +
          log(gsl_matrix_get(substModelF, seq1Int[i-1], seq2Int[j-1])) +
          temp +
          log(A*exp(gsl_matrix_get(L1F, i-1, j-1) - temp) +
              pf*B*(exp(gsl_matrix_get(L1, i-1, j-1) - temp) +
                    exp(gsl_matrix_get(L2, i-1, j-1) - temp) +
                    exp(gsl_matrix_get(L3, i-1, j-1) - temp) +
                    exp(gsl_matrix_get(L4, i-1, j-1) - temp) +
                    exp(gsl_matrix_get(L5, i-1, j-1) - temp) +
                    exp(gsl_matrix_get(L6, i-1, j-1) - temp))
            )
          );

      temp = GSL_MAX(gsl_matrix_get(L1S, i, j), gsl_matrix_get(L1F, i, j));
      gsl_matrix_set(L1, i, j,
         temp + log(exp(gsl_matrix_get(L1S, i, j) - temp) +
                    exp(gsl_matrix_get(L1F, i, j) - temp))); 


      A = r * Km;
      B = exp(lP01t) * (1 - r)  * lambda / mu;
      temp = GSL_MAX(GSL_MAX(GSL_MAX(GSL_MAX(GSL_MAX(
                  gsl_matrix_get(L1, i-1, j),
                  gsl_matrix_get(L2, i-1, j)),
                  gsl_matrix_get(L3, i-1, j)),
                  gsl_matrix_get(L4, i-1, j)),
                  gsl_matrix_get(L5, i-1, j)),
                  gsl_matrix_get(L6, i-1, j)
                  );
      gsl_matrix_set(L2, i, j,
          log(gsl_vector_get(eqFrequencies, seq1Int[i-1])) +
          temp +
          log(A*exp(gsl_matrix_get(L2, i-1, j) - temp) +
              B*(exp(gsl_matrix_get(L1, i-1, j) - temp) +
                 exp(gsl_matrix_get(L2, i-1, j) - temp) +
                 exp(gsl_matrix_get(L3, i-1, j) - temp) +
                 exp(gsl_matrix_get(L4, i-1, j) - temp) + 
                 exp(gsl_matrix_get(L5, i-1, j) - temp) + 
                 exp(gsl_matrix_get(L6, i-1, j) - temp))
            )
          );


      if(isfinite(gsl_matrix_get(L1, i, j-1)) || isfinite(gsl_matrix_get(L3, i, j-1)) || isfinite(gsl_matrix_get(L4, i, j-1)) || isfinite(gsl_matrix_get(L5, i, j-1)) || isfinite(gsl_matrix_get(L6, i, j-1))){
        A = r * Kn;
        B = lambda * beta * (1 - r);
        temp = GSL_MAX(GSL_MAX(GSL_MAX(GSL_MAX(
                    gsl_matrix_get(L1, i, j-1),
                    gsl_matrix_get(L3, i, j-1)),
                    gsl_matrix_get(L4, i, j-1)),
                    gsl_matrix_get(L5, i, j-1)),
                    gsl_matrix_get(L6, i, j-1)
              );
        gsl_matrix_set(L3, i, j,
            log(gsl_vector_get(eqFrequencies, seq2Int[j-1])) + 
            temp +
            log(A*exp(gsl_matrix_get(L3, i, j-1) - temp) +
                B*(exp(gsl_matrix_get(L1, i, j-1) - temp) +
                   exp(gsl_matrix_get(L3, i, j-1) - temp) +
                   exp(gsl_matrix_get(L4, i, j-1) - temp) + 
                   exp(gsl_matrix_get(L5, i, j-1) - temp) +
                   exp(gsl_matrix_get(L6, i, j-1) - temp))
              )
            );
      }else{
        gsl_matrix_set(L3, i, j, -INFINITY);
      }

      A = r * r * Km * Kn;
      B = P11t * pow((1 - r), 2) * lambda / mu;
      temp = GSL_MAX(GSL_MAX(GSL_MAX(GSL_MAX(GSL_MAX(
                  gsl_matrix_get(L1, i-1, j-1),
                  gsl_matrix_get(L2, i-1, j-1)),
                  gsl_matrix_get(L3, i-1, j-1)),
                  gsl_matrix_get(L4, i-1, j-1)),
                  gsl_matrix_get(L5, i-1, j-1)),
                  gsl_matrix_get(L6, i-1, j-1)
                  );
      gsl_matrix_set(L4, i, j,
          log(gsl_vector_get(eqFrequencies, seq1Int[i-1])) + 
          log(gsl_vector_get(eqFrequencies, seq2Int[j-1])) +
          temp +
          log(A*exp(gsl_matrix_get(L4, i-1, j-1) - temp) +
              B*(
                exp(gsl_matrix_get(L1, i-1, j-1) - temp) +
                exp(gsl_matrix_get(L2, i-1, j-1) - temp) +
                exp(gsl_matrix_get(L3, i-1, j-1) - temp) +
                exp(gsl_matrix_get(L4, i-1, j-1) - temp) +
                exp(gsl_matrix_get(L5, i-1, j-1) - temp) +
                exp(gsl_matrix_get(L6, i-1, j-1) - temp))
            )
          );

      if(isfinite(gsl_matrix_get(L4, i-1, j)) || isfinite(gsl_matrix_get(L5, i-1, j))){
        temp = GSL_MAX(gsl_matrix_get(L4, i-1, j),
                       gsl_matrix_get(L5, i-1, j)
                       );
        gsl_matrix_set(L5, i, j,
            log(gsl_vector_get(eqFrequencies, seq1Int[i-1])) +
            log(Km * Kn) + log(r) + temp +
            log(exp(gsl_matrix_get(L4, i-1, j) - temp) + 
                exp(gsl_matrix_get(L5, i-1, j) - temp))
            );
      }else{
        gsl_matrix_set(L5, i, j, -INFINITY);
      }

      if(isfinite(gsl_matrix_get(L4, i, j-1)) || isfinite(gsl_matrix_get(L6, i,j-1))){
        temp = GSL_MAX(gsl_matrix_get(L4, i, j-1), gsl_matrix_get(L6, i,j-1));
        gsl_matrix_set(L6, i, j,
            log(gsl_vector_get(eqFrequencies, seq2Int[j-1])) +
            log(Km * Kn) + log(r) + temp +
            log(exp(gsl_matrix_get(L4, i, j-1) - temp) +
                exp(gsl_matrix_get(L6, i, j-1) - temp))
            );
      }else{
        gsl_matrix_set(L6, i, j, -INFINITY);
      }
    }
  }

  temp = GSL_MAX(GSL_MAX(GSL_MAX(GSL_MAX(GSL_MAX(gsl_matrix_get(L1, SA, SB), gsl_matrix_get(L2, SA, SB)), gsl_matrix_get(L3, SA, SB)), gsl_matrix_get(L4, SA, SB)), gsl_matrix_get(L5, SA, SB)), gsl_matrix_get(L6, SA, SB));
  double likelihood;
  likelihood = -(temp + log(exp(gsl_matrix_get(L1, SA, SB) - temp) + 
                            exp(gsl_matrix_get(L2, SA, SB) - temp) +
                            exp(gsl_matrix_get(L3, SA, SB) - temp) + 
                            exp(gsl_matrix_get(L4, SA, SB) - temp) +
                            exp(gsl_matrix_get(L5, SA, SB) - temp) + 
                            exp(gsl_matrix_get(L6, SA, SB) - temp))
                );

  // free the allocated matrix
  gsl_matrix_free(L1);
  gsl_matrix_free(L1S);
  gsl_matrix_free(L1F);
  gsl_matrix_free(L2);
  gsl_matrix_free(L3);
  gsl_matrix_free(L4);
  gsl_matrix_free(L5);
  gsl_matrix_free(L6);

//  Rprintf("The distance, mu, r and likelihood are: \n");
//  Rprintf("%f\t%f\t%f\t%f\n", distance, mu, r, likelihood);

  return likelihood;
}


double TKF92HGLikelihoodFunction1D(double distance, void *params){
  struct TKF92HGLikelihoodFunction1D_params *p = (struct TKF92HGLikelihoodFunction1D_params *) params;
  double len = p->len;
  double mu = p->mu;
  double r = p->r;
  double ps = p->ps;
  double kf = p->kf;

  gsl_matrix *substModelS = gsl_matrix_alloc(p->substModel->size1, p->substModel->size2);
  gsl_matrix *substModelF = gsl_matrix_alloc(p->substModel->size1, p->substModel->size2);
  PAMn(p->substModel, distance, substModelS);
  PAMn(p->substModel, distance * kf, substModelF);
  gsl_vector *eqFrequencies = p->eqFrequencies;
  int *seq1Int = p->seq1Int;
  int *seq2Int = p->seq2Int;
  int SA = p->SA;
  int SB = p->SB;
  double likelihood;
  likelihood = TKF92HGLikelihoodFunction(seq1Int, seq2Int, len, 
      mu, r, ps, kf, distance, 
      substModelS, substModelF, eqFrequencies, SA, SB);
 
  // free the allocated matrix
  gsl_matrix_free(substModelS);
  gsl_matrix_free(substModelF);
  return likelihood;
}

SEXP TKF92HGLikelihoodFunctionWrapper(SEXP seq1IntR, SEXP seq2IntR, SEXP distanceR, SEXP muR, SEXP rR, SEXP psR, SEXP kfR, SEXP expectedLength, SEXP probMatR, SEXP eqFrequenciesR){
  double len = REAL(expectedLength)[0];
  double mu = REAL(muR)[0];
  double r = REAL(rR)[0];
  double ps = REAL(psR)[0];
  double kf = REAL(kfR)[0];

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
  gsl_matrix *substModelS = gsl_matrix_alloc(probMat->size1, probMat->size2);
  gsl_matrix *substModelF = gsl_matrix_alloc(probMat->size1, probMat->size2);
  PAMn(probMat, distance, substModelS);
  PAMn(probMat, distance * kf, substModelF);

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
  likelihood = TKF92HGLikelihoodFunction(seq1Int, seq2Int, len, 
      mu, r, ps, kf, distance,
      substModelS, substModelF, eqFrequencies, SA, SB);

  SEXP ans, ansNames;
  PROTECT(ans = NEW_NUMERIC(6)); // a vector of distance, mu, r, ps, kf and the negative log likelihood
  PROTECT(ansNames = NEW_CHARACTER(6));
  REAL(ans)[0] = REAL(distanceR)[0];
  REAL(ans)[1] = REAL(muR)[0];
  REAL(ans)[2] = REAL(rR)[0];
  REAL(ans)[3] = REAL(psR)[0];
  REAL(ans)[4] = REAL(kfR)[0];
  REAL(ans)[5] = likelihood;

  SET_STRING_ELT(ansNames, 0, mkChar("PAM"));
  SET_STRING_ELT(ansNames, 1, mkChar("Mu"));
  SET_STRING_ELT(ansNames, 2, mkChar("r"));
  SET_STRING_ELT(ansNames, 3, mkChar("Ps"));
  SET_STRING_ELT(ansNames, 4, mkChar("Kf")); 
  SET_STRING_ELT(ansNames, 5, mkChar("negLogLikelihood"));
  SET_NAMES(ans, ansNames);

  // free the allocated matrix
  gsl_matrix_free(substModelS);
  gsl_matrix_free(substModelF);
  gsl_matrix_free(probMat);
  gsl_vector_free(eqFrequencies);
  
  UNPROTECT(2);
  return ans;
}

double TKF92HGLikelihoodFunction5D(const gsl_vector *v,  void *params){
  if(gsl_vector_ispos(v) != 1){
    return GSL_POSINF;
  }
  double distance, mu, r, ps, kf;
  struct TKF92HGLikelihoodFunction5D_params *p = (struct TKF92HGLikelihoodFunction5D_params *) params;
  double len = p->len;
  distance = gsl_vector_get(v, 0);
  mu = gsl_vector_get(v, 1);
  r = gsl_vector_get(v, 2);
  ps = gsl_vector_get(v, 3);
  kf = gsl_vector_get(v, 4);
  gsl_matrix *substModelS = gsl_matrix_alloc(p->substModel->size1, p->substModel->size2);
  gsl_matrix *substModelF = gsl_matrix_alloc(p->substModel->size1, p->substModel->size2);
  PAMn(p->substModel, distance, substModelS);
  PAMn(p->substModel, distance * kf, substModelF);
  gsl_vector *eqFrequencies = p->eqFrequencies;
  int *seq1Int = p->seq1Int;
  int *seq2Int = p->seq2Int;
  int SA = p->SA;
  int SB = p->SB;
  double likelihood;
  likelihood = TKF92HGLikelihoodFunction(seq1Int, seq2Int, len, 
      mu, r, ps, kf, distance,
      substModelS, substModelF, eqFrequencies, SA, SB);
  // free the allocated matrix
  gsl_matrix_free(substModelS);
  gsl_matrix_free(substModelF);
  //Rprintf("The mu, distance and 3D likelihood is \n");
  //Rprintf("%f\t%f\t%f\n", mu, distance, likelihood);
  return likelihood;
}


SEXP TKF92HGLikelihoodFunction1DMain(SEXP seq1IntR, SEXP seq2IntR, SEXP muR,
    SEXP rR, SEXP psR, SEXP kfR, SEXP expectedLength, SEXP probMatR, 
    SEXP eqFrequenciesR){
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
  struct TKF92HGLikelihoodFunction1D_params params;
  params.len = REAL(expectedLength)[0];
  params.mu = REAL(muR)[0];
  params.r = REAL(rR)[0];
  params.ps = REAL(psR)[0];
  params.kf = REAL(kfR)[0];
  params.substModel = probMat;
  params.eqFrequencies = eqFrequencies;
  params.seq1Int = seq1Int;
  params.seq2Int = seq2Int;
  params.SA = GET_LENGTH(seq1IntR);
  params.SB = GET_LENGTH(seq2IntR);
  F.function = &TKF92HGLikelihoodFunction1D;
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
  PROTECT(ans = NEW_NUMERIC(6)); // a vector of distance, mu and the negative log likelihood
  PROTECT(ansNames = NEW_CHARACTER(6));
  REAL(ans)[0] = x;
  REAL(ans)[1] = REAL(muR)[0];
  REAL(ans)[2] = REAL(rR)[0];
  REAL(ans)[3] = REAL(psR)[0];
  REAL(ans)[4] = REAL(kfR)[0];
  REAL(ans)[5] = gsl_min_fminimizer_f_minimum(s); 

  SET_STRING_ELT(ansNames, 0, mkChar("PAM"));
  SET_STRING_ELT(ansNames, 1, mkChar("Mu"));
  SET_STRING_ELT(ansNames, 2, mkChar("r"));
  SET_STRING_ELT(ansNames, 3, mkChar("Ps"));
  SET_STRING_ELT(ansNames, 4, mkChar("Kf"));
  SET_STRING_ELT(ansNames, 5, mkChar("negLogLikelihood"));
  SET_NAMES(ans, ansNames);
  
  // free the allocation
  gsl_matrix_free(probMat);
  gsl_min_fminimizer_free (s);
  gsl_vector_free(eqFrequencies);

  UNPROTECT(2);
  return ans;
}


SEXP TKF92HGLikelihoodFunction5DMainNM(SEXP seq1IntR, SEXP seq2IntR,
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
  struct TKF92HGLikelihoodFunction5D_params params;
  params.len = REAL(expectedLength)[0];
  params.substModel = probMat;
  params.eqFrequencies = eqFrequencies;
  params.seq1Int = seq1Int;
  params.seq2Int = seq2Int;
  params.SA = GET_LENGTH(seq1IntR);
  params.SB = GET_LENGTH(seq2IntR);
  
  // starting points
  gsl_vector *x;
  x = gsl_vector_alloc(5);
  gsl_vector_set(x, 0, 100); // distance, from Tools/aligner.cpp,  Vector2d(-3.0, 2.0); in Exp scale.
  gsl_vector_set(x, 1, exp(-3)); // mu, same
  gsl_vector_set(x, 2, 0.5); // r
  gsl_vector_set(x, 3, 0.5); // ps
  gsl_vector_set(x, 4, 1.2); // kf

  // Set initial step sizes 
  gsl_vector *ss;
  ss = gsl_vector_alloc(5);
  gsl_vector_set(ss, 0, 1);
  gsl_vector_set(ss, 1, 0.01);
  gsl_vector_set(ss, 2, 0.1);
  gsl_vector_set(ss, 3, 0.1);
  gsl_vector_set(ss, 4, 0.1);

  // Initialize method and iterate
  F.n = 5;
  F.f = &TKF92HGLikelihoodFunction5D;
  F.params = &params;
  T = gsl_multimin_fminimizer_nmsimplex2;
  s = gsl_multimin_fminimizer_alloc(T, 5);
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
      Rprintf("converged to minimu at distance  mu  r  ps, kf neg-log-likelihood  size\n");
    }
    Rprintf("%5d %10.3e %10.3e %10.3e %10.3e %10.3e f() = %7.3f size = %.5f\n", 
        iter,
        gsl_vector_get (s->x, 0), 
        gsl_vector_get (s->x, 1), 
        gsl_vector_get (s->x, 2),
        gsl_vector_get (s->x, 3),
        gsl_vector_get (s->x, 4),
        s->fval, size);
  }
  while(status == GSL_CONTINUE && iter < max_iter);

  SEXP ans, ansNames;
  PROTECT(ans = NEW_NUMERIC(6)); // a vector of distance, mu and the negative log likelihood
  PROTECT(ansNames = NEW_CHARACTER(6));
  REAL(ans)[0] = gsl_vector_get (s->x, 0);
  REAL(ans)[1] = gsl_vector_get (s->x, 1);
  REAL(ans)[2] = gsl_vector_get(s->x, 2);
  REAL(ans)[3] = gsl_vector_get(s->x, 3);
  REAL(ans)[4] = gsl_vector_get(s->x, 4);
  REAL(ans)[5] = s->fval;
  SET_STRING_ELT(ansNames, 0, mkChar("PAM"));
  SET_STRING_ELT(ansNames, 1, mkChar("Mu"));
  SET_STRING_ELT(ansNames, 2, mkChar("r"));
  SET_STRING_ELT(ansNames, 3, mkChar("Ps"));
  SET_STRING_ELT(ansNames, 4, mkChar("Kf"));
  SET_STRING_ELT(ansNames, 5, mkChar("negLogLikelihood"));
  SET_NAMES(ans, ansNames);

  // free everything
  gsl_vector_free(x);
  gsl_vector_free(ss);
  gsl_multimin_fminimizer_free(s);

  UNPROTECT(2);
  return ans;
}


