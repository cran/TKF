#include <R.h>
#include <Rmath.h>
#include <stdio.h>
#include <math.h>

#ifndef max
    #define max( a, b ) ( ((a) > (b)) ? (a) : (b) )
#endif

#ifndef min
    #define min( a, b ) ( ((a) < (b)) ? (a) : (b) )
#endif

void dynprog(double *Mu, double *distance, double *r, int *SA, int *SB, double *AvgLen, int *sequenceA, int *sequenceB, double *SubMat, double *equi, double *result)
{
    double expLen = *AvgLen;
    double Lambda = *Mu * (expLen * (1.0 - *r) - 2.0 * *r + sqrt(expLen * expLen * (1.0 - *r) * (1.0 - *r) + 4.0 * *r)) / (2.0 * (expLen + 1.0) * (1.0 - *r));
    double Beta = (1.0 - exp((Lambda - *Mu) * *distance)) / (*Mu - Lambda * exp((Lambda - *Mu) * *distance));
    double P1t = exp(- *Mu * *distance) * (1.0 - Lambda * Beta);
    double P11t = (1.0 - exp(- *Mu * *distance) - *Mu * Beta) * (1.0 - Lambda * Beta);
    double P12t = 1.0 - Lambda * Beta;
    double P01t = *Mu * Beta;
    // construct the Substitution Matrix from the *SubMat
    double SubMatrix[20][20];
    int k = 0;
    int i = 0;
    int j = 0;
    for(j=0;j<20;j++){
        for(i=0;i<20;i++){
            SubMatrix[i][j] = SubMat[k];
            k++;
        }
    }
    //check the parameters
    //Rprintf("expLen:%f\n", expLen);
    //Rprintf("Lambda:%f\n", Lambda);
    //Rprintf("Beta:%f\n", Beta);
    //Rprintf("P1t:%f\n", P1t);
    //Rprintf("P11t:%f\n", P11t);
    //Rprintf("P12t:%f\n", P12t);
    // initialize the entries tables, only log-likelihood is stored in thie table
    double L1[*SA+1][*SB+1];
    double L2[*SA+1][*SB+1];
    double L3[*SA+1][*SB+1];
    double L4[*SA+1][*SB+1];
    double L5[*SA+1][*SB+1];
    double L6[*SA+1][*SB+1];

    // initialize the boundary conditions
    double temp;
    L1[0][0] = log(1.0 - Lambda / *Mu) + log(P12t);
    L2[0][0] = -INFINITY;
    L3[0][0] = -INFINITY;
    L4[0][0] = -INFINITY;
    L5[0][0] = -INFINITY;
    L6[0][0] = -INFINITY;
    
    L2[1][0] = log((1.0 - Lambda / *Mu) * Lambda / *Mu * (1.0 - *r)) + log(P12t) + log(equi[sequenceA[0]-1]) + log(P01t);
    L1[1][0] = -INFINITY;
    L3[1][0] = -INFINITY;
    L4[1][0] = -INFINITY;
    L5[1][0] = -INFINITY;
    L6[1][0] = -INFINITY;
    for(i=2;i<(*SA+1);i++){
        temp = 0;
        for(j=2;j<=i;j++){
            temp = temp + log(equi[sequenceA[j-1]-1]) + log(*r + P01t * (1.0 - *r) * Lambda / *Mu);
        }
        L2[i][0] = log((1.0 - Lambda / *Mu) * Lambda / *Mu * (1.0 - *r)) + log(P12t) + log(equi[sequenceA[0]-1]) + log(P01t) + temp;
        L1[i][0] = -INFINITY;
        L3[i][0] = -INFINITY;
        L4[i][0] = -INFINITY;
        L5[i][0] = -INFINITY;
        L6[i][0] = -INFINITY;
    }
    // check the boundary conditions
    //Rprintf("L1[0][0]:%f\n", L1[0][0]);
    //for(i=0;i<= *SA;i++){
    //    Rprintf("%f\n",L6[i][0]);
    //}
    

    L3[0][1] = log(1.0 - Lambda / *Mu) + log(1.0 - Lambda * Beta) + log(Lambda * Beta) + log(1.0 - *r) + log(equi[sequenceB[0]-1]);
    L1[0][1] = -INFINITY;
    L2[0][1] = -INFINITY;
    L4[0][1] = -INFINITY;
    L5[0][1] = -INFINITY;
    L6[0][1] = -INFINITY;
    for(j=2;j<(*SB+1);j++){
        temp = 0;
        for(i=2;i<=j;i++){
            temp = temp + log(equi[sequenceB[i-1]-1]) + log(*r + Lambda * Beta * (1.0 - *r));
        }
        L3[0][j] = log(1.0 - Lambda / *Mu) + log(1.0 - Lambda * Beta) + log(Lambda * Beta) + log(1.0 - *r) + log(equi[sequenceB[0]-1]) + temp;
        L1[0][j] = -INFINITY;
        L2[0][j] = -INFINITY;
        L4[0][j] = -INFINITY;
        L5[0][j] = -INFINITY;
        L6[0][j] = -INFINITY;
    }
    // check the L4
    //for(j=0;j<= *SB;j++){
    //    Rprintf("%f\n",L2[0][j]);
    //}
    //double test = exp(L1[0][0] - (-INFINITY));
    //Rprintf("test: %f\n",test);
    double Km;
    double Kn;
    for(i=1;i<(*SA+1);i++){
        for(j=1;j<(*SB+1);j++){
    //for(i=1;i<2;i++){
    //    for(j=1;j<2;j++){
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
            temp = max(max(max(max(max(L1[i-1][j-1], L2[i-1][j-1]), L3[i-1][j-1]), L4[i-1][j-1]), L5[i-1][j-1]), L6[i-1][j-1]);
            L1[i][j] = log(equi[sequenceA[i-1]-1]) + log(SubMatrix[sequenceA[i-1]-1][sequenceB[j-1]-1]) + temp + log(*r * Km * Kn * exp(L1[i-1][j-1] - temp) + P1t * (1.0 - *r) * Lambda / *Mu * (exp(L1[i-1][j-1] - temp) + exp(L2[i-1][j-1] - temp) + exp(L3[i-1][j-1] - temp) + exp(L4[i-1][j-1] - temp) + exp(L5[i-1][j-1] - temp) + exp(L6[i-1][j-1] - temp)));
            //Rprintf("L1: %f\n",L1[i][j]);
            
            temp = max(max(max(max(max(L1[i-1][j], L2[i-1][j]), L3[i-1][j]), L4[i-1][j]), L5[i-1][j]), L6[i-1][j]);
            L2[i][j] = log(equi[sequenceA[i-1]-1]) + temp + log(*r * Km * exp(L2[i-1][j] - temp) + P01t * (1.0 - *r) * Lambda / *Mu * (exp(L1[i-1][j] - temp) + exp(L2[i-1][j] - temp) + exp(L3[i-1][j] - temp) +exp(L4[i-1][j] - temp) + exp(L5[i-1][j] - temp) + exp(L6[i-1][j] - temp)));
            //Rprintf("L2: %f\n",L2[i][j]);
            
            if(isfinite(L1[i][j-1]) || isfinite(L3[i][j-1]) || isfinite(L4[i][j-1]) || isfinite(L5[i][j-1]) || isfinite(L6[i][j-1])){
                temp = max(max(max(max(L1[i][j-1], L3[i][j-1]), L4[i][j-1]), L5[i][j-1]), L6[i][j-1]);
                L3[i][j] = log(equi[sequenceB[j-1]-1]) + temp + log(*r * Kn * exp(L3[i][j-1] - temp) + Lambda * Beta * (1.0 - *r) * (exp(L1[i][j-1] - temp) + exp(L3[i][j-1] - temp) + exp(L4[i][j-1] - temp) + exp(L5[i][j-1] - temp) + exp(L6[i][j-1] - temp)));
            }else{
                L3[i][j] = -INFINITY;
            }
            //Rprintf("L3: %f\n",L3[i][j]);
            
            temp = max(max(max(max(max(L1[i-1][j-1], L2[i-1][j-1]), L3[i-1][j-1]), L4[i-1][j-1]), L5[i-1][j-1]), L6[i-1][j-1]);
            L4[i][j] = log(equi[sequenceA[i-1]-1]) + log(equi[sequenceB[j-1]-1]) + temp + log(*r * *r * Km * Kn * exp(L4[i-1][j-1] - temp) + P11t * (1.0 - *r) * (1.0 - *r) * Lambda / *Mu * (exp(L1[i-1][j-1] - temp) + exp(L2[i-1][j-1] - temp) + exp(L3[i-1][j-1] - temp) + exp(L4[i-1][j-1] - temp) + exp(L5[i-1][j-1] - temp) + exp(L6[i-1][j-1] - temp)));
            //Rprintf("L4: %f\n",L4[i][j]);
            
            if(isfinite(L4[i-1][j]) || isfinite(L5[i-1][j])){
                temp = max(L4[i-1][j], L5[i-1][j]);
                L5[i][j] = log(equi[sequenceA[i-1]-1]) + log(Km * Kn) + log(*r) + temp + log(exp(L4[i-1][j] - temp) + exp(L5[i-1][j] - temp));
            }else{
                L5[i][j] = -INFINITY;
            }
            //Rprintf("L5: %f\n",L5[i][j]);
            
            if(isfinite(L4[i][j-1]) || isfinite(L6[i][j-1])){
                temp = max(L4[i][j-1], L6[i][j-1]);
                L6[i][j] = log(equi[sequenceB[j-1]-1]) + log(Km * Kn) + log(*r) + temp + log(exp(L4[i][j-1] - temp) + exp(L6[i][j-1] - temp));
            }else{
                L6[i][j] = -INFINITY;
            }
            //Rprintf("L6: %f\n",L6[i][j]);
        }
    }
    // test the output
    //for(i=0;i<=*SB;i++){
    //    Rprintf("L3: %f\n",L3[1][i]);
    //}

    temp = max(max(max(max(max(L6[*SA][*SB], L1[*SA][*SB]), L2[*SA][*SB]), L3[*SA][*SB]), L4[*SA][*SB]), L5[*SA][*SB]);
    *result = -(temp + log(exp(L1[*SA][*SB] - temp) + exp(L2[*SA][*SB] - temp) + exp(L3[*SA][*SB] - temp) + exp(L4[*SA][*SB] - temp) + exp(L5[*SA][*SB] - temp) + exp(L6[*SA][*SB] - temp)));

}

