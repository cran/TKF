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

void dynprog(double *Mu, double *distance, int *SA, int *SB, double *AvgLen, int *sequenceA, int *sequenceB, double *SubMat, double *equi, double *result)
{
    //double Lambda = ((double) *SA + (double) *SB)/((double) *SA + (double) *SB + 2.0) * *Mu;
    double Lambda = *AvgLen / (*AvgLen + 1.0) * *Mu;
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
    // test the output of submatrix
    //for(i=0;i<20;i++){
    //    for(j=0;j<20;j++){
    //        Rprintf("%f ", SubMatrix[i][j]);
    //    }
    //    Rprintf("\n");
    //}
    // test the PI
    //for(i=0;i<20;i++){
    //    Rprintf("%f ", equi[i]);
    //}
    // initialize the entries tables, only log-likelihood is stored in thie table
    double L0[*SA+1][*SB+1];
    double L1[*SA+1][*SB+1];
    double L2[*SA+1][*SB+1];

    // initialize the boundary conditions
    double temp;
    L0[0][0] = -INFINITY;
    L2[0][0] = -INFINITY;
    L1[0][0] = log(1.0 - Lambda / *Mu) + log(P12t);
    
    for(i=1;i<(*SA+1);i++){
        temp = 0;
        for(j=1;j<=i;j++){
            temp = temp + log(equi[sequenceA[j-1] - 1]) + log(P01t);
        }
        L0[i][0] = log(1.0 - Lambda / *Mu) + i * log(Lambda / *Mu) + log(P12t) + temp;
        L1[i][0] = -INFINITY;
        L2[i][0] = -INFINITY;
    }

    for(j=1;j<(*SB+1);j++){
        temp = 0;
        for(i=1;i<=j;i++){
            temp = temp + log(equi[sequenceB[i-1]-1]);
        }
        L2[0][j] = log(1.0 - Lambda / *Mu) + log(P12t) + j * (log(Lambda) + log(Beta)) + temp;
        L1[0][j] = -INFINITY;
        L0[0][j] = -INFINITY;
    }
    //Rprintf("%f", max(4,max(L0[0][0], L1[0][0])));
    // recursive iteration
    for(i=1;i<(*SA+1);i++){
        for(j=1;j<(*SB+1);j++){
            temp = max(max(L0[i-1][j], L1[i-1][j]), L2[i-1][j]);
            //Rprintf("%f\n", temp);
            L0[i][j] = log(Lambda) - log(*Mu) + log(equi[sequenceA[i-1]-1]) + log(P01t) + temp + log(exp(L0[i-1][j] - temp) + exp(L1[i-1][j] - temp) + exp(L2[i-1][j] - temp));
            temp = max(max(L0[i-1][j-1], L1[i-1][j-1]), L2[i-1][j-1]);
            L1[i][j] = log(Lambda) - log(*Mu) + log(equi[sequenceA[i-1]-1]) + log(SubMatrix[sequenceA[i-1]-1][sequenceB[j-1]-1] * P1t + equi[sequenceB[j-1]-1] * P11t) + temp + log(exp(L0[i-1][j-1] - temp) + exp(L1[i-1][j-1] - temp) + exp(L2[i-1][j-1] - temp));
            if(isfinite(L1[i][j-1]) || isfinite(L2[i][j-1])){
                temp = max(L1[i][j-1], L2[i][j-1]);
                L2[i][j] = log(equi[sequenceB[j-1]-1]) + log(Lambda) + log(Beta) + temp + log(exp(L1[i][j-1] - temp) + exp(L2[i][j-1] - temp));
            }
            else
                L2[i][j] = -INFINITY;
        }
    }
    temp = max(max(L0[*SA][*SB], L1[*SA][*SB]), L2[*SA][*SB]);
    *result = -(temp + log(exp(L0[*SA][*SB] - temp) + exp(L1[*SA][*SB] - temp) + exp(L2[*SA][*SB] - temp)));

}

