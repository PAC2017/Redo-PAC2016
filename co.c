#include <stdlib.h>
#include <omp.h>
#include <math.h>

#define TRUE 1
#define FALSE 0
#define IsMissingPheno(T) !(!isnan(T) && (T < 64.999 || T > 65.001))

void findcov_(int *NVAR, int *NROW, int *NV, double MType[*NROW][*NVAR], double COV[*NV][*NVAR]){
    //In Fortran is MType[NVAR][NROW]
    //In Fortran is COV[NVAR][NV]
    int i, j;
    int I2;
    double Covariance(int *NObs, int *nvar, int c1, int c2, double MType[*NObs][*nvar]);
    
    #pragma omp parallel for private(i, j, I2) shared(NVAR, COV, MType)
    for(j = 0; j < *NVAR; ++ j){
        for(i = 0; i < *NV; ++ i){
            I2 = (i+j+1)%(*NVAR);
            COV[i][j] = Covariance(NROW, NVAR, j, I2, MType);
        }
    }

    return;
}

double Covariance(int *NObs, int* nvar, int c1, int c2, double MType[*NObs][*nvar]){
    int IObs, NumMissing = 0;
    double result = 0.0;
    double MeanX1 = 0.0, MeanX2 = 0.0;
    //char IsMissingPheno(double*);
/*
    for(IObs = 0; IObs < *NObs; ++ IObs){
        if(IsMissingPheno((MType[IObs][c1])) == TRUE || IsMissingPheno((MType[IObs][c2])) == TRUE)
            ++ NumMissing;
        else{
            MeanX1 += MType[IObs][c1];
            MeanX2 += MType[IObs][c2];
            result += MType[IObs][c1]*MType[IObs][c2];
        }
    }
*/

    #pragma omp parallel shared(MType, NumMissing, result, MeanX1, MeanX2)
    {
        double t1, t2;
        #pragma omp for reduction(+:MeanX1) reduction(+:MeanX2) reduction(+:result)
        for(IObs = 0; IObs < *NObs; ++ IObs){
            t1 = (MType[IObs][c1]);
            t2 = (MType[IObs][c2]);
            if(IsMissingPheno(t1) == TRUE || IsMissingPheno(t2) == TRUE)
                ++ NumMissing;
            else{
                MeanX1 += t1;
                MeanX2 += t2;
                result += t1*t2;
            }
        }
    }

    MeanX1 = MeanX1/(*NObs-NumMissing);
    MeanX2 = MeanX2/(*NObs-NumMissing);

    if((*NObs <= NumMissing))
        result = 0.0;
    else
        result = result/(*NObs-NumMissing) - MeanX1*MeanX2;
    
    return result;
}


// char IsMissingPheno(double *Phenotype){
//     static double MissingPheno = 65.0;
//     char result = FALSE;

//     if(isnan(*Phenotype))result = TRUE;
//     else
//         if(fabs(*Phenotype - MissingPheno) <= 0.001)result = TRUE;
    
//     return result;
// }