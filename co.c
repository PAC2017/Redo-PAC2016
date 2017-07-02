#include <stdlib.h>
#include <omp.h>
#include <math.h>

#define Malloc(type,n)  (type *)malloc((n)*sizeof(type)) 
#define TRUE 1
#define FALSE 0

void findcov_(int *NVAR, int *NROW, int *NV, double MType[*NROW][*NVAR], double COV[*NV][*NVAR]){
    //In Fortran is MType[NVAR][NROW]
    //In Fortran is COV[NVAR][NV]
    int i, j;
    int I2;
    double Covariance(int *NObs, int *nvar, int c1, int c2, double MType[*NObs][*nvar]);
    
    printf("IN FINDCOV #1\n");
    
    #pragma omp parallel for private(i, j, I2) shared(NVAR, COV, MType)
    for(j = 0; j < *NVAR; ++ j){
        for(i = 0; i < *NV; ++ i){
            I2 = (i+j+1)%(*NVAR);
            COV[i][j] = Covariance(NROW, NVAR, j, I2, MType);

            if(i == 10 && j == 221){
                printf("COV[222][11] == %f\n", COV[10][221]);
                printf("I2 == %d\n", I2);
            }
        }
    }

    printf("IN FINDCOV #2\n");

    return;
}

double Covariance(int *NObs, int* nvar, int c1, int c2, double MType[*NObs][*nvar]){
    int IObs, NumMissing = 0;
    double result = 0.0;
    double MeanX1 = 0.0, MeanX2 = 0.0;
    char IsMissingPheno(double*);

    for(IObs = 0; IObs < *NObs; ++ IObs){
        if(IsMissingPheno(&(MType[IObs][c1])) == TRUE || IsMissingPheno(&(MType[IObs][c2])) == TRUE)
            ++ NumMissing;
        else{
            MeanX1 += MType[IObs][c1];
            MeanX2 += MType[IObs][c2];
            result += MType[IObs][c1]*MType[IObs][c2];
        }
    }

    MeanX1 = MeanX1/(*NObs-NumMissing);
    MeanX2 = MeanX2/(*NObs-NumMissing);

    if((*NObs - NumMissing) <= 0)
        result = 0.0;
    else
        result = result/(*NObs-NumMissing) - MeanX1*MeanX2;
    
    return result;
}

char IsMissingPheno(double *Phenotype){
    static double MissingPheno = 65.0;
    char result = FALSE;

    if(isnan(*Phenotype))result = TRUE;
    else
        if(fabs(*Phenotype - MissingPheno) <= 0.001)result = TRUE;
    
    return result;
}