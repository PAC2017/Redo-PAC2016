#include <stdlib.h>
#include <omp.h>
#include <math.h>

#define Malloc(type,n)  (type *)malloc((n)*sizeof(type)) 
#define TRUE 1
#define FALSE 0

void findcov_(int *NVAR, int *NROW, int *NV, double MType[*NROW][*NVAR], double ** COV){
    //In Fortran is MType[NVAR][NROW]
    int i, j;
    int I2;
    double Covariance(int *NObs, int *nvar, int c1, int c2, double MType[*NObs][*nvar]);
    /*
    double *elements;
    elements = Malloc(double, NVAR*NROW);
    MType = Malloc(double *, NROW);
    for(i = 0, j = 0; i < NROW; ++ i, j += NVAR){MType[i] = &elements[j];}
    */

    printf("IN FINDCOV #1\n");

    //TODO: BETTER METHOD
    // double tmp;
    // int nrow_half = NROW/2, nvar_half = NVAR/2;
    // for(i = 0; i <= nrow_half; ++ i){
    //     for(j = 0; j <= nvar_half; ++ j){
    //         tmp = MType[i][j];
    //         MType[i][j] = MType[NROW-i-1][NVAR-j-1];
    //         MType[NROW-i-1][NVAR-j-1] = tmp;
    //     }
    // }

    double *elements2;
    elements2 = Malloc(double, (*NVAR)*(*NV));
    COV = Malloc(double *, *NV);
    for(i = 0, j = 0; i < *NV; ++ i, j += *NVAR){COV[i] = &elements2[j];}

    printf("IN FINDCOV #3\n");

    double tmp_cov;
    #pragma omp parallel
    {
        #pragma omp for
        for(i = 0; i < *NV; ++ i){
            for(j = 0; j < *NVAR; ++ j){
                I2 = (i+j-1)%(*NVAR) + 1;
                tmp_cov = Covariance(NROW, NVAR, i, I2, MType);
                COV[i][j] = tmp_cov;
            }
        }
    }

    printf("IN FINDCOV #4\n");

    free(elements2);
    free(COV);

    printf("IN FINDCOV #5\n");

    return;
}

double Covariance(int *NObs, int* nvar, int c1, int c2, double MType[*NObs][*nvar]){
    //X1[?] == MType[?][c1]
    //X2[?] == MType[?][c2]
    int IObs, NumMissing = 0;
    double result = 0.0;
    double MeanX1 = 0.0, MeanX2 = 0.0;
    char IsMissingPheno(double*);

    for(IObs = 0; IObs < *NObs; ++ IObs){
        //printf("IObs == %d\n", IObs);
        //printf("X1[IObs] == %f\n", X1[IObs]);
        if(IsMissingPheno(&(MType[IObs][c1])) || IsMissingPheno(&(MType[IObs][c2])))
            ++ NumMissing;
        else{
            MeanX1 += MType[IObs][c1];
            MeanX2 += MType[IObs][c2];
            result += MType[IObs][c1]*MType[IObs][c2];
        }
    }

    MeanX1 /= (*NObs-NumMissing);
    MeanX2 /= (*NObs-NumMissing);

    if(*NObs - NumMissing <= 0)
        result = 0.0;
    else
        result =result/(*NObs-NumMissing) - MeanX1*MeanX2;
    
    return result;
}

char IsMissingPheno(double *Phenotype){
    static double MissingPheno = 65.0;
    char result = FALSE;

    //printf("$$$$ Phenotype == %f\n", *Phenotype);

    if(isnan(*Phenotype))result = TRUE;
    else
        if(abs(*Phenotype - MissingPheno) <= 0.001)result = TRUE;
    
    return result;
}