#include <stdlib.h>
#include <omp.h>
#include <math.h>

#define TRUE 1
#define FALSE 0
//#define isnan(T) (T == 0.0/0.0)?(1):(0)
#define IsMissingPheno(T) !(!isnan(T) && (T < 64.999 || T > 65.001))
#define Malloc(type,n)  (type *)malloc((n)*sizeof(type)) 
#define THREAD_NUM 32

struct M{
    int n;
    double m1;
    double m2;
    double m3;
};

void findcov_(int *NVAR, int *NROW, int *NV, double MType[*NROW][*NVAR], double COV[*NV][*NVAR]){
    //In Fortran is MType[NVAR][NROW]
    //In Fortran is COV[NVAR][NV]
    int i, j, k;

    //Malloc space for Mean
    struct M **Mean;
    struct M *elements;
    elements = Malloc(struct M, (*NVAR)*(*NV));
    Mean = Malloc(struct M*, (*NVAR));
    omp_set_num_threads(THREAD_NUM);
    #pragma omp parallel for
    for(i = 0, j = 0; i < *NVAR; ++i, j += (*NV)){Mean[i] = &elements[j];}
    //Mean[j][I2]

    printf("INNER TAG #1\n");

    //init Mean[:][:]
    omp_set_num_threads(THREAD_NUM);
    #pragma omp parallel for
    for(i = 0; i < *NVAR; ++ i){
        for(j = 0; j < *NV; ++ j){
            Mean[i][j].n = 0;
            Mean[i][j].m1 = 0.0;
            Mean[i][j].m2 = 0.0;
            Mean[i][j].m3 = 0.0;
        }
    }

    printf("INNER TAG #2\n");

    //cal raw cov and Mean(mean and missing num)
    double t1, t2;
    int I2;
    omp_set_num_threads(THREAD_NUM);
    #pragma omp parallel for private(t1, t2, I2, i, j, k) schedule(static)
    for(j = 0; j < (*NVAR); ++ j){
        for(k = 0; k < (*NROW); ++ k){
            t1 = MType[k][j];
            for(i = 0; i < (*NV); ++ i){
                I2 = (i+j+1)%(*NVAR);
                t2 = MType[k][I2];
            
                if(IsMissingPheno(t1) || IsMissingPheno(t2)){
                    (Mean[j][i].n) += 1;
                }else{
                    (Mean[j][i].m1) += t1;
                    (Mean[j][i].m2) += t2;
                    (Mean[j][i].m3) += t1*t2;
                }
            }
        }
    }

    printf("INNER TAG #3\n");

    //cal cov
    double mm1, mm2;
    int unmissing;
    omp_set_num_threads(THREAD_NUM);
    #pragma omp parallel for private(mm1, mm2, unmissing, i, j)
    for(i = 0; i < *NV; ++ i){
        for(j = 0; j < *NVAR; ++ j){
            unmissing = *NROW - (Mean[j][i].n);
            if(unmissing == 0){
                COV[i][j] = 0.0;
            }else{
                mm1 = (Mean[j][i].m1) / unmissing;
                mm2 = (Mean[j][i].m2) / unmissing;
                COV[i][j] = (Mean[j][i].m3) / unmissing - mm1*mm2;
            }
        }
    }

    printf("########### COV[222][11] == %f\n", COV[10][221]);

    return;
}
