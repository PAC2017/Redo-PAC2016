#include <stdlib.h>
#include <omp.h>

#define TRUE 1
#define FALSE 0
#define isnan(T) (T == 0.0/0.0)?(1):(0)
#define IsMissingPheno(T) !(!isnan(T) && (T < 64.999 || T > 65.001))
#define Malloc(type,n)  (type *)malloc((n)*sizeof(type)) 

void findcov_(int *NVAR, int *NROW, int *NV, double MType[*NROW][*NVAR], double COV[*NV][*NVAR]){
    //In Fortran is MType[NVAR][NROW]
    //In Fortran is COV[NVAR][NV]
    int i, j, I2, k;

    //Malloc space for Mean, 8 is for the size of cache(flase sharing)
    double **Mean;
    double *elements;
    elements = Malloc(double, (*NVAR)*8);
    Mean = Malloc(double*, (*NVAR));
    #pragma omp parallel for
    for(i = 0, j = 0; i < *NVAR; ++i, j += 8){Mean[i] = &elements[j]; Mean[i][0] = 0.0; Mean[i][1] = 0.0; Mean[i][2] = 0.0;}
    //Mean[c][0] intend if have mean or not
    //Mean[c][1] is for mean of MType[:][c]
    //Mean[c][2] is for missing num of MType[:][c]

    //Malloc a cov buffer
    double *T;
    int s_T = (*NV <= 8)?(8):(*NV);
    T = Malloc(double, s_T);
    //init cov buffer
    for(i = 0; i < s_T; ++ i)T[i] = 0.0;

    //cal raw cov and mean and missing num
    double t1, t2;
    int itr, v_num;
    #pragma omp parallel for
    for(j = 0; j < *NVAR; ++ j){
        for(k = 0; k < *NROW; ++ k){
            t1 = MType[k][j];
            if(IsMissingPheno(t1)){
                #pragma omp atomic
                Mean[j][2] += 1;
                continue;
            }else{
                #pragma omp atomic
                Mean[j][1] += t1;
            }

            for(i = 0; i < *NV; ++ i){
                I2 = (i+j+1)%(*NVAR);
                if(I2 >= *NVAR)printf("##################### %d\n", I2);
                if(k >= *NROW)printf("$$$$$$$$$$$$$$$$$$$$$ %d\n", k);
                t2 = MType[k][I2];
                
                if(IsMissingPheno(t2)){
                    continue;
                }
                #pragma omp atomic
                T[i] += t1*t2;
            }
        }
        //raw cov here
        for(itr = 0; itr < *NV; ++ itr){
            COV[itr][j] = T[itr];
            T[itr] = 0.0;
        }
    } 

    //cal cov
    double m1, m2;
    int n1, n2;
    int I3, missing;
    #pragma omp parallel for
    for(i = 0; i < *NV; ++ i){
        for(j = 0; j < *NVAR; ++ j){
            I3 = (i+j+1)%(*NVAR); 
            n1 = (int)Mean[j][2]; n2 = (int)Mean[I3][2];
            missing = (n1 >= n2)?(n1):(n2);
            
            if(missing == *NROW){
                COV[i][j] = 0.0;
                continue;
            }else{
                m1 = Mean[j][1] / missing;
                m2 = Mean[I3][1] / missing;
                #pragma omp atomic
                COV[i][j] = COV[i][j]/missing - m1*m2;
            }
        }
    }

    return;
}
