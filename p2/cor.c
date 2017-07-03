#include <stdlib.h>
#include <omp.h>

#define TRUE 1
#define FALSE 0
#define isnan(T) (T == 0/0)?(1):(0)
#define IsMissingPheno(T) !(!isnan(T) && (T < 64.999 || T > 65.001))
#define Malloc(type,n)  (type *)malloc((n)*sizeof(type)) 

void findcov_(int *NVAR, int *NROW, int *NV, double MType[*NROW][*NVAR], double COV[*NV][*NVAR]){
    //In Fortran is MType[NVAR][NROW]
    //In Fortran is COV[NVAR][NV]
    int i, j;

    //Malloc space for Mean, 64 is for the size of cache(flase sharing)
    double **Mean;
    double *elements;
    elements = Malloc(double, (*NVAR)*8);
    Mean = Malloc(double*, (*NVAR))
    for(i = 0, j = 0; i < *NVAR; ++i, j += 8){Mean[i] = &elements[j]; Mean[i][0] = 0.0;}
    //Mean[c][0] intend if have mean or not
    //Mean[c][1] is for mean of MType[:][c]
    //Mean[c][2] is for missing num of MType[:][c]

    //Calculate the Covariance
    #pragma omp paralle shared(MType, Mean, COV)
    {
        double t1, t2;
        double t3, t4;
        int I2;
        int i_, j_, k_;

        //i_ == 0 can we calculate all mean
        #pragma omp for shared(MType, Mean, COV, NVAR, NROW) private(i_, j_, k_)
        for(j_= 0; j_ < *NVAR; ++ j_){
                I2 = (j_+1)%(*NVAR);
                t3 = Mean[j_][0]; t4 = Mean[I2][0];
                if(t3 == 0.0 && t4 == 0.0){
                    #pragma omp for shared(MType, Mean, NROW) private(i_, j_, k_)
                    for(k_ = 0; k_ < *NROW; ++ k_){
                        t1 = MType[k_][j_];
                        t2 = MType[k_][I2];
                        if(IsMissingPheno(t1)){
                            ++ Mean[j_][2];
                        }else if(IsMissingPheno(t2)){
                            ++ Mena[I2][2];
                        }else{
                            Mean[j_][1] += t1;
                            Mean[I2][1] += t2;
                            COV[0][j_] += t1*t2;
                        }
                    }
                    Mean[j_][0] = 1;
                    Mean[I2][0] = 1;
                }else if(t3 == 0.0){
                    #pragma omp for shared(MType, Mean, NROW) private(i_, j_, k_)
                    for(k_ = 0; k_ < *NROW; ++ k_){
                        t1 = MType[k_][j_];
                        t2 = MType[k_][I2];
                        if(IsMissingPheno(t1)){
                            ++ Mean[j_][2];
                        }else if(IsMissingPheno(t2)){
                            //++ Mena[I2][2];
                        }else{
                            Mean[j_][1] += t1;
                            //Mean[I2][1] += t2;
                            COV[0][j_] += t1*t2;
                        }
                    }
                    Mean[j_][0] = 1;
                }else if(t4 == 0.0){
                    #pragma omp for shared(MType, Mean, NROW) private(i_, j_, k_)
                    for(k_ = 0; k_ < *NROW; ++ k_){
                        t1 = MType[k_][j_];
                        t2 = MType[k_][I2];
                        if(IsMissingPheno(t1)){
                            //++ Mean[j_][2];
                        }else if(IsMissingPheno(t2)){
                            ++ Mena[I2][2];
                        }else{
                            //Mean[j_][1] += t1;
                            Mean[I2][1] += t2;
                            COV[0][j_] += t1*t2;
                        }
                    }
                    Mean[I2][0] = 1;
                }else{//both are not 0.0
                    #pragma omp for shared(MType, Mean, NROW) private(i_, j_, k_) reduction(+:COV[i_][j_])
                    for(k_ = 0; k_ < *NROW; ++ k_){
                        t1 = MType[k_][j_];
                        t2 = MType[k_][I2];
                        if(IsMissingPheno(t1) || IsMissingPheno(t2)){
                            
                        }else{
                            COV[0][j_] += t1*t2;
                        }
                    }
                }
        }
        //i > 0
        #pragma omp for shared(MType, Mean, COV, NV, NVAR, NROW) private(i_, j_, k_) reduction(+:COV[i_][j_])
        for(i_ = 1; i_ < *NV; ++ i_){
            for(j_= 0; j_ < *NVAR; ++ j_){
                I2 = (i_+j_+1)%(*NVAR);
                for(k_ = 0; k_ < *NROW; ++ k_){
                        t1 = MType[k_][j_];
                        t2 = MType[k_][I2];
                        if(IsMissingPheno(t1) || IsMissingPheno(t2)){

                        }else{
                            COV[i][j_] += t1*t2;
                        }
                }
            }
        }

        int mn;
        int m1, m2;
        double n1, n2;

        //process finnal cov
        #pragma omp for private(mn, m1, m2, n1, n2) shared(Mean, COV, NROW, NV, NVAR)
        for(i_ = 0; i_ < *NV; ++ i_){
            for(j_ = 0; j_ < *NVAR; ++ j_){
                I2 = (i_+j_+1)%(*NVAR);
                m1 = Mean[j_][2];
                m2 = Mean[I2][2];
                mn = (m1 >= m2)?(m1):(m2);
                if(mn == *NROW)
                    COV[i_][j_] = 0.0;
                else{
                    n1 = Mean[j_][1] / (*NROW-m1);
                    n2 = Mean[I2][1] / (*NROW-m2);
                    COV[i_][j_] = COV[i_][j_]/(*NROW - mn) - n1*n2;
                }
            }
        }

    }

    return;
}
