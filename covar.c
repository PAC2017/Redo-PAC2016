#include <omp.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#define SSIZE 21
#define num_threads 32
#define IsMissingPheno(T) (isnan(T) || (T >= 64.999 && T <= 65.001))
#define IsMissingPheno2(p1, p2) (isnan(p1) || isnan(p2) || p1 <= 65.001 || p2 <= 65.001)
void findcov_(int *NVAR,int *NROW,int *NV, double MType[*NROW][*NVAR], double COV[*NV][*NVAR])
{
    double * tmp[num_threads];
    int * missnum[num_threads];
    omp_set_num_threads(num_threads);
    #pragma omp parallel shared(NVAR,NROW,COV,MType,tmp,missnum)
    {
        int ID = omp_get_thread_num();
        tmp[ID] = calloc(*NV * *NVAR * 3,sizeof(double));
        missnum[ID] = calloc(*NV * *NVAR,sizeof(int));
        int per = (*NROW + num_threads - 1) / num_threads;
        int k = ID * per,i,j,I2;
        double x,y;
        int end = k + per < *NROW ? k + per : *NROW ;

        int itr;
        double sbf[SSIZE];
        char csbf[SSIZE];
        double tmp_read;
        int w1, w2;
        int nvar = *NVAR;
        int nrow = *NROW;
        int nv = *NV;
        char x_, y_;
        int I2_;

        for(; k < end ; k++)
        {
            for(j = 0; j < *NV ; j++)
            {
                if(i == 0){
                    for(itr = 0; itr < SSIZE; ++ itr){
                        tmp_read = MType[k][itr];
                        sbf[itr] = tmp_read;
                        csbf[itr] = IsMissingPheno(tmp_read);
                        w1 = 0;
                        w2 = SSIZE;
                    }
                }else{
                    tmp_read = MType[k][w2];
                    sbf[w1] = tmp_read;
                    csbf[w1] = IsMissingPheno(tmp_read);
                    ++w1;
                    ++w2;
                    if(w1 == SSIZE)w1 = 0;
                    if(w2 == nvar)w2 = 0;
                }

                x = sbf[w1];
                x_ = csbf[w1];

                for(i = 0 ; i < *NVAR ; i++)
                {
                    //I2 = (i + j + 1) % *NVAR;
                    //I2 = (i + j + 1) >= *NVAR ? i + j +1 - *NVAR : i + j + 1;
                    I2_ = w1+j+1;
                    I2 = (I2_ >= SSIZE)?(I2_-SSIZE):(I2_);
                    //x = MType[k][i];
                    //y = MType[k][I2];
                    y = sbf[I2];
                    y_ = csbf[I2];

                    if( x_ || y_ )
                    {
                        (*(missnum[ID] + j * (*NVAR) + i) )++;
                    }
                    else
                    {
                        (*(tmp[ID] + j * (*NVAR) * 3 + i * 3) ) += x;
                        (*(tmp[ID] + j * (*NVAR) * 3 + i * 3 + 1) ) += y;
                        (*(tmp[ID] + j * (*NVAR) * 3 + i * 3 + 2) ) += x*y;
                    }
                }
            }
        }
    }
    int j;
    for(j = 0; j < *NV ; j++)
    {
        #pragma omp parallel shared(j,COV,tmp,missnum,NVAR,NROW)
        {
            int ID = omp_get_thread_num();
            int per = (*NVAR + num_threads - 1)/ num_threads;
            int i = ID * per,k,leavenum;
            int end = i + per < *NVAR ? k=i + per : *NVAR ;
            double x,y,z;
            for( ; i < end ; i++)
            {
                for(k = 1; k < num_threads ; k++)
                {
                    (*(tmp[0] + j * (*NVAR) * 3 + i * 3) ) += (*(tmp[k] + j * (*NVAR) * 3 + i * 3) );
                    (*(tmp[0] + j * (*NVAR) * 3 + i * 3 + 1) ) += (*(tmp[k] + j * (*NVAR) * 3 + i * 3 + 1) );
                    (*(tmp[0] + j * (*NVAR) * 3 + i * 3 + 2) ) += (*(tmp[k] + j * (*NVAR) * 3 + i * 3 + 2) );
                    (*(missnum[0] + j * (*NVAR) + i) ) += (*(missnum[k] + j * (*NVAR) + i) );
                }
                leavenum = *NROW - (*(missnum[0] + j * (*NVAR) + i) );
                x = (*(tmp[0] + j * (*NVAR) * 3 + i * 3) );
                y = (*(tmp[0] + j * (*NVAR) * 3 + i * 3 + 1) );
                z = (*(tmp[0] + j * (*NVAR) * 3 + i * 3 + 2) );
                if(leavenum <= 0)
                {
                    COV[j][i] = 0;
                }
                else
                {
                    COV[j][i] = z/leavenum - (x * y)/(leavenum * leavenum);
                }
            }
        }
    }
    int i;
    #pragma omp parallel for
    for(i = 0; i< num_threads; i++ )
    {
        free(missnum[i]);
        free(tmp[i]);
    }

    printf("=========== COV[222][11] == %f\n", COV[10][221]);
}