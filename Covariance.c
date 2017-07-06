#include <omp.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#define num_threads 272
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
        for(; k < end ; k++)
        {
            for(j = 0; j < *NV ; j++)
            {
                for(i = 0 ; i < *NVAR ; i++)
                {
                    //I2 = (i + j + 1) % *NVAR;
                    I2 = (i + j + 1) >= *NVAR ? i + j +1 - *NVAR : i + j + 1;
                    x = MType[k][i];
                    y = MType[k][I2];
                    if( IsMissingPheno2(x,y) )
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
    /*
    int I2;
    double meanxi,meanxI2,Covar,x,y;
    #pragma omp parallel for
    for(i = 0 ; i < *NVAR ; i++ )
    {
        for(j = 0 ; j < *NV ; j++)
        {
            I2 = (i + j + 1) % *NVAR;
            Covar = 0.0;
            meanxi = 0.0;
            meanxI2 = 0.0;
            leavenum = *NROW;
            for(k = 0;k < *NROW ; k++)
            {
                x = MType[k][i];
                y = MType[k][I2];
                if( IsMissingPheno2(x,y) )
                {
                    leavenum --;
                }
                else
                {
                    meanxi += x;
                    meanxI2 += y;
                    Covar += x*y;
                }
            }
            if(leavenum <= 0)
            {
                Covar = 0;
            }
            else
            {
                if(i==0&j==0)
                {
                    printf("%d;%lf;%lf;%lf\n",leavenum,meanxi,meanxI2,Covar);
                }
                Covar = Covar/leavenum - (meanxi *  meanxI2)/ (leavenum * leavenum);
            }
            COV[j][i] = Covar;
        }
    }*/

}
