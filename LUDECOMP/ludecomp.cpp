#include "ludecomp.hpp"
namespace LinearAlgebra{

    void LU::ludecomp(double *matA, double *xvec, double *bvec,int M, int N){
        //double *mat_aux = new double[M*N];
        double l_ik;
        for(int k=0; k<N-1; k++){
            for(int i=k+1; i<N; i++){
                l_ik = matA[i + k*N]/matA[k*N + k];
                for(int j=k; j<N; j++){
                    matA[i*N + j] = matA[i*N + j] - l_ik*matA[k*N + j];
                }
                bvec[i] = bvec[i] - l_ik*bvec[k];
            }
        }
        //Backsust
        bvec[N-1] = bvec[N-1]/matA[N*M-1];
        for(int k=N-2; k>=0; k--){
            for(int j=k+1; j<N; j++){
                bvec[k] = bvec[k] - matA[k*N + j]*bvec[j];
            }
            bvec[k] = bvec[k]/matA[k*N+k];
        }


        //delete [] mat_aux;
    }



}