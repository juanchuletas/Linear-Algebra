#include "ludecomp_gauss.hpp"

namespace LinearAlgebra{

    int GaussLU::getMaxValueAt(int index, double *matIn, int N){

        double maxValue = matIn[index*N + index]; 
        int maxIndx = index; 
        for(int j=index+1; j<N; j++){
            double target = matIn[j*N + index];
            if(maxValue<target){
                maxValue = target;
                maxIndx = j; //Gets the max index value
            } 
        }

        return maxIndx;
    }
    void GaussLU::swapRows(int row1, int row2, double *matrix ,int N){

        double *aux = new double[N]; 
        for(int i=0; i<N; i++){
            aux[i] = matrix[row1*N + i]; 
            matrix[row1*N + i] = matrix[row2*N + i]; 
            matrix[row2*N + i] = aux[i]; 
        }

        delete [] aux; 
    }
    void GaussLU::decomp(double *matIn, int matSize){


        for(int i=0; i<matSize-1; i++){//Main loop 

            //First step: Find the highest value of this column "i"
            int target = getMaxValueAt(i, matIn, matSize);
            if(target!=i){
                swapRows(i,target, matIn, matSize);
            } 

            for(int j=i+1; j<matSize; j++){
                double pivot = matIn[j*matSize + i]/matIn[i*matSize+i];
                matIn[j*matSize + i] = pivot; 

            }




        }

    }



}
