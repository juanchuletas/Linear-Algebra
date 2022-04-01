#include<iostream>
#include"ludecomp.hpp"

int main(){

   double matA[25]={2.,1.,1.,3.,2.,
                    1.,2.,2.,1.,1.,
                    1.,2.,9.,1.,5.,
                    3.,1.,1.,7.,1.,
                    2.,1.,5.,1.,8.};

   double LU[25];

   double bvec[5]={-2.,4.,3.,-5.,1.};
   double *x;
   int rows = 5;


    LinearAlgebra::LU::ludecomp(matA,x,bvec,rows,rows);



    /*lucmp.dcmp();
    lucmp.print();
    lucmp.solve(bvec,x);*/
    for(int i=0; i<rows; i++){
        for(int j=0; j<rows; j++){
                //printf("\t%.10lf ",matA[i*rows+ j]);
        }
        printf("\n");
        printf("%.10lf\n",bvec[i]);
    }




   return 0; 
}