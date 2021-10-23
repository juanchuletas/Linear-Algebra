#if !defined(_LU_DECOMP)
#define _LU_DECOMP

#include<iostream>
#include<cmath>

namespace LinearAlgebra{

    class LU{
        public:
            static void ludecomp( double *MatA, double *xvec, double *bvec,int M, int N);
    };



}


#endif // _LU_DECOMP
