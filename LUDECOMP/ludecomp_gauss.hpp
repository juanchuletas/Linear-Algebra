#if !defined(_LU_DECOMP_GAUSS_)
#define _LU_DECOMP_GAUSS_

namespace LinearAlgebra{


    class GaussLU{
        public:
            static int getMaxValueAt(int target, double *matIn,int N);
            static void decomp(double *matIn,int N);
            static void swapRows(int row1, int row2, double *matrix ,int N);  

        


    };



}

#endif // _LU_DECOMP_GAUSS_
