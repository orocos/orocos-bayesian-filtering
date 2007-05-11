#include "matrix_wrapper.h"
#include <math.h>

namespace MatrixWrapper
{

    bool
    SymmetricMatrix_Wrapper::cholesky_semidefinite(MyMatrix& a) const
    {
        // apply cholesky to *this
        // put result in a
        // result is lower triangular matrix
        a = (*(MySymmetricMatrix*)this);
        int sz = a.rows(); 
          for (int k=1; k<sz+1; ++k) {
           if (a(k,k) < std::numeric_limits<double>::epsilon()) {
                 std::cout<< "Warning: matrix of which cholesky decomposition is asked, is negative definite!: returning zero matrix" << std::endl;
                 a = 0.0; return false;//matrix is negative definite
           } 
           else  a(k,k)=sqrt(a(k,k));
           for (int i=k+1; i<sz+1; ++i) {
               if (a(k,k)< std::numeric_limits<double>::epsilon()){
                   a(i,k) = 0.0; //set to zero, matrix is semidefinite
               }
               else a(i,k)/=a(k,k);
           }
           for (int j=k+1; j<sz+1; ++j) {
             for (int i=j; i<sz+1; ++i)  a(i,j)-=a(i,k)*a(j,k);
           }
          }
          //delete upper  triangle
          for (int i=1; i<sz+1; i++) {
            for (int j=i+1; j<sz+1; j++) a(i,j)=0.0;
          }
        return true;
    }
}
