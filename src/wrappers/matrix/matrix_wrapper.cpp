#include "matrix_wrapper.h"
#include <math.h>
#include <vector>

using namespace std;

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
           // check if close to zero => put to zero
           if (a(k,k)< 100.0* std::numeric_limits<double>::epsilon() && a(k,k)> -100.0*std::numeric_limits<double>::epsilon()){
                   a(k,k) = 0.0; //set to zero, matrix is semidefinite
           }
           if (a(k,k) < 0.0) {
                 std::cout<< "Warning: matrix of which cholesky decomposition is asked, is negative definite!: returning zero matrix" << std::endl;
                 std::cout<< "a(" << k << "," << k << ")=" << a(k,k)  << std::endl;
                 std::cout<< "std::numeric_limits<double>::epsilon()=" << std::numeric_limits<double>::epsilon()  << std::endl;
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

    bool
    Matrix_Wrapper::SVD(MyColumnVector& w, MyMatrix& U, MyMatrix& V) const 
    {  
        //get the rows of the matrix
        const int rows=this->rows();
        //get the columns of the matrix
        const int cols=this->columns();

        U = *((MyMatrix*)(this)); 
        
        static const int maxIter=150;
    
        w.resize(cols);
        V.resize(cols,cols,false,true);
        int i(-1),its(-1),j(-1),jj(-1),k(-1),nm(0);
        int ppi(0);
        bool flag;
		double maxarg1, maxarg2;
        double  anorm(0),c(0),f(0),g(0),h(0),s(0),scale(0),x(0),y(0),z(0);
    
        // Householder reduction to bidiagonal form
        std::vector<double> rv1(cols,0.0);
    
        g=scale=anorm=0.0; 
        
        for (i=1; i <= cols; i++){
          ppi=i+1;
          rv1.at(i-1)=scale*g;
          g=s=scale=0.0; 
          if (i <= rows ) {
            // compute the sum of the i-th column, starting from the i-th row
            for(k=i;k<=rows;k++) scale += fabs(U(k,i));
            if (scale) {
              // multiply the i-th column by 1.0/scale, start from the i-th element
              // sum of squares of column i, start from the i-th element
              for (k=i;k<=rows;k++){
                U(k,i)/= scale;
                s += U(k,i) * U(k,i);
              }
              f=U(i,i);  // f is the diag elem
              g=-SIGN(sqrt(s),f);
              h=f*g-s;
              U(i,i)=f-g;
              for (j=ppi; j <= cols; j++) { 
                // dot product of columns i and j, starting from the i-th row
                for (s=0.0,k=i;k<=rows;k++) s += U(k,i) * U(k,j);
                f=s/h;
                // copy the scaled i-th column into the j-th column
                for (k=i;k<=rows;k++) U(k,j) += f*U(k,i);
              }
                for (k=i;k<=rows;k++) U(k,i) *= scale;
            }
          }
          // save singular value
          w(i) = scale*g;
          g=s=scale=0.0;
          if ((i <= rows) && (i != cols)) {
            // sum of row i, start from columns i+1
            for(k=ppi;k<=cols;k++) scale += fabs(U(i,k));
            if (scale) {
              for(k=ppi;k<=cols;k++){
                U(i,k) /= scale;
                s += U(i,k)*U(i,k); 
              }
              f=U(i,ppi);
              g=-SIGN(sqrt(s),f); //<---- do something
              h=f*g-s;
              U(i,ppi)=f-g;
              for ( k=ppi; k <= cols; k++)   rv1.at(k-1)=U(i,k)/h;
              for ( j=ppi; j <= rows; j++) {
                for( s=0.0,k=ppi;k<=cols;k++) s += U(j,k) * U(i,k);
                for ( k=ppi; k <= cols; k++)  U(j,k) += s*rv1.at(k-1);
              }
                for( k=ppi;k<=cols;k++) U(i,k) *= scale;
            }
          }
          maxarg1=anorm;
          maxarg2=(fabs(w(i))+fabs(rv1.at(i-1)));
          anorm = maxarg1 > maxarg2 ? maxarg1 : maxarg2;

        }
    
        // Accumulation of right-hand transformation
        for (i= cols ; i>=1; i--) {
          if ( i < cols ) {
            if (g) {
              for ( j=ppi; j <= cols; j++)  V(j,i)=(U(i,j)/U(i,ppi))/g;
              for ( j=ppi; j <= cols; j++) {
                for (s=0.0,k=ppi;k <= cols;k++)  s += U(i,k)*V(k,j);
                for ( k=ppi; k <= cols;k++) V(k,j) += s*V(k,i);
              }
            }
            for( j=ppi; j<=cols ; j++) V(i,j)=V(j,i)=0.0;
          }
          V(i,i)=1.0;
          g=rv1.at(i-1);
          ppi=i; 
        }
    
        // Accumulation of left-hand transformation
        for (i= cols < rows ? cols: rows; i>=1; i--) { 
          ppi=i+1;
          g=w(i);
          for( j=ppi; j<=cols ; j++) U(i,j)=0.0;
          if (g) {
            g=1.0/g;
            for ( j=ppi; j <= cols; j++) {
              for( s=0.0, k=ppi; k<=rows ; k++) s += U(k,i)*U(k,j);
              f=(s/U(i,i))*g;
              for ( k=i; k <= rows;k++)  U(k,j) += f*U(k,i);
            }
            for( j=i; j<=rows ; j++) U(j,i) *=g;
          } else {
            for( j=i; j<=rows ; j++) U(j,i) = 0.0;
          }
          ++U(i,i);
        }
    
        // Diagonalization of the bidiagonal form:
        // Loop over singular values,and over allowed iterations
        for ( k=cols; k >= 1; k--) {
          for (its=1; its <= maxIter; its++) {
            flag=true;
            //Test for splitting. Note that rv1[i] is always 0
            for (ppi=k; ppi >= 1; ppi--) {
              nm=ppi-1;
              if ((fabs(rv1.at(ppi-1))+anorm) == anorm) {
                flag=false;
                break;
              }
              if ((fabs(w(nm)+anorm) == anorm)) {
                break;
              }
            }
            //Cancellation of rv1[ppi],if ppi>1.
            if (flag) {
              c=0.0;
              s=1.0;
              for ( i=ppi; i <= k ;i++) {
                f=s*rv1.at(i-1);
                rv1.at(i-1)=c*rv1.at(i-1);
                if ((fabs(f)+anorm) == anorm) {
                  break;
                }
                g=w(i);
                h=PYTHAG(f,g);
                w(i)=h;
                h=1.0/h;
                c=g*h;
                s=-f*h;
                for ( j=1;j <= rows; j++) {
                  y=U(j,nm);
                  z=U(j,i);
                  U(j,nm)=y*c+z*s;
                  U(j,i)=z*c-y*s;
                }
              }
            }
            z=w(k);
    
            // Convergence. Singular value is made nonnegative.
            if (ppi == k) {
              if (z < 0.0) {
                w(k)=-z;
                for (j=1; j <= cols; j++) V(j,k)=-V(j,k);
              }
              break;
            }
    
            if (its == maxIter) {
              //char x[80];
              std::cout << "SVD did not converge after " <<  maxIter << " iterations!" << std::endl;
              // make all singular values zero -> rank 0
              w = 0.0; 
              return false;
            }
            // shift from bottom 2-by-2 minor
            x=w(ppi);
            nm=k-1;
            y=w(nm);
            g=rv1.at(nm-1);
            h=rv1.at(k-1);
            f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
    
            g = PYTHAG(f,1.0);
            f = ((x-z)*(x+z)+h*((y/(f+SIGN(g,f)))-h))/x;
    
            //Next QR transformation
            c=s=1.0;
            for ( j=ppi; j<=nm ;j++){
              i=j+1;
              g=rv1.at(i-1);
              y=w(i);
              h=s*g;
              g=c*g;
              z=PYTHAG(f,h);
              rv1.at(j-1)=z;
              c=f/z;
              s=h/z;
              f=x*c+g*s;
              g=g*c-x*s;
              h=y*s;
              y*=c;
              for (jj=1; jj<=cols ;jj++){
                x=V(jj,j);
                z=V(jj,i);
                V(jj,j)=x*c+z*s;
                V(jj,i)=z*c-x*s;
              }
              z=PYTHAG(f,h);
              w(j)=z;
    
              if (z) {
                z=1.0/z;
                c=f*z;
                s=h*z;
              }
              f=c*g+s*y;
              x=c*y-s*g;
              for (jj=1; jj<=rows; jj++){
                y=U(jj,j);
                z=U(jj,i);
                U(jj,j)=y*c+z*s;
                U(jj,i)=z*c-y*s;
              }
            }
            rv1.at(ppi-1)=0.0;
            rv1.at(k-1)=f;
            w(k)=x;
    
          }
        }
        return true;
    }


  double
  Matrix_Wrapper::PYTHAG(double a,double b) const
 {
     double at,bt,ct;
     at = fabs(a);
     bt = fabs(b);
     if (at > bt ) {
         ct=bt/at;
         return at*sqrt(1.0+ct*ct);
     } else {
         if (bt==0)
             return 0.0;
         else {
             ct=at/bt;
             return bt*sqrt(1.0+ct*ct);
         }
     }
 }
 
 
 double
 Matrix_Wrapper::SIGN(double a,double b) const
 {
     return ((b) >= 0.0 ? fabs(a) : -fabs(a));
 }
 
 // See <http://dsp.ee.sun.ac.za/~schwardt/dsp813/lecture10/node7.html>
 MyMatrix
 Matrix_Wrapper::pseudoinverse(double epsilon) const
 {
   int rows;
   rows = this->rows();
   int cols = this->columns();
   // calculate SVD decomposition
   MyMatrix U,V;
   MyColumnVector D;
   
   bool res;
   res = SVD(D,U,V);  // U=mxn  D=n  V=nxn
   assert(res);
   
   Matrix Dinv(cols,cols);
   Dinv = 0;
   for (unsigned int i=0; i<D.rows(); i++)
     if ( D(i+1) < epsilon )
       Dinv(i+1,i+1) = 0;
     else
       Dinv(i+1,i+1) = 1/D(i+1);
 
   #ifdef __DEBUG__
     std::cout << "MATRIX::pseudoinverse() Dinv =\n" << Dinv << std::endl;
   #endif //__DEBUG__
 
   return V * Dinv * U.transpose();
 }
 
} //namespace
