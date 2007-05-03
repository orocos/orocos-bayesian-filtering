// $id: matrix_LTI.cpp 6578 2006-03-28 13:50:58Z wmeeusse $
// Copyright (C) 2002 Klaas Gadeyne <first dot last at gmail dot com>

//  
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation; either version 2.1 of the License, or
// (at your option) any later version.
//  
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
//  

#include "../config.h"
#ifdef __MATRIXWRAPPER_LTI__

#include <iostream>
#include "matrix_LTI.h"

#include <ltilib/ltiMatrixInversion.h>
#include <ltilib/ltiSVD.h>
#include <ltilib/ltiMatrixDecomposition.h>
#include <ltilib/ltiSymmetricMatrixInversion.h>
#include <ltilib/ltiCholeskyDecomposition.h>



///////////////////
// CLASS MATRIX  //
///////////////////
// Passing the constructor arguments...
MyMatrix::Matrix() : ltiMatrix() {}
MyMatrix::Matrix(int num_rows, int num_cols) : ltiMatrix(num_rows,
							 num_cols,
							 0.0){}
// Destructor
MyMatrix::~Matrix(){}

// Copy constructor
MyMatrix::Matrix(const MyMatrix& a) : ltiMatrix(a){}

// This is a bad solution, but necessary if the base class is
// ill-designed
MyMatrix::Matrix(const ltiMatrix & a) : ltiMatrix(a){}

// Number of Rows/Cols
unsigned int MyMatrix::rows() const { return ((ltiMatrix)(*this)).rows();}
unsigned int MyMatrix::columns() const { return ((ltiMatrix)(*this)).columns();}


MyRowVector MyMatrix::rowCopy(unsigned int r) const
{
  ltiMatrix temp = (ltiMatrix) *this;
  return (MyRowVector) temp.getRowCopy(r-1);
}

MyColumnVector MyMatrix::columnCopy(unsigned int c) const
{
  ltiMatrix temp = (ltiMatrix) *this;
  return (MyColumnVector) temp.getColumnCopy(c-1);
}

double& MyMatrix::operator()(unsigned int a, unsigned int b) 
{
  //ltiMatrix & op1 = (*this);
  return this->at(a-1,b-1);
}

const double MyMatrix::operator()(unsigned int a, unsigned int b) const
{
  //ltiMatrix  op1(*this);
  return this->at(a-1,b-1);
}


// MATRIX - SCALAR operators
MyMatrix& MyMatrix::operator+= (double a)
{
  ltiMatrix & op1 = (*this);
  op1 += a;
  return (MyMatrix&) op1;
}

MyMatrix& MyMatrix::operator-= (double a)
{
  ltiMatrix & op1 = (*this);
  op1 -= a;
  return (MyMatrix&) op1;
}

MyMatrix& MyMatrix::operator*= (double a)
{
  ltiMatrix & op1 = (*this);
  op1 *= a;
  return (MyMatrix&) op1;
}

MyMatrix MyMatrix::operator+ (double a) const
{
  ltiMatrix op1(*this);
  op1 += a;
  return (MyMatrix) op1;
}

MyMatrix MyMatrix::operator- (double a) const
{
  ltiMatrix op1(*this);
  op1 -= a;
  return (MyMatrix) op1;
}

MyMatrix MyMatrix::operator* (double a) const
{
  ltiMatrix op1(*this);
  op1 *= a;
  return (MyMatrix) op1;
}


MyMatrix& MyMatrix::operator/= (double a)
{
  ltiMatrix & op1 = (*this);
  op1 /= a;
  return (MyMatrix&) op1;
}


// MATRIX - MATRIX Operators
MyMatrix MyMatrix::operator- (const MyMatrix& a) const
{
  ltiMatrix op1 = (*this);
  ltiMatrix op2 = a;
  ltiMatrix result = (ltiMatrix) (op1.subtract(op2));
  return (MyMatrix) result;
}

MyMatrix MyMatrix::operator+ (const MyMatrix& a) const
{
  ltiMatrix op1 = (*this);
  ltiMatrix op2 = a;
  ltiMatrix result = (ltiMatrix) (op1.add(op2));
  return (MyMatrix) result;
}

MyMatrix MyMatrix::operator* (const MyMatrix& a) const
{
  ltiMatrix op1 = (*this);
  ltiMatrix op2 = a;
  ltiMatrix result = (ltiMatrix) (op1.multiply(op2));
  return (MyMatrix) result;
}

MyMatrix MyMatrix::operator/ (double b) const
{
  ltiMatrix op1(*this);
  op1 /= b;
  return (MyMatrix) op1;
}


// THIS SHOULD BE THE IDEAL SOLUTION (IMPOSSIBLE DUE TO BASE CLASS??)
// const MyMatrix & MyMatrix::operator* (const MyMatrix& a)
// {
//   ltiMatrix & op1 = (*this);
//   ltiMatrix & op2 = a;
//   ltiMatrix result = (ltiMatrix) (op1 * op2);
//   return (MyMatrix &) result;
// }


MyMatrix & MyMatrix::operator+= (const MyMatrix& a)
{
  ltiMatrix & op1 = (*this);
  const ltiMatrix & op2 = a;
  op1 += op2;
  return (MyMatrix &) op1;
}

MyMatrix & MyMatrix::operator-= (const MyMatrix& a)
{
  ltiMatrix & op1 = (*this);
  const ltiMatrix & op2 = a;
  op1 -= op2;
  return (MyMatrix &) op1;
}


// MATRIX - VECTOR Operators
MyColumnVector MyMatrix::operator* (const MyColumnVector &b) const
{
  ltiMatrix op1(*this);
  ltiColumnVector op2(b);
  ltiColumnVector result = op1.multiply(op2);
  return (MyColumnVector) result;
}

// Set all elements equal to a
MyMatrix &MyMatrix::operator=(const double a)
{
  ltiMatrix & op1 = (*this);
  op1.fill(a,0,0,(*this).rows(),(*this).columns());
  return (MyMatrix&) op1;
}



MyMatrix MyMatrix::transpose() const
{
  ltiMatrix base(*this);
  base.transpose();
  return (MyMatrix) base;
}

double MyMatrix::determinant() const
{
  ltiMatrix & base = (ltiMatrix &) *this;
  lti::matrix<double> tmp;
  tmp.resize(base.size());
  for (int i=0;i<tmp.rows();i++)
    for (int j=0;j<tmp.columns();j++) {
      tmp.at(i,j)=(double)(base.at(i,j));
    }
  lti::luDecomposition<double> lu;
  double result = lu.det(tmp);
  return result;
}


MyMatrix MyMatrix::inverse() const
{
  lti::matrixInversion<double> inv;
  ltiMatrix base(*this);
  inv.apply(base);
  return (MyMatrix) base;
}

// See <http://dsp.ee.sun.ac.za/~schwardt/dsp813/lecture10/node7.html>
MyMatrix MyMatrix::pseudoinverse(double epsilon) const
{
  MyMatrix result;
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

  result = V * Dinv * U.transpose();
  return result;
}

int 
MyMatrix::convertToSymmetricMatrix(MySymmetricMatrix& sym)
{
  // test if matrix is square matrix
  assert(this->rows() == this->columns() );
  
  // if necessairy, resize sym
  // only check cols or rows. Symmetric matrix is square.
  if ( sym.rows() != this->rows() )
    sym.resize(this->rows());
  
  // copy elements 
  for ( unsigned int i=0; i<this->rows(); i++ )
    for ( unsigned int j=0; j<=i; j++ )
    {
	    sym[i][j] = (*this)[i][j];
      	    sym[j][i] = (*this)[i][j];
    }
  return 0;
}

// get sub matrix
MyMatrix MyMatrix::sub(int i_start, int i_end, int j_start , int j_end) const
{
  ltiMatrix m(*this,i_start-1,i_end-1, j_start-1,j_end-1);
  return (MyMatrix) m;
}


void
MyMatrix::resize(unsigned int i, unsigned int j, bool copy, bool initialize)
{
  ltiMatrix& base = (ltiMatrix &) *this;
  base.resize(i,j, copy, initialize);
}



bool 
MyMatrix::SVD(MyColumnVector& D, MyMatrix& U, MyMatrix& V) const
{
  // BIG FAT WARNING: THIS ALGORITHM HAS PROBLEMS WHEN VALUES OF MATRIX ARE TOO SMALL
  // BIG FAT WARNING: THIS ALGORITHM HAS PROBLEMS WHEN VALUES OF MATRIX ARE TOO SMALL
  // BIG FAT WARNING: THIS ALGORITHM HAS PROBLEMS WHEN VALUES OF MATRIX ARE TOO SMALL
  // BIG FAT WARNING: THIS ALGORITHM HAS PROBLEMS WHEN VALUES OF MATRIX ARE TOO SMALL
  // BIG FAT WARNING: THIS ALGORITHM HAS PROBLEMS WHEN VALUES OF MATRIX ARE TOO SMALL
  // BIG FAT WARNING: THIS ALGORITHM HAS PROBLEMS WHEN VALUES OF MATRIX ARE TOO SMALL
  // BIG FAT WARNING: THIS ALGORITHM HAS PROBLEMS WHEN VALUES OF MATRIX ARE TOO SMALL
  // BIG FAT WARNING: THIS ALGORITHM HAS PROBLEMS WHEN VALUES OF MATRIX ARE TOO SMALL
  // BIG FAT WARNING: THIS ALGORITHM HAS PROBLEMS WHEN VALUES OF MATRIX ARE TOO SMALL
  // BIG FAT WARNING: THIS ALGORITHM HAS PROBLEMS WHEN VALUES OF MATRIX ARE TOO SMALL
  // BIG FAT WARNING: THIS ALGORITHM HAS PROBLEMS WHEN VALUES OF MATRIX ARE TOO SMALL

  lti::singularValueDecomp<double> mysvd(true);  //true = sort singular values
  return mysvd.apply(*this,U,D,V);

  // BIG FAT WARNING: THIS ALGORITHM HAS PROBLEMS WHEN VALUES OF MATRIX ARE TOO SMALL
  // BIG FAT WARNING: THIS ALGORITHM HAS PROBLEMS WHEN VALUES OF MATRIX ARE TOO SMALL
  // BIG FAT WARNING: THIS ALGORITHM HAS PROBLEMS WHEN VALUES OF MATRIX ARE TOO SMALL
  // BIG FAT WARNING: THIS ALGORITHM HAS PROBLEMS WHEN VALUES OF MATRIX ARE TOO SMALL
  // BIG FAT WARNING: THIS ALGORITHM HAS PROBLEMS WHEN VALUES OF MATRIX ARE TOO SMALL
  // BIG FAT WARNING: THIS ALGORITHM HAS PROBLEMS WHEN VALUES OF MATRIX ARE TOO SMALL
  // BIG FAT WARNING: THIS ALGORITHM HAS PROBLEMS WHEN VALUES OF MATRIX ARE TOO SMALL
  // BIG FAT WARNING: THIS ALGORITHM HAS PROBLEMS WHEN VALUES OF MATRIX ARE TOO SMALL
  // BIG FAT WARNING: THIS ALGORITHM HAS PROBLEMS WHEN VALUES OF MATRIX ARE TOO SMALL
  // BIG FAT WARNING: THIS ALGORITHM HAS PROBLEMS WHEN VALUES OF MATRIX ARE TOO SMALL
  // BIG FAT WARNING: THIS ALGORITHM HAS PROBLEMS WHEN VALUES OF MATRIX ARE TOO SMALL
}



/////////////////////////////
// CLASS SYMMETRIC MATRIX  //
/////////////////////////////

MySymmetricMatrix::SymmetricMatrix() : ltiSymmetricMatrix() {}
MySymmetricMatrix::SymmetricMatrix(int n) : ltiSymmetricMatrix(n,n) {}

// Copy constructor
MySymmetricMatrix::SymmetricMatrix(const SymmetricMatrix& a) : ltiSymmetricMatrix(a){}
MySymmetricMatrix::SymmetricMatrix(const ltiSymmetricMatrix & a) : ltiSymmetricMatrix(a){}

// Destructor
MySymmetricMatrix::~SymmetricMatrix(){}

// Ask Number of Rows and Columns
unsigned int MySymmetricMatrix::rows() const { return (((ltiSymmetricMatrix)(*this)).rows());}
unsigned int MySymmetricMatrix::columns() const { return (((ltiSymmetricMatrix)(*this)).rows());}

MySymmetricMatrix MySymmetricMatrix::transpose() const 
{return (*this);}

MySymmetricMatrix MySymmetricMatrix::inverse() const
{
  lti::matrixInversion<double> inv;
  ltiSymmetricMatrix base(*this);
  inv.apply(base);
  return (MySymmetricMatrix) base;
}

double MySymmetricMatrix::determinant() const
{
  ltiSymmetricMatrix & base = (ltiSymmetricMatrix &) *this;
  lti::matrix<double> tmp;
  tmp.resize(base.size());
  for (int i=0;i<tmp.rows();i++)
    for (int j=0;j<tmp.columns();j++) {
      tmp.at(i,j)=(double)(base.at(i,j));
    }
  lti::luDecomposition<double> lu;
  double result = lu.det(tmp);
  return result;
}


double& MySymmetricMatrix::operator()(unsigned int a, unsigned int b) 
{
  ltiSymmetricMatrix & op1 = (*this);
  return op1.at(a-1,b-1);
}
const double MySymmetricMatrix::operator()(unsigned int a, unsigned int b) const
{
  ltiSymmetricMatrix op1(*this);
  return op1.at(a-1,b-1);
}


// MATRIX - SCALAR operators
MySymmetricMatrix& MySymmetricMatrix::operator=(const double a)
{
  ltiSymmetricMatrix temp = (ltiSymmetricMatrix) *this;
  temp.fill(a,0,0,temp.rows(),temp.columns());
  *this = (MySymmetricMatrix) temp;

  return *this;
}

// MATRIX - SCALAR operators
MySymmetricMatrix& 
MySymmetricMatrix::operator +=(double a)
{
  ltiSymmetricMatrix & op1 = (*this);
  op1 += a;
  return (MySymmetricMatrix&) op1;
}

MySymmetricMatrix& 
MySymmetricMatrix::operator -=(double a)
{
  ltiSymmetricMatrix & op1 = (*this);
  op1 -= a;
  return (MySymmetricMatrix&) op1;
}

MySymmetricMatrix& 
MySymmetricMatrix::operator *=(double b)
{
  ltiSymmetricMatrix & op1 = (*this);
  op1 *= b;
  return (MySymmetricMatrix&) op1;
}
 
MySymmetricMatrix& 
MySymmetricMatrix::operator /=(double b)
{
  ltiSymmetricMatrix & op1 = (*this);
  op1 /= b;
  return (MySymmetricMatrix&) op1;
}

MySymmetricMatrix 
MySymmetricMatrix::operator+ (double a) const
{
  ltiSymmetricMatrix op1(*this);
  op1 += a;
  return (MySymmetricMatrix) op1;
}

MySymmetricMatrix 
MySymmetricMatrix::operator- (double a) const
{
  ltiSymmetricMatrix op1(*this);
  op1 -= a;
  return (MySymmetricMatrix) op1;
}

MySymmetricMatrix MySymmetricMatrix::operator* (double a) const
{
  ltiSymmetricMatrix op1(*this);
  op1 *= a;
  return (MySymmetricMatrix) op1;
}

MySymmetricMatrix MySymmetricMatrix::operator/ (double b) const
{
  ltiSymmetricMatrix op1(*this);
  op1 /= b;
  return (MySymmetricMatrix) op1;
}





// SYMMETRICMATRIX - MATRIX operators
MyMatrix& MySymmetricMatrix::operator +=(const MyMatrix& a)
{
  ltiSymmetricMatrix & op1 = (*this);
  const ltiMatrix & op2 = a;
  op1 += op2;
  return (MyMatrix &) op1;
}

MyMatrix& 
MySymmetricMatrix::operator -=(const MyMatrix& a)
{
  ltiSymmetricMatrix & op1 = (*this);
  const ltiMatrix & op2 = a;
  op1 -= op2;
  return (MyMatrix &) op1;
}


MyMatrix 
MySymmetricMatrix::operator+ (const MyMatrix &a) const
{
  ltiMatrix op1(*this);
  ltiMatrix op2 = a;
  ltiMatrix result = (ltiMatrix) (op1.add(op2));
  return (MyMatrix) result;
}

MyMatrix 
MySymmetricMatrix::operator- (const MyMatrix &a) const
{

  ltiMatrix op1(*this);
  ltiMatrix op2 = a;
  ltiMatrix result = (ltiMatrix) (op1.subtract(op2));
  return (MyMatrix) result;
}

MyMatrix 
MySymmetricMatrix::operator* (const MyMatrix &a) const
{
  ltiMatrix op1(*this);
  ltiMatrix op2 = a;
  ltiMatrix result = (ltiMatrix) (op1.multiply(op2));
  return (MyMatrix) result;
}

MyMatrix& 
MyMatrix::operator =(const MySymmetricMatrix& a)
{
  *this =(MyMatrix) a;

  return *this;
}

MySymmetricMatrix& 
MySymmetricMatrix::operator +=(const MySymmetricMatrix& a)
{
  ltiSymmetricMatrix & op1 = (*this);
  const ltiSymmetricMatrix & op2 = a;
  op1 += op2;
  return (MySymmetricMatrix &) op1;
}

MySymmetricMatrix& 
MySymmetricMatrix::operator -=(const MySymmetricMatrix& a)
{
  ltiSymmetricMatrix & op1 = (*this);
  const ltiSymmetricMatrix & op2 = a;
  op1 -= op2;
  return (MySymmetricMatrix &) op1;
}


MySymmetricMatrix
MySymmetricMatrix::operator+ (const MySymmetricMatrix &a) const
{
  ltiSymmetricMatrix op1 = (*this);
  const ltiSymmetricMatrix & op2 = a;
  op1 += op2;
  return (MySymmetricMatrix &) op1;
}

MySymmetricMatrix
MySymmetricMatrix::operator- (const MySymmetricMatrix &a) const
{
  ltiSymmetricMatrix op1 = (*this);
  const ltiSymmetricMatrix & op2 = a;
  op1 -= op2;
  return (MySymmetricMatrix &) op1;
}

MySymmetricMatrix 
MySymmetricMatrix::operator* (const MySymmetricMatrix &a) const
{
  ltiSymmetricMatrix op1 = (*this);
  ltiSymmetricMatrix op2 = a;
  ltiSymmetricMatrix result = (ltiSymmetricMatrix) (op1.multiply(op2));
  return (MySymmetricMatrix) result;
}






MyColumnVector MySymmetricMatrix::operator* (const MyColumnVector &b) const
{
  ltiSymmetricMatrix op1 = (ltiSymmetricMatrix) *this;
  ltiColumnVector op2(b);
  ltiColumnVector result = op1.multiply(op2);
  return (MyColumnVector) result;
}

MyMatrix MySymmetricMatrix::sub(int i_start, int i_end, 
				int j_start , int j_end) const
{
  ltiMatrix m(*this,i_start-1,i_end-1, j_start-1,j_end-1);
  return (MyMatrix) m;
}


void
MySymmetricMatrix::resize(unsigned int i, bool copy, bool initialize)
{
  ltiSymmetricMatrix& base = (ltiSymmetricMatrix &) *this;
  base.resize(i, i, copy, initialize);
}


bool
MySymmetricMatrix::cholesky(MyMatrix& m) const
{
  lti::choleskyDecomposition<double> Cholesky;
  return Cholesky.apply(*this,m);
}
#endif
