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
#include <ltilib/ltiMatrixDecomposition.h>
#include <ltilib/ltiSymmetricMatrixInversion.h>



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

const bool MyMatrix::operator==(const MyMatrix& a) const
{
  const ltiMatrix& op1 = *this;
  const ltiMatrix& op2 = a;
  return (op1 == op2);
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
  return (MyMatrix) (op1.subtract(a));
}

MyMatrix MyMatrix::operator+ (const MyMatrix& a) const
{
  ltiMatrix op1 = (*this);
  return (MyMatrix) (op1.add(a));
}

MyMatrix MyMatrix::operator* (const MyMatrix& a) const
{
  ltiMatrix op1 = (*this);
  return (MyMatrix) (op1.multiply(a));
}

MyMatrix MyMatrix::operator/ (double b) const
{
  ltiMatrix op1(*this);
  op1 /= b;
  return (MyMatrix) op1;
}


MyMatrix & MyMatrix::operator+= (const MyMatrix& a)
{
  ltiMatrix & op1 = (*this);
  op1 += a;
  return (MyMatrix &) op1;
}

MyMatrix & MyMatrix::operator-= (const MyMatrix& a)
{
  ltiMatrix & op1 = (*this);
  op1 -= a;
  return (MyMatrix &) op1;
}


// MATRIX - VECTOR Operators
MyColumnVector MyMatrix::operator* (const MyColumnVector &b) const
{
  const ltiMatrix& op1 = *this;
  ltiColumnVector op2(b);
  return (MyColumnVector) op1.multiply(op2);
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
  return lu.det(tmp);
}


MyMatrix MyMatrix::inverse() const
{
  lti::matrixInversion<double> inv;
  ltiMatrix base(*this);
  inv.apply(base);
  return (MyMatrix) base;
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
  return lu.det(tmp);
}


double& MySymmetricMatrix::operator()(unsigned int a, unsigned int b)
{
  ltiSymmetricMatrix & op1 = (*this);
  // only fill in lower triangle
  if (a < b)
    return op1.at(b-1,a-1);
  else
    return op1.at(a-1,b-1);
}
const double MySymmetricMatrix::operator()(unsigned int a, unsigned int b) const
{
  ltiSymmetricMatrix op1(*this);
  // only fill in lower triangle
  if (a < b)
    return op1.at(b-1,a-1);
  else
    return op1.at(a-1,b-1);
}


const bool MySymmetricMatrix::operator==(const MySymmetricMatrix& a) const
{
  const ltiSymmetricMatrix& op1 = *this;
  const ltiSymmetricMatrix& op2 = a;
  return (op1 == op2);
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
  return (MyMatrix) (op1.add(a));
}

MyMatrix
MySymmetricMatrix::operator- (const MyMatrix &a) const
{

  ltiMatrix op1(*this);
  return (MyMatrix) (op1.subtract(a));
}

MyMatrix
MySymmetricMatrix::operator* (const MyMatrix &a) const
{
  ltiMatrix op1(*this);
  return (MyMatrix) (op1.multiply(a));
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
  op1 += a;
  return (MySymmetricMatrix &) op1;
}

MySymmetricMatrix
MySymmetricMatrix::operator- (const MySymmetricMatrix &a) const
{
  ltiSymmetricMatrix op1 = (*this);
  op1 -= a;
  return (MySymmetricMatrix &) op1;
}

MyMatrix
MySymmetricMatrix::operator* (const MySymmetricMatrix &a) const
{
  ltiSymmetricMatrix op1 = (*this);
  return (MyMatrix) (op1.multiply(a));
}






MyColumnVector MySymmetricMatrix::operator* (const MyColumnVector &b) const
{
  const ltiSymmetricMatrix& op1 = (const ltiSymmetricMatrix&) *this;
  ltiColumnVector op2 = b;
  return (MyColumnVector) op1.multiply(op2);
}

MyMatrix MySymmetricMatrix::sub(int i_start, int i_end,
				int j_start , int j_end) const
{
  // first copy all elements from lower triangle to upper triangle
  unsigned int r = this->rows();
  unsigned int c = this->columns();
  ltiMatrix copy = *this;
  for (unsigned int i=0; i<r; i++)
    for (unsigned int j=0; j<=i; j++)
      copy.at(j,i) = copy.at(i,j);
  ltiMatrix m(copy,i_start-1,i_end-1, j_start-1,j_end-1);
  return (MyMatrix) m;
}


void
MySymmetricMatrix::resize(unsigned int i, bool copy, bool initialize)
{
  ltiSymmetricMatrix& base = (ltiSymmetricMatrix &) *this;
  base.resize(i, i, copy, initialize);
}


#endif
