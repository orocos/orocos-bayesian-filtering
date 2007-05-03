// $Id: matrix_BOOST.cpp 27906 2007-04-27 11:50:53Z wmeeusse $
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
#ifdef __MATRIXWRAPPER_BOOST__

#include "matrix_BOOST.h"


// Passing the constructor arguments...
MyMatrix::Matrix() : BoostMatrix() {}
MyMatrix::Matrix(int num_rows, int num_cols) : BoostMatrix(num_rows,
							   num_cols){}

// Destructor
MyMatrix::~Matrix(){}

// Copy constructor
MyMatrix::Matrix(const MyMatrix& a) : BoostMatrix(a){}

// ill-designed
MyMatrix::Matrix(const BoostMatrix & a) : BoostMatrix(a){}

// Number of Rows/Cols
unsigned int MyMatrix::rows() const { return this->size1();}
unsigned int MyMatrix::columns() const { return this->size2();}

// MATRIX - SCALAR operators
MyMatrix& MyMatrix::operator+= (double a)
{
  BoostMatrix & op1 = *this;
  op1 += boost::numeric::ublas::scalar_matrix<double>(rows(),columns(),a);
  return (MyMatrix&)op1;
}

MyMatrix& MyMatrix::operator-= (double a)
{
  BoostMatrix & op1 = (*this);
  op1 -= boost::numeric::ublas::scalar_matrix<double>(rows(),columns(),a);
  return (MyMatrix&) op1;
}

MyMatrix& MyMatrix::operator*= (double a)
{
  BoostMatrix & op1 = (*this);
  op1 *= a;
  return (MyMatrix&) op1;
}

MyMatrix& MyMatrix::operator/= (double a)
{
  BoostMatrix & op1 = (*this);
  op1 /= a;
  return (MyMatrix&) op1;
}

MyMatrix MyMatrix::operator+ (double a) const
{
  return (BoostMatrix)(((BoostMatrix)(*this)) + boost::numeric::ublas::scalar_matrix<double>(rows(),columns(),a));
}

MyMatrix MyMatrix::operator- (double a) const
{
  return (BoostMatrix)(((BoostMatrix)(*this)) - boost::numeric::ublas::scalar_matrix<double>(rows(),columns(),a));
}

MyMatrix MyMatrix::operator* (double a) const
{
  BoostMatrix op1 = (*this);
  op1 *= a;
  return (MyMatrix&) op1;
}

MyMatrix MyMatrix::operator/ (double a) const
{
  BoostMatrix op1 = (*this);
  op1 /= a;
  return (MyMatrix&) op1;
}

MyMatrix& 
MyMatrix::operator =(const MySymmetricMatrix& a)
{
  *this =(MyMatrix) a;

  return *this;
}

// MATRIX - MATRIX Operators
MyMatrix MyMatrix::operator- (const MyMatrix& a) const
{
  return (MyMatrix)((BoostMatrix&)(*this) - (BoostMatrix&)(a));
}

MyMatrix MyMatrix::operator+ (const MyMatrix& a) const
{
  BoostMatrix op2 = a;
  BoostMatrix result = (BoostMatrix) ((BoostMatrix)(*this) + op2);
  return (MyMatrix) result;
}

MyMatrix MyMatrix::operator* (const MyMatrix& a) const
{
  const BoostMatrix & op1 = (*this);
  return (MyMatrix) prod(op1, (const BoostMatrix&)a);
}

MyMatrix & MyMatrix::operator+= (const MyMatrix& a)
{
  BoostMatrix & op1 = (*this);
  const BoostMatrix & op2 = a;
  op1 += op2;
  return (MyMatrix &) op1;
}

MyMatrix & MyMatrix::operator-= (const MyMatrix& a)
{
  BoostMatrix & op1 = (*this);
  const BoostMatrix & op2 = a;
  op1 -= op2;
  return (MyMatrix &) op1;
}


// MATRIX - VECTOR Operators
MyColumnVector MyMatrix::operator* (const MyColumnVector &b) const
{
  const BoostMatrix op1 = (BoostMatrix) *this;
  BoostColumnVector result = prod(op1, ((const BoostColumnVector&)b));
  return (MyColumnVector) result;
}


// Set all elements equal to a
MyMatrix&
 MyMatrix::operator=(double a)
{
  for (unsigned int i=0; i<this->rows(); i++)
    for (unsigned int j=0; j<this->columns(); j++)
      (*this)(i,j) = a;

  return *this;
}


MyRowVector MyMatrix::rowCopy(unsigned int r) const
{
  unsigned int cols = columns();
  BoostRowVector temp(cols);
  for (unsigned int i=0; i<cols; i++)
    temp(i) = (*this)(r,i);
  return (MyRowVector) temp;
}

MyColumnVector MyMatrix::columnCopy(unsigned int c) const
{
  unsigned int ro = rows();
  BoostColumnVector temp(ro);
  for (unsigned int i=0; i<ro; i++)
    temp(i) = (*this)(ro,i);
  return (MyColumnVector) temp;
}




MyMatrix MyMatrix::transpose() const
{
  return (MyMatrix) trans((BoostMatrix &) *this);
}

double MyMatrix::determinant() const
{
  assert(rows() == columns());
 
  // create a working copy of the input 
  BoostMatrix mat(*this);
  boost::numeric::ublas::permutation_matrix<std::size_t> pivots(rows());
  lu_factorize(mLu, pivots);
  double det = 1.0;
  for (std::size_t i=0; i < pivots.size(); ++i) {
    if (pivots(i) != i) 
      det *= -1.0;
    det *= mLu(i,i);
  }
  return det;
}


MyMatrix MyMatrix::inverse() const
{
  BoostMatrix & base = (BoostMatrix &) *this;
  BoostMatrix inverted = base.i();
  return (MyMatrix) inverted;
}

MyMatrix MyMatrix::pseudoinverse(double epsilon) const
{
  MyMatrix result, U,V;  MyColumnVector D;
  int rows = this->rows();
  int cols = this->columns();

  if (cols > rows)
    this->transpose().SVD(D,U,V);
  else
    this->SVD(D,U,V);  // U=rxc  D=c   V=cxc
  
  int D_size = min(rows, cols);
  MyMatrix Dinv(D_size,D_size);
  Dinv = 0;
  for (unsigned int i=0; i<(unsigned int)D_size; i++)
    if ( D(i+1) < epsilon )
      Dinv(i+1,i+1) = 0;
    else
      Dinv(i+1,i+1) = 1/D(i+1);

  if (cols > rows)
    result = U * Dinv * V.transpose();
  else
    result = V * Dinv * U.transpose();

  return result;
}

int 
MyMatrix::convertToSymmetricMatrix(MySymmetricMatrix& sym)
{
  // test if matrix is square matrix
  assert(this->rows() == this->columns());
  
  // if necessairy, resize sym
  // only check cols or rows. Symmetric matrix is square.
  if ( sym.rows() != this->rows() )
    sym.ReSize(this->rows());
  
  
  // copy elements 
  for ( unsigned int i=0; i<this->rows(); i++ )
    for ( unsigned int j=0; j<=i; j++ )
      sym(i+1,j+1) = (*this)(i+1,j+1);
  return 0;
}

void
MyMatrix::resize(unsigned int i, unsigned int j, bool copy, bool initialize)
{
  BoostMatrix & temp = (BoostMatrix &) (*this);
  temp.ReSize(i,j);
}

// get sub matrix
MyMatrix MyMatrix::sub(int i_start, int i_end, int j_start , int j_end) const
{
  return (MyMatrix)(this->SubMatrix(i_start, i_end, j_start, j_end));
}


bool 
MyMatrix::SVD(MyColumnVector& D, MyMatrix& U, MyMatrix& V) const
{
  if (columns() > rows()){
    cout << endl << "ERROR" << endl;    
    cout << "Boost doesn't support svd for columns > rows" << endl;
    cout << "you can avoid this problem by using the LTI library" << endl;
    return false;
  }

  int Acolumns = columns();
  BOOST::DiagonalMatrix d(Acolumns);
  BoostMatrix u(rows(),Acolumns); 
  BoostMatrix v(Acolumns, Acolumns);
  BOOST::SVD(*this,d,u,v);
  D.ReSize(Acolumns);
  for ( int row = 0; row < Acolumns; row++) 
      D(row+1) = d(row+1);

  U = (MyMatrix &) u;
  V = (MyMatrix &) v;

  return true;
}


/////////////////////////////
// CLASS SYMMETRIC MATRIX  //
/////////////////////////////

MySymmetricMatrix::SymmetricMatrix() : BoostSymmetricMatrix() {}
MySymmetricMatrix::SymmetricMatrix(int n) : BoostSymmetricMatrix(n) {}

// Copy constructor
MySymmetricMatrix::SymmetricMatrix(const SymmetricMatrix& a) : BoostSymmetricMatrix(a){}
MySymmetricMatrix::SymmetricMatrix(const BoostSymmetricMatrix & a) : BoostSymmetricMatrix(a){}

// Destructor
MySymmetricMatrix::~SymmetricMatrix(){}

// Ask Number of Rows and Columns
unsigned int MySymmetricMatrix::rows() const { return this->Nrows();}
unsigned int MySymmetricMatrix::columns() const { return this->Ncols();}


MySymmetricMatrix MySymmetricMatrix::transpose() const {return (*this);}

MySymmetricMatrix MySymmetricMatrix::inverse() const
{
  BoostSymmetricMatrix & base = (BoostSymmetricMatrix &) *this;
  BoostSymmetricMatrix inverted = base.i();
  return (MySymmetricMatrix) inverted;
}

double MySymmetricMatrix::determinant() const
{
  BoostSymmetricMatrix & base = (BoostSymmetricMatrix &) *this;
  BOOST::LogAndSign temp = base.LogDeterminant();
  double result = temp.Value();
  return result;
}


double& MyMatrix::operator()(unsigned int a, unsigned int b) 
{
  BoostMatrix & op1 = (*this);
  return op1(a,b);
}

const double MyMatrix::operator()(unsigned int a, unsigned int b) const
{
  BoostMatrix  op1(*this);
  return op1(a,b);
}


// void MySymmetricMatrix::operator=(const MySymmetricMatrix &a)
// {
//   BoostSymmetricMatrix temp = (BoostSymmetricMatrix) *this;
//   const BoostSymmetricMatrix temp2 = (const BoostSymmetricMatrix &) a; 
//   temp = temp2;
//   *this = (MySymmetricMatrix) temp;
// }

// Set all elements equal to a
MySymmetricMatrix& MySymmetricMatrix::operator=(const double a)
{
  BoostSymmetricMatrix temp = (BoostSymmetricMatrix) *this;
  temp = a;
  *this = (MySymmetricMatrix) temp;

  return *this;
}


// SYMMETRICMATRIX - SCALAR operators
MySymmetricMatrix& MySymmetricMatrix::operator +=(double a)
{
  BoostSymmetricMatrix & op1 = (*this);
  op1 += a;
  return (MySymmetricMatrix&) op1;
}

MySymmetricMatrix& MySymmetricMatrix::operator -=(double a)
{
  BoostSymmetricMatrix & op1 = (*this);
  op1 -= a;
  return (MySymmetricMatrix&) op1;
}

MySymmetricMatrix& MySymmetricMatrix::operator *=(double b)
{
  BoostSymmetricMatrix & op1 = (*this);
  op1 *= b;
  return (MySymmetricMatrix&) op1;
}
  
MySymmetricMatrix& MySymmetricMatrix::operator /=(double b)
{
  BoostSymmetricMatrix & op1 = (*this);
  op1 /= b;
  return (MySymmetricMatrix&) op1;
}

MySymmetricMatrix MySymmetricMatrix::operator +(double a) const
{
  BoostSymmetricMatrix op1 = (*this);
  op1 += a;
  return (MySymmetricMatrix) op1;
}

MySymmetricMatrix MySymmetricMatrix::operator -(double a) const
{
  BoostSymmetricMatrix op1 = (*this);
  op1 -= a;
  return (MySymmetricMatrix) op1;
}

MySymmetricMatrix MySymmetricMatrix::operator *(double b) const
{
  BoostSymmetricMatrix op1 = (*this);
  op1 *= b;
  return (MySymmetricMatrix) op1;
}
  
MySymmetricMatrix MySymmetricMatrix::operator /(double b) const
{
  BoostSymmetricMatrix op1 = (*this);
  op1 /= b;
  return (MySymmetricMatrix) op1;
}




// SYMMETRICMATRIX - MATRIX operators
MyMatrix& MySymmetricMatrix::operator +=(const MyMatrix& a)
{
  BoostSymmetricMatrix & op1 = (*this);
  const BoostMatrix & op2 = a;
  op1 += op2;
  return (MyMatrix &) op1;
}

MyMatrix& MySymmetricMatrix::operator -=(const MyMatrix& a)
{
  BoostSymmetricMatrix & op1 = (*this);
  const BoostMatrix & op2 = a;
  op1 -= op2;
  return (MyMatrix &) op1;
}


MyMatrix MySymmetricMatrix::operator+ (const MyMatrix &a) const
{
  const BoostMatrix op1 = (*this);
  BoostMatrix op2 = a;
  BoostMatrix result = (BoostMatrix) (op1 + op2);
  return (MyMatrix) result;
}

MyMatrix MySymmetricMatrix::operator- (const MyMatrix &a) const
{
  const BoostMatrix op1 = (*this);
  BoostMatrix op2 = a;
  BoostMatrix result = (BoostMatrix) (op1 - op2);
  return (MyMatrix) result;
}

MyMatrix MySymmetricMatrix::operator* (const MyMatrix &a) const
{
  const BoostMatrix op1 = (*this);
  BoostMatrix op2 = a;
  BoostMatrix result = (BoostMatrix) (op1 * op2);
  return (MyMatrix) result;
}



// SYMMETRICMATRIX - SYMMETRICMATRIX operators
MySymmetricMatrix& MySymmetricMatrix::operator +=(const MySymmetricMatrix& a)
{
  BoostSymmetricMatrix & op1 = (*this);
  const BoostSymmetricMatrix & op2 = a;
  op1 += op2;
  return (MySymmetricMatrix &) op1;
}

MySymmetricMatrix& MySymmetricMatrix::operator -=(const MySymmetricMatrix& a)
{
  BoostSymmetricMatrix & op1 = (*this);
  const BoostSymmetricMatrix & op2 = a;
  op1 -= op2;
  return (MySymmetricMatrix &) op1;
}

MySymmetricMatrix MySymmetricMatrix::operator+ (const MySymmetricMatrix &a) const
{
  BoostSymmetricMatrix op1 = (*this);
  BoostSymmetricMatrix op2 = a;
  BoostSymmetricMatrix result = (BoostSymmetricMatrix) (op1 + op2);
  return (MySymmetricMatrix) result;
}

MySymmetricMatrix MySymmetricMatrix::operator- (const MySymmetricMatrix &a) const
{
  BoostSymmetricMatrix op1 = (*this);
  BoostSymmetricMatrix op2 = a;
  BoostSymmetricMatrix result = (BoostSymmetricMatrix) (op1 - op2);
  return (MySymmetricMatrix) result;
}

MySymmetricMatrix MySymmetricMatrix::operator* (const MySymmetricMatrix &a) const
{
  BoostSymmetricMatrix op1 = (*this);
  BoostSymmetricMatrix op2 = a;
  BoostSymmetricMatrix result = (BoostSymmetricMatrix) (op1 * op2);
  return (MySymmetricMatrix) result;
}




MyColumnVector MySymmetricMatrix::operator* (const MyColumnVector &b) const
{
  const BoostSymmetricMatrix op1 = (BoostSymmetricMatrix) *this;
  BoostColumnVector op2 = (BoostMatrix) b;
  BoostColumnVector result = op1 * op2;
  return (MyColumnVector) result;
}

MyMatrix MySymmetricMatrix::sub(int i_start, int i_end, 
				int j_start , int j_end) const
{
  return (MyMatrix)(this->SubMatrix(i_start, i_end, j_start, j_end));
}


void
MySymmetricMatrix::resize(unsigned int i, bool copy, bool initialize)
{
  BoostSymmetricMatrix & temp = (BoostSymmetricMatrix &) (*this);
  temp.ReSize(i);
}



bool
MySymmetricMatrix::cholesky(MyMatrix& m) const
{
  BoostMatrix & op1 = m;
  try{op1 = Cholesky(*this);}
  catch(BOOST::NPDException t){ return false;}
  return true;
}


double& MySymmetricMatrix::operator()(unsigned int a, unsigned int b) 
{
  BoostSymmetricMatrix & op1 = (*this);
  return op1(a,b);
}

const double MySymmetricMatrix::operator()(unsigned int a, unsigned int b) const
{
  BoostSymmetricMatrix op1(*this);
  return op1(a,b);
}


#endif
