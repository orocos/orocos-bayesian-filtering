// $Id$
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
#ifdef __MATRIXWRAPPER_NEWMAT__

#include <iostream>
#include "matrix_NEWMAT.h"
#include <assert.h>


// Passing the constructor arguments...
MyMatrix::Matrix() : NewMatMatrix() {}
MyMatrix::Matrix(int num_rows, int num_cols) : NewMatMatrix(num_rows,
							    num_cols){}

// Destructor
MyMatrix::~Matrix(){}

// Copy constructor
MyMatrix::Matrix(const MyMatrix& a) : NewMatMatrix(a){}

// ill-designed
MyMatrix::Matrix(const NewMatMatrix & a) : NewMatMatrix(a){}

// Number of Rows/Cols
unsigned int MyMatrix::rows() const { return this->Nrows();}
unsigned int MyMatrix::columns() const { return this->Ncols();}

// MATRIX - SCALAR operators
MyMatrix& MyMatrix::operator+= (double a)
{
  NewMatMatrix & op1 = (*this);
  op1 += a;
  return (MyMatrix&) op1;
}

MyMatrix& MyMatrix::operator-= (double a)
{
  NewMatMatrix & op1 = (*this);
  op1 -= a;
  return (MyMatrix&) op1;
}

MyMatrix& MyMatrix::operator*= (double a)
{
  NewMatMatrix & op1 = (*this);
  op1 *= a;
  return (MyMatrix&) op1;
}

MyMatrix& MyMatrix::operator/= (double a)
{
  NewMatMatrix & op1 = (*this);
  op1 /= a;
  return (MyMatrix&) op1;
}

MyMatrix MyMatrix::operator+ (double a) const
{
  // make copy
  NewMatMatrix op1 = (*this);
  op1 += a;
  return (MyMatrix&) op1;
}

MyMatrix MyMatrix::operator- (double a) const
{
  // make copy
  NewMatMatrix op1 = (*this);
  op1 -= a;
  return (MyMatrix&) op1;
}

MyMatrix MyMatrix::operator* (double a) const
{
  // make copy
  NewMatMatrix op1 = (*this);
  op1 *= a;
  return (MyMatrix&) op1;
}

MyMatrix MyMatrix::operator/ (double a) const
{
  // make copy
  NewMatMatrix op1 = (*this);
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
  const NewMatMatrix& op1 = (*this);
  const NewMatMatrix& op2 = a;
  NewMatMatrix result = (NewMatMatrix) (op1 - op2);
  return (MyMatrix) result;
}

MyMatrix MyMatrix::operator+ (const MyMatrix& a) const
{
  const NewMatMatrix& op1 = (*this);
  const NewMatMatrix& op2 = a;
  NewMatMatrix result = (NewMatMatrix) (op1 + op2);
  return (MyMatrix) result;
}

MyMatrix MyMatrix::operator* (const MyMatrix& a) const
{
  const NewMatMatrix& op1 = (*this);
  const NewMatMatrix& op2 = a;
  NewMatMatrix result = (NewMatMatrix) (op1 * op2);
  return (MyMatrix) result;
}

MyMatrix & MyMatrix::operator+= (const MyMatrix& a)
{
  NewMatMatrix & op1 = (*this);
  const NewMatMatrix & op2 = a;
  op1 += op2;
  return (MyMatrix &) op1;
}

MyMatrix & MyMatrix::operator-= (const MyMatrix& a)
{
  NewMatMatrix & op1 = (*this);
  const NewMatMatrix & op2 = a;
  op1 -= op2;
  return (MyMatrix &) op1;
}


// MATRIX - VECTOR Operators
MyColumnVector MyMatrix::operator* (const MyColumnVector &b) const
{
  const NewMatMatrix& op1 = (*this);
  const NewMatColumnVector& op2 = b;
  return (MyColumnVector) (op1 * op2);
}


// Set all elements equal to a
MyMatrix& MyMatrix::operator=(double a)
{
  NewMatMatrix temp = (NewMatMatrix) *this;
  temp = a;
  *this = (MyMatrix) temp;

  return *this;
}


MyRowVector MyMatrix::rowCopy(unsigned int r) const
{
  const NewMatMatrix& temp = (*this);
  return (MyRowVector) temp.Row(r);
}

MyColumnVector MyMatrix::columnCopy(unsigned int c) const
{
  const NewMatMatrix& temp = (*this);
  return (MyColumnVector) temp.Column(c);
}




MyMatrix MyMatrix::transpose() const
{
  const NewMatMatrix & base = (*this);
  NewMatMatrix transposedbase = base.t();
  return (MyMatrix) transposedbase;
}

double MyMatrix::determinant() const
{
  const NewMatMatrix& base = (*this);
  NEWMAT::LogAndSign temp = base.LogDeterminant();
  double result = temp.Value();
  return result;
}


MyMatrix MyMatrix::inverse() const
{
  const NewMatMatrix & base = (*this);
  NewMatMatrix inverted = base.i();
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
  assert( this->rows() == this->columns() );
  
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
  NewMatMatrix & temp = (NewMatMatrix &) (*this);
  temp.ReSize(i,j);
}

// get sub matrix
MyMatrix MyMatrix::sub(int i_start, int i_end, int j_start , int j_end) const
{
  assert(i_start >= 0 && i_start <= rows());
  assert(j_start >= 0 && j_start <= columns());
  assert(i_end   >= 0 && i_end   <= rows());
  assert(j_end   >= 0 && j_end   <= columns());
  assert(i_start <= i_end);
  assert(j_start <= j_end);

  return (MyMatrix)(this->SubMatrix(i_start, i_end, j_start, j_end));
}


bool 
MyMatrix::SVD(MyColumnVector& D, MyMatrix& U, MyMatrix& V) const
{
  if (columns() > rows()){
    cout << endl << "ERROR" << endl;    
    cout << "Newmat doesn't support svd for columns > rows" << endl;
    cout << "you can avoid this problem by using the LTI library" << endl;
    return false;
  }

  int Acolumns = columns();
  NEWMAT::DiagonalMatrix d(Acolumns);
  NewMatMatrix u(rows(),Acolumns); 
  NewMatMatrix v(Acolumns, Acolumns);
  NEWMAT::SVD(*this,d,u,v);
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

MySymmetricMatrix::SymmetricMatrix() : NewMatSymmetricMatrix() {}
MySymmetricMatrix::SymmetricMatrix(int n) : NewMatSymmetricMatrix(n) {}

// Copy constructor
MySymmetricMatrix::SymmetricMatrix(const SymmetricMatrix& a) : NewMatSymmetricMatrix(a){}
MySymmetricMatrix::SymmetricMatrix(const NewMatSymmetricMatrix & a) : NewMatSymmetricMatrix(a){}

// Destructor
MySymmetricMatrix::~SymmetricMatrix(){}

// Ask Number of Rows and Columns
unsigned int MySymmetricMatrix::rows() const { return this->Nrows();}
unsigned int MySymmetricMatrix::columns() const { return this->Ncols();}


MySymmetricMatrix MySymmetricMatrix::transpose() const {return (*this);}

MySymmetricMatrix MySymmetricMatrix::inverse() const
{
  const NewMatSymmetricMatrix & base = (NewMatSymmetricMatrix &) *this;
  NewMatSymmetricMatrix inverted = base.i();
  return (MySymmetricMatrix) inverted;
}

double MySymmetricMatrix::determinant() const
{
  const NewMatSymmetricMatrix & base = (NewMatSymmetricMatrix &) *this;
  NEWMAT::LogAndSign temp = base.LogDeterminant();
  double result = temp.Value();
  return result;
}


double& MyMatrix::operator()(unsigned int a, unsigned int b) 
{
  NewMatMatrix & op1(*this);
  return op1(a,b);
}

const double MyMatrix::operator()(unsigned int a, unsigned int b) const
{
  const NewMatMatrix& op1(*this);
  return op1(a,b);
}

const bool MyMatrix::operator==(const MyMatrix& a) const
{
  const NewMatMatrix& op1 = *this;
  const NewMatMatrix& op2 = a;
  return (op1 == op2);
}

// Set all elements equal to a
MySymmetricMatrix& MySymmetricMatrix::operator=(const double a)
{
  NewMatSymmetricMatrix temp = (NewMatSymmetricMatrix) *this;
  temp = a;
  *this = (MySymmetricMatrix) temp;

  return *this;
}


// SYMMETRICMATRIX - SCALAR operators
MySymmetricMatrix& MySymmetricMatrix::operator +=(double a)
{
  NewMatSymmetricMatrix & op1 = (*this);
  op1 += a;
  return (MySymmetricMatrix&) op1;
}

MySymmetricMatrix& MySymmetricMatrix::operator -=(double a)
{
  NewMatSymmetricMatrix & op1 = (*this);
  op1 -= a;
  return (MySymmetricMatrix&) op1;
}

MySymmetricMatrix& MySymmetricMatrix::operator *=(double b)
{
  NewMatSymmetricMatrix & op1 = (*this);
  op1 *= b;
  return (MySymmetricMatrix&) op1;
}
  
MySymmetricMatrix& MySymmetricMatrix::operator /=(double b)
{
  NewMatSymmetricMatrix & op1 = (*this);
  op1 /= b;
  return (MySymmetricMatrix&) op1;
}

MySymmetricMatrix MySymmetricMatrix::operator +(double a) const
{
  // make copy
  NewMatSymmetricMatrix op1 = (*this);
  op1 += a;
  return (MySymmetricMatrix) op1;
}

MySymmetricMatrix MySymmetricMatrix::operator -(double a) const
{
  // make copy
  NewMatSymmetricMatrix op1 = (*this);
  op1 -= a;
  return (MySymmetricMatrix) op1;
}

MySymmetricMatrix MySymmetricMatrix::operator *(double b) const
{
  // make copy
  NewMatSymmetricMatrix op1 = (*this);
  op1 *= b;
  return (MySymmetricMatrix) op1;
}
  
MySymmetricMatrix MySymmetricMatrix::operator /(double b) const
{
  // make copy
  NewMatSymmetricMatrix op1 = (*this);
  op1 /= b;
  return (MySymmetricMatrix) op1;
}




// SYMMETRICMATRIX - MATRIX operators
MyMatrix& MySymmetricMatrix::operator +=(const MyMatrix& a)
{
  NewMatSymmetricMatrix & op1 = (*this);
  const NewMatMatrix & op2 = a;
  op1 += op2;
  return (MyMatrix &) op1;
}

MyMatrix& MySymmetricMatrix::operator -=(const MyMatrix& a)
{
  NewMatSymmetricMatrix & op1 = (*this);
  const NewMatMatrix & op2 = a;
  op1 -= op2;
  return (MyMatrix &) op1;
}


MyMatrix MySymmetricMatrix::operator+ (const MyMatrix &a) const
{
  const NewMatMatrix& op1 = (*this);
  const NewMatMatrix& op2 = a;
  NewMatMatrix result = (NewMatMatrix) (op1 + op2);
  return (MyMatrix) result;
}

MyMatrix MySymmetricMatrix::operator- (const MyMatrix &a) const
{
  const NewMatMatrix& op1 = (*this);
  const NewMatMatrix& op2 = a;
  NewMatMatrix result = (NewMatMatrix) (op1 - op2);
  return (MyMatrix) result;
}

MyMatrix MySymmetricMatrix::operator* (const MyMatrix &a) const
{
  const NewMatMatrix& op1 = (*this);
  const NewMatMatrix& op2 = a;
  NewMatMatrix result = (NewMatMatrix) (op1 * op2);
  return (MyMatrix) result;
}



// SYMMETRICMATRIX - SYMMETRICMATRIX operators
MySymmetricMatrix& MySymmetricMatrix::operator +=(const MySymmetricMatrix& a)
{
  NewMatSymmetricMatrix & op1 = (*this);
  const NewMatSymmetricMatrix & op2 = a;
  op1 += op2;
  return (MySymmetricMatrix &) op1;
}

MySymmetricMatrix& MySymmetricMatrix::operator -=(const MySymmetricMatrix& a)
{
  NewMatSymmetricMatrix & op1 = (*this);
  const NewMatSymmetricMatrix & op2 = a;
  op1 -= op2;
  return (MySymmetricMatrix &) op1;
}

MySymmetricMatrix MySymmetricMatrix::operator+ (const MySymmetricMatrix &a) const
{
  const NewMatSymmetricMatrix& op1 = (*this);
  const NewMatSymmetricMatrix& op2 = a;
  NewMatSymmetricMatrix result = (NewMatSymmetricMatrix) (op1 + op2);
  return (MySymmetricMatrix) result;
}

MySymmetricMatrix MySymmetricMatrix::operator- (const MySymmetricMatrix &a) const
{
  const NewMatSymmetricMatrix& op1 = (*this);
  const NewMatSymmetricMatrix& op2 = a;
  NewMatSymmetricMatrix result = (NewMatSymmetricMatrix) (op1 - op2);
  return (MySymmetricMatrix) result;
}

MySymmetricMatrix MySymmetricMatrix::operator* (const MySymmetricMatrix &a) const
{
  const NewMatSymmetricMatrix& op1 = (*this);
  const NewMatSymmetricMatrix& op2 = a;
  return (MySymmetricMatrix) (op1 * op2);
}




MyColumnVector MySymmetricMatrix::operator* (const MyColumnVector &b) const
{
  const NewMatSymmetricMatrix& op1 = *this;
  const NewMatColumnVector& op2 = b;
  return (MyColumnVector) (op1 * op2);
}

MyMatrix MySymmetricMatrix::sub(int i_start, int i_end, 
				int j_start , int j_end) const
{
  return (MyMatrix)(this->SubMatrix(i_start, i_end, j_start, j_end));
}


void
MySymmetricMatrix::resize(unsigned int i, bool copy, bool initialize)
{
  NewMatSymmetricMatrix & temp = (NewMatSymmetricMatrix &) (*this);
  temp.ReSize(i);
}



double& MySymmetricMatrix::operator()(unsigned int a, unsigned int b) 
{
  NewMatSymmetricMatrix & op1 = (*this);
  return op1(a,b);
}

const double MySymmetricMatrix::operator()(unsigned int a, unsigned int b) const
{
  const NewMatSymmetricMatrix op1(*this);
  return op1(a,b);
}

const bool MySymmetricMatrix::operator==(const MySymmetricMatrix& a) const
{
  const NewMatSymmetricMatrix& op1 = *this;
  const NewMatSymmetricMatrix& op2 = a;
  return (op1 == op2);
}


#endif
