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

using namespace std;

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
  return *this;
}

MyMatrix& MyMatrix::operator/= (double a)
{
  BoostMatrix & op1 = (*this);
  op1 /= a;
  return (MyMatrix&) op1;
}

MyMatrix MyMatrix::operator+ (double a) const
{
  return (MyMatrix)(((BoostMatrix)(*this)) + boost::numeric::ublas::scalar_matrix<double>(rows(),columns(),a));
}

MyMatrix MyMatrix::operator- (double a) const
{
  return (MyMatrix)(((BoostMatrix)(*this)) - boost::numeric::ublas::scalar_matrix<double>(rows(),columns(),a));
}

MyMatrix MyMatrix::operator* (double a) const
{
  const BoostMatrix& op1 = (*this);
  return (MyMatrix) (op1 *  a);
}

MyMatrix MyMatrix::operator/ (double a) const
{
  const BoostMatrix& op1 = (*this);
  return (MyMatrix) (op1 /  a);
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
  const BoostMatrix& op1 = *this;
  const BoostMatrix& op2 = a;

  return (MyMatrix)(op1 - op2);
}

MyMatrix MyMatrix::operator+ (const MyMatrix& a) const
{
  const BoostMatrix& op1 = *this;
  const BoostMatrix& op2 = a;

  return (MyMatrix)(op1 + op2);
}

MyMatrix MyMatrix::operator* (const MyMatrix& a) const
{
  const BoostMatrix& op1 = *this;
  const BoostMatrix& op2 = a;

  return (MyMatrix) prod(op1,op2);
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
  const BoostMatrix& op1 = (*this);
  return (MyColumnVector) prod(op1, ((const BoostColumnVector&)b));
}



double& MyMatrix::operator()(unsigned int a, unsigned int b)
{
  BoostMatrix & op1 = (*this);
  return op1(a-1,b-1);
}

const double MyMatrix::operator()(unsigned int a, unsigned int b) const
{
  BoostMatrix  op1(*this);
  return op1(a-1,b-1);
}

const bool MyMatrix::operator==(const MyMatrix& a) const
{
  if (this->rows() != a.rows()) return false;
  if (this->columns() != a.columns()) return false;
  return(norm_inf((BoostMatrix)(*this)-(BoostMatrix)a) == 0);
}


// Set all elements equal to a
MyMatrix&
 MyMatrix::operator=(double a)
{
  *this = (MyMatrix)boost::numeric::ublas::scalar_matrix<double>(rows(),columns(),a);

  return *this;
}


MyRowVector MyMatrix::rowCopy(unsigned int r) const
{
  unsigned int cols = columns();
  BoostRowVector temp(cols);
  for (unsigned int i=0; i<cols; i++)
    temp(i) = (*this)(r,i+1);
  return (MyRowVector) temp;
}

MyColumnVector MyMatrix::columnCopy(unsigned int c) const
{
  unsigned int ro = rows();
  BoostColumnVector temp(ro);
  for (unsigned int i=0; i<ro; i++)
    temp(i) = (*this)(i+1,c);
  return (MyColumnVector) temp;
}




MyMatrix MyMatrix::transpose() const
{
  const BoostMatrix &op1 = (*this);
  return (MyMatrix) trans(op1);
}

double MyMatrix::determinant() const
{
  unsigned int r = this->rows();
  assert(r == this->columns());
  double result = 1.0;
  const BoostMatrix& A = (*this);
  switch (r)
  {
        case 1:
            return A(0,0);
        case 2: 
            return ( ( A(0,0) * A(1,1)) - ( A(1,0) * A(0,1)) );
        default: 
            BoostMatrix LU(r,r);
            boost::numeric::ublas::permutation_matrix<> ndx(r);
            noalias(LU) = A;
            int res = lu_factorize(LU,ndx);
            assert(res == 0);

            int s = 1;
            for (boost::numeric::ublas::matrix<double>::size_type i=0;i<LU.size1();++i) {
              result *= LU(i,i);
              if (ndx(i)!=i) s = -s;
            }
            return result*s;
  }
}


MyMatrix MyMatrix::inverse() const
{
  unsigned int r = this->rows();
  assert(r == this->columns());
  const BoostMatrix& A = (*this);
  BoostMatrix Ai(r,r);
  switch (r) 
  {
     case 1:
     {
        Ai(0,0) = 1/A(0,0);
        break;
     }
     case 2:
     {
       double det = A(0,0)*A(1,1)-A(0,1)*A(1,0);
       Ai(0,0) = A(1,1)/det;
       Ai(1,1) = A(0,0)/det;
       Ai(0,1) = -A(0,1)/det;
       Ai(1,0) = -A(1,0)/det;
       break;
     }
     default:
     {
       BoostMatrix LU(r,r);
       boost::numeric::ublas::permutation_matrix<> ndx(r);
       noalias(LU) = A;
       int res = lu_factorize(LU,ndx);
       assert(res == 0);
       noalias(Ai) = boost::numeric::ublas::identity_matrix<double>(r);
       lu_substitute(LU,ndx,Ai);
       break;
     }
  }
  return Ai;
}


int
MyMatrix::convertToSymmetricMatrix(MySymmetricMatrix& sym)
{
  // test if matrix is square matrix
  assert(this->rows() == this->columns());

  // if necessairy, resize sym
  // only check cols or rows. Symmetric matrix is square.
  if ( sym.rows() != this->rows() )
    sym = MySymmetricMatrix(this->rows());

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
  temp.resize(i,j,copy);
}

// get sub matrix
MyMatrix MyMatrix::sub(int i_start, int i_end, int j_start , int j_end) const
{
  MyMatrix submatrix(i_end-i_start+1, j_end-j_start+1);
  for (int i=i_start; i<=i_end; i++)
    for (int j=j_start; j<=j_end; j++)
      submatrix(i-i_start+1,j-j_start+1) = (*this)(i,j);

  return submatrix;
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
unsigned int MySymmetricMatrix::rows() const { return this->size1();}
unsigned int MySymmetricMatrix::columns() const { return this->size2();}


MySymmetricMatrix MySymmetricMatrix::transpose() const {return (*this);}

MySymmetricMatrix MySymmetricMatrix::inverse() const
{
  unsigned int r = this->rows();
  assert(r == this->columns());
  const BoostMatrix& A = (*this);
  BoostSymmetricMatrix Ai(r,r);
  switch (r) 
  {
     case 1:
     {
        Ai(0,0) = 1/A(0,0);
        break;
     }
     case 2:
     {
       double det = A(0,0)*A(1,1)-A(0,1)*A(1,0);
       Ai(0,0) = A(1,1)/det;
       Ai(1,1) = A(0,0)/det;
       Ai(0,1) = -A(0,1)/det;
       Ai(1,0) = -A(1,0)/det;
       break;
     }
     default:
     {
       BoostSymmetricMatrix LU(r,r);
       boost::numeric::ublas::permutation_matrix<> ndx(r);
       noalias(LU) = A;
       int res = lu_factorize(LU,ndx);
       assert(res == 0);
       noalias(Ai) = boost::numeric::ublas::identity_matrix<double>(r);
       lu_substitute(LU,ndx,Ai);
       break;
     }
  }
  return Ai;
}

double MySymmetricMatrix::determinant() const
{
  unsigned int r = this->rows();
  assert(r == this->columns());
  const BoostMatrix& A = (*this);
  switch (r) 
  {
     case 1:
     {
        return A(0,0);
     }
     case 2:
     {
       return A(0,0)*A(1,1)-A(0,1)*A(1,0);
     }
     default:
     {
        BoostSymmetricMatrix LU(r,r);
        boost::numeric::ublas::permutation_matrix<> ndx(r);
        noalias(LU) = A;
        int res = lu_factorize(LU,ndx);
        assert(res == 0);

        double result = 1.0;
        int s = 1;
        for (boost::numeric::ublas::matrix<double>::size_type i=0;i<LU.size1();++i) {
          result *= LU(i,i);
          if (ndx(i)!=i) s = -s;
        }
        return result*s;
     }
  }
}


// Set all elements equal to a
MySymmetricMatrix& MySymmetricMatrix::operator=(const double a)
{
  *this = (MySymmetricMatrix)boost::numeric::ublas::scalar_matrix<double>(rows(),columns(),a);

  return *this;
}


// SYMMETRICMATRIX - SCALAR operators
MySymmetricMatrix& MySymmetricMatrix::operator +=(double a)
{
  BoostSymmetricMatrix & op1 = *this;
  op1 += boost::numeric::ublas::scalar_matrix<double>(rows(),columns(),a);
  return (MySymmetricMatrix&)op1;
}

MySymmetricMatrix& MySymmetricMatrix::operator -=(double a)
{
  BoostSymmetricMatrix & op1 = *this;
  op1 -= boost::numeric::ublas::scalar_matrix<double>(rows(),columns(),a);
  return (MySymmetricMatrix&)op1;
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
  return (MySymmetricMatrix)(((BoostSymmetricMatrix)(*this)) + boost::numeric::ublas::scalar_matrix<double>(rows(),columns(),a));
}

MySymmetricMatrix MySymmetricMatrix::operator -(double a) const
{
  return (MySymmetricMatrix)(((BoostSymmetricMatrix)(*this)) - boost::numeric::ublas::scalar_matrix<double>(rows(),columns(),a));
}

MySymmetricMatrix MySymmetricMatrix::operator *(double b) const
{
 const BoostSymmetricMatrix& op1 = (*this);
  return (MySymmetricMatrix) (op1 *  b);
}

MySymmetricMatrix MySymmetricMatrix::operator /(double b) const
{
  const BoostSymmetricMatrix& op1 = (*this);
  return (MySymmetricMatrix) (op1 /  b);
}




// SYMMETRICMATRIX - MATRIX operators
MyMatrix& MySymmetricMatrix::operator +=(const MyMatrix& a)
{
  BoostSymmetricMatrix & op1 = (*this);
  op1 += a;
  return (MyMatrix &) op1;
}

MyMatrix& MySymmetricMatrix::operator -=(const MyMatrix& a)
{
  BoostSymmetricMatrix & op1 = (*this);
  op1 -= a;
  return (MyMatrix &) op1;
}


MyMatrix MySymmetricMatrix::operator+ (const MyMatrix &a) const
{
  const BoostSymmetricMatrix& op1 = *this;
  const BoostMatrix& op2 = a;

  return (MyMatrix) (op1 + op2);
}

MyMatrix MySymmetricMatrix::operator- (const MyMatrix &a) const
{
  const BoostSymmetricMatrix& op1 = *this;
  const BoostMatrix& op2 = a;

  return (MyMatrix) (op1 - op2);
}

MyMatrix MySymmetricMatrix::operator* (const MyMatrix &a) const
{
  const BoostSymmetricMatrix& op1 = *this;
  const BoostMatrix& op2 = a;

  return (MyMatrix) prod(op1, op2);
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
  const BoostSymmetricMatrix& op1 = *this;
  const BoostSymmetricMatrix& op2 = a;

  return (MySymmetricMatrix) (op1 + op2);
}

MySymmetricMatrix MySymmetricMatrix::operator- (const MySymmetricMatrix &a) const
{
  const BoostSymmetricMatrix& op1 = *this;
  const BoostSymmetricMatrix& op2 = a;

  return (MySymmetricMatrix) (op1 - op2);
}

MyMatrix MySymmetricMatrix::operator* (const MySymmetricMatrix &a) const
{
  const BoostSymmetricMatrix& op1 = *this;
  const BoostSymmetricMatrix& op2 = a;

  return (MyMatrix) prod(op1, op2);
}




MyColumnVector MySymmetricMatrix::operator* (const MyColumnVector &b) const
{
  const BoostSymmetricMatrix& op1 = (BoostSymmetricMatrix) *this;
  return (MyColumnVector) prod(op1, ((const BoostColumnVector&)b));
}

void MySymmetricMatrix::multiply (const MyColumnVector &b, MyColumnVector &result) const
{
  const BoostSymmetricMatrix& op1 = (BoostSymmetricMatrix) *this;
  result = (MyColumnVector) prod(op1, ((const BoostColumnVector&)b));
}

MyMatrix MySymmetricMatrix::sub(int i_start, int i_end, int j_start , int j_end) const
{
  MyMatrix submatrix(i_end-i_start+1, j_end-j_start+1);
  for (int i=i_start; i<=i_end; i++)
    for (int j=j_start; j<=j_end; j++)
      submatrix(i-i_start+1,j-j_start+1) = (*this)(i,j);

  return submatrix;
}



double& MySymmetricMatrix::operator()(unsigned int a, unsigned int b)
{
  BoostSymmetricMatrix & op1 = (*this);
  return op1(a-1,b-1);
}

const double MySymmetricMatrix::operator()(unsigned int a, unsigned int b) const
{
  BoostSymmetricMatrix op1(*this);
  return op1(a-1,b-1);
}

const bool MySymmetricMatrix::operator==(const MySymmetricMatrix& a) const
{
  if (this->rows() != a.rows()) return false;
  if (this->columns() != a.columns()) return false;
  return(norm_inf((BoostSymmetricMatrix)(*this)-(BoostSymmetricMatrix)a) == 0);
}

void
MySymmetricMatrix::resize(unsigned int i, bool copy, bool initialize)
{
  BoostSymmetricMatrix & temp = (BoostSymmetricMatrix &) (*this);
  temp.resize(i, copy);
}


#endif
