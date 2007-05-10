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
#ifdef __MATRIXWRAPPER_LTI__

#include "vector_LTI.h"


using std::cout;
using std::endl;


// Constructors
MyColumnVector::ColumnVector() : ltiColumnVector() {}
MyColumnVector::ColumnVector(int num_rows) : ltiColumnVector(num_rows){}
MyColumnVector::ColumnVector(const MyColumnVector& a, 
			     const MyColumnVector& b) 
  : ltiColumnVector(a.rows() + b.rows()) 
{
  ltiColumnVector& opl = (*this);
  
  unsigned int i;
  
  // copy elements of a to opl
  for (i=0; i< a.rows(); i++)
    opl.at(i) = a(i+1);

  // copy elements of b to opl
  for (i=0; i< b.rows(); i++)
    opl.at(a.rows() + i) = b(i+1);
}

// Create from matrix
MyColumnVector::ColumnVector(const MyMatrix &a) 
: ltiColumnVector(a.rows())
{
  assert(a.columns() == 1);

  // copy elements of a to opl
  ltiColumnVector& opl = (*this);
  for (unsigned int i=0; i< a.rows(); i++)
    {
      opl.at(i) = a(i+1,1);
    }
}

// Destructor
MyColumnVector::~ColumnVector(){}

// Copy constructor
MyColumnVector::ColumnVector(const MyColumnVector& a) : 
  ltiColumnVector(a){}
MyColumnVector::ColumnVector(const ltiColumnVector & a) : 
  ltiColumnVector(a){}

// Resizing
void MyColumnVector::resize(int num_rows)
{
  ltiColumnVector & op1 = (*this);
  op1.resize(num_rows);
}

// Number of Rows / Cols
unsigned int MyColumnVector::rows() const { return this->size();}
unsigned int MyColumnVector::columns() const { return 1;}

MyColumnVector 
MyColumnVector::vectorAdd(const MyColumnVector& v2) const
{
  const MyColumnVector& v1 = *this;
  MyColumnVector res(v1.rows() + v2.rows());
  
  for (unsigned int i=0; i<v1.rows(); i++)
    res(i) = v1(i);

  for (unsigned int i=0; i<v2.rows(); i++)
    res(v1.rows()+i) = v2(i);
 
  return res;
}

const bool MyColumnVector::operator==(const MyColumnVector& a) const
{
  const ltiColumnVector& op1 = *this;
  const ltiColumnVector& op2 = a;
  return (op1 == op2);
}


// Operators
MyColumnVector & MyColumnVector::operator+= (const MyColumnVector& a)
{
  ltiColumnVector & op1 = (*this);
  const ltiColumnVector & op2 = a;
  op1 += op2;
  return (MyColumnVector &) op1 ;
}

MyColumnVector & MyColumnVector::operator-= (const MyColumnVector& a)
{
  ltiColumnVector & op1 = (*this);
  const ltiColumnVector & op2 = a;
  op1 -= op2;
  return (MyColumnVector &) op1 ;
}

MyColumnVector MyColumnVector::operator+ (const MyColumnVector &a) const
{
  ltiColumnVector op1(*this);
  const ltiColumnVector & op2 = a;
  ltiColumnVector result = (ltiColumnVector) (op1.add(op2));
  return (MyColumnVector) result;
}

MyColumnVector MyColumnVector::operator- (const MyColumnVector &a) const
{
  ltiColumnVector op1(*this);
  const ltiColumnVector & op2 = a;
  ltiColumnVector result = (ltiColumnVector) (op1.subtract(op2));
  return (MyColumnVector) result;
}



MyColumnVector& MyColumnVector::operator+= (double a)
{
  ltiColumnVector& op1 = *this;
  op1 += a;
  return (MyColumnVector&) op1;
}

MyColumnVector& MyColumnVector::operator-= (double a)
{
  ltiColumnVector& op1 = *this;
  op1 -= a;
  return (MyColumnVector&) op1;
}


MyColumnVector& MyColumnVector::operator*= (double a)
{
  ltiColumnVector& op1 = *this;
  op1 *= a;
  return (MyColumnVector&) op1;
}

MyColumnVector& MyColumnVector::operator/= (double a)
{
  ltiColumnVector& op1 = *this;
  op1 /= a;
  return (MyColumnVector&) op1;
}


MyColumnVector MyColumnVector::operator+ (double a) const
{
  ltiColumnVector op1(*this);
  op1 += a;
  return (MyColumnVector&) op1;
}


MyColumnVector MyColumnVector::operator- (double a) const
{
  ltiColumnVector op1(*this);
  op1 -= a;
  return (MyColumnVector&) op1;
}

MyColumnVector MyColumnVector::operator* (double a) const
{
  ltiColumnVector op1(*this);
  return (MyColumnVector) (op1.multiply(a));
}


MyColumnVector MyColumnVector::operator/ (double a) const
{
  ltiColumnVector op1(*this);
  return (MyColumnVector) (op1.divide(a));
}







MyRowVector MyColumnVector::transpose() const
{
  ltiRowVector transposedbase((*this));
  return (MyRowVector) transposedbase;
}

MyMatrix MyColumnVector::operator* (const MyRowVector &a) const
{
  const ltiColumnVector& op1 = *this;
  const ltiRowVector & op2 = a;
  ltiMatrix result(op1.size(),op2.size());
  result.outerProduct(op1,op2);
  return (MyMatrix) result;
}

MyColumnVector& MyColumnVector::operator=(const MyColumnVector &a)
{
  // Both these implementations result in the same
  ltiColumnVector * op1; const ltiColumnVector * op2;
  op1 = this; op2 = &a;
  *op1 = *op2;

  return *this;
}

MyColumnVector& MyColumnVector::operator=(const double a)
{
  // Both these implementations result in the same
  ltiColumnVector * op1;
  op1 = this;
  (*op1).fill(a,0,(*op1).size());

  return *this;
}

double &MyColumnVector::operator()(unsigned int i)
{
  assert (i != 0);
  //ltiColumnVector * op1;
  //op1 = this;
  return this->at(i-1);
}

const double MyColumnVector::operator()(unsigned int i) const 
{
  assert (i != 0);
  //ltiColumnVector op1(*this);
  return this->at(i-1);
}


MyColumnVector MyColumnVector::sub(int j_start , int j_end) const
{
  ltiColumnVector result(*this,j_start-1,j_end-1);
  return (MyColumnVector) result;
}

//////////////////////////////////////////////////////////////////////
////////////////////////////// ROWVECTOR /////////////////////////////
//////////////////////////////////////////////////////////////////////

// Constructors
MyRowVector::RowVector() : ltiRowVector() {}
MyRowVector::RowVector(int num_rows) : ltiRowVector(num_rows){}

// Destructor
MyRowVector::~RowVector(){}

// Copy constructor
MyRowVector::RowVector(const MyRowVector& a) : ltiRowVector(a){}
MyRowVector::RowVector(const ltiRowVector & a) : ltiRowVector(a){}

// Number of Rows / Cols
unsigned int MyRowVector::rows() const { return 1;}
unsigned int MyRowVector::columns() const { return this->size();}

const bool MyRowVector::operator==(const MyRowVector& a) const
{
  const ltiRowVector& op1 = *this;
  const ltiRowVector& op2 = a;
  return (op1 == op2);
}


// Operators
MyRowVector & MyRowVector::operator+= (const MyRowVector& a)
{
  ltiRowVector & op1 = (*this);
  const ltiRowVector & op2 = a;
  op1 += op2;
  return (MyRowVector &) op1 ;
}

MyRowVector & MyRowVector::operator-= (const MyRowVector& a)
{
  ltiRowVector & op1 = (*this);
  const ltiRowVector & op2 = a;
  op1 -= op2;
  return (MyRowVector &) op1 ;
}

MyRowVector MyRowVector::operator+ (const MyRowVector &a) const
{
  ltiRowVector op1(*this);
  const ltiRowVector & op2 = a;
  ltiRowVector result = (ltiRowVector) (op1.add(op2));
  return (MyRowVector) result;
}

MyRowVector MyRowVector::operator- (const MyRowVector &a) const
{
  ltiRowVector op1(*this);
  const ltiRowVector & op2 = a;
  ltiRowVector result = (ltiRowVector) (op1.subtract(op2));
  return (MyRowVector) result;
}






MyRowVector& MyRowVector::operator+= (double a)
{
  ltiRowVector& op1 = *this;
  op1 += a;
  return (MyRowVector&) op1;
}

MyRowVector& MyRowVector::operator-= (double a)
{
  ltiRowVector& op1 = *this;
  op1 -= a;
  return (MyRowVector&) op1;
}


MyRowVector& MyRowVector::operator*= (double a)
{
  ltiRowVector& op1 = *this;
  op1 *= a;
  return (MyRowVector&) op1;
}

MyRowVector& MyRowVector::operator/= (double a)
{
  ltiRowVector& op1 = *this;
  op1 /= a;
  return (MyRowVector&) op1;
}


MyRowVector MyRowVector::operator+ (double a) const
{
  ltiRowVector op1(*this);
  op1 += a;
  return (MyRowVector&) op1;
}


MyRowVector MyRowVector::operator- (double a) const
{
  ltiRowVector op1(*this);
  op1 -= a;
  return (MyRowVector&) op1;
}

MyRowVector MyRowVector::operator* (double a) const
{
  ltiRowVector op1(*this);
  return (MyRowVector) (op1.multiply(a));
}


MyRowVector MyRowVector::operator/ (double a) const
{
  ltiRowVector op1(*this);
  return (MyRowVector) (op1.divide(a));
}




MyRowVector MyRowVector::operator*(const MyMatrix& a)
{
  assert(this->columns() == a.rows());

  ltiRowVector & op1 = (*this);
  const ltiMatrix & op2 = a;
  ltiRowVector result = op2.leftMultiply(op1);
  return (MyRowVector) result;
}


double MyRowVector::operator*(const MyColumnVector& a) const
{
  assert(this->columns() == a.rows());

  const ltiRowVector & op1 = (*this);
  const ltiColumnVector & op2 = a;
  double result=op1.dot(op2);
  return result;
}

MyRowVector& MyRowVector::operator=(const MyRowVector &a)
{
  // Both these implementations result in the same
  ltiRowVector * op1; 
  const ltiRowVector * op2;
  op1 = this; 
  op2 = &a;
  *op1 = *op2;

  return *this;
}


MyRowVector& MyRowVector::operator=(const double a)
{
  // Both these implementations result in the same
  ltiRowVector * op1;
  op1 = this;
  (*op1).fill(a,0,(*op1).size());

  return *this;
}


MyColumnVector MyRowVector::transpose() const
{
  ltiColumnVector transposedbase((*this));
  return (MyColumnVector) transposedbase;
}


double &MyRowVector::operator()(unsigned int i)
{
  assert (i != 0);
  ltiRowVector * op1;
  op1 = this;
  return op1->at(i-1);
}

const double MyRowVector::operator()(unsigned int i) const 
{
  assert (i != 0);
  ltiRowVector op1(*this);
  return ((op1.at(i-1)));
}



MyRowVector MyRowVector::sub(int j_start , int j_end) const
{
  ltiRowVector result(*this,j_start-1,j_end-1);
  return (MyRowVector) result;
}

// Resizing
void MyRowVector::resize(int num_cols)
{
  ltiRowVector & op1 = (*this);
  op1.resize(num_cols);
}



MyRowVector 
MyRowVector::vectorAdd(const MyRowVector& v2) const
{
  const MyRowVector& v1 = *this;
  MyRowVector res(v1.columns() + v2.columns());
  
  for (unsigned int i=0; i<v1.columns(); i++)
    res(i) = v1(i);

  for (unsigned int i=0; i<v2.columns(); i++)
    res(v1.columns()+i) = v2(i);
 
  return res;
}

#endif
