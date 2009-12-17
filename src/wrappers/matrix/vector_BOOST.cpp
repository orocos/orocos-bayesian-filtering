// $Id: vector_BOOST.cpp 27906 2007-04-27 11:50:53Z wmeeusse $
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

#include "vector_BOOST.h"
#include <iostream>


// Constructors
MyColumnVector::ColumnVector() : BoostColumnVector() {}
MyColumnVector::ColumnVector(int num_rows) : BoostColumnVector(num_rows){}
MyColumnVector::ColumnVector(int num_rows,double value) : BoostColumnVector(num_rows){
  this->assign(boost::numeric::ublas::scalar_vector<double>(num_rows,value));
}
MyColumnVector::ColumnVector(const MyColumnVector& a, const MyColumnVector& b) : BoostColumnVector(a.rows() + b.rows())
{
  BoostColumnVector& opl = (*this);

  unsigned int i;

  // copy elements of a to opl
  for (i=0; i< a.rows(); i++)
    opl(i) = a(i+1);

  // copy elements of b to opl
  for (i=0; i< b.rows(); i++)
    opl(a.rows() + i) = b(i+1);
}

// Destructor
MyColumnVector::~ColumnVector(){}

// Copy constructor
MyColumnVector::ColumnVector(const MyColumnVector& a) :
  BoostColumnVector(a){}
MyColumnVector::ColumnVector(const BoostColumnVector & a) :
  BoostColumnVector(a){}

// Resizing
void MyColumnVector::resize(int num_rows)
{
  BoostColumnVector & op1 = (*this);
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
    res(i+1) = v1(i+1);

  for (unsigned int i=0; i<v2.rows(); i++)
    res(v1.rows()+i+1) = v2(i+1);

  return res;
}

double& MyColumnVector::operator()(unsigned int i)
{
  //std::cout << "(BOOSTVECTOR) operator() called" << std::endl;
  BoostColumnVector& op1 = *(this);
  return op1(i-1);
}

double MyColumnVector::operator()(unsigned int i) const
{
  //std::cout << "(BOOSTVECTOR) operator() called" << std::endl;
  const BoostColumnVector op1 = (*this);
  return op1(i-1);
}


bool MyColumnVector::operator==(const MyColumnVector& a) const
{
  if (this->rows() != a.rows()) return false;
  return(norm_inf((BoostColumnVector)(*this)-(BoostColumnVector)a) == 0);
}

// Operators
MyColumnVector & MyColumnVector::operator+= (const MyColumnVector& a)
{
  BoostColumnVector & op1 = (*this);
  const BoostColumnVector & op2 = a;
  op1 += op2;
  return (MyColumnVector &) op1;
}

MyColumnVector & MyColumnVector::operator-= (const MyColumnVector& a)
{
  BoostColumnVector & op1 = (*this);
  const BoostColumnVector & op2 = a;
  op1 -= op2;
  return (MyColumnVector &) op1;
}

MyColumnVector MyColumnVector::operator+ (const MyColumnVector &a) const
{
  return (MyColumnVector) ((BoostColumnVector)(*this) + (BoostColumnVector)a);
}

MyColumnVector MyColumnVector::operator- (const MyColumnVector &a) const
{
  return (MyColumnVector) ((BoostColumnVector)(*this) - (BoostColumnVector)a);
}



MyColumnVector& MyColumnVector::operator+= (double a)
{
  BoostColumnVector & op1 = *this;
  op1 += boost::numeric::ublas::scalar_vector<double>(rows(),a);
  return (MyColumnVector&)op1;
}

MyColumnVector& MyColumnVector::operator-= (double a)
{
  BoostColumnVector & op1 = *this;
  op1 -= boost::numeric::ublas::scalar_vector<double>(rows(),a);
  return (MyColumnVector&)op1;
}

MyColumnVector& MyColumnVector::operator*= (double a)
{
  BoostColumnVector& op1 = *this;
  op1 *= a;
  return (MyColumnVector&) op1;
}

MyColumnVector& MyColumnVector::operator/= (double a)
{
  BoostColumnVector& op1 = *this;
  op1 /= a;
  return (MyColumnVector&) op1;
}


MyColumnVector MyColumnVector::operator+ (double a) const
{
  return (MyColumnVector)(((BoostColumnVector)(*this)) + boost::numeric::ublas::scalar_vector<double>(rows(),a));
}

MyColumnVector MyColumnVector::operator- (double a) const
{
  return (MyColumnVector)(((BoostColumnVector)(*this)) - boost::numeric::ublas::scalar_vector<double>(rows(),a));
}

MyColumnVector MyColumnVector::operator* (double a) const
{
  const BoostColumnVector & op1 = (*this);
  return (MyColumnVector) (op1 * a);
}

MyColumnVector MyColumnVector::operator/ (double a) const
{
  const BoostColumnVector & op1 = (*this);
  return (MyColumnVector) (op1 / a);
}



MyRowVector MyColumnVector::transpose() const
{
  unsigned int r = this->rows();
  MyRowVector result(r);
  for (unsigned int i=0; i<r; i++)
    result(i+1) = (*this)(i+1);
  return result;
}

MyMatrix MyColumnVector::operator* (const MyRowVector &a) const
{
  unsigned int r = this->rows();
  unsigned int c = a.columns();

  MyMatrix result(r,c);
  for (unsigned int i=0; i<r; i++)
    for (unsigned int j=0; j<c; j++)
      result(i+1,j+1) = (*this)(i+1) * a(j+1);
  return result;
}

MyColumnVector&
MyColumnVector::operator=(const MyColumnVector &a)
{
  BoostColumnVector& op1 = *this;
  op1 = (BoostColumnVector)a;
  return *this;
}

MyColumnVector&
MyColumnVector::operator=(double a)
{
  BoostColumnVector& op1 = *this;
  op1 = boost::numeric::ublas::scalar_vector<double>(this->rows(),a);
  return *this;
}

MyColumnVector MyColumnVector::sub(int j_start , int j_end) const
{
  MyColumnVector subvector(j_end-j_start+1);
  for (int j=j_start; j<=j_end; j++)
    subvector(j-j_start+1) = (*this)(j);

  return subvector;
}




//////////////////////////////////////////////////////////////////////
////////////////////////////// ROWVECTOR /////////////////////////////
//////////////////////////////////////////////////////////////////////

// Constructors
MyRowVector::RowVector() : BoostRowVector() {}
MyRowVector::RowVector(int num_cols) : BoostRowVector(num_cols){}
MyRowVector::RowVector(int num_cols,double value) : BoostRowVector(num_cols){
  this->assign(boost::numeric::ublas::scalar_vector<double>(num_cols,value));
}

// Destructor
MyRowVector::~RowVector(){}

// Copy constructor
MyRowVector::RowVector(const MyRowVector& a) :
  BoostRowVector(a){}
MyRowVector::RowVector(const BoostRowVector & a) :
  BoostRowVector(a){}

// Resizing
void MyRowVector::resize(int num_columns)
{
  BoostRowVector & op1 = (*this);
  op1.resize(num_columns);
}

// Number of Rows / Cols
unsigned int MyRowVector::rows() const { return 1;}
unsigned int MyRowVector::columns() const { return this->size();}

MyRowVector
MyRowVector::vectorAdd(const MyRowVector& v2) const
{
  const MyRowVector& v1 = *this;
  MyRowVector res(v1.columns() + v2.columns());

  for (unsigned int i=0; i<v1.columns(); i++)
    res(i+1) = v1(i+1);

  for (unsigned int i=0; i<v2.columns(); i++)
    res(v1.columns()+i+1) = v2(i+1);

  return res;
}

double& MyRowVector::operator()(unsigned int i)
{
  //std::cout << "(BOOSTVECTOR) operator() called" << std::endl;
  BoostRowVector& op1 = *(this);
  return op1(i-1);
}

double MyRowVector::operator()(unsigned int i) const
{
  //std::cout << "(BOOSTVECTOR) operator() called" << std::endl;
  BoostRowVector op1 = (*this);
  return op1(i-1);
}

bool MyRowVector::operator==(const MyRowVector& a) const
{
  if (this->columns() != a.columns()) return false;
  return(norm_inf((BoostRowVector)(*this)-(BoostRowVector)a) == 0);
}

// Operators
MyRowVector & MyRowVector::operator+= (const MyRowVector& a)
{
  BoostRowVector & op1 = (*this);
  const BoostRowVector & op2 = a;
  op1 += op2;
  return (MyRowVector &) op1;
}

MyRowVector & MyRowVector::operator-= (const MyRowVector& a)
{
  BoostRowVector & op1 = (*this);
  const BoostRowVector & op2 = a;
  op1 -= op2;
  return (MyRowVector &) op1;
}

MyRowVector MyRowVector::operator+ (const MyRowVector &a) const
{
  return (MyRowVector) ((BoostRowVector)(*this) + (BoostRowVector)a);
}

MyRowVector MyRowVector::operator- (const MyRowVector &a) const
{
  return (MyRowVector) ((BoostRowVector)(*this) - (BoostRowVector)a);
}



MyRowVector& MyRowVector::operator+= (double a)
{
  BoostRowVector & op1 = *this;
  op1 += boost::numeric::ublas::scalar_vector<double>(columns(),a);
  return (MyRowVector&)op1;
}

MyRowVector& MyRowVector::operator-= (double a)
{
  BoostRowVector & op1 = *this;
  op1 -= boost::numeric::ublas::scalar_vector<double>(columns(),a);
  return (MyRowVector&)op1;
}

MyRowVector& MyRowVector::operator*= (double a)
{
  BoostRowVector& op1 = *this;
  op1 *= a;
  return (MyRowVector&) op1;
}

MyRowVector& MyRowVector::operator/= (double a)
{
  BoostRowVector& op1 = *this;
  op1 /= a;
  return (MyRowVector&) op1;
}


MyRowVector MyRowVector::operator+ (double a) const
{
  return (MyRowVector)(((BoostRowVector)(*this)) + boost::numeric::ublas::scalar_vector<double>(columns(),a));
}

MyRowVector MyRowVector::operator- (double a) const
{
  return (MyRowVector)(((BoostRowVector)(*this)) - boost::numeric::ublas::scalar_vector<double>(columns(),a));
}

MyRowVector MyRowVector::operator* (double a) const
{
  const BoostRowVector & op1 = (*this);
  return (MyRowVector) (op1 * a);
}

MyRowVector MyRowVector::operator/ (double a) const
{
  const BoostRowVector & op1 = (*this);
  return (MyRowVector) (op1 / a);
}



MyColumnVector MyRowVector::transpose() const
{
  unsigned int c = this->columns();
  MyColumnVector result(c);
  for (unsigned int i=0; i<c; i++)
    result(i+1) = (*this)(i+1);
  return result;
}

double MyRowVector::operator* (const MyColumnVector &a) const
{
  unsigned int r = a.rows();
  unsigned int c = this->columns();

  assert(c == r);

  double result = 0;
  for (unsigned int i=0; i<r; i++)
      result += (*this)(i+1) * a(i+1);
  return result;
}

MyRowVector&
MyRowVector::operator=(const MyRowVector &a)
{
  BoostRowVector& op1 = *this;
  op1 = (BoostRowVector)a;
  return *this;
}

MyRowVector&
MyRowVector::operator=(double a)
{
  BoostRowVector& op1 = *this;
  op1 = boost::numeric::ublas::scalar_vector<double>(this->columns(),a);
  return *this;
}

MyRowVector MyRowVector::sub(int j_start , int j_end) const
{
  MyRowVector subvector(j_end-j_start+1);
  for (int j=j_start; j<=j_end; j++)
    subvector(j-j_start+1) = (*this)(j);

  return subvector;
}

#endif
