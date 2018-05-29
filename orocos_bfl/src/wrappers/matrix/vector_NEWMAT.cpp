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
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
//

#include "../config.h"
#ifdef __MATRIXWRAPPER_NEWMAT__

#include "vector_NEWMAT.h"

// Constructors
MyColumnVector::ColumnVector() : NewMatColumnVector() {}
MyColumnVector::ColumnVector(int num_rows) : NewMatColumnVector(num_rows){}
MyColumnVector::ColumnVector(const MyColumnVector& a, const MyColumnVector& b) : NewMatColumnVector(a.rows() + b.rows())
{
  NewMatColumnVector& opl = (*this);

  unsigned int i;

  // copy elements of a to opl
  for (i=0; i< a.rows(); i++)
    opl(i+1) = a(i+1);

  // copy elements of b to opl
  for (i=0; i< b.rows(); i++)
    opl(a.rows() + i+1) = b(i+1);
}

// Destructor
MyColumnVector::~ColumnVector(){}

// Copy constructor
MyColumnVector::ColumnVector(const MyColumnVector& a) :
  NewMatColumnVector(a){}
MyColumnVector::ColumnVector(const NewMatColumnVector & a) :
  NewMatColumnVector(a){}

// Resizing
void MyColumnVector::resize(int num_rows)
{
  NewMatColumnVector & op1 = (*this);
  op1.ReSize(num_rows);
}

// Assign
void MyColumnVector::assign(int num_rows, double value)
{
  NewMatColumnVector & op1 = (*this);
  op1.resize(num_rows);
  for (unsigned int i=0; i<num_rows; i++)
    op1(i+1) = value;
}

// Number of Rows / Cols
unsigned int MyColumnVector::rows() const { return this->Nrows();}
unsigned int MyColumnVector::columns() const { return this->Ncols();}
unsigned int MyColumnVector::capacity() const { return this->Nrows();}

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
  NewMatColumnVector& op1 = *(this);
  return op1(i);
}

const double MyColumnVector::operator()(unsigned int i) const
{
  NewMatColumnVector op1 = (*this);
  return op1(i);
}


// Operators
MyColumnVector & MyColumnVector::operator+= (const MyColumnVector& a)
{
  NewMatColumnVector & op1 = (*this);
  const NewMatColumnVector & op2 = a;
  op1 += op2;
  return (MyColumnVector &) op1 ;
}

MyColumnVector & MyColumnVector::operator-= (const MyColumnVector& a)
{
  NewMatColumnVector & op1 = (*this);
  const NewMatColumnVector & op2 = a;
  op1 -= op2;
  return (MyColumnVector &) op1 ;
}

MyColumnVector MyColumnVector::operator+ (const MyColumnVector &a) const
{
  const NewMatColumnVector& op1 = (*this);
  const NewMatColumnVector & op2 = a;
  NewMatColumnVector result = (NewMatColumnVector) (op1 + op2);
  return (MyColumnVector) result;
}

MyColumnVector MyColumnVector::operator- (const MyColumnVector &a) const
{
  const NewMatColumnVector & op1 = (*this);
  const NewMatColumnVector & op2 = a;
  NewMatColumnVector result = (NewMatColumnVector) (op1 - op2);
  return (MyColumnVector) result;
}



MyColumnVector& MyColumnVector::operator+= (double a)
{
  NewMatColumnVector& op1 = *this;
  op1 += a;
  return (MyColumnVector&) op1;
}

MyColumnVector& MyColumnVector::operator-= (double a)
{
  NewMatColumnVector& op1 = *this;
  op1 -= a;
  return (MyColumnVector&) op1;
}

MyColumnVector& MyColumnVector::operator*= (double a)
{
  NewMatColumnVector& op1 = *this;
  op1 *= a;
  return (MyColumnVector&) op1;
}

MyColumnVector& MyColumnVector::operator/= (double a)
{
  NewMatColumnVector& op1 = *this;
  op1 /= a;
  return (MyColumnVector&) op1;
}


MyColumnVector MyColumnVector::operator+ (double a) const
{
  const NewMatColumnVector & op1 = (*this);
  return (MyColumnVector) (op1 + a);
}

MyColumnVector MyColumnVector::operator- (double a) const
{
  const NewMatColumnVector & op1 = (*this);
  return (MyColumnVector) (op1 - a);
}

MyColumnVector MyColumnVector::operator* (double a) const
{
  const NewMatColumnVector & op1 = (*this);
  return (MyColumnVector) (op1 * a);
}

MyColumnVector MyColumnVector::operator/ (double a) const
{
  const NewMatColumnVector & op1 = (*this);
  return (MyColumnVector) (op1 / a);
}

const bool MyColumnVector::operator==(const MyColumnVector& a) const
{
  const NewMatColumnVector& op1 = *this;
  const NewMatColumnVector& op2 = a;
  return (op1 == op2);
}

MyRowVector MyColumnVector::transpose() const
{
  NewMatColumnVector & base = (NewMatColumnVector &) *this;
  NewMatRowVector transposedbase = base.t();
  return (MyRowVector) transposedbase;
}

MyMatrix MyColumnVector::operator* (const MyRowVector &a) const
{
  const NewMatColumnVector & op1 = (*this);
  const NewMatRowVector & op2 = a;
  NewMatMatrix result = (NewMatMatrix) (op1 * op2);
  return (MyMatrix) result;
}

MyColumnVector&
MyColumnVector::operator=(const MyColumnVector &a)
{
  // Both these implementations result in the same
  NewMatColumnVector * op1; const NewMatColumnVector * op2;
  op1 = this; op2 = &a;
  *op1 = *op2;
  return *(this);
}

MyColumnVector&
MyColumnVector::operator=(double a)
{
  // Both these implementations result in the same
  NewMatColumnVector * op1;
  op1 = this;
  *op1 = a;
  return *this;
}

// KG temporary: to read in files fast
istream& operator >> (istream& is, MyColumnVector& a)
{
  int nr = a.rows();
  int nc = 1;

  if (nr < 1 || nc < 1)
    is.clear (ios::badbit);
  else
    {
      double tmp;
      for (int i = 0; i < nr; i++)
	{
	  is >> tmp;
	  if (is)
	    a(i+1) = tmp; // check if this is ok
	  else
	    goto done;
	}
    }
done:
  return is;
}


MyColumnVector MyColumnVector::sub(int j_start , int j_end) const
{
  return (MyColumnVector)(this->SubMatrix(j_start, j_end, 1, 1));
}

//////////////////////////////////////////////////////////////////////
////////////////////////////// ROWVECTOR /////////////////////////////
//////////////////////////////////////////////////////////////////////

// Constructors
MyRowVector::RowVector() : NewMatRowVector() {}
MyRowVector::RowVector(int num_rows) : NewMatRowVector(num_rows){}

// Destructor
MyRowVector::~RowVector(){}

// Copy constructor
MyRowVector::RowVector(const MyRowVector& a) : NewMatRowVector(a){}
MyRowVector::RowVector(const NewMatRowVector & a) : NewMatRowVector(a){}

// Number of Rows / Cols
unsigned int MyRowVector::rows() const { return this->Nrows();}
unsigned int MyRowVector::columns() const { return this->Ncols();}
unsigned int MyRowVector::capacity() const { return this->Ncols();}

// Operators
MyRowVector & MyRowVector::operator+= (const MyRowVector& a)
{
  NewMatRowVector & op1 = (*this);
  const NewMatRowVector & op2 = a;
  op1 += op2;
  return (MyRowVector &) op1 ;
}

MyRowVector & MyRowVector::operator-= (const MyRowVector& a)
{
  NewMatRowVector & op1 = (*this);
  const NewMatRowVector & op2 = a;
  op1 -= op2;
  return (MyRowVector &) op1 ;
}

MyRowVector MyRowVector::operator+ (const MyRowVector &a) const
{
  const NewMatRowVector & op1 = (*this);
  const NewMatRowVector & op2 = a;
  NewMatRowVector result = (NewMatRowVector) (op1 + op2);
  return (MyRowVector) result;
}

MyRowVector MyRowVector::operator- (const MyRowVector &a) const
{
  const NewMatRowVector & op1 = (*this);
  const NewMatRowVector & op2 = a;
  NewMatRowVector result = (NewMatRowVector) (op1 - op2);
  return (MyRowVector) result;
}




double
MyRowVector::operator*(const MyColumnVector& a) const
{
  assert (this->columns() == a.rows());

  const NewMatRowVector& op1 = (*this);
  const NewMatColumnVector & op2 = a;
  NewMatMatrix matrixresult = op1 * op2;
  double result = matrixresult.AsScalar();
  return result;
}


MyRowVector&
MyRowVector::operator=(const MyRowVector &a)
{
  // Both these implementations result in the same
  NewMatRowVector * op1; const NewMatRowVector * op2;
  op1 = this; op2 = &a;
  *op1 = *op2;

  return *this;
}


MyColumnVector MyRowVector::transpose() const
{
  NewMatRowVector & base = (NewMatRowVector &) *this;
  NewMatColumnVector transposedbase = base.t();
  return (MyColumnVector) transposedbase;
}


MyRowVector MyRowVector::sub(int j_start , int j_end) const
{
  return (MyRowVector)(this->SubMatrix(1, 1, j_start, j_end));
}



MyRowVector& MyRowVector::operator+= (double a)
{
  NewMatRowVector& op1 = *this;
  op1 += a;
  return (MyRowVector&) op1;
}

MyRowVector& MyRowVector::operator-= (double a)
{
  NewMatRowVector& op1 = *this;
  op1 -= a;
  return (MyRowVector&) op1;
}

MyRowVector& MyRowVector::operator*= (double a)
{
  NewMatRowVector& op1 = *this;
  op1 *= a;
  return (MyRowVector&) op1;
}

MyRowVector& MyRowVector::operator/= (double a)
{
  NewMatRowVector& op1 = *this;
  op1 /= a;
  return (MyRowVector&) op1;
}


MyRowVector MyRowVector::operator+ (double a) const
{
  const NewMatRowVector & op1 = (*this);
  return (MyRowVector) (op1 + a);
}

MyRowVector MyRowVector::operator- (double a) const
{
  const NewMatRowVector & op1 = (*this);
  return (MyRowVector) (op1 - a);
}

MyRowVector MyRowVector::operator* (double a) const
{
  const NewMatRowVector & op1 = (*this);
  return (MyRowVector) (op1 * a);
}

MyRowVector MyRowVector::operator/ (double a) const
{
  const NewMatRowVector & op1 = (*this);
  return (MyRowVector) (op1 / a);
}



MyRowVector&
MyRowVector::operator=(double a)
{
  // Both these implementations result in the same
  NewMatRowVector * op1;
  op1 = this;
  *op1 = a;
  return *this;
}


double& MyRowVector::operator()(unsigned int i)
{
  NewMatRowVector& op1 = *(this);
  return  op1(i);
}

const double MyRowVector::operator()(unsigned int i) const
{
  NewMatRowVector op1 = (*this);
  return  op1(i);
}

const bool MyRowVector::operator==(const MyRowVector& a) const
{
  const NewMatRowVector& op1 = *this;
  const NewMatRowVector& op2 = a;
  return (op1 == op2);
}


MyRowVector
MyRowVector::vectorAdd(const MyRowVector& v2) const
{
  const MyRowVector& v1 = *this;
  MyRowVector res(v1.rows() + v2.rows());

  for (unsigned int i=0; i<v1.rows(); i++)
    res(i+1) = v1(i+1);

  for (unsigned int i=0; i<v2.rows(); i++)
    res(v1.rows()+i+1) = v2(i+1);

  return res;
}


// Resizing
void MyRowVector::resize(int num_cols)
{
  NewMatRowVector & op1 = (*this);
  op1.ReSize(num_cols);
}

// Assign
void MyRowVector::assign(int num_columns, double value)
{
  NewMatRowVector & op1 = (*this);
  op1.resize(num_columns);
  for (unsigned int i=0; i<num_columns; i++)
    op1(i+1) = value;
}


#endif
