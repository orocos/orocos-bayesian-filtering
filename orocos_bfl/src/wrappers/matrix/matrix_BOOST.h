// $Id: matrix_BOOST.h 27912 2007-04-27 14:25:33Z wmeeusse $
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

#ifndef __MATRIX_BOOST__
#define __MATRIX_BOOST__
#include "../../bfl_constants.h"
#include "matrix_wrapper.h"
#include "vector_wrapper.h"
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/symmetric.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/storage.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/matrix_expression.hpp>
#include <assert.h>


typedef boost::numeric::ublas::matrix<double> BoostMatrix;
typedef boost::numeric::ublas::symmetric_matrix<double, boost::numeric::ublas::lower> BoostSymmetricMatrix;

namespace MatrixWrapper
{

/// Implementation of Matrixwrapper using Boost
class Matrix : public BoostMatrix, public Matrix_Wrapper
{
 private: // No private members:  We don't add anything.

 public: // Public Members

  // Constructors
  Matrix();
  Matrix(int m, int n);

  // Destructor
  virtual ~Matrix();

  // Copy constructor
  Matrix (const MyMatrix& a);
  Matrix(const BoostMatrix & a);

  Matrix(int num_rows,const RowVector& v);

  
  virtual unsigned int rows() const;
  virtual unsigned int columns() const;
  virtual double& operator()(unsigned int,unsigned int);
  virtual double operator()(unsigned int,unsigned int) const;
  virtual RowVector operator[](unsigned int)const;

  virtual bool operator==(const MyMatrix& a) const;

  virtual MyMatrix& operator =(double a);

  virtual MyMatrix& operator +=(double a);
  virtual MyMatrix& operator -=(double a);
  virtual MyMatrix& operator *=(double b);
  virtual MyMatrix& operator /=(double b);
  virtual MyMatrix operator+ (double b) const;
  virtual MyMatrix operator- (double b) const;
  virtual MyMatrix operator* (double b) const;
  virtual MyMatrix operator/ (double b) const;

  virtual MyMatrix& operator =(const MySymmetricMatrix& a);
  virtual MyMatrix& operator +=(const MyMatrix& a);
  virtual MyMatrix& operator -=(const MyMatrix& a);
  virtual MyMatrix operator+ (const MyMatrix &a) const;
  virtual MyMatrix operator- (const MyMatrix &a) const;
  virtual MyMatrix operator* (const MyMatrix &a) const;

  virtual MyColumnVector operator* ( const MyColumnVector &b) const;

  virtual MyRowVector rowCopy(unsigned int r) const;
  virtual MyColumnVector columnCopy(unsigned int c) const;

  virtual void resize(unsigned int i, unsigned int j,
		      bool copy=true, bool initialize=true);
  virtual MyMatrix inverse() const;
  virtual MyMatrix transpose() const;
  virtual double determinant() const;
  virtual int convertToSymmetricMatrix(MySymmetricMatrix& sym);
  virtual MyMatrix sub(int i_start, int i_end, int j_start , int j_end) const;

};

class SymmetricMatrix : public BoostSymmetricMatrix, public SymmetricMatrix_Wrapper
{
 private: //

 public: //
  // Constructors
  SymmetricMatrix();
  SymmetricMatrix(int n);

  // Copy constructors
  SymmetricMatrix(const MySymmetricMatrix& a);
  SymmetricMatrix(const BoostSymmetricMatrix & a);

  SymmetricMatrix(int num_rows,const RowVector& v);

  // Destructor
  virtual ~SymmetricMatrix();

  virtual unsigned int rows() const;
  virtual unsigned int columns() const;
  virtual MySymmetricMatrix inverse() const;
  virtual MySymmetricMatrix transpose() const;
  virtual double determinant() const;

  virtual double& operator()(unsigned int,unsigned int);
  virtual double operator()(unsigned int,unsigned int) const;
  virtual RowVector operator[](unsigned int)const;
  virtual bool operator==(const MySymmetricMatrix& a) const;

  virtual MySymmetricMatrix& operator=(double a);

  virtual MySymmetricMatrix& operator +=(double a);
  virtual MySymmetricMatrix& operator -=(double a);
  virtual MySymmetricMatrix& operator *=(double b);
  virtual MySymmetricMatrix& operator /=(double b);
  virtual MySymmetricMatrix  operator + (double b) const;
  virtual MySymmetricMatrix  operator - (double b) const;
  virtual MySymmetricMatrix  operator * (double b) const;
  virtual MySymmetricMatrix  operator / (double b) const;

  virtual MyRowVector rowCopy(unsigned int r) const;

  virtual MyMatrix& operator +=(const MyMatrix& a);
  virtual MyMatrix& operator -=(const MyMatrix& a);
  virtual MyMatrix operator  + (const MyMatrix &a) const;
  virtual MyMatrix operator  - (const MyMatrix &a) const;
  virtual MyMatrix operator  * (const MyMatrix &a) const;

  virtual MySymmetricMatrix& operator +=(const MySymmetricMatrix& a);
  virtual MySymmetricMatrix& operator -=(const MySymmetricMatrix& a);
  virtual MySymmetricMatrix  operator + (const MySymmetricMatrix &a) const;
  virtual MySymmetricMatrix  operator - (const MySymmetricMatrix &a) const;
  virtual MyMatrix  operator * (const MySymmetricMatrix& a) const;

  virtual MyColumnVector operator* (const MyColumnVector &b) const;
  virtual void multiply (const MyColumnVector &b, MyColumnVector &result) const;

  virtual void resize(unsigned int i, bool copy=true, bool initialize=true);
  virtual MyMatrix sub(int i_start, int i_end, int j_start , int j_end) const;

};

}

#endif

#endif
