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
#ifndef __OROVECTOR__
#define __OROVECTOR__

#define use_namespace

// We define here the oro_math namespace:  This should be done
// explicitly because of the nesting of the namespaces and the fact
// that the 2 classes can have the same name...
// A vector is nothing more then a matrix with a restriction to the
// number of cols/rows, but it facilitates the readings...

#define MyColumnVector MatrixWrapper::ColumnVector
#define MyRowVector MatrixWrapper::RowVector
#define MyMatrix MatrixWrapper::Matrix

namespace MatrixWrapper{

class Matrix;
class ColumnVector;
class RowVector;

/// Class ColumnVectorWrapper
class ColumnVector_Wrapper
{
public:

  /// Constructor
  ColumnVector_Wrapper() {};

  /// Destructor
  virtual ~ColumnVector_Wrapper() {};

  /// resize
  virtual void resize(int num_rows) = 0;

  /// Ask number of rows
  virtual unsigned int rows() const = 0;

  /// Ask numbers of columns (=1)
  virtual unsigned int columns() const = 0;

  /// join two vectors
  virtual MyColumnVector vectorAdd(const MyColumnVector& v2) const = 0;

  /// element indexing
  virtual double operator()(unsigned int i) const = 0;

  /// element indexing
  virtual double& operator()(unsigned int i) = 0;

  /// element indexing STARTING FROM 0
  virtual double operator[](unsigned int i) const
    { return (*this)(i+1);}

  /// element indexing STARTING FROM 0
  virtual double& operator[](unsigned int i) 
    { return (*this)(i+1);}

  /// Operator ==
  virtual bool operator==(const MyColumnVector& a) const = 0;

  /// operator =
  virtual MyColumnVector& operator =(const MyColumnVector& a) = 0;

  /// Initialise all elements to a
  virtual MyColumnVector& operator =(double a) = 0;



  /// Operators
  virtual MyColumnVector& operator+= (const MyColumnVector& a) = 0;

  /// Operators
  virtual MyColumnVector& operator-= (const MyColumnVector& a) = 0;

  /// Operators
  virtual MyColumnVector operator+ (const MyColumnVector &a) const = 0;

  /// Operators
  virtual MyColumnVector operator- (const MyColumnVector &a) const = 0;



  /// Operators
  virtual MyColumnVector& operator+= (double b) = 0;

  /// Operators
  virtual MyColumnVector& operator-= (double b) = 0;

  /// Operators
  virtual MyColumnVector& operator*= (double b) = 0;

  /// Operators
  virtual MyColumnVector& operator/= (double b) = 0;

  /// Operators
  virtual MyColumnVector operator+ (double b) const = 0;

  /// Operators
  virtual MyColumnVector operator- (double b) const = 0;

  /// Operators
  virtual MyColumnVector operator* (double b) const = 0;

  /// Operators
  virtual MyColumnVector operator/ (double b) const = 0;



  /// Operators
  virtual MyMatrix operator* (const MyRowVector &a) const = 0;

  /// get sub matrix
  virtual MyColumnVector sub(int j_start , int j_end) const = 0;

  /// get transpose
  virtual MyRowVector transpose() const = 0;

}; // class


/// Class RowVectorWrapper
class RowVector_Wrapper
{
public:

  /// Constructor
  RowVector_Wrapper() {};

  /// Destructor
  virtual ~RowVector_Wrapper() {};

  /// resize
  virtual void resize(int num_cols) = 0;

  /// Ask number of rows
  virtual unsigned int rows() const = 0;

  /// Ask numbers of columns (=1)
  virtual unsigned int columns() const = 0;

  /// join two vectors
  virtual MyRowVector vectorAdd(const MyRowVector& v2) const = 0;

  /// element indexing
  virtual double operator()(unsigned int) const = 0;

  /// element indexing
  virtual double& operator()(unsigned int) = 0;

  /// Operator ==
  virtual bool operator==(const MyRowVector& a) const = 0;

  /// operator =
  virtual MyRowVector& operator =(const MyRowVector& a) = 0;

  /// Initialise all elements to a
  virtual MyRowVector& operator =(double a) = 0;



  /// Operators
  virtual MyRowVector& operator+= (const MyRowVector& a) = 0;

  /// Operators
  virtual MyRowVector& operator-= (const MyRowVector& a) = 0;

  /// Operators
  virtual MyRowVector operator+ (const MyRowVector &a) const = 0;

  /// Operators
  virtual MyRowVector operator- (const MyRowVector &a) const = 0;


  /// Operators
  virtual MyRowVector& operator+= (double b) = 0;

  /// Operators
  virtual MyRowVector& operator-= (double b) = 0;

  /// Operators
  virtual MyRowVector& operator*= (double b) = 0;

  /// Operators
  virtual MyRowVector& operator/= (double b) = 0;

  /// Operators
  virtual MyRowVector operator+(double b) const = 0;

  /// Operators
  virtual RowVector operator- (double b) const = 0;

  /// Operators
  virtual MyRowVector operator* (double b) const = 0;

  /// Operators
  virtual RowVector operator/ (double b) const = 0;

  /// Operators
  virtual double operator* (const MyColumnVector &a) const = 0;

  /// get sub matrix
  virtual MyRowVector sub(int j_start , int j_end) const = 0;

  /// get transpose
  virtual MyColumnVector transpose() const = 0;

}; // class



} // namespace

#include "vector_NEWMAT.h"
#include "vector_LTI.h"
#include "vector_BOOST.h"

#endif // __OROVECTOR__
