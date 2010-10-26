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

#ifndef __VECTOR_NEWMAT__
#define __VECTOR_NEWMAT__

#include "matrix_wrapper.h"
#include "vector_wrapper.h"
#include <newmat/newmatio.h>
#include <assert.h>


#define NewMatColumnVector NEWMAT::ColumnVector
#define NewMatRowVector    NEWMAT::RowVector

namespace MatrixWrapper
{

/// Wrapper class for ColumnVectors (Newmat implementation)
class ColumnVector : public NewMatColumnVector, public ColumnVector_Wrapper
{
public:

  /// Constructor
  ColumnVector();

  /// Constructor
  ColumnVector(int nrows);

  /// Constructor
  ColumnVector(const MyColumnVector& a, const MyColumnVector& b);

  /// Destructor
  virtual ~ColumnVector();

  /// Copy constructor
  ColumnVector (const MyColumnVector& a);

  /// Copy constructor
  ColumnVector (const NewMatColumnVector& a);

  virtual void resize(int num_rows);
  virtual void assign(int size, double value) ;
  virtual unsigned int rows() const;
  virtual unsigned int columns() const;
  virtual unsigned int capacity() const;
  virtual ColumnVector vectorAdd(const MyColumnVector& v2) const;
  virtual ColumnVector& operator =(const MyColumnVector& a);
  virtual ColumnVector& operator =(double a);
  virtual const bool operator==(const MyColumnVector& a) const;

  virtual MyColumnVector & operator+= (const MyColumnVector& a);
  virtual MyColumnVector & operator-= (const MyColumnVector& a);
  virtual MyColumnVector operator+ (const MyColumnVector &a) const;
  virtual MyColumnVector operator- (const MyColumnVector &a) const;

  virtual MyColumnVector& operator+= (double b);
  virtual MyColumnVector& operator-= (double b);
  virtual MyColumnVector& operator*= (double b);
  virtual MyColumnVector& operator/= (double b);
  virtual MyColumnVector operator+ (double b) const;
  virtual MyColumnVector operator- (double b) const;
  virtual MyColumnVector operator* (double b) const;
  virtual MyColumnVector operator/ (double b) const;

  virtual const double operator()(unsigned int) const;
  virtual double& operator()(unsigned int);
  virtual MyMatrix operator* (const MyRowVector &a) const;
  virtual MyColumnVector sub(int j_start , int j_end) const;
  virtual MyRowVector transpose() const;


};

/// Wrapper class for RowVectors (Newmat implementation)
class RowVector : public NewMatRowVector, public RowVector_Wrapper
{
  // No private member:  We don't add anything.

  // Public Members
 public:
  RowVector();
  RowVector(int ncols);
  // If you have another constructor in the matrix library you
  // want to use, you'll have to redefine it yourself

  // Copy constructor
  RowVector (const MyRowVector& a);
  // Copy constructor for newmat
  RowVector (const NewMatRowVector& a);

  virtual ~RowVector();

  virtual void resize(int num_cols);
  virtual void assign(int size, double value) ;
  virtual RowVector vectorAdd(const MyRowVector& v2) const;
  virtual unsigned int rows() const;
  virtual unsigned int columns() const;
  virtual unsigned int capacity() const;
  virtual RowVector& operator =(double a);
  virtual RowVector& operator =(const MyRowVector& a);
  virtual const bool operator==(const MyRowVector& a) const;

  virtual MyRowVector & operator+= (const MyRowVector& a);
  virtual MyRowVector & operator-= (const MyRowVector& a);
  virtual MyRowVector operator+ (const MyRowVector &a) const;
  virtual MyRowVector operator- (const MyRowVector &a) const;

  virtual MyRowVector& operator+= (double b);
  virtual MyRowVector& operator-= (double b);
  virtual MyRowVector& operator*= (double b);
  virtual MyRowVector& operator/= (double b);
  virtual MyRowVector operator+ (double b) const;
  virtual MyRowVector operator- (double b) const;
  virtual MyRowVector operator* (double b) const;
  virtual MyRowVector operator/ (double b) const;

  virtual const double operator()(unsigned int) const;
  virtual double& operator()(unsigned int);
  virtual MyRowVector sub(int j_start , int j_end) const;
  virtual MyColumnVector transpose() const;
  virtual double operator*(const MyColumnVector& a) const;

};

}

#endif

#endif
