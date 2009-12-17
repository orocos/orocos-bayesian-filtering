// $Id: vector_BOOST.h 27906 2007-04-27 11:50:53Z wmeeusse $
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

#ifndef __VECTOR_BOOST__
#define __VECTOR_BOOST__

#include "matrix_wrapper.h"
#include "vector_wrapper.h"
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>


#define BoostColumnVector boost::numeric::ublas::vector<double>
#define BoostRowVector    boost::numeric::ublas::vector<double>


namespace MatrixWrapper
{

/// Wrapper class for ColumnVectors (Boost implementation)
class ColumnVector : public BoostColumnVector, public ColumnVector_Wrapper
{
public:

  /// Constructor
  ColumnVector();

  /// Constructor
  ColumnVector(int nrows);
  ColumnVector(int nrows,double value);

  /// Constructor
  ColumnVector(const MyColumnVector& a, const MyColumnVector& b);

  /// Destructor
  virtual ~ColumnVector();

  /// Copy constructor
  ColumnVector (const MyColumnVector& a);

  /// Copy constructor
  ColumnVector (const BoostColumnVector& a);

  virtual void resize(int num_rows);
  virtual unsigned int rows() const;
  virtual unsigned int columns() const;
  virtual ColumnVector vectorAdd(const MyColumnVector& v2) const;
  virtual ColumnVector& operator =(const MyColumnVector& a);
  virtual ColumnVector& operator =(double a);

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

  /// element indexing STARTING FROM 0
  virtual double operator[](unsigned int i) const
  { 
  //std::cout << "(BOOSTVECTOR) operator[] called " << i << std::endl;
   // if (i==0)
   //     std::cout << "(BOOSTVECTOR) operator[0]" << std::endl;
    
   return (*this)(i+1);
  }

  /// element indexing STARTING FROM 0
  virtual double& operator[](unsigned int i) 
  { 
  //std::cout << "(BOOSTVECTOR) operator[] called " << i << std::endl;
  //  if (i==0)
  //      std::cout << "(BOOSTVECTOR) operator[0]" << std::endl;
     return (*this)(i+1);
  }

  virtual double operator()(unsigned int) const;
  virtual bool operator==(const MyColumnVector& a) const;
  virtual double& operator()(unsigned int);
  virtual MyMatrix operator* (const MyRowVector &a) const;
  virtual MyColumnVector sub(int j_start , int j_end) const;
  virtual MyRowVector transpose() const;


};

/// Wrapper class for RowVectors (Boost implementation)
class RowVector : public BoostRowVector, public RowVector_Wrapper
{
  // No private member:  We don't add anything.

  // Public Members
 public:
  RowVector();
  RowVector(int ncols);
  RowVector(int ncols,double value);
  // If you have another constructor in the matrix library you
  // want to use, you'll have to redefine it yourself

  // Copy constructor
  RowVector (const MyRowVector& a);
  // Copy constructor for boost
  RowVector (const BoostRowVector& a);

  virtual ~RowVector();

  virtual void resize(int num_cols);
  virtual RowVector vectorAdd(const MyRowVector& v2) const;
  virtual unsigned int rows() const;
  virtual unsigned int columns() const;
  virtual RowVector& operator =(double a);
  virtual RowVector& operator =(const MyRowVector& a);

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

  virtual double operator()(unsigned int) const;
  virtual bool operator==(const MyRowVector& a) const;
  virtual double& operator()(unsigned int);
  virtual MyRowVector sub(int j_start , int j_end) const;
  virtual MyColumnVector transpose() const;
  virtual double operator*(const MyColumnVector& a) const;

};

}

#endif

#endif
