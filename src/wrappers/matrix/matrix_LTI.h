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

#ifndef __MATRIX_LTI__
#define __MATRIX_LTI__

#include "matrix_wrapper.h"
#include "vector_wrapper.h"
#include <ltilib/ltiMatrix.h> 
#include <assert.h>



//  To avoid problems and clarify things 
//  with namespaces, we use some macros
#define ltiMatrix          lti::matrix<double>
#define ltiSymmetricMatrix lti::matrix<double>
#define ltiColumnVector    lti::vector<double>

/// Wrapper class for Matrix using LTI
class MyMatrix : public ltiMatrix, public MatrixWrapper::Matrix_Wrapper
{
 private: // No private members:  We don't add anything.
      
 public: // Public Members

  /// Constructor
  Matrix();

  /// Constructor
  Matrix(int m, int n);

  /// Destructor
  virtual ~Matrix();

  /// Copy constructor
  Matrix (const MyMatrix& a);

  /// Copy constructor
  Matrix(const ltiMatrix & a);

  virtual unsigned int rows() const;
  virtual unsigned int columns() const;
  virtual double& operator()(unsigned int,unsigned int);
  virtual const double operator()(unsigned int,unsigned int) const;

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
  virtual MyMatrix pseudoinverse(double epsilon=0.01) const;
  virtual MyMatrix inverse() const;
  virtual MyMatrix transpose() const;
  virtual double determinant() const;
  virtual int convertToSymmetricMatrix(MySymmetricMatrix& sym);
  virtual MyMatrix sub(int i_start, int i_end, int j_start , int j_end) const;
  virtual bool SVD(MyColumnVector& D, MyMatrix& U, MyMatrix& V) const;

};



class MySymmetricMatrix : public ltiSymmetricMatrix, public SymmetricMatrix_Wrapper
{
    public: 
        // Constructors
        SymmetricMatrix();
        SymmetricMatrix(int n);

        // Copy constructors
        SymmetricMatrix(const SymmetricMatrix& a);
        SymmetricMatrix(const ltiSymmetricMatrix & a);

        // Destructor
        virtual ~SymmetricMatrix();


        virtual unsigned int rows() const;
        virtual unsigned int columns() const;
        virtual MySymmetricMatrix inverse() const;
        virtual MySymmetricMatrix transpose() const;
        virtual double determinant() const;

        virtual double& operator()(unsigned int,unsigned int);
        virtual const double operator()(unsigned int,unsigned int) const;

        virtual MySymmetricMatrix& operator=(double a);

        virtual MySymmetricMatrix& operator +=(double a);
        virtual MySymmetricMatrix& operator -=(double a);
        virtual MySymmetricMatrix& operator *=(double b);
        virtual MySymmetricMatrix& operator /=(double b);
        virtual MySymmetricMatrix  operator + (double b) const;
        virtual MySymmetricMatrix  operator - (double b) const;
        virtual MySymmetricMatrix  operator * (double b) const;
        virtual MySymmetricMatrix  operator / (double b) const;

        virtual MyMatrix& operator +=(const MyMatrix& a);
        virtual MyMatrix& operator -=(const MyMatrix& a);
        virtual MyMatrix operator  + (const MyMatrix &a) const;
        virtual MyMatrix operator  - (const MyMatrix &a) const;
        virtual MyMatrix operator  * (const MyMatrix &a) const;

        virtual MySymmetricMatrix& operator +=(const MySymmetricMatrix& a);
        virtual MySymmetricMatrix& operator -=(const MySymmetricMatrix& a);
        virtual MySymmetricMatrix  operator + (const MySymmetricMatrix &a) const;
        virtual MySymmetricMatrix  operator - (const MySymmetricMatrix &a) const;
        virtual MySymmetricMatrix  operator * (const MySymmetricMatrix& a) const;
      
        virtual MyColumnVector operator* (const MyColumnVector &b) const;
      
        virtual void resize(unsigned int i, bool copy=true, bool initialize=true);
        virtual MyMatrix sub(int i_start, int i_end, int j_start , int j_end) const;
        virtual bool cholesky(MyMatrix& m) const;

};


#endif

#endif
