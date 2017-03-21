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
#ifndef __MATRIX_WRAPPER__
#define __MATRIX_WRAPPER__



// define
#define use_namespace
#define MyMatrix          MatrixWrapper::Matrix
#define MyColumnVector    MatrixWrapper::ColumnVector
#define MyRowVector       MatrixWrapper::RowVector
#define MySymmetricMatrix MatrixWrapper::SymmetricMatrix

namespace MatrixWrapper{

class Matrix;
class ColumnVector;
class RowVector;
class SymmetricMatrix;

/// Class Matrixwrapper
class Matrix_Wrapper
{
public:

  /// Constructor
  Matrix_Wrapper() {};

  /// Destructor
  virtual ~Matrix_Wrapper() {};

  /// Ask Number of Rows
  virtual unsigned int size() const = 0;

  /// Ask Number of Rows
  virtual unsigned int capacity() const = 0;

  /// Ask Number of Rows
  virtual unsigned int rows() const = 0;

  /// Ask Number of Columns
  virtual unsigned int columns() const = 0;

  /// Operator ()
  virtual double& operator()(unsigned int,unsigned int) = 0;

  /// Operator ()
  virtual double operator()(unsigned int,unsigned int) const = 0;

  /// Operator ==
  virtual bool operator==(const MyMatrix& a) const = 0;


  /// Set all elements of the Matrix equal to a
  virtual MyMatrix& operator =(double a) = 0;



  /// MATRIX - SCALAR operator
  virtual MyMatrix& operator +=(double a) = 0;

  /// MATRIX - SCALAR operator
  virtual MyMatrix& operator -=(double a) = 0;

  /// MATRIX - SCALAR operator
  virtual MyMatrix& operator *=(double b) = 0;

  /// MATRIX - SCALAR operator
  virtual MyMatrix& operator /=(double b) = 0;

  /// MATRIX - SCALAR operator
  virtual MyMatrix operator+ (double b) const = 0;

  /// MATRIX - SCALAR operator
  virtual MyMatrix operator- (double b) const = 0;

  /// MATRIX - SCALAR operator
  virtual MyMatrix operator* (double b) const = 0;

  /// MATRIX - SCALAR operator
  virtual MyMatrix operator/ (double b) const = 0;

  /// MATRIX - SYMMETRICMATRIX operators
  virtual MyMatrix& operator =(const MySymmetricMatrix& a) = 0;

  /// MATRIX - MATRIX operator
  virtual MyMatrix& operator +=(const MyMatrix& a) = 0;

  /// MATRIX - MATRIX operator
  virtual MyMatrix& operator -=(const MyMatrix& a) = 0;

  /// MATRIX - MATRIX operator
  virtual MyMatrix operator+ (const MyMatrix &a) const = 0;

  /// MATRIX - MATRIX operator
  virtual MyMatrix operator- (const MyMatrix &a) const = 0;

  /// MATRIX - MATRIX operator
  virtual MyMatrix operator* (const MyMatrix &a) const = 0;



  /// MATRIX - VECTOR operator
  virtual MyColumnVector operator* ( const MyColumnVector &b) const = 0;


  /// Get row from matrix
  virtual MyRowVector rowCopy(unsigned int r) const = 0;

  /// Get column from matrix
  virtual MyColumnVector columnCopy(unsigned int c) const = 0;

  /// resize matrix
  virtual void resize(unsigned int i, unsigned int j,
		      bool copy=true, bool initialize=true) = 0;

  /// get pseudoinverse
  virtual MyMatrix pseudoinverse(double epsilon = 0.01 ) const;

  /// get inverse
  virtual MyMatrix inverse() const = 0;

  /// get transpose
  virtual MyMatrix transpose() const = 0;

  /// get determinant
  virtual double determinant() const = 0;

  /// Turn matrix into Symmetric one
  /** Convert Matrix to SymmetricMatrix
      Elements of matrix are copied to lower triangle of new symmetric matrix
  */
  virtual int convertToSymmetricMatrix(MySymmetricMatrix& sym) = 0;

  /// get sub matrix
  virtual MyMatrix sub(int i_start, int i_end, int j_start , int j_end) const = 0;

  /// SVD Decomposition (for pseudo-inverse properties)
  virtual bool SVD(MyColumnVector& D, MyMatrix& U, MyMatrix& V) const ;

  double PYTHAG(double a,double b) const;

  double SIGN(double a,double b) const;

}; // class Matrix_Wrapper


/// Class SymmetricMatrixWrapper
class SymmetricMatrix_Wrapper
{
public:
  /// Constructor
  SymmetricMatrix_Wrapper() {};

  /// Destructor
  virtual ~SymmetricMatrix_Wrapper() {};

  /// Ask Number of Rows
  virtual unsigned int size() const = 0;

  /// Ask Number of Rows
  virtual unsigned int capacity() const = 0;


  /// Ask Number of Rows
  virtual unsigned int rows() const = 0;

  /// Ask Number of Columns
  virtual unsigned int columns() const = 0;

  /// Operator ()
  virtual double& operator()(unsigned int,unsigned int) = 0;

  /// Operator ()
  virtual double operator()(unsigned int,unsigned int) const = 0;

  /// Operator ==
  virtual bool operator==(const MySymmetricMatrix& a) const = 0;

  /// Set all elements of the Matrix equal to a
  virtual MySymmetricMatrix& operator =(double a) = 0;



  /// SYMMETRICMATRIX - SCALAR operator
  virtual MySymmetricMatrix& operator +=(double a) = 0;

  /// SYMMETRICMATRIX - SCALAR operator
  virtual MySymmetricMatrix& operator -=(double a) = 0;

  /// SYMMETRICMATRIX - SCALAR operator
  virtual MySymmetricMatrix& operator *=(double b) = 0;

  /// SYMMETRICMATRIX - SCALAR operator
  virtual MySymmetricMatrix& operator /=(double b) = 0;

  /// SYMMETRICMATRIX - SCALAR operator
  virtual MySymmetricMatrix operator+ (double b) const = 0;

  /// SYMMETRICMATRIX - SCALAR operator
  virtual MySymmetricMatrix operator- (double b) const = 0;

  /// SYMMETRICMATRIX - SCALAR operator
  virtual MySymmetricMatrix operator* (double b) const = 0;

  /// SYMMETRICMATRIX - SCALAR operator
  virtual MySymmetricMatrix operator/ (double b) const = 0;


  /// SYMMETRICMATRIX - MATRIX operator
  virtual MyMatrix& operator +=(const MyMatrix& a) = 0;

  /// SYMMETRICMATRIX - MATRIX operator
  virtual MyMatrix& operator -=(const MyMatrix& a) = 0;

  /// SYMMETRICMATRIX - MATRIX operator
  virtual MyMatrix operator+ (const MyMatrix &a) const = 0;

  /// SYMMETRICMATRIX - MATRIX operator
  virtual MyMatrix operator- (const MyMatrix &a) const = 0;

  /// SYMMETRICMATRIX - MATRIX operator
  virtual MyMatrix operator* (const MyMatrix &a) const = 0;

  /// SYMMETRICMATRIX - SYMMETRICMATRIX operators
  virtual MySymmetricMatrix& operator +=(const MySymmetricMatrix& a) = 0;

  /// SYMMETRICMATRIX - SYMMETRICMATRIX operators
  virtual MySymmetricMatrix& operator -=(const MySymmetricMatrix& a) = 0;

  /// SYMMETRICMATRIX - SYMMETRICMATRIX operators
  virtual MySymmetricMatrix operator+ (const MySymmetricMatrix &a) const = 0;

  /// SYMMETRICMATRIX - SYMMETRICMATRIX operators
  virtual MySymmetricMatrix operator- (const MySymmetricMatrix &a) const= 0;

  /// SYMMETRICMATRIX - SYMMETRICMATRIX operators
  virtual MyMatrix operator* (const MySymmetricMatrix &a) const = 0;



  /// SYMMETRICMATRIX - VECTOR operator
  virtual ColumnVector operator* ( const MyColumnVector &b) const = 0;

  /// SYMMETRICMATRIX - VECTOR operator
  virtual void multiply( const MyColumnVector &b, MyColumnVector &result) const = 0;

  /// resize symmetric matrix
  virtual void resize(unsigned int i, bool copy=true, bool initialize=true) = 0;

  /// get inverse
  virtual MySymmetricMatrix inverse() const = 0;

  /// get transpose
  virtual MySymmetricMatrix transpose() const = 0;

  /// get determinant
  virtual double determinant() const = 0;

  /// get sub matrix
  virtual MyMatrix sub(int i_start, int i_end, int j_start , int j_end) const = 0;

  /// Cholesky Decomposition for semidefinite matrices
  virtual bool cholesky_semidefinite(MyMatrix& m) const ;

}; //class SymmetricMatrix_Wrapper


} // namespace



// include
#include "matrix_NEWMAT.h"
#include "matrix_LTI.h"
#include "matrix_BOOST.h"
#include "matrix_EIGEN.h"


#endif // __MATRIX_WRAPPER__
