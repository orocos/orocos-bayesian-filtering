#include "../config.h"
#ifdef __MATRIXWRAPPER_EIGEN__

#ifndef __MATRIX_EIGEN__
#define __MATRIX_EIGEN__
#include "../../bfl_constants.h"
#include "matrix_wrapper.h"
#include "vector_wrapper.h"

#include <Eigen/Core>
#include <assert.h>

typedef Eigen::MatrixXd EigenMatrix;
typedef Eigen::MatrixXd EigenSymmetricMatrix;
typedef Eigen::MatrixXd::ConstSelfAdjointViewReturnType<Eigen::Upper>::Type EigenSymmetricView;

namespace MatrixWrapper
{

/// Implementation of Matrixwrapper using Eigen
class Matrix : public EigenMatrix, public Matrix_Wrapper
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
  Matrix(const EigenMatrix & a);

  Matrix(int num_rows,const RowVector& v);

  
  virtual unsigned int size() const;
  virtual unsigned int capacity() const;
  virtual unsigned int rows() const;
  virtual unsigned int columns() const;
  virtual double& operator()(unsigned int,unsigned int);
  virtual double operator()(unsigned int,unsigned int) const;
  virtual RowVector operator[](unsigned int)const;

  using EigenMatrix::operator ==;
  using EigenMatrix::operator =;
  using EigenMatrix::operator +=;
  using EigenMatrix::operator -=;
  using EigenMatrix::operator +;
  using EigenMatrix::operator -;

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

class SymmetricMatrix : public EigenSymmetricMatrix, public SymmetricMatrix_Wrapper
{
 private: //

 public: //
  // Constructors
  SymmetricMatrix();
  SymmetricMatrix(int n);

  // Copy constructors
  SymmetricMatrix(const MySymmetricMatrix& a);
  SymmetricMatrix(const EigenSymmetricMatrix& a);
  SymmetricMatrix(const EigenSymmetricView & a);

  SymmetricMatrix(int num_rows,const RowVector& v);

  // Destructor
  virtual ~SymmetricMatrix();

  virtual unsigned int size() const;
  virtual unsigned int capacity() const;
  virtual unsigned int rows() const;
  virtual unsigned int columns() const;
  virtual MySymmetricMatrix inverse() const;
  virtual MySymmetricMatrix transpose() const;
  virtual double determinant() const;

  virtual double& operator()(unsigned int,unsigned int);
  virtual double operator()(unsigned int,unsigned int) const;
  virtual RowVector operator[](unsigned int)const;

  using EigenSymmetricMatrix::operator ==;
  using EigenSymmetricMatrix::operator =;
  using EigenSymmetricMatrix::operator +=;
  using EigenSymmetricMatrix::operator -=;
  using EigenSymmetricMatrix::operator +;
  using EigenSymmetricMatrix::operator -;

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
