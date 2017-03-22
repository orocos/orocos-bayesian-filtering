#include "../config.h"
#ifdef __MATRIXWRAPPER_EIGEN__

#include "vector_EIGEN.h"
#include <iostream>


// Constructors
MyColumnVector::ColumnVector() : EigenColumnVector() {}
MyColumnVector::ColumnVector(int num_rows) : EigenColumnVector(num_rows){}
MyColumnVector::ColumnVector(int num_rows,double value) : EigenColumnVector(num_rows){
  ((EigenColumnVector*)this)->setConstant(value);
}
MyColumnVector::ColumnVector(const MyColumnVector& a, const MyColumnVector& b) : EigenColumnVector(a.rows() + b.rows())
{
  EigenColumnVector& opl = (*this);
  opl.head(a.rows()) = (const EigenColumnVector &)(a);
  opl.tail(b.rows()) = (const EigenColumnVector &)(b);
}

// Destructor
MyColumnVector::~ColumnVector(){}

// Copy constructor
MyColumnVector::ColumnVector(const MyColumnVector& a) :
  EigenColumnVector(a){}
MyColumnVector::ColumnVector(const EigenColumnVector & a) :
  EigenColumnVector(a){}

// Resizing
void MyColumnVector::resize(int num_rows)
{
  EigenColumnVector & op1 = (*this);
  op1.resize(num_rows);
}

// Assign
void MyColumnVector::assign(int num_rows, double value)
{
  EigenColumnVector & op1 = (*this);
  op1.resize(num_rows);
  op1.setConstant(value);
}

// Number of Rows / Cols
unsigned int MyColumnVector::rows() const { return ((const EigenColumnVector *)this)->rows();}
unsigned int MyColumnVector::columns() const { return ((const EigenColumnVector *)this)->cols();}
unsigned int MyColumnVector::capacity() const { return ((const EigenColumnVector *)this)->size();}

MyColumnVector
MyColumnVector::vectorAdd(const MyColumnVector& v2) const
{
  const MyColumnVector& v1 = *this;
  MyColumnVector res(v1.rows() + v2.rows());
  EigenColumnVector& opl = res;
  opl.head(v1.rows()) = (const EigenColumnVector &)(v1);
  opl.tail(v2.rows()) = (const EigenColumnVector &)(v2);

  return res;
}

double& MyColumnVector::operator()(unsigned int i)
{
  //std::cout << "(BOOSTVECTOR) operator() called" << std::endl;
  EigenColumnVector& op1 = *(this);
  return op1(i-1);
}

double MyColumnVector::operator()(unsigned int i) const
{
  //std::cout << "(BOOSTVECTOR) operator() called" << std::endl;
  const EigenColumnVector op1 = (*this);
  return op1(i-1);
}


bool MyColumnVector::operator==(const MyColumnVector& a) const
{
  if (this->rows() != a.rows()) return false;
  return(((EigenColumnVector)(*this)-(EigenColumnVector)a).isApproxToConstant(0.0));
}

// Operators
MyColumnVector & MyColumnVector::operator+= (const MyColumnVector& a)
{
  EigenColumnVector & op1 = (*this);
  const EigenColumnVector & op2 = a;
  op1 += op2;
  return (MyColumnVector &) op1;
}

MyColumnVector & MyColumnVector::operator-= (const MyColumnVector& a)
{
  EigenColumnVector & op1 = (*this);
  const EigenColumnVector & op2 = a;
  op1 -= op2;
  return (MyColumnVector &) op1;
}

MyColumnVector MyColumnVector::operator+ (const MyColumnVector &a) const
{
  return (MyColumnVector) ((EigenColumnVector)(*this) + (EigenColumnVector)a);
}

MyColumnVector MyColumnVector::operator- (const MyColumnVector &a) const
{
  return (MyColumnVector) ((EigenColumnVector)(*this) - (EigenColumnVector)a);
}



MyColumnVector& MyColumnVector::operator+= (double a)
{
  EigenColumnVector & op1 = *this;
  op1 += EigenColumnVector::Constant(rows(), a);
  return (MyColumnVector&)op1;
}

MyColumnVector& MyColumnVector::operator-= (double a)
{
  EigenColumnVector & op1 = *this;
  op1 -= EigenColumnVector::Constant(rows(), a);
  return (MyColumnVector&)op1;
}

MyColumnVector& MyColumnVector::operator*= (double a)
{
  EigenColumnVector& op1 = *this;
  op1 *= a;
  return (MyColumnVector&) op1;
}

MyColumnVector& MyColumnVector::operator/= (double a)
{
  EigenColumnVector& op1 = *this;
  op1 /= a;
  return (MyColumnVector&) op1;
}


MyColumnVector MyColumnVector::operator+ (double a) const
{
  return (MyColumnVector)(((EigenColumnVector)(*this)) + EigenColumnVector::Constant(rows(), a));
}

MyColumnVector MyColumnVector::operator- (double a) const
{
  return (MyColumnVector)(((EigenColumnVector)(*this)) - EigenColumnVector::Constant(rows(), a));
}

MyColumnVector MyColumnVector::operator* (double a) const
{
  const EigenColumnVector & op1 = (*this);
  return (MyColumnVector) (op1 * a);
}

MyColumnVector MyColumnVector::operator/ (double a) const
{
  const EigenColumnVector & op1 = (*this);
  return (MyColumnVector) (op1 / a);
}



MyRowVector MyColumnVector::transpose() const
{
  const EigenColumnVector & op1 = (*this);
  return MyRowVector(op1.transpose());
}

MyMatrix MyColumnVector::operator* (const MyRowVector &a) const
{
  const EigenColumnVector & op1 = (*this);
  const EigenRowVector & op2 = a;

  return MyMatrix(op1 * op2);
}

MyColumnVector&
MyColumnVector::operator=(const MyColumnVector &a)
{
  EigenColumnVector& op1 = *this;
  op1 = (EigenColumnVector)a;
  return *this;
}

MyColumnVector&
MyColumnVector::operator=(double a)
{
  EigenColumnVector& op1 = *this;
  op1.setConstant(a);
  return *this;
}

MyColumnVector MyColumnVector::sub(int j_start , int j_end) const
{
  const EigenColumnVector& op1 = *this;
  return MyColumnVector(op1.segment(j_start-1,j_end-j_start+1));
}



//////////////////////////////////////////////////////////////////////
////////////////////////////// ROWVECTOR /////////////////////////////
//////////////////////////////////////////////////////////////////////

// Constructors
MyRowVector::RowVector() : EigenRowVector() {}
MyRowVector::RowVector(int num_cols) : EigenRowVector(num_cols){}
MyRowVector::RowVector(int num_cols,double value) : EigenRowVector(num_cols){
  ((EigenRowVector*)this)->setConstant(value);
}

// Destructor
MyRowVector::~RowVector(){}

// Copy constructor
MyRowVector::RowVector(const MyRowVector& a) :
  EigenRowVector(a){}
MyRowVector::RowVector(const EigenRowVector & a) :
  EigenRowVector(a){}

// Resizing
void MyRowVector::resize(int num_columns)
{
  EigenRowVector & op1 = (*this);
  op1.resize(num_columns);
}

// Assign
void MyRowVector::assign(int num_columns, double value)
{
  EigenRowVector & op1 = (*this);
  op1.resize(num_columns);
  op1.setConstant(value);
}

// Number of Rows / Cols
unsigned int MyRowVector::rows() const { return ((const EigenRowVector *)this)->rows();}
unsigned int MyRowVector::columns() const { return ((const EigenRowVector *)this)->cols();}
unsigned int MyRowVector::capacity() const { return ((const EigenRowVector *)this)->size();}

MyRowVector
MyRowVector::vectorAdd(const MyRowVector& v2) const
{
  const MyRowVector& v1 = *this;
  MyRowVector res(v1.rows() + v2.rows());
  EigenRowVector& opl = res;
  opl.head(v1.rows()) = (const EigenRowVector &)(v1);
  opl.tail(v2.rows()) = (const EigenRowVector &)(v2);
  return res;
}

double& MyRowVector::operator()(unsigned int i)
{
  EigenRowVector& op1 = *(this);
  return op1(i-1);
}

double MyRowVector::operator()(unsigned int i) const
{
  const EigenRowVector& op1 = (*this);
  return op1(i-1);
}

bool MyRowVector::operator==(const MyRowVector& a) const
{
  if (this->columns() != a.columns()) return false;
  return(((EigenRowVector)(*this)-(EigenRowVector)a).isApproxToConstant(0.0));
}

// Operators
MyRowVector & MyRowVector::operator+= (const MyRowVector& a)
{
  EigenRowVector & op1 = (*this);
  const EigenRowVector & op2 = a;
  op1 += op2;
  return (MyRowVector &) op1;
}

MyRowVector & MyRowVector::operator-= (const MyRowVector& a)
{
  EigenRowVector & op1 = (*this);
  const EigenRowVector & op2 = a;
  op1 -= op2;
  return (MyRowVector &) op1;
}

MyRowVector MyRowVector::operator+ (const MyRowVector &a) const
{
  return (MyRowVector) ((EigenRowVector)(*this) + (EigenRowVector)a);
}

MyRowVector MyRowVector::operator- (const MyRowVector &a) const
{
  return (MyRowVector) ((EigenRowVector)(*this) - (EigenRowVector)a);
}



MyRowVector& MyRowVector::operator+= (double a)
{
  EigenRowVector & op1 = *this;
  op1 += EigenRowVector::Constant(columns(),a);
  return (MyRowVector&)op1;
}

MyRowVector& MyRowVector::operator-= (double a)
{
  EigenRowVector & op1 = *this;
  op1 -= EigenRowVector::Constant(columns(),a);
  return (MyRowVector&)op1;
}

MyRowVector& MyRowVector::operator*= (double a)
{
  EigenRowVector& op1 = *this;
  op1 *= a;
  return (MyRowVector&) op1;
}

MyRowVector& MyRowVector::operator/= (double a)
{
  EigenRowVector& op1 = *this;
  op1 /= a;
  return (MyRowVector&) op1;
}


MyRowVector MyRowVector::operator+ (double a) const
{
  return (MyRowVector)(((EigenRowVector)(*this)) + EigenRowVector::Constant(columns(),a));
}

MyRowVector MyRowVector::operator- (double a) const
{
  return (MyRowVector)(((EigenRowVector)(*this)) - EigenRowVector::Constant(columns(),a));
}

MyRowVector MyRowVector::operator* (double a) const
{
  const EigenRowVector & op1 = (*this);
  return (MyRowVector) (op1 * a);
}

MyRowVector MyRowVector::operator/ (double a) const
{
  const EigenRowVector & op1 = (*this);
  return (MyRowVector) (op1 / a);
}



MyColumnVector MyRowVector::transpose() const
{
  const EigenRowVector & op1 = (*this);
  return MyColumnVector(op1.transpose());
}

double MyRowVector::operator* (const MyColumnVector &a) const
{
  const EigenRowVector & op1 = (*this);
  const EigenColumnVector & op2 = a;
  return (op1 * op2)(0,0);
}

MyRowVector&
MyRowVector::operator=(const MyRowVector &a)
{
  EigenRowVector& op1 = *this;
  op1 = (EigenRowVector)a;
  return *this;
}

MyRowVector&
MyRowVector::operator=(double a)
{
  EigenRowVector& op1 = *this;
  op1.setConstant(a);
  return *this;
}

MyRowVector MyRowVector::sub(int j_start , int j_end) const
{
  const EigenRowVector& op1 = *this;
  return MyRowVector(op1.segment(j_start-1,j_end-j_start+1));
}

#endif
