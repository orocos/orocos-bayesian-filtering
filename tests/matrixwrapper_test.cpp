// Copyright (C) 2007 Wim Meeussen <wim.meeussen@mech.kuleuven.be>
//  
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//  
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//  
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
//  
 

#include "matrixwrapper_test.hpp"
#include "approxEqual.hpp"


// Registers the fixture into the 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION( MatrixwrapperTest );

using namespace MatrixWrapper;


void 
MatrixwrapperTest::setUp()
{
}


void 
MatrixwrapperTest::tearDown()
{
}

void 
MatrixwrapperTest::testMatrixwrapperValue()
{
  double epsilon = 0.00001;

  unsigned int r = 4;
  unsigned int c = 3;

  // test dimensions
  Matrix A(r,c);
  CPPUNIT_ASSERT_EQUAL(A.rows(), r);
  CPPUNIT_ASSERT_EQUAL(A.columns(), c);

  // test operator = double
  double v = 3.5;
  A = v;
  for (unsigned int i=0; i<r; i++)
    for (unsigned int j=0; j<c; j++)
      CPPUNIT_ASSERT_EQUAL(A(i+1,j+1), v);

  // test operator ()
  Matrix B(r,c);
  for (unsigned int i=0; i<r; i++){
    for (unsigned int j=0; j<c; j++){
      B(i+1,j+1) = r*i+j;
      CPPUNIT_ASSERT_EQUAL(B(i+1,j+1), (double)(r*i+j));
    }
  }

  // test operator = matrix
  A = B;
  for (unsigned int i=0; i<r; i++){
    for (unsigned int j=0; j<c; j++){
      CPPUNIT_ASSERT_EQUAL(A(i+1,j+1), (double)(r*i+j));
      CPPUNIT_ASSERT_EQUAL(B(i+1,j+1), (double)(r*i+j));
    }
  }

  // test resize
  Matrix K(r+2,c+2); K = v;
  Matrix Ktest(r,c); Ktest = v;
  CPPUNIT_ASSERT_EQUAL(K.rows(), r+2);
  CPPUNIT_ASSERT_EQUAL(K.columns(), c+2);
  K.resize(r,c);
  CPPUNIT_ASSERT_EQUAL(K.rows(), r);
  CPPUNIT_ASSERT_EQUAL(K.columns(), c);
  CPPUNIT_ASSERT_EQUAL(K, Ktest);

  // test operator ==
  Matrix Beq; Beq = B;
  CPPUNIT_ASSERT_EQUAL(Beq == B, true);
  B(1,1) = B(1,1) + 1;
  CPPUNIT_ASSERT_EQUAL(B == Beq, false);
  Matrix Bres(B.rows(), B.columns()+1);
  CPPUNIT_ASSERT_EQUAL(B == Bres, false);
  Bres = B;
  CPPUNIT_ASSERT_EQUAL(B == Bres, true);
  CPPUNIT_ASSERT_EQUAL(B.rows(), Bres.rows());
  CPPUNIT_ASSERT_EQUAL(B.columns(), Bres.columns());

  // test transpose
  Matrix At = A.transpose();
  for (unsigned int i=0; i<r; i++){
    for (unsigned int j=0; j<c; j++){
      CPPUNIT_ASSERT_EQUAL(A(i+1,j+1),  (double)(r*i+j));
      CPPUNIT_ASSERT_EQUAL(At(j+1,i+1), (double)(r*i+j));
    }
  }

  // test sub matrix
  Matrix Asub = A.sub(1,c,1,c);
  CPPUNIT_ASSERT_EQUAL(Asub.rows(), c);
  CPPUNIT_ASSERT_EQUAL(Asub.columns(), c);
  for (unsigned int i=0; i<c; i++){
    for (unsigned int j=0; j<c; j++){
      CPPUNIT_ASSERT_EQUAL(Asub(i+1,j+1), A(i+1,j+1));
      CPPUNIT_ASSERT_EQUAL(Asub(i+1,j+1), (double)(r*i+j));
    }
  }
  
  // test operator * matrix
  Matrix C = A * At;
  for (unsigned int i=0; i<r; i++){
    for (unsigned int j=0; j<c; j++){
      CPPUNIT_ASSERT_EQUAL(A(i+1,j+1),  (double)(r*i+j));
      CPPUNIT_ASSERT_EQUAL(At(j+1,i+1), (double)(r*i+j));
    }
  }

  // test inverse
  Matrix R(c,c);
  R(1,1) = 3;   R(1,2) = 3;    R(1,3) = 3;
  R(2,1) = 5;   R(2,2) = 2;    R(2,3) = 9;
  R(3,1) = 9;   R(3,2) = 7;    R(3,3) = 0;
  Matrix Rinv = R.inverse();
  Matrix I(c,c); I = 0;
  for (unsigned int i=0; i<c; i++)
    I(i+1,i+1) = 1;
  Matrix Itest = R * Rinv;
  CPPUNIT_ASSERT_EQUAL(approxEqual(Itest, I, epsilon),true);

  // test operator - -=
  Matrix Rbak; Rbak = R;
  Matrix D(c,c); D = v;
  Matrix E = R - D;
  Matrix F = R - v;
  Matrix G(R);
  G -= D;
  Matrix H; H = R;
  H -= v;
  CPPUNIT_ASSERT_EQUAL(R, Rbak);
  CPPUNIT_ASSERT_EQUAL( E, F);
  CPPUNIT_ASSERT_EQUAL( F, G);
  CPPUNIT_ASSERT_EQUAL( G, H);
  for (unsigned int i=0; i<c; i++)
    for (unsigned int j=0; j<c; j++)
      CPPUNIT_ASSERT_EQUAL( D(i+1,j+1), v);
  
  // test pseudoinverse
  Matrix Bbak; Bbak = B;
  Matrix Bpinv = B.pseudoinverse();
  CPPUNIT_ASSERT_EQUAL(B, Bbak);
  Itest = Bpinv * B;
  CPPUNIT_ASSERT_EQUAL(approxEqual(Itest, I, epsilon),true);




}



