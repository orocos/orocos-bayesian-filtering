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

  unsigned int one = 1;
  unsigned int r = 4;
  unsigned int c = 3;

  // reference numbers
  vector< vector<double> > REF;
  for (unsigned int i=0; i<r+c; i++){
    vector<double> row;
    for (unsigned int j=0; j<r+c; j++)
      row.push_back( (r+c)*i+j );
    REF.push_back(row);
  }

  // test dimensions
  Matrix Am(r,c);
  CPPUNIT_ASSERT_EQUAL(Am.rows(), r);
  CPPUNIT_ASSERT_EQUAL(Am.columns(), c);
  SymmetricMatrix As(r);
  CPPUNIT_ASSERT_EQUAL(As.rows(), r);
  CPPUNIT_ASSERT_EQUAL(As.columns(), r);
  ColumnVector Ac(r);
  CPPUNIT_ASSERT_EQUAL(Ac.rows(), r);
  CPPUNIT_ASSERT_EQUAL(Ac.columns(), one);
  RowVector Ar(c);
  CPPUNIT_ASSERT_EQUAL(Ar.rows(), one);
  CPPUNIT_ASSERT_EQUAL(Ar.columns(), c);

  // test operator = double
  double v = 3.5;
  Am = v;  As = v;  Ac = v;  Ar = v;
  for (unsigned int i=0; i<r; i++)
    for (unsigned int j=0; j<c; j++)
      CPPUNIT_ASSERT_EQUAL(Am(i+1,j+1), v);
  for (unsigned int i=0; i<r; i++)
    for (unsigned int j=0; j<r; j++)
      CPPUNIT_ASSERT_EQUAL(As(i+1,j+1), v);
  for (unsigned int i=0; i<r; i++)
    CPPUNIT_ASSERT_EQUAL(Ac(i+1), v);
  for (unsigned int i=0; i<c; i++)
    CPPUNIT_ASSERT_EQUAL(Ar(i+1), v);

  // test operator ()
  Matrix Bm(r,c);
  for (unsigned int i=0; i<r; i++){
    for (unsigned int j=0; j<c; j++){
      Bm(i+1,j+1) = REF[i][j];
      CPPUNIT_ASSERT_EQUAL(Bm(i+1,j+1), REF[i][j]);
    }
  }
  SymmetricMatrix Bs(r); 
  for (unsigned int i=0; i<r; i++){ // fill in upper triangle
    for (unsigned int j=0; j<=i; j++){
      Bs(j+1,i+1) = REF[i][j];
      CPPUNIT_ASSERT_EQUAL(Bs(i+1,j+1), REF[i][j]);
      CPPUNIT_ASSERT_EQUAL(Bs(j+1,i+1), REF[i][j]);
    }
  }
  for (unsigned int i=0; i<r; i++){   // fill in lower triangle
    for (unsigned int j=0; j<=i; j++){
      Bs(i+1,j+1) = REF[i][j];
      CPPUNIT_ASSERT_EQUAL(Bs(i+1,j+1), REF[i][j]);
      CPPUNIT_ASSERT_EQUAL(Bs(j+1,i+1), REF[i][j]);
    }
  }
  ColumnVector Bc(r);
  for (unsigned int i=0; i<r; i++){
    Bc(i+1) = REF[0][i];
    CPPUNIT_ASSERT_EQUAL(Bc(i+1), REF[0][i]);
  }
  RowVector Br(c);
  for (unsigned int i=0; i<c; i++){
    Br(i+1) = REF[0][i];
    CPPUNIT_ASSERT_EQUAL(Br(i+1), REF[0][i]);
  }

  // test operator = 
  Am = Bm;
  for (unsigned int i=0; i<r; i++){
    for (unsigned int j=0; j<c; j++){
      CPPUNIT_ASSERT_EQUAL(Am(i+1,j+1), REF[i][j]);
      CPPUNIT_ASSERT_EQUAL(Bm(i+1,j+1), REF[i][j]);
    }
  }
  As = Bs;
  for (unsigned int i=0; i<r; i++){
    for (unsigned int j=0; j<=i; j++){
      CPPUNIT_ASSERT_EQUAL(As(i+1,j+1), REF[i][j]);
      CPPUNIT_ASSERT_EQUAL(Bs(i+1,j+1), REF[i][j]);
      CPPUNIT_ASSERT_EQUAL(As(j+1,i+1), REF[i][j]);
      CPPUNIT_ASSERT_EQUAL(Bs(j+1,i+1), REF[i][j]);
    }
  }
  Ac = Bc;
  for (unsigned int i=0; i<r; i++){
    CPPUNIT_ASSERT_EQUAL(Ac(i+1), REF[0][i]);
    CPPUNIT_ASSERT_EQUAL(Bc(i+1), REF[0][i]);
  }
  Ar = Br;
  for (unsigned int i=0; i<c; i++){
    CPPUNIT_ASSERT_EQUAL(Ar(i+1), REF[0][i]);
    CPPUNIT_ASSERT_EQUAL(Br(i+1), REF[0][i]);
  }


  // test resize
  Matrix Km(r+2,c+2); Km = v;
  Matrix Km_resize(r,c); Km_resize = v;
  CPPUNIT_ASSERT_EQUAL(Km.rows(), r+2);
  CPPUNIT_ASSERT_EQUAL(Km.columns(), c+2);
  Km.resize(r,c);
  CPPUNIT_ASSERT_EQUAL(Km.rows(), r);
  CPPUNIT_ASSERT_EQUAL(Km.columns(), c);
  //CPPUNIT_ASSERT_EQUAL(Km, Km_resize);

  // test operator ==
  Matrix Bm_eq; Bm_eq = Bm;
  CPPUNIT_ASSERT_EQUAL(Bm_eq == Bm, true);
  Bm(1,1) = Bm(1,1) + v;
  CPPUNIT_ASSERT_EQUAL(Bm == Bm_eq, false);
  Matrix Bm_res(Bm.rows(), Bm.columns()+1);
  CPPUNIT_ASSERT_EQUAL(Bm == Bm_res, false);
  Bm_res = Bm;
  CPPUNIT_ASSERT_EQUAL(Bm == Bm_res, true);

  SymmetricMatrix Bs_eq; Bs_eq = Bs;
  CPPUNIT_ASSERT_EQUAL(Bs_eq == Bs, true);
  Bs(1,1) = Bs(1,1) + v;
  CPPUNIT_ASSERT_EQUAL(Bs == Bs_eq, false);
  SymmetricMatrix Bs_res(Bs.columns()+1);
  CPPUNIT_ASSERT_EQUAL(Bs == Bs_res, false);
  Bs_res = Bs;
  CPPUNIT_ASSERT_EQUAL(Bs == Bs_res, true);

  ColumnVector Bc_eq; Bc_eq = Bc;
  CPPUNIT_ASSERT_EQUAL(Bc_eq == Bc, true);
  Bc(1) = Bc(1) + v;
  CPPUNIT_ASSERT_EQUAL(Bc == Bc_eq, false);
  ColumnVector Bc_res(Bc.rows()+1);
  CPPUNIT_ASSERT_EQUAL(Bc == Bc_res, false);
  Bc_res = Bc;
  CPPUNIT_ASSERT_EQUAL(Bc == Bc_res, true);

  RowVector Br_eq; Br_eq = Br;
  CPPUNIT_ASSERT_EQUAL(Br_eq == Br, true);
  Br(1) = Br(1) + v;
  CPPUNIT_ASSERT_EQUAL(Br == Br_eq, false);
  RowVector Br_res(Br.columns()+1);
  CPPUNIT_ASSERT_EQUAL(Br == Br_res, false);
  Br_res = Br;
  CPPUNIT_ASSERT_EQUAL(Br == Br_res, true);


  // test transpose
  Matrix Am_trans = Am.transpose();
  CPPUNIT_ASSERT_EQUAL(Am.rows(),  Am_trans.columns());
  CPPUNIT_ASSERT_EQUAL(Am.columns(),  Am_trans.rows());
  for (unsigned int i=0; i<r; i++){
    for (unsigned int j=0; j<c; j++){
      CPPUNIT_ASSERT_EQUAL(Am(i+1,j+1),  REF[i][j]);
      CPPUNIT_ASSERT_EQUAL(Am_trans(j+1,i+1), REF[i][j]);
    }
  }
  SymmetricMatrix As_trans = As.transpose();
  CPPUNIT_ASSERT_EQUAL(As.rows(),  As_trans.columns());
  CPPUNIT_ASSERT_EQUAL(As.columns(),  As_trans.rows());
  for (unsigned int i=0; i<r; i++){
    for (unsigned int j=0; j<=i; j++){
      CPPUNIT_ASSERT_EQUAL(As(i+1,j+1),  REF[i][j]);
      CPPUNIT_ASSERT_EQUAL(As_trans(j+1,i+1), REF[i][j]);
      CPPUNIT_ASSERT_EQUAL(As(j+1,i+1),  REF[i][j]);
      CPPUNIT_ASSERT_EQUAL(As_trans(i+1,j+1), REF[i][j]);
    }
  }
  RowVector Ac_trans = Ac.transpose();
  CPPUNIT_ASSERT_EQUAL(Ac.rows(),  Ac_trans.columns());
  CPPUNIT_ASSERT_EQUAL(Ac.columns(),  Ac_trans.rows());
  for (unsigned int i=0; i<r; i++){
    CPPUNIT_ASSERT_EQUAL(Ac(i+1),  REF[0][i]);
    CPPUNIT_ASSERT_EQUAL(Ac_trans(i+1), REF[0][i]);
  }
  ColumnVector Ar_trans = Ar.transpose();
  CPPUNIT_ASSERT_EQUAL(Ar.rows(),  Ar_trans.columns());
  CPPUNIT_ASSERT_EQUAL(Ar.columns(),  Ar_trans.rows());
  for (unsigned int i=0; i<c; i++){
    CPPUNIT_ASSERT_EQUAL(Ar(i+1),  REF[0][i]);
    CPPUNIT_ASSERT_EQUAL(Ar_trans(i+1), REF[0][i]);
  }

  // test sub matrix
  Matrix Am_sub = Am.sub(1,c,1,c);
  CPPUNIT_ASSERT_EQUAL(Am_sub.rows(), c);
  CPPUNIT_ASSERT_EQUAL(Am_sub.columns(), c);
  for (unsigned int i=0; i<c; i++){
    for (unsigned int j=0; j<c; j++){
      CPPUNIT_ASSERT_EQUAL(Am_sub(i+1,j+1), Am(i+1,j+1));
      CPPUNIT_ASSERT_EQUAL(Am_sub(i+1,j+1), REF[i][j]);
    }
  }
  Matrix As_sub = As.sub(1,c,1,c);
  CPPUNIT_ASSERT_EQUAL(As_sub.rows(), c);
  CPPUNIT_ASSERT_EQUAL(As_sub.columns(), c);
  for (unsigned int i=0; i<c; i++){
    for (unsigned int j=0; j<=i; j++){ 
      CPPUNIT_ASSERT_EQUAL(As_sub(i+1,j+1), As(i+1,j+1));
      CPPUNIT_ASSERT_EQUAL(As_sub(i+1,j+1), REF[i][j]);
      CPPUNIT_ASSERT_EQUAL(As_sub(j+1,i+1), As(i+1,j+1));
      CPPUNIT_ASSERT_EQUAL(As_sub(j+1,i+1), REF[i][j]);
    }
  }
  ColumnVector Ac_sub = Ac.sub(1,c);
  CPPUNIT_ASSERT_EQUAL(Ac_sub.rows(), c);
  CPPUNIT_ASSERT_EQUAL(Ac_sub.columns(), one);
  for (unsigned int i=0; i<c; i++){
    CPPUNIT_ASSERT_EQUAL(Ac_sub(i+1), Ac(i+1));
    CPPUNIT_ASSERT_EQUAL(Ac_sub(i+1), REF[0][i]);
  }
  RowVector Ar_sub = Ar.sub(1,r-1);
  CPPUNIT_ASSERT_EQUAL(Ar_sub.rows(), one);
  CPPUNIT_ASSERT_EQUAL(Ar_sub.columns(), r-1);
  for (unsigned int i=0; i<r-1; i++){
    CPPUNIT_ASSERT_EQUAL(Ar_sub(i+1), Ar(i+1));
    CPPUNIT_ASSERT_EQUAL(Ar_sub(i+1), REF[0][i]);
  }
  
  // test operator *
  Matrix Cm = Am * Am_trans;
  for (unsigned int i=0; i<r; i++){
    for (unsigned int j=0; j<c; j++){
      CPPUNIT_ASSERT_EQUAL(Am(i+1,j+1),  REF[i][j]);
      CPPUNIT_ASSERT_EQUAL(Am_trans(j+1,i+1), REF[i][j]);
    }
  }

  Cm = Am * v;
  for (unsigned int i=0; i<r; i++){
    for (unsigned int j=0; j<c; j++){
      CPPUNIT_ASSERT_EQUAL(Cm(i+1,j+1),  REF[i][j] * v);
    }
  }

  SymmetricMatrix Cs = As * As_trans;
  for (unsigned int i=0; i<r; i++){
    for (unsigned int j=0; j<=i; j++){
      CPPUNIT_ASSERT_EQUAL(As(i+1,j+1),  REF[i][j]);
      CPPUNIT_ASSERT_EQUAL(As_trans(j+1,i+1), REF[i][j]);
      CPPUNIT_ASSERT_EQUAL(As(j+1,i+1),  REF[i][j]);
      CPPUNIT_ASSERT_EQUAL(As_trans(i+1,j+1), REF[i][j]);
    }
  }

  Cs = As * v;
  for (unsigned int i=0; i<r; i++){
    for (unsigned int j=0; j<=i; j++){
      CPPUNIT_ASSERT_EQUAL(Cs(i+1,j+1),  REF[i][j] * v);
      CPPUNIT_ASSERT_EQUAL(Cs(j+1,i+1),  REF[i][j] * v);
    }
  }
  ColumnVector Cc = Ac * v;
  for (unsigned int i=0; i<r; i++)
    CPPUNIT_ASSERT_EQUAL(Cc(i+1),  REF[0][i] * v);
  RowVector Cr = Ar * v;
  for (unsigned int i=0; i<c; i++)
    CPPUNIT_ASSERT_EQUAL(Cr(i+1),  REF[0][i] * v);


  // test operator *
  Cm = Am / v;
  for (unsigned int i=0; i<r; i++){
    for (unsigned int j=0; j<c; j++){
      CPPUNIT_ASSERT_EQUAL(Cm(i+1,j+1),  REF[i][j] / v);
    }
  }
  Cs = As / v;
  for (unsigned int i=0; i<r; i++){
    for (unsigned int j=0; j<=i; j++){
      CPPUNIT_ASSERT_EQUAL(Cs(i+1,j+1),  REF[i][j] / v);
      CPPUNIT_ASSERT_EQUAL(Cs(j+1,i+1),  REF[i][j] / v);
    }
  }
  Cc = Ac / v;
  for (unsigned int i=0; i<r; i++)
    CPPUNIT_ASSERT_EQUAL(Cc(i+1),  REF[0][i] / v);
  Cr = Ar / v;
  for (unsigned int i=0; i<c; i++)
    CPPUNIT_ASSERT_EQUAL(Cr(i+1),  REF[0][i] / v);

  // test inverse
  Matrix Rm(c,c);
  Rm(1,1) = 3;   Rm(1,2) = 3;    Rm(1,3) = 3;
  Rm(2,1) = 5;   Rm(2,2) = 2;    Rm(2,3) = 9;
  Rm(3,1) = 9;   Rm(3,2) = 7;    Rm(3,3) = 0;
  Matrix Rm_inv = Rm.inverse();
  Matrix Im(c,c); Im = 0;
  for (unsigned int i=0; i<c; i++)
    Im(i+1,i+1) = 1;
  Matrix Im_test = Rm * Rm_inv;
  CPPUNIT_ASSERT_EQUAL(approxEqual(Im_test, Im, epsilon),true);

  Matrix Rs(c,c);
  Rs(1,1) = 3;   Rs(1,2) = 5;    Rs(1,3) = 3;
  Rs(2,1) = 5;   Rs(2,2) = 2;    Rs(2,3) = 7;
  Rs(3,1) = 3;   Rs(3,2) = 7;    Rs(3,3) = 0;
  Matrix Rs_inv = Rs.inverse();
  Matrix Is(c,c); Is = 0;
  for (unsigned int i=0; i<c; i++)
    Is(i+1,i+1) = 1;
  Matrix Is_test = Rs * Rs_inv;
  CPPUNIT_ASSERT_EQUAL(approxEqual(Is_test, Is, epsilon),true);

  // test determinant
  CPPUNIT_ASSERT_EQUAL(approxEqual(Rm.determinant(), 105, epsilon),true);
  CPPUNIT_ASSERT_EQUAL(approxEqual(Rs.determinant(), 45, epsilon),true);

  // test cholesky
  SymmetricMatrix Ps(c);
  Ps(1,1) = 3;   Ps(1,2) = 2;    Ps(1,3) = 1;
  Ps(2,1) = 2;   Ps(2,2) = 2;    Ps(2,3) = 1;
  Ps(3,1) = 1;   Ps(3,2) = 1;    Ps(3,3) = 1;
  Matrix CHs;
  Matrix CHs_check(c,c);
  CHs_check(1,1) = 1.73205;   CHs_check(1,2) = 0.00000;    CHs_check(1,3) = 0.00000;
  CHs_check(2,1) = 1.15470;   CHs_check(2,2) = 0.81650;    CHs_check(2,3) = 0.00000;
  CHs_check(3,1) = 0.57735;   CHs_check(3,2) = 0.40825;    CHs_check(3,3) = 0.70711;
  Ps.cholesky_semidefinite(CHs);
  CPPUNIT_ASSERT_EQUAL(approxEqual(CHs, CHs_check, epsilon),true);  

  // test operator - -=
  Matrix Rm_bak; Rm_bak = Rm;
  Matrix Dm(c,c); Dm = v;
  Matrix Em = Rm - Dm;
  Matrix Fm = Rm - v;
  Matrix Gm(Rm);
  Gm -= Dm;
  Matrix Hm; Hm = Rm;
  Hm -= v;
  CPPUNIT_ASSERT_EQUAL(Rm, Rm_bak);
  CPPUNIT_ASSERT_EQUAL( Em, Fm);
  CPPUNIT_ASSERT_EQUAL( Fm, Gm);
  CPPUNIT_ASSERT_EQUAL( Gm, Hm);
  for (unsigned int i=0; i<c; i++)
    for (unsigned int j=0; j<c; j++)
      CPPUNIT_ASSERT_EQUAL( Dm(i+1,j+1), v);

  // test pseudoinverse
  Matrix Bm_bak; Bm_bak = Bm;
  Matrix Bm_pinv = Bm.pseudoinverse();
  CPPUNIT_ASSERT_EQUAL(Bm, Bm_bak);
  Im_test = Bm_pinv * Bm;
  CPPUNIT_ASSERT_EQUAL(approxEqual(Im_test, Im, epsilon),true);




}



