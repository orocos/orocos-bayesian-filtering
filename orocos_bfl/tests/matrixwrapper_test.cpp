// Copyright (C) 2007 Wim Meeussen <wim.meeussen@mech.kuleuven.be>
//                    Tinne De Laet<first DOT last AT mech.kuleuven.be>
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
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
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

  // TEST DIMENSIONS
  //
  // TEST DIMENSIONS  MATRIX
  Matrix Am(r,c);
  // rows()
  CPPUNIT_ASSERT_EQUAL(Am.rows(), r);
  // columns()
  CPPUNIT_ASSERT_EQUAL(Am.columns(), c);
  // TEST DIMENSIONS SYMMETRICMATRIX
  SymmetricMatrix As(r);
  // rows()
  CPPUNIT_ASSERT_EQUAL(As.rows(), r);
  // columns()
  CPPUNIT_ASSERT_EQUAL(As.columns(), r);
  // TEST DIMENSIONS COLUMNVECTOR
  ColumnVector Ac(r);
  // rows()
  CPPUNIT_ASSERT_EQUAL(Ac.rows(), r);
  // columns()
  CPPUNIT_ASSERT_EQUAL(Ac.columns(), one);
  // TEST DIMENSIONS ROWVECTOR
  RowVector Ar(c);
  // rows()
  CPPUNIT_ASSERT_EQUAL(Ar.rows(), one);
  // columns()
  CPPUNIT_ASSERT_EQUAL(Ar.columns(), c);

  // test operator = double
  double v = 3.5;
  Am = v;  As = v;  Ac = v;  Ar = v;
  // MATRIX
  for (unsigned int i=0; i<r; i++)
    for (unsigned int j=0; j<c; j++)
      CPPUNIT_ASSERT_EQUAL(Am(i+1,j+1), v);
  // SYMMETRICMATRIX
  for (unsigned int i=0; i<r; i++)
    for (unsigned int j=0; j<r; j++)
      CPPUNIT_ASSERT_EQUAL(As(i+1,j+1), v);
  // COLUMNVECTOR
  for (unsigned int i=0; i<r; i++)
    CPPUNIT_ASSERT_EQUAL(Ac(i+1), v);
  // ROWVECTOR
  for (unsigned int i=0; i<c; i++)
    CPPUNIT_ASSERT_EQUAL(Ar(i+1), v);

  // test operator ()
  // MATRIX
  Matrix Bm(r,c);
  for (unsigned int i=0; i<r; i++){
    for (unsigned int j=0; j<c; j++){
      Bm(i+1,j+1) = REF[i][j];
      CPPUNIT_ASSERT_EQUAL(Bm(i+1,j+1), REF[i][j]);
    }
  }

  // SYMMETRICMATRIX
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
  // COLUMNVECTOR
  ColumnVector Bc(r);
  for (unsigned int i=0; i<r; i++){
    Bc(i+1) = REF[0][i];
    CPPUNIT_ASSERT_EQUAL(Bc(i+1), REF[0][i]);
  }
  // ROWVECTOR
  RowVector Br(c);
  for (unsigned int i=0; i<c; i++){
    Br(i+1) = REF[0][i];
    CPPUNIT_ASSERT_EQUAL(Br(i+1), REF[0][i]);
  }

  // test operator =
  // MATRIX
  Am = Bm;
  for (unsigned int i=0; i<r; i++){
    for (unsigned int j=0; j<c; j++){
      CPPUNIT_ASSERT_EQUAL(Am(i+1,j+1), REF[i][j]);
      CPPUNIT_ASSERT_EQUAL(Bm(i+1,j+1), REF[i][j]);
    }
  }
  // SYMMETRICMATRIX
  As = Bs;
  for (unsigned int i=0; i<r; i++){
    for (unsigned int j=0; j<=i; j++){
      CPPUNIT_ASSERT_EQUAL(As(i+1,j+1), REF[i][j]);
      CPPUNIT_ASSERT_EQUAL(Bs(i+1,j+1), REF[i][j]);
      CPPUNIT_ASSERT_EQUAL(As(j+1,i+1), REF[i][j]);
      CPPUNIT_ASSERT_EQUAL(Bs(j+1,i+1), REF[i][j]);
    }
  }
  // COLUMNVECTOR
  Ac = Bc;
  for (unsigned int i=0; i<r; i++){
    CPPUNIT_ASSERT_EQUAL(Ac(i+1), REF[0][i]);
    CPPUNIT_ASSERT_EQUAL(Bc(i+1), REF[0][i]);
  }
  // ROWVECTOR
  Ar = Br;
  for (unsigned int i=0; i<c; i++){
    CPPUNIT_ASSERT_EQUAL(Ar(i+1), REF[0][i]);
    CPPUNIT_ASSERT_EQUAL(Br(i+1), REF[0][i]);
  }



  // test resize
  // MATRIX
  Matrix Km(r+2,c+2); Km = v;
  Matrix Km_resize(r,c); Km_resize = v;
  CPPUNIT_ASSERT_EQUAL(Km.rows(), r+2);
  CPPUNIT_ASSERT_EQUAL(Km.columns(), c+2);
  Km.resize(r,c);
  CPPUNIT_ASSERT_EQUAL(Km.rows(), r);
  CPPUNIT_ASSERT_EQUAL(Km.columns(), c);
  CPPUNIT_ASSERT_EQUAL(Km, Km_resize);
  // SYMMETRICMATRIX
  SymmetricMatrix Ks(r+2); Ks = v;
  SymmetricMatrix Ks_resize(r); Ks_resize = v;
  CPPUNIT_ASSERT_EQUAL(Ks.rows(), r+2);
  CPPUNIT_ASSERT_EQUAL(Ks.columns(), r+2);
//  Ks.resize(r);
//  CPPUNIT_ASSERT_EQUAL(Ks.rows(), r);
//  CPPUNIT_ASSERT_EQUAL(Ks.columns(), r);
//  CPPUNIT_ASSERT_EQUAL(Ks, Ks_resize);
  // COLUMNVECTOR
  ColumnVector Kc(r+2); Kc = v;
  ColumnVector Kc_resize(r); Kc_resize = v;
  CPPUNIT_ASSERT_EQUAL(Kc.rows(), r+2);
  Kc.resize(r);
  CPPUNIT_ASSERT_EQUAL(Kc.rows(), r);
  CPPUNIT_ASSERT_EQUAL(Kc, Kc_resize);
  // ROWVECTOR
  RowVector Kr(c+2); Kr = v;
  RowVector Kr_resize(c); Kr_resize = v;
  CPPUNIT_ASSERT_EQUAL(Kr.columns(), c+2);
  Kr.resize(c);
  CPPUNIT_ASSERT_EQUAL(Kr.columns(), c);
  CPPUNIT_ASSERT_EQUAL(Kr, Kr_resize);

  // test operator ==
  // MATRIX
  Matrix Bm_eq; Bm_eq = Bm;
  CPPUNIT_ASSERT_EQUAL(Bm_eq == Bm, true);
  Bm(1,1) = Bm(1,1) + v;
  CPPUNIT_ASSERT_EQUAL(Bm == Bm_eq, false);
  Matrix Bm_res(Bm.rows(), Bm.columns()+1);
  CPPUNIT_ASSERT_EQUAL(Bm == Bm_res, false);
  Bm_res = Bm;
  CPPUNIT_ASSERT_EQUAL(Bm == Bm_res, true);
  // SYMMETRICMATRIX
  SymmetricMatrix Bs_eq; Bs_eq = Bs;
  CPPUNIT_ASSERT_EQUAL(Bs_eq == Bs, true);
  Bs(1,1) = Bs(1,1) + v;
  CPPUNIT_ASSERT_EQUAL(Bs == Bs_eq, false);
  SymmetricMatrix Bs_res(Bs.columns()+1);
  CPPUNIT_ASSERT_EQUAL(Bs == Bs_res, false);
  Bs_res = Bs;
  CPPUNIT_ASSERT_EQUAL(Bs == Bs_res, true);
  // COLUMNVECTOR
  ColumnVector Bc_eq; Bc_eq = Bc;
  CPPUNIT_ASSERT_EQUAL(Bc_eq == Bc, true);
  Bc(1) = Bc(1) + v;
  CPPUNIT_ASSERT_EQUAL(Bc == Bc_eq, false);
  ColumnVector Bc_res(Bc.rows()+1);
  CPPUNIT_ASSERT_EQUAL(Bc == Bc_res, false);
  Bc_res = Bc;
  CPPUNIT_ASSERT_EQUAL(Bc == Bc_res, true);
  // ROWVECTOR
  RowVector Br_eq; Br_eq = Br;
  CPPUNIT_ASSERT_EQUAL(Br_eq == Br, true);
  Br(1) = Br(1) + v;
  CPPUNIT_ASSERT_EQUAL(Br == Br_eq, false);
  RowVector Br_res(Br.columns()+1);
  CPPUNIT_ASSERT_EQUAL(Br == Br_res, false);
  Br_res = Br;
  CPPUNIT_ASSERT_EQUAL(Br == Br_res, true);


  // test transpose
  // MATRIX
  Matrix Am_trans = Am.transpose();
  CPPUNIT_ASSERT_EQUAL(Am.rows(),  Am_trans.columns());
  CPPUNIT_ASSERT_EQUAL(Am.columns(),  Am_trans.rows());
  for (unsigned int i=0; i<r; i++){
    for (unsigned int j=0; j<c; j++){
      CPPUNIT_ASSERT_EQUAL(Am(i+1,j+1),  REF[i][j]);
      CPPUNIT_ASSERT_EQUAL(Am_trans(j+1,i+1), REF[i][j]);
    }
  }
  // SYMMETRICMATRIX
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
  // COLUMNVECTOR
  RowVector Ac_trans = Ac.transpose();
  CPPUNIT_ASSERT_EQUAL(Ac.rows(),  Ac_trans.columns());
  CPPUNIT_ASSERT_EQUAL(Ac.columns(),  Ac_trans.rows());
  for (unsigned int i=0; i<r; i++){
    CPPUNIT_ASSERT_EQUAL(Ac(i+1),  REF[0][i]);
    CPPUNIT_ASSERT_EQUAL(Ac_trans(i+1), REF[0][i]);
  }
  // ROWVECTOR
  ColumnVector Ar_trans = Ar.transpose();
  CPPUNIT_ASSERT_EQUAL(Ar.rows(),  Ar_trans.columns());
  CPPUNIT_ASSERT_EQUAL(Ar.columns(),  Ar_trans.rows());
  for (unsigned int i=0; i<c; i++){
    CPPUNIT_ASSERT_EQUAL(Ar(i+1),  REF[0][i]);
    CPPUNIT_ASSERT_EQUAL(Ar_trans(i+1), REF[0][i]);
  }

  // test sub matrix
  // MATRIX
  Matrix Am_sub = Am.sub(2,r,2,c);
  CPPUNIT_ASSERT_EQUAL(Am_sub.rows(), r-1);
  CPPUNIT_ASSERT_EQUAL(Am_sub.columns(), c-1);
  for (unsigned int i=0; i<c-1; i++){
    for (unsigned int j=0; j<c-1; j++){
      CPPUNIT_ASSERT_EQUAL(Am_sub(i+1,j+1), Am(i+2,j+2));
      CPPUNIT_ASSERT_EQUAL(Am_sub(i+1,j+1), REF[i+1][j+1]);
    }
  }
  // SYMMETRICMATRIX
  Matrix As_sub = As.sub(2,c,2,c);
  CPPUNIT_ASSERT_EQUAL(As_sub.rows(), c-1);
  CPPUNIT_ASSERT_EQUAL(As_sub.columns(), c-1);
  for (unsigned int i=0; i<c-1; i++){
    for (unsigned int j=0; j<=i; j++){
      CPPUNIT_ASSERT_EQUAL(As_sub(i+1,j+1), As(i+2,j+2));
      CPPUNIT_ASSERT_EQUAL(As_sub(i+1,j+1), REF[i+1][j+1]);
      CPPUNIT_ASSERT_EQUAL(As_sub(j+1,i+1), As(i+2,j+2));
      CPPUNIT_ASSERT_EQUAL(As_sub(j+1,i+1), REF[i+1][j+1]);
    }
  }
  // COLUMNVECTOR
  ColumnVector Ac_sub = Ac.sub(2,c);
  CPPUNIT_ASSERT_EQUAL(Ac_sub.rows(), c-1);
  CPPUNIT_ASSERT_EQUAL(Ac_sub.columns(), one);
  for (unsigned int i=0; i<c-1; i++){
    CPPUNIT_ASSERT_EQUAL(Ac_sub(i+1), Ac(i+2));
    CPPUNIT_ASSERT_EQUAL(Ac_sub(i+1), REF[0][i+1]);
  }
  // ROWVECTOR
  RowVector Ar_sub = Ar.sub(2,r-1);
  CPPUNIT_ASSERT_EQUAL(Ar_sub.rows(), one);
  CPPUNIT_ASSERT_EQUAL(Ar_sub.columns(), r-2);
  for (unsigned int i=0; i<r-2; i++){
    CPPUNIT_ASSERT_EQUAL(Ar_sub(i+1), Ar(i+2));
    CPPUNIT_ASSERT_EQUAL(Ar_sub(i+1), REF[0][i+1]);
  }

  // test operator *
  // MATRIX * MATRIX
  Matrix Cm = Am * Am_trans;
  Matrix Cm_check(r,c);
  for (unsigned int i=0; i<r; i++){
    for (unsigned int j=0; j<c; j++){
      // test if original elements were maintained
      CPPUNIT_ASSERT_EQUAL(Am(i+1,j+1),  REF[i][j]);
      CPPUNIT_ASSERT_EQUAL(Am_trans(j+1,i+1), REF[i][j]);
      // test correct multiplication
      double sum = 0.0;
      for (unsigned int t=0; t<c; t++){
          sum += Am(i+1,t+1) * Am_trans(t+1,j+1);
      }
      Cm_check(i+1,j+1) = sum;
      CPPUNIT_ASSERT_EQUAL(Cm(i+1,j+1),  Cm_check(i+1,j+1));
    }
  }

  // MATRIX * DOUBLE
  Cm = Am * v;
  for (unsigned int i=0; i<r; i++){
    for (unsigned int j=0; j<c; j++){
      // test if original elements were maintained
      CPPUNIT_ASSERT_EQUAL(Am(i+1,j+1),  REF[i][j]);
      CPPUNIT_ASSERT_EQUAL(Am_trans(j+1,i+1), REF[i][j]);
      // test correct multiplication
      CPPUNIT_ASSERT_EQUAL(Cm(i+1,j+1),  REF[i][j] * v);
    }
  }

  // SYMMETRICMATRIX * SYMMETRICMATRIX
  Matrix Cs_check(r,r);
  Matrix Cs = As * As_trans;
  for (unsigned int i=0; i<r; i++){
    for (unsigned int j=0; j<=i; j++){
      // test if original elements were maintained
      CPPUNIT_ASSERT_EQUAL(As(i+1,j+1),  REF[i][j]);
      CPPUNIT_ASSERT_EQUAL(As_trans(j+1,i+1), REF[i][j]);
      CPPUNIT_ASSERT_EQUAL(As(j+1,i+1),  REF[i][j]);
      CPPUNIT_ASSERT_EQUAL(As_trans(i+1,j+1), REF[i][j]);
      // test correct multiplication
      double sum = 0.0;
      for (unsigned int t=0; t<r; t++){
          sum += As(i+1,t+1) * As_trans(t+1,j+1);
      }
      Cs_check(i+1,j+1) = sum;
      CPPUNIT_ASSERT_EQUAL(Cs(i+1,j+1),  Cs_check(i+1,j+1));
    }
  }

  // SYMMETRICMATRIX * DOUBLE
  Cs = As * v;
  for (unsigned int i=0; i<r; i++){
    for (unsigned int j=0; j<=i; j++){
      // test if original elements were maintained
      CPPUNIT_ASSERT_EQUAL(As(i+1,j+1),  REF[i][j]);
      CPPUNIT_ASSERT_EQUAL(As(j+1,i+1),  REF[i][j]);
      // test correct multiplication
      CPPUNIT_ASSERT_EQUAL(Cs(i+1,j+1),  REF[i][j] * v);
      CPPUNIT_ASSERT_EQUAL(Cs(j+1,i+1),  REF[i][j] * v);
    }
  }

  // MATRIX * SYMMETRICMATRIX
  // TODO: not implemented?

  // SYMMETRICMATRIX * MATRIX
  Matrix Csm_check(r,c);
  Matrix Csm = As * Am;
  for (unsigned int i=0; i<r; i++){
    for (unsigned int j=0; j<c; j++){
      // test if original elements were maintained
      if (j<=i){
         CPPUNIT_ASSERT_EQUAL(As(i+1,j+1),  REF[i][j]);
         CPPUNIT_ASSERT_EQUAL(As(j+1,i+1),  REF[i][j]);
      }
      CPPUNIT_ASSERT_EQUAL(Am(i+1,j+1), REF[i][j]);
      // test correct multiplication
      double sum = 0.0;
      for (unsigned int t=0; t<r; t++){
          sum += As(i+1,t+1) * Am(t+1,j+1);
      }
      Csm_check(i+1,j+1) = sum;
      CPPUNIT_ASSERT_EQUAL(Csm(i+1,j+1),  Csm_check(i+1,j+1));
    }
  }

  // COLUMNVECTOR * DOUBLE
  ColumnVector Cc = Ac * v;
  for (unsigned int i=0; i<r; i++){
      // test if original elements were maintained
      CPPUNIT_ASSERT_EQUAL(Ac(i+1),  REF[0][i]);
      // test correct multiplication
      CPPUNIT_ASSERT_EQUAL(Cc(i+1),  REF[0][i] * v);
  }

  // ROWVECTOR * DOUBLE
  RowVector Cr = Ar * v;
  for (unsigned int i=0; i<c; i++){
      // test if original elements were maintained
      CPPUNIT_ASSERT_EQUAL(Ar(i+1),  REF[0][i]);
      // test correct multiplication
      CPPUNIT_ASSERT_EQUAL(Cr(i+1),  REF[0][i] * v);
  }

  // COLUMNVECTOR * ROWVECTOR
  Matrix Ccr = Ac * Ar;
  Matrix Ccr_check(r,c);
  for (unsigned int j=0; j<c; j++)  CPPUNIT_ASSERT_EQUAL(Ar(j+1), REF[0][j]);
  for (unsigned int i=0; i<r; i++){
     // test if original elements were maintained
     CPPUNIT_ASSERT_EQUAL(Ac(i+1),  REF[0][i]);
    for (unsigned int j=0; j<c; j++){
      // test correct multiplication
      Ccr_check(i+1,j+1) = Ac(i+1) * Ar(j+1);
      CPPUNIT_ASSERT_EQUAL(Ccr(i+1,j+1),  Ccr_check(i+1,j+1));
    }
  }

  // ROWVECTOR * COLUMNVECTOR
  double rc = Ac_trans * Ac;
  double rc_check;
  // test if original elements were maintained
  for (unsigned int j=0; j<c; j++)  CPPUNIT_ASSERT_EQUAL(Ac_trans(j+1), REF[0][j]);
  for (unsigned int i=0; i<r; i++)  CPPUNIT_ASSERT_EQUAL(Ac(i+1),  REF[0][i]);
  // test correct multiplication
  double sum = 0.0;
  for (unsigned int t=0; t<r; t++){
      sum += Ac_trans(t+1) * Ac(t+1);
  }
  rc_check = sum;
  CPPUNIT_ASSERT_EQUAL(rc,  rc_check);

  // ROWVECTOR * MATRIX
  // TODO: only implemented for lti
  //RowVector Cr2= Ar * Am;
  //Matrix Cr2_check(r);
  //for (unsigned int j=0; j<c; j++){
  //  // test if original elements were maintained
  //  CPPUNIT_ASSERT_EQUAL(Ar(j+1),  REF[0][j]);
  //  for (unsigned int i=0; i<r; i++){
  //    // test if original elements were maintained
  //    CPPUNIT_ASSERT_EQUAL(Am(i+1,j+1),  REF[i][j]);
  //  }
  //}
  //for (unsigned int i=0; i<r; i++){
  //  // test correct multiplication
  //  double sum = 0.0;
  //  for (unsigned int t=0; t<c; t++){
  //      sum += Ar(t+1) * Am(t+1,j+1);
  //  }
  //  Cr2_check(i+1) = sum;
  //  CPPUNIT_ASSERT_EQUAL(Cr2(i+1),  Cr2_check(i+1));
  //}

  // MATRIX * COLUMNVECTOR
  ColumnVector Cc2= Am_trans * Ac;
  ColumnVector Cc2_check(c);
  for (unsigned int j=0; j<r; j++){
    // test if original elements were maintained
    CPPUNIT_ASSERT_EQUAL(Ac(j+1),  REF[0][j]);
    for (unsigned int i=0; i<c; i++){
      // test if original elements were maintained
      CPPUNIT_ASSERT_EQUAL(Am_trans(i+1,j+1),  REF[j][i]);
    }
  }
  for (unsigned int i=0; i<c; i++){
    // test correct multiplication
    double sum = 0.0;
    for (unsigned int t=0; t<r; t++){
        sum += Am_trans(i+1,t+1) * Ac(t+1);
    }
    Cc2_check(i+1) = sum;
    CPPUNIT_ASSERT_EQUAL(Cc2(i+1),  Cc2_check(i+1));
  }


  // test operator /
  Cm = Am / v;
  for (unsigned int i=0; i<r; i++){
    for (unsigned int j=0; j<c; j++){
    //  CPPUNIT_ASSERT_EQUAL(Cm(i+1,j+1),  REF[i][j] / v);
      CPPUNIT_ASSERT_EQUAL(approxEqual(Cm(i+1,j+1),  REF[i][j] / v,epsilon),true);
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

  // SymmetricMatrix Rs(c);
  Matrix Rs(c,c);
  Rs(1,1) = 3;   Rs(1,2) = 5;    Rs(1,3) = 3;
  Rs(2,1) = 5;   Rs(2,2) = 2;    Rs(2,3) = 7;
  Rs(3,1) = 3;   Rs(3,2) = 7;    Rs(3,3) = 0;
  // SymmetricMatrix Rs_inv = Rs.inverse();
  Matrix Rs_inv = Rs.inverse();
  Matrix Is(c,c); Is = 0;
  for (unsigned int i=0; i<c; i++)
    Is(i+1,i+1) = 1;
  Matrix Is_test = Rs * Rs_inv;
  CPPUNIT_ASSERT_EQUAL(approxEqual(Is_test, Is, epsilon),true);

  // Issue #35
  SymmetricMatrix MI35(c);
  MI35(1,1) = 3; MI35(1,2) = 2; MI35(1,3) = 0;
  MI35(2,1) = 2; MI35(2,2) = 2; MI35(2,3) = 0;
  MI35(3,1) = 0; MI35(3,2) = 0; MI35(3,3) = 0.5;
  
  // test determinant
  CPPUNIT_ASSERT_EQUAL(approxEqual(Rm.determinant(), 105, epsilon),true);
  CPPUNIT_ASSERT_EQUAL(approxEqual(Rs.determinant(), 45, epsilon),true);
  CPPUNIT_ASSERT_EQUAL(approxEqual(MI35.determinant(), 1, epsilon),true); // Issue #35

  // test symmetric inverse
  SymmetricMatrix MI35_inv = MI35.inverse();

  SymmetricMatrix MI35_inv_test(c);
  MI35_inv_test(1,1) = 1.;  MI35_inv_test(1,2) = -1.; MI35_inv_test(1,3) = 0.;
  MI35_inv_test(2,1) = -1.; MI35_inv_test(2,2) = 1.5; MI35_inv_test(2,3) = 0.;
  MI35_inv_test(3,1) = 0.;  MI35_inv_test(3,2) = 0.;  MI35_inv_test(3,3) = 2.;

  CPPUNIT_ASSERT_EQUAL(approxEqual(MI35_inv_test, MI35_inv, epsilon), true);
 

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

  // test svd
  int rows = 4;
  int cols = 3;
  Matrix A_svd(rows,cols);
  Matrix W_svd(cols,cols);
  Matrix U_svd,V_svd;
  W_svd = 0.0;

  ColumnVector w_svd;
  A_svd(1,1)=1; A_svd(2,2)=2; A_svd(3,3)=3;
  A_svd(1,2)=-0.5; A_svd(1,3)=-0.8;
  A_svd(2,1)=-1.5; A_svd(2,3)=-2.8;
  A_svd(3,1)=2.5;  A_svd(3,2)=0.8;
  A_svd(4,1)=0.5;  A_svd(4,2)=1.8; A_svd(4,3)=1.6 ;

  A_svd.SVD(w_svd,U_svd,V_svd);
  for (int i=1; i<=A_svd.columns() ; i++)  W_svd(i,i) = w_svd(i);
  CPPUNIT_ASSERT_EQUAL(approxEqual(A_svd, U_svd * W_svd * V_svd.transpose(), epsilon),true);

  int rows2 = 3;
  int cols2 = 4;
  Matrix A2_svd(rows2,cols2);
  Matrix W2_svd(cols2,cols2);
  Matrix U2_svd,V2_svd;
  W2_svd = 0.0;

  ColumnVector w2_svd;
  A2_svd(1,1)=1; A2_svd(2,2)=2; A2_svd(3,3)=3; //A(4,4)=4;
  A2_svd(1,2)=-0.5; A2_svd(1,3)=-0.8; A2_svd(1,4)=-0.1 ;
  A2_svd(2,1)=-1.5; A2_svd(2,3)=-2.8; A2_svd(2,4)=3.1 ;
  A2_svd(3,1)=2.5;  A2_svd(3,2)=-0.8; A2_svd(3,4)=1.1 ;

  A2_svd.SVD(w2_svd,U2_svd,V2_svd);
  for (int i=1; i<=A2_svd.columns() ; i++)  W2_svd(i,i) = w2_svd(i);
  CPPUNIT_ASSERT_EQUAL(approxEqual(A2_svd, U2_svd * W2_svd * V2_svd.transpose(), epsilon),true);

  // TEST SPECIAL CASES
  // Inverse for 1x1 Matrix
  Matrix M1(1,1);
  M1(1,1)= 1.4;
  Matrix M1_inv = M1.inverse();
  Matrix I1(1,1);
  I1(1,1)= 1.0;
  CPPUNIT_ASSERT_EQUAL(M1_inv * M1, I1);
  // Inverse for 2x2 Matrix
  Matrix M2(2,2);
  M2(1,1)= 1.4;
  M2(2,2)= 0.4;
  M2(1,2)= 2.1;
  M2(2,1)= -0.8;
  Matrix M2_inv = M2.inverse();
  Matrix I2(2,2);
  I2=0.0;
  I2(1,1)= 1.0;
  I2(2,2)= 1.0;
  CPPUNIT_ASSERT_EQUAL(approxEqual(M2_inv * M2, I2,epsilon),true);
  // Determinant for 1x1 Matrix
  CPPUNIT_ASSERT_EQUAL(M1.determinant(), M1(1,1));
  // Determinant for 2x2 Matrix
  CPPUNIT_ASSERT_EQUAL(M2.determinant(), M2(1,1)*M2(2,2)-M2(1,2)*M2(2,1));
  // Inverse for 1x1 SymmetricMatrix
  SymmetricMatrix SM1(1);
  SM1(1,1)= 1.4;
  SymmetricMatrix SM1_inv = SM1.inverse();
  CPPUNIT_ASSERT_EQUAL(Matrix(SM1_inv * SM1), I1);
  // Inverse for 2x2 Matrix
  SymmetricMatrix SM2(2);
  SM2(1,1)= 1.4;
  SM2(2,2)= 0.4;
  SM2(1,2)= 2.1;
  SM2(2,1)= -0.8;
  SymmetricMatrix SM2_inv = SM2.inverse();
  CPPUNIT_ASSERT_EQUAL(approxEqual(Matrix(SM2_inv * SM2), I2,epsilon),true);
  // Determinant for 1x1 Matrix
  CPPUNIT_ASSERT_EQUAL(SM1.determinant(), SM1(1,1));
  // Determinant for 2x2 Matrix
  CPPUNIT_ASSERT_EQUAL(SM2.determinant(), SM2(1,1)*SM2(2,2)-SM2(1,2)*SM2(2,1));

  Matrix M3(3,3);
  M3(1,1)=1;
  M3(1,2)=2;
  M3(1,3)=3;
  M3(2,1)=4;
  M3(2,2)=5;
  M3(2,3)=6;
  M3(3,1)=7;
  M3(3,2)=8;
  M3(3,3)=9;
  // test rowCopy()
  RowVector rcopy1 = M3.rowCopy(1);
  RowVector rcopy1test(3);
  rcopy1test(1) = M3(1,1);
  rcopy1test(2) = M3(1,2);
  rcopy1test(3) = M3(1,3);
  CPPUNIT_ASSERT_EQUAL(rcopy1,rcopy1test);
  RowVector rcopy2 = M3.rowCopy(2);
  RowVector rcopy2test(3);
  rcopy2test(1) = M3(2,1);
  rcopy2test(2) = M3(2,2);
  rcopy2test(3) = M3(2,3);
  CPPUNIT_ASSERT_EQUAL(rcopy2,rcopy2test);
  RowVector rcopy3 = M3.rowCopy(3);
  RowVector rcopy3test(3);
  rcopy3test(1) = M3(3,1);
  rcopy3test(2) = M3(3,2);
  rcopy3test(3) = M3(3,3);
  CPPUNIT_ASSERT_EQUAL(rcopy3,rcopy3test);

  // test columnCopy()
  ColumnVector copy1 = M3.columnCopy(1);
  ColumnVector copy1test(3);
  copy1test(1) = M3(1,1);
  copy1test(2) = M3(2,1);
  copy1test(3) = M3(3,1);
  CPPUNIT_ASSERT_EQUAL(copy1,copy1test);
  ColumnVector copy2 = M3.columnCopy(2);
  ColumnVector copy2test(3);
  copy2test(1) = M3(1,2);
  copy2test(2) = M3(2,2);
  copy2test(3) = M3(3,2);
  CPPUNIT_ASSERT_EQUAL(copy2,copy2test);
  ColumnVector copy3 = M3.columnCopy(3);
  ColumnVector copy3test(3);
  copy3test(1) = M3(1,3);
  copy3test(2) = M3(2,3);
  copy3test(3) = M3(3,3);
  CPPUNIT_ASSERT_EQUAL(copy3,copy3test);
}



