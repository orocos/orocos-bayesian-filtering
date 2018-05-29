// Copyright (C) 2007 Tinne De Laet <first dot last at mech dot kuleuven dot be>
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

#include "model_test_ginac.hpp"
#include <cmath> // For sinus
#include <ginac/ginac.h>

#define MU 0.0
#define SIGMA 0.5

// Registers the fixture into the 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION( ModelTestGinac );

using namespace BFL;

void
ModelTestGinac::setUp()
{
}

void
ModelTestGinac::tearDown()
{
}

void
ModelTestGinac::testNonLinearAnalyticSystemModelGaussianUncertaintyGinac()
{
  GiNaC::symbol x0("x0"), x1("x1"), x2("x2"), x3("x3"), x4("x4"), x5("x5");
  GiNaC::symbol u0("u0"), u1("u1");

  // *** FUNCTION ***
  GiNaC::matrix f_sym(6,1);
  f_sym(0,0) = x0*x1+sin(u0);
  f_sym(1,0) = u1*x3*x2+u0*x0*x3;
  f_sym(2,0) = 2;
  f_sym(3,0) = x3*x3;
  f_sym(4,0) = x1;
  f_sym(5,0) = u0*x1;

  // *** STATE ***
  vector<GiNaC::symbol> x_sym(6);
  ColumnVector x_num(6);
  x_num = 3;
  x_sym[0] = x0;
  x_sym[1] = x1;
  x_sym[2] = x2;
  x_sym[3] = x3;
  x_sym[4] = x4;
  x_sym[5] = x5;

  // *** INPUT ***
  vector<GiNaC::symbol> u_sym(2);
  ColumnVector u_num(2);
  u_num = 3.4;
  u_sym[0] = u0;
  u_sym[1] = u1;

  // *** NOISE ON SYSTEM ***
  ColumnVector mu(6);  SymmetricMatrix sigma(6);
  mu = MU;
  for (int index_sigma_rows=0; index_sigma_rows < 6; index_sigma_rows++)
  {
      for (int index_sigma_cols=0; index_sigma_cols < 6; index_sigma_cols++)
       	{
    	  if (index_sigma_cols == index_sigma_rows)
    	    sigma(index_sigma_rows+1,index_sigma_cols+1)=SIGMA;
    	}
  }
  Gaussian My_Noise(mu,sigma);

 // *** SYSTEM MODEL ***

  NonLinearAnalyticConditionalGaussian_Ginac pdf(f_sym, u_sym, x_sym, My_Noise);
  NonLinearAnalyticSystemModelGaussianUncertainty_Ginac a_nonLinSysModel(&pdf);

  /* FunctionGet */
// a_nonLinSysModel.FunctionGet();
 //CPPUNIT_ASSERT_EQUAL(f_sym , a_nonLinSysModel.FunctionGet());
 //vector<GiNaC::symbol> x_sym_val = a_nonLinSysModel.StateGet();
 //CPPUNIT_ASSERT_EQUAL( x_sym.size() , x_sym_val.size());
 // /* StateGet */
 // for(int i = 0 ;i < x_sym.size() ; i++)
 // {
 //    CPPUNIT_ASSERT_EQUAL(x_sym[i] ,x_sym_val[i]);
 // }

  /* InputGet */
  //CPPUNIT_ASSERT_EQUAL(u_sym , a_nonLinSysModel.InputGet());

  /* dfdxGet */
//  CPPUNIT_ASSERT_EQUAL(pdf.dfGet(0) , a_nonLinSysModel.df_dxGet());


  // test ExpectedValueGet()
  // test CovarianceGet()
  // test HGet()
  //cout << "J = \n" << a_nonLinSysModel.PredictionGet(u_num, x_num) << endl;
  //cout << "Q = \n" << a_nonLinSysModel.CovarianceGet(u_num, x_num) << endl;
  //cout << "F = \n" << a_nonLinSysModel.df_dxGet(u_num, x_num);
}

void
ModelTestGinac::testNonLinearAnalyticMeasurementModelGaussianUncertaintyGinac()
{
}
