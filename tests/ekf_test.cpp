// Copyright (C) 2008 Willow Garage inc.
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

#include "ekf_test.hpp"
#include "approxEqual.hpp"

using namespace BFL;
using namespace MatrixWrapper;


// Registers the fixture into the 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION( EKFTest );

static const unsigned int state_size = 6;
static const unsigned int meas_size = 4;
static const double epsilon = 0.0001;




void EKFTest::setUp()
{
  assert(meas_size <= state_size);

  // create MEASUREMENT MODEL
  ColumnVector measNoise_Mu(meas_size);      measNoise_Mu = 0;
  SymmetricMatrix measNoise_Cov(meas_size);  measNoise_Cov = 0;
  Matrix meas_matrix(meas_size, state_size); meas_matrix = 0;
  for (unsigned int i=1; i<=meas_size; i++){
    measNoise_Mu(i) = 9.2-i;
    measNoise_Cov(i,i) = pow(2.8, 2);
    meas_matrix(i,i) = 2.8-3*i;
  }
  Gaussian measurement_Uncertainty(measNoise_Mu, measNoise_Cov);
  meas_pdf_   = new LinearAnalyticConditionalGaussian(meas_matrix, measurement_Uncertainty);
  meas_model_ = new LinearAnalyticMeasurementModelGaussianUncertainty(meas_pdf_);

  // create EKF filter
  ColumnVector prior_Mu(state_size);      prior_Mu = 0;
  SymmetricMatrix prior_Cov(state_size);  prior_Cov = 0;
  for (unsigned int i=1; i<=state_size; i++) {
    prior_Mu(i) = 29+i*i*i;
    prior_Cov(i,i) = pow(29.34*i+23.23, 2);
  }
  prior_  = new Gaussian(prior_Mu,prior_Cov);
  filter_ = new ExtendedKalmanFilter(prior_);
}



void EKFTest::tearDown()
{
  delete filter_;
  delete meas_model_;
  delete meas_pdf_;
  delete prior_;
}



void EKFTest::testMeasUpdate()
{
  ColumnVector meas(meas_size);
  for (unsigned int i=1; i<= meas_size; i++)
    meas(i) = 2.4*i;
  filter_->Update(meas_model_, meas);

  
  // expected value should be like this
  ColumnVector expectedvalue(state_size);
  expectedvalue(1) = 29.066225;
  expectedvalue(2) = 0.754136;
  expectedvalue(3) = -0.160365;
  expectedvalue(4) = -0.477823;
  expectedvalue(5) = 154.000000;
  expectedvalue(6) = 245.000000;
  CPPUNIT_ASSERT_EQUAL(approxEqual(expectedvalue, filter_->PostGet()->ExpectedValueGet(), epsilon), true);

  // covariance should be like this
  SymmetricMatrix covariance(state_size);  covariance = 0;
  covariance(1,1) = 183.019889;
  covariance(2,2) = 0.765538;
  covariance(3,3) = 0.203951;
  covariance(4,4) = 0.092627;
  covariance(5,5) = 28876.204900;
  covariance(6,6) = 39708.532900;
  CPPUNIT_ASSERT_EQUAL(approxEqual(covariance, filter_->PostGet()->CovarianceGet(), epsilon), true);
}
