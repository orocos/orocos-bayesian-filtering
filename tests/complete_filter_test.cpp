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
 

#include "complete_filter_test.hpp"
#include "approxEqual.hpp"




// Registers the fixture into the 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION( Complete_FilterTest );

using namespace MatrixWrapper;
using namespace BFL;


void 
Complete_FilterTest::setUp()
{
}


void 
Complete_FilterTest::tearDown()
{
}

void 
Complete_FilterTest::testComplete_FilterValue()
{
  double epsilon       = 0.01;
  double epsilon_large = 0.5;
  double epsilon_huge  = 2.0;

  /****************************
   * Initialise system model *
   ***************************/
  ColumnVector SysNoise_Mu(STATE_SIZE);
  SymmetricMatrix SysNoise_Cov(STATE_SIZE);
  SysNoise_Mu = 0.0;
  SysNoise_Cov = 0.0;
  // Uncertainty or Noice (Additive) and Matrix A
  SysNoise_Cov(1,1) = SIGMA_SYSTEM_NOISE_X;
  SysNoise_Cov(2,2) = SIGMA_SYSTEM_NOISE_Y;
  SysNoise_Cov(3,3) = SIGMA_SYSTEM_NOISE_THETA;

  Gaussian System_Uncertainty(SysNoise_Mu, SysNoise_Cov);
  NonLinearAnalyticConditionalGaussianMobile sys_pdf(System_Uncertainty);
  AnalyticSystemModelGaussianUncertainty sys_model(&sys_pdf);


  /*********************************
   * Initialise measurement model *
   ********************************/
  // Fill up H
  double wall_ct = 2/(sqrt(pow(RICO_WALL,2.0) + 1));
  Matrix H(MEAS_SIZE,STATE_SIZE);
  H = 0.0;
  H(1,1) = wall_ct * RICO_WALL;
  H(1,2) = 0 - wall_ct;
  
  // Construct the measurement noise (a scalar in this case)
  ColumnVector MeasNoise_Mu(MEAS_SIZE);
  SymmetricMatrix MeasNoise_Cov(MEAS_SIZE);
  MeasNoise_Mu(1) = MU_MEAS_NOISE;
  MeasNoise_Cov(1,1) = SIGMA_MEAS_NOISE;
  
  Gaussian Measurement_Uncertainty(MeasNoise_Mu,MeasNoise_Cov);
  LinearAnalyticConditionalGaussian meas_pdf(H,Measurement_Uncertainty);
  LinearAnalyticMeasurementModelGaussianUncertainty meas_model(&meas_pdf);

  /****************************
   * Initialise prior DENSITY *
   ***************************/
  // Continuous Gaussian prior (for Kalman filters)
  ColumnVector prior_mu(STATE_SIZE);  
  SymmetricMatrix prior_sigma(STATE_SIZE);
  prior_mu(1) = PRIOR_MU_X;
  prior_mu(2) = PRIOR_MU_Y;
  prior_mu(STATE_SIZE) = PRIOR_MU_THETA;
  prior_sigma = 0.0;
  prior_sigma(1,1) = PRIOR_COV_X;
  prior_sigma(2,2) = PRIOR_COV_Y;
  prior_sigma(3,3) = PRIOR_COV_THETA;
  Gaussian prior_cont(prior_mu,prior_sigma); 
  
  // Discrete prior for Particle filter (using the continuous Gaussian prior)
  vector<Sample<ColumnVector> > prior_samples(NUM_SAMPLES);
  MCPdf<ColumnVector> prior_discr(NUM_SAMPLES,STATE_SIZE);
  prior_cont.SampleFrom(prior_samples,NUM_SAMPLES,CHOLESKY,NULL);
  prior_discr.ListOfSamplesSet(prior_samples);

  // check
  ColumnVector mean_check(STATE_SIZE);
  mean_check(1) = PRIOR_MU_X; mean_check(2) = PRIOR_MU_Y; mean_check(3) = PRIOR_MU_THETA;
  SymmetricMatrix cov_check(STATE_SIZE);
  cov_check(1,1) = PRIOR_COV_X; cov_check(1,2) = 0; cov_check(1,3) = 0;
  cov_check(2,1) = 0; cov_check(2,2) = PRIOR_COV_Y; cov_check(2,3) = 0;
  cov_check(3,1) = 0; cov_check(3,2) = 0; cov_check(3,3) = PRIOR_COV_THETA;
  CPPUNIT_ASSERT_EQUAL(approxEqual(prior_cont.ExpectedValueGet(), mean_check, epsilon),true);
  CPPUNIT_ASSERT_EQUAL(approxEqual(prior_cont.CovarianceGet(), cov_check, epsilon),true);
  
  
  /***************************
   * initialise MOBILE ROBOT *
   **************************/
  // Model of mobile robot in world with one wall
  // The model is used to simultate the distance measurements.
  MobileRobot mobile_robot;
  ColumnVector input(INPUT_SIZE);
  input(1) = LIN_SPEED * DELTA_T;
  input(2)  = ROT_SPEED * DELTA_T;


  /******************************
   * Construction of the Filter *
   ******************************/
  Filter<ColumnVector,ColumnVector> *my_filter_extendedkalman, *my_filter_iteratedextendedkalman, *my_filter_bootstrap, *my_filter_ekparticle;
  my_filter_extendedkalman = new ExtendedKalmanFilter(&prior_cont);
  my_filter_iteratedextendedkalman = new IteratedExtendedKalmanFilter(&prior_cont,NUM_ITERATIONS);
  my_filter_bootstrap = new BootstrapFilter<ColumnVector,ColumnVector> (&prior_discr, RESAMPLE_PERIOD, RESAMPLE_THRESHOLD);
  my_filter_ekparticle = new EKParticleFilter(&prior_discr, 0, RESAMPLE_THRESHOLD);

  /*******************
   * ESTIMATION LOOP *
   *******************/
  cout << "Running 4 different filters. This may take a few minutes... " << endl; 
  unsigned int time_step;
  for (time_step = 0; time_step < NUM_TIME_STEPS-1; time_step++)
    {
      // DO ONE STEP WITH MOBILE ROBOT
      mobile_robot.Move(input);

      // DO ONE MEASUREMENT
      ColumnVector measurement = mobile_robot.Measure();
     
      // UPDATE FILTER                                      
      my_filter_extendedkalman->Update(&sys_model,input,&meas_model, measurement);
      my_filter_iteratedextendedkalman->Update(&sys_model,input,&meas_model, measurement);
      my_filter_bootstrap->Update(&sys_model,input,&meas_model, measurement);
      //my_filter_ekparticle->Update(&sys_model,input,&meas_model, measurement);	  
    }


  // ek_check
  Pdf<ColumnVector> * posterior_extendedkalman = my_filter_extendedkalman->PostGet();
  ColumnVector mean_ek_check(STATE_SIZE);
  mean_ek_check(1) = 7.18713; mean_ek_check(2) = -7.15689; mean_ek_check(3) = -0.783556;
  SymmetricMatrix cov_ek_check(STATE_SIZE);
  cov_ek_check(1,1) = 0.0599729;   cov_ek_check(1,2) = 0.000291386; cov_ek_check(1,3) = 0.00223255;
  cov_ek_check(2,1) = 0.000291386; cov_ek_check(2,2) = 0.000277528; cov_ek_check(2,3) = 0.000644136;
  cov_ek_check(3,1) = 0.00223255;  cov_ek_check(3,2) = 0.000644136; cov_ek_check(3,3) = 0.00766009;
  CPPUNIT_ASSERT_EQUAL(approxEqual(posterior_extendedkalman->ExpectedValueGet(), mean_ek_check, epsilon_large),true);
  CPPUNIT_ASSERT_EQUAL(approxEqual(posterior_extendedkalman->CovarianceGet(), cov_ek_check, epsilon),true);

  // it_check
  Pdf<ColumnVector> * posterior_iteratedextendedkalman = my_filter_iteratedextendedkalman->PostGet();
  ColumnVector mean_it_check(STATE_SIZE);
  mean_it_check(1) = 7.00657; mean_it_check(2) = -7.28003; mean_it_check(3) = -0.773119;
  SymmetricMatrix cov_it_check(STATE_SIZE);
  cov_it_check(1,1) = 0.0611143;   cov_it_check(1,2) = 0.000315923; cov_it_check(1,3) = 0.00238938;
  cov_it_check(2,1) = 0.000315923; cov_it_check(2,2) = 0.000280736; cov_it_check(2,3) = 0.000665735;
  cov_it_check(3,1) = 0.00238938;  cov_it_check(3,2) = 0.000665735; cov_it_check(3,3) = 0.00775776;
  CPPUNIT_ASSERT_EQUAL(approxEqual(posterior_iteratedextendedkalman->ExpectedValueGet(), mean_it_check, epsilon_large),true);
  CPPUNIT_ASSERT_EQUAL(approxEqual(posterior_iteratedextendedkalman->CovarianceGet(), cov_it_check, epsilon),true);

  // bs_check
  Pdf<ColumnVector> * posterior_bootstrap = my_filter_bootstrap->PostGet();
  ColumnVector mean_bs_check(STATE_SIZE);
  mean_bs_check(1) = 6.64581; mean_bs_check(2) = -7.05499; mean_bs_check(3) = -0.76974;
  SymmetricMatrix cov_bs_check(STATE_SIZE);
  cov_bs_check(1,1) = 0.0160492;   cov_bs_check(1,2) = 0.000193798; cov_bs_check(1,3) = 0.0013101;
  cov_bs_check(2,1) = 0.000193798; cov_bs_check(2,2) = 0.000289425; cov_bs_check(2,3) = 0.000701263;
  cov_bs_check(3,1) = 0.0013101;   cov_bs_check(3,2) = 0.000701263; cov_bs_check(3,3) = 0.00682061;
  CPPUNIT_ASSERT_EQUAL(approxEqual(posterior_bootstrap->ExpectedValueGet(), mean_bs_check, epsilon_huge),true);
  CPPUNIT_ASSERT_EQUAL(approxEqual(posterior_bootstrap->CovarianceGet(), cov_bs_check, epsilon),true);

  // ep_check
  /*
  Pdf<ColumnVector> * posterior_ekparticle = my_filter_ekparticle->PostGet();
  cout << " Posterior Mean = " << endl << posterior_ekparticle->ExpectedValueGet() << endl
       << " Covariance = " << endl << posterior_ekparticle->CovarianceGet() << "" << endl;
  ColumnVector mean_ep_check(STATE_SIZE);
  mean_ep_check(1) = 6.64581; mean_ep_check(2) = -7.05499; mean_ep_check(3) = -0.76974;
  SymmetricMatrix cov_ep_check(STATE_SIZE);
  cov_ep_check(1,1) = 0.0160492;   cov_ep_check(1,2) = 0.000193798; cov_ep_check(1,3) = 0.0013101;
  cov_ep_check(2,1) = 0.000193798; cov_ep_check(2,2) = 0.000289425; cov_ep_check(2,3) = 0.000701263;
  cov_ep_check(3,1) = 0.0013101;   cov_ep_check(3,2) = 0.000701263; cov_ep_check(3,3) = 0.00682061;
  cout << "mean_ep_check " << mean_ep_check << endl;
  cout << "cov_ep_check " << cov_ep_check << endl;
  CPPUNIT_ASSERT_EQUAL(approxEqual(posterior_ekparticle->ExpectedValueGet(), mean_ep_check, epsilon_huge),true);
  CPPUNIT_ASSERT_EQUAL(approxEqual(posterior_ekparticle->CovarianceGet(), cov_ep_check, epsilon_large),true);
  */

  // delete the filters
  delete my_filter_extendedkalman;
  delete my_filter_iteratedextendedkalman;
  delete my_filter_bootstrap;
  delete my_filter_ekparticle;
}
