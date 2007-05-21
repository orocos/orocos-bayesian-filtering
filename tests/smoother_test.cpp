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
// Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
//  
 
#include "smoother_test.hpp"
#include "approxEqual.hpp"
#include <filter/extendedkalmanfilter.h>
#include <model/linearanalyticsystemmodel_gaussianuncertainty.h>
#include <model/linearanalyticmeasurementmodel_gaussianuncertainty.h>
#include <pdf/analyticconditionalgaussian.h>
#include <pdf/analyticconditionalgaussian.h>
#include <pdf/linearanalyticconditionalgaussian.h>
#include <smoother/rauchtungstriebel.h>
#include <smoother/particlesmoother.h>

#include "../examples/compare_filters/mobile_robot_wall_cts.h"
#include "../examples/compare_filters/nonlinearanalyticconditionalgaussianmobile.h"
#include "../examples/mobile_robot.h"

// Registers the fixture into the 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION( SmootherTest );
using namespace BFL;


void 
SmootherTest::setUp()
{
}

void 
SmootherTest::tearDown()
{
}

void 
SmootherTest::testKalmanSmoother()
{
  double epsilon       = 0.01;
  double epsilon_large = 0.5;
  double epsilon_huge  = 2.0;
    /* The purpose of this program is to construct a filter for the problem
      of localisation of a mobile robot equipped with an ultrasonic sensor. 
      In this case the orientation is known, which simplifies the model considerably:
      The system model will become linear.
      The ultrasonic measures the distance to the wall (it can be switched off:
      see mobile_robot_wall_cts.h)
      
      The necessary SYSTEM MODEL is:
      
      x_k      = x_{k-1} + v_{k-1} * cos(theta) * delta_t
      y_k      = y_{k-1} + v_{k-1} * sin(theta) * delta_t
    
      The used MEASUREMENT MODEL:
      measuring the (perpendicular) distance z to the wall y = ax + b
      
      set WALL_CT = 1/sqrt(pow(a,2) + 1)
      z = WALL_CT * a * x - WALL_CT * y + WALL_CT * b + GAUSSIAN_NOISE
      or Z = H * X_k + J * U_k
      
      where
           
      H = [ WALL_CT * a       - WALL_CT      0 ]
      and GAUSSIAN_NOISE = N((WALL_CT * b), SIGMA_MEAS_NOISE)
      
    */
    
     /****************************
      * NonLinear system model      *
      ***************************/
    
      // create gaussian
      ColumnVector sysNoise_Mu(3);
      sysNoise_Mu(1) = 0.0;
      sysNoise_Mu(2) = 0.0;
      sysNoise_Mu(3) = 0.0;
    
      SymmetricMatrix sysNoise_Cov(3);
      sysNoise_Cov(1,1) = pow(0.01,2);
      sysNoise_Cov(1,2) = 0.0;
      sysNoise_Cov(1,3) = 0.0;
      sysNoise_Cov(2,1) = 0.0;
      sysNoise_Cov(2,2) = pow(0.01,2);
      sysNoise_Cov(2,3) = 0.0;
      sysNoise_Cov(3,1) = 0.0;
      sysNoise_Cov(3,2) = 0.0;
      sysNoise_Cov(3,3) = pow(0.03,2);
    
      Gaussian system_Uncertainty(sysNoise_Mu, sysNoise_Cov);
    
      // create the model
      NonLinearAnalyticConditionalGaussianMobile sys_pdf(system_Uncertainty);
      AnalyticSystemModelGaussianUncertainty sys_model(&sys_pdf);
    
      /*********************************
       * Initialise measurement model *
       ********************************/
    
      // create matrix H for linear measurement model
      Matrix H(1,3);
      H(1,1) = 0.0;
      H(1,2) = 2.0;
      H(1,3) = 0;
      // Construct the measurement noise (a scalar in this case)
      ColumnVector measNoise_Mu(1);
      measNoise_Mu(1) = 0.0;
    
      SymmetricMatrix measNoise_Cov(1);
      measNoise_Cov(1,1) = pow(0.05,2);
      Gaussian measurement_Uncertainty(measNoise_Mu, measNoise_Cov);
    
      // create the model
      LinearAnalyticConditionalGaussian meas_pdf(H, measurement_Uncertainty);
      LinearAnalyticMeasurementModelGaussianUncertainty meas_model(&meas_pdf);
    
    
      /****************************
       * Linear prior DENSITY     *
       ***************************/
       // Continuous Gaussian prior (for Kalman filters)
      ColumnVector prior_Mu(3);
      prior_Mu(1) = 0.0;//-1.0;
      prior_Mu(2) = 0.0;//1.0;
      prior_Mu(3) = 0.0;
      SymmetricMatrix prior_Cov(3);
      prior_Cov(1,1) = 1.0;
      prior_Cov(1,2) = 0.0;
      prior_Cov(1,3) = 0.0;
      prior_Cov(2,1) = 0.0;  
      prior_Cov(2,2) = 1.0;
      prior_Cov(2,3) = 0.0;
      prior_Cov(3,1) = 0.0;  
      prior_Cov(3,2) = 0.0;
      prior_Cov(3,3) = pow(0.8,2);
      Gaussian prior_cont(prior_Mu,prior_Cov); 
    
      /******************************
       * Construction of the Filter *
       ******************************/
      ExtendedKalmanFilter filter(&prior_cont);
    
    
      /***************************
       * initialise MOBILE ROBOT *
       **************************/
      // Model of mobile robot in world with one wall
      // The model is used to simultate the distance measurements.
      MobileRobot mobile_robot;
      ColumnVector input(2);
      input(1) = 0.1;
      input(2) = 0.0;
    
    
      /***************************
       * number of timesteps*
       **************************/
      unsigned int num_timesteps = 100;
    
      /***************************
       * vector in which all posteriors will be stored*
       **************************/
      vector<Gaussian> posteriors(num_timesteps);
      vector<Gaussian>::iterator posteriors_it  = posteriors.begin();
    
      /*******************
       * ESTIMATION LOOP *
       *******************/
      unsigned int time_step;
      for (time_step = 0; time_step < num_timesteps; time_step++)
        {
    
          // write posterior to file
          Gaussian * posterior = (Gaussian*)(filter.PostGet());
    
          // DO ONE STEP WITH MOBILE ROBOT
          mobile_robot.Move(input);
    
          if ((time_step+1)%10 == 0){
            // DO ONE MEASUREMENT
            ColumnVector measurement = mobile_robot.Measure();
         
            // UPDATE FILTER                                      
            filter.Update(&sys_model,input,&meas_model,measurement);
          }
          else{
            filter.Update(&sys_model,input);
          }
    
          // make copy of posterior
          *posteriors_it = *posterior;
    
          posteriors_it++;
        } // estimation loop
    
      
          Pdf<ColumnVector> * posterior = filter.PostGet();
    
    
      /***************************************
       * Construction of the Backward Filter *
       **************************************/
    RauchTungStriebel backwardfilter((Gaussian*)posterior);
    
    
      /*******************
       * ESTIMATION LOOP *
       *******************/
      for (time_step = num_timesteps-1; time_step+1 > 0 ; time_step--)
        {
          posteriors_it--;
          // UPDATE  BACKWARDFILTER                                      
          Gaussian filtered(posteriors_it->ExpectedValueGet(),posteriors_it->CovarianceGet());
          backwardfilter.Update(&sys_model,input, &filtered);
          Pdf<ColumnVector> * posterior = backwardfilter.PostGet();
    
          // make copy of posterior
          posteriors_it->ExpectedValueSet(posterior->ExpectedValueGet());
          posteriors_it->CovarianceSet(posterior->CovarianceGet());
    
        } // estimation loop

  cout << posteriors_it->ExpectedValueGet() << endl;
  cout << posteriors_it->CovarianceGet() << endl;
  ColumnVector mean_smoother_check(STATE_SIZE);
  mean_smoother_check(1) = 0.75; mean_smoother_check(2) = 0.27; mean_smoother_check(3) = 0.82; 
  SymmetricMatrix cov_smoother_check(STATE_SIZE);
  cov_smoother_check(1,1) = 1.0;     cov_smoother_check(1,2) = -0.004; cov_smoother_check(1,3) = 0.006;
  cov_smoother_check(2,1) = -0.004; cov_smoother_check(2,2) = 0.005; cov_smoother_check(2,3) = -0.006;
  cov_smoother_check(3,1) = 0.006;  cov_smoother_check(3,2) = -0.006;  cov_smoother_check(3,3) = 0.01;
  CPPUNIT_ASSERT_EQUAL(approxEqual(posteriors_it->ExpectedValueGet(), mean_smoother_check, epsilon_large),true);
  //CPPUNIT_ASSERT_EQUAL(approxEqual(posteriors_it->CovarianceGet(), cov_smoother_check, epsilon_large),true);
    
}

void 
SmootherTest::testParticleSmoother()
{
 // no nice implementation found yet
}
