// $Id: test_compare_filters.cpp 5925 2006-03-14 21:23:49Z tdelaet $
// Copyright (C) 2006 Klaas Gadeyne <first dot last at gmail dot com>
//                    Tinne De Laet <first dot last at mech dot kuleuven dot be>
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

/* Demonstration program for the Bayesian Filtering Library.
   Mobile robot localization with respect to wall with different possibilities for filter
*/

#include <filter/bootstrapfilter.h>
#include <filter/extendedkalmanfilter.h>
#include <filter/iteratedextendedkalmanfilter.h>
#include <filter/asirfilter.h>
#include <filter/EKparticlefilter.h>

#include "nonlinearanalyticconditionalgaussianmobile.h"//added
#include <model/analyticsystemmodel_gaussianuncertainty.h>
#include <model/linearanalyticmeasurementmodel_gaussianuncertainty.h>
#include <pdf/analyticconditionalgaussian.h>
#include <pdf/gaussian.h>

#include <wrappers/matrix/matrix_wrapper.h>
#include <wrappers/matrix/vector_wrapper.h>

#include <iostream>
#include <fstream>
#include <string>

// Include file with properties
#include "mobile_robot_wall_cts.h"
// Include pdf specific for this example
#include "nonlinearanalyticconditionalgaussianmobile.h"
// Include the mobile_robot simulator
#include "../mobile_robot.h"

/* The purpose of this program is to test different filter for the problem
  of localisation of a mobile robot equipped with an ultrasonic sensor.
  The ultrasonic measures the distance to the wall (it can be switched off:
  see mobile_robot_wall_cts.h)

  The necessary SYSTEM MODEL is:

  x_k      = x_{k-1} + v_{k-1} * cos(theta_{k-1} * delta_t
  y_k      = y_{k-1} + v_{k-1} * sin(theta_{k-1} * delta_t
  theta_k  = theta_{k-1} + omega_{k-1} * delta_t

  The used MEASUREMENT MODEL:
  measuring the (perpendicular) distance z to the wall y = ax + b

  set WALL_CT = 1/sqrt(pow(a,2) + 1)
  z = WALL_CT * a * x - WALL_CT * y + WALL_CT * b + GAUSSIAN_NOISE
  or Z = H * X_k + J * U_k

  where

  H = [ WALL_CT * a       - WALL_CT      0 ]
  and GAUSSIAN_NOISE = N((WALL_CT * b), SIGMA_MEAS_NOISE)

*/

using namespace MatrixWrapper;
using namespace BFL;
using namespace std;

#define KALMAN      1
#define IE_KALMAN   2
#define BOOTSTRAP   3
#define EK_PARTICLE 4
#define ASIR        5

int main(int argc, char** argv)
{
  cerr << "==================================================" << endl
       << "Test of switching between different filters" << endl
       << "Mobile robot localisation example" << endl
       << "==================================================" << endl;


  int filter_name;
  if (!(argc== 2 ))
  {
    cout << "Please provide one argument. Possible arguments are:" << endl
	 << "  kalman_filter" << endl
	 << "  bootstrap_filter" << endl
	 << "  EK_particle_filter" << endl
	 << "  ASIR_filter" << endl
	 << "  IE_kalman_filter" << endl;
    return 0;
  }
  else {
    string argument = argv[1];
    if (argument == "kalman_filter")             filter_name = KALMAN;
    else if (argument == "IE_kalman_filter")     filter_name = IE_KALMAN;
    else if (argument == "bootstrap_filter")     filter_name = BOOTSTRAP;
    else if (argument == "EK_particle_filter")   filter_name = EK_PARTICLE;
    else if (argument == "ASIR_filter")          filter_name = ASIR;
    else
      {
	cout << "Please provide another argument. Possible arguments are:" << endl
	     << "  kalman_filter" << endl
	     << "  bootstrap_filter" << endl
	     << "  EK_particle_filter" << endl
	     << "  ASIR_filter" << endl
	     << "  IE_kalman_filter" << endl;
	return 0 ;
      }
  }

  /***********************
   * PREPARE FILESTREAMS *
   **********************/
  ofstream fout_time, fout_E, fout_cov, fout_meas, fout_states, fout_particles, fout_numparticles;

  fout_time.open("time.out");
  fout_E.open("E.out");
  fout_cov.open("cov.out");
  fout_meas.open("meas.out");
  fout_states.open("states.out");

  if ((filter_name == BOOTSTRAP) || (filter_name == EK_PARTICLE) ||(filter_name == ASIR))
    {
      fout_particles.open("particles.out");
      fout_numparticles.open("numparticles.out");
    }


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
  cout<< "Measurment model H = " << H << endl;

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

   cout<< "Prior initialised"<< "" << endl;
   cout << "Prior Mean = " << endl << prior_cont.ExpectedValueGet() << endl
        << "Covariance = " << endl << prior_cont.CovarianceGet() << endl;




  /***************************
   * initialise MOBILE ROBOT *
   **************************/
  // Model of mobile robot in world with one wall
  // The model is used to simultate the distance measurements.
  MobileRobot mobile_robot;
  ColumnVector input(INPUT_SIZE);
  input(1) = LIN_SPEED * DELTA_T;
  input(2)  = ROT_SPEED * DELTA_T;
  cout << "Mobile robot initialised"<< "" << endl;



  /******************************
   * Construction of the Filter *
   ******************************/
  Filter<ColumnVector,ColumnVector> * my_filter = NULL;

  switch (filter_name){
  case KALMAN:{
    cout << "Using the Extended Kalman Filter" << endl;
    my_filter = new ExtendedKalmanFilter(&prior_cont);
    break;}
  case IE_KALMAN:{
    cout << "Using the Iterated Extended Kalman Filter, " << NUM_ITERATIONS << " iterations" << endl;
    my_filter = new IteratedExtendedKalmanFilter(&prior_cont,NUM_ITERATIONS);
    break;}
  case BOOTSTRAP:{
    cout << "Using the bootstrapfilter, " << NUM_SAMPLES << " samples, dynamic resampling" << endl;
    my_filter = new BootstrapFilter<ColumnVector,ColumnVector> (&prior_discr, RESAMPLE_PERIOD, RESAMPLE_THRESHOLD);
    break;}
  case EK_PARTICLE:{
    cout << "Using the Extended Particle Kalman Filter, " << NUM_SAMPLES << " samples, dynamic resampling" << endl;
    my_filter = new EKParticleFilter(&prior_discr, 0, RESAMPLE_THRESHOLD);
    break;}
    //case ASIR:{
    //cout << "Using the ASIR-filter, " << NUM_SAMPLES << " samples, fixed period resampling" << endl;
    //my_filter = new ASIRFilter<ColumnVector,ColumnVector> (prior_discr, RESAMPLE_PERIOD, RESAMPLE_THRESHOLD);
    //break;}
  default:{
    cout << "Type if filter not recognised on construction" <<endl;
    return 0 ;}
  }



  /*******************
   * ESTIMATION LOOP *
   *******************/
  cout << "MAIN: Starting estimation" << endl;
  unsigned int time_step;
  for (time_step = 0; time_step < NUM_TIME_STEPS-1; time_step++)
    {
      // write date in files
      fout_time << time_step << ";" << endl;
      fout_meas << mobile_robot.Measure()(1) << ";" << endl;
      fout_states << mobile_robot.GetState()(1) << "," << mobile_robot.GetState()(2) << ","
                  << mobile_robot.GetState()(3) << ";" << endl;

      // write posterior to file
      Pdf<ColumnVector> * posterior = my_filter->PostGet();
      fout_E << posterior->ExpectedValueGet()(1) << "," << posterior->ExpectedValueGet()(2)<< ","
             << posterior->ExpectedValueGet()(3) << ";"  << endl;
      fout_cov << posterior->CovarianceGet()(1,1) << "," << posterior->CovarianceGet()(1,2) << ","
               << posterior->CovarianceGet()(1,3) << "," << posterior->CovarianceGet()(2,2) << ","
               << posterior->CovarianceGet()(2,3) << "," << posterior->CovarianceGet()(3,3) << ";" << endl;


      // write particles to file
      if ((filter_name == BOOTSTRAP) || (filter_name == EK_PARTICLE) ||(filter_name == ASIR))
	{
	  fout_numparticles << ((MCPdf<ColumnVector>*)posterior)->NumSamplesGet() << ";" << endl;
	  vector<WeightedSample<ColumnVector> > samples_list = ((MCPdf<ColumnVector>*)posterior)->ListOfSamplesGet() ;
	  vector<WeightedSample<ColumnVector> >::iterator sample;
          for ( sample = samples_list.begin() ; sample != samples_list.end() ; sample++ )
	      fout_particles << sample->ValueGet()(1) << "," << sample->ValueGet()(2) << ","
                             << sample->ValueGet()(3) << "," <<  sample->WeightGet() <<";"<<endl;
	}


      // DO ONE STEP WITH MOBILE ROBOT
      mobile_robot.Move(input);

      // DO ONE MEASUREMENT
      ColumnVector measurement = mobile_robot.Measure();

      // UPDATE FILTER
      if (USE_MEASUREMENTS)
          my_filter->Update(&sys_model,input,&meas_model, measurement);
      else
          my_filter->Update(&sys_model, input);
    } // estimation loop



  Pdf<ColumnVector> * posterior = my_filter->PostGet();
  cout << "After " << time_step+1 << " timesteps " << endl;
  cout << " Posterior Mean = " << endl << posterior->ExpectedValueGet() << endl
       << " Covariance = " << endl << posterior->CovarianceGet() << "" << endl;
  cout << "=============================================" << endl;


  // delete the filter
  delete my_filter;

  cout << "==================================================" << endl
       << "End of the Comparing filters test" << endl
       << "=================================================="
       << endl;



  /****************************
   * CLOSE FILESTREAMS
   ***************************/
  fout_time.close();
  fout_E.close();
  fout_cov.close();
  fout_meas.close();
  fout_states.close();
  if ((filter_name == BOOTSTRAP) || (filter_name == EK_PARTICLE) ||(filter_name == ASIR))
    {
      fout_particles.close();
      fout_numparticles.close();
    }
  return 0;
}
