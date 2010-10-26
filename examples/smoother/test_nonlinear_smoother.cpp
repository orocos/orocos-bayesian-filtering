// $Id: test_nonlinear_smoother.cpp 5925 2006-03-14 21:23:49Z tdelaet $
// Copyright (C) 2006 Tinne De Laet <first dot last at mech dot kuleuven dot be>
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


#include <filter/extendedkalmanfilter.h>
#include <filter/bootstrapfilter.h>

#include <model/linearanalyticsystemmodel_gaussianuncertainty.h>
#include <model/linearanalyticmeasurementmodel_gaussianuncertainty.h>

#include <pdf/analyticconditionalgaussian.h>
#include <pdf/analyticconditionalgaussian.h>
#include <pdf/linearanalyticconditionalgaussian.h>
#include <pdf/mcpdf.h>
#include <pdf/pdf.h>

#include "nonlinearanalyticconditionalgaussianmobile.h"
#include "mobile_robot.h"

#include <smoother/rauchtungstriebel.h>
#include <smoother/particlesmoother.h>

#include <iostream>
#include <fstream>

using namespace MatrixWrapper;
using namespace BFL;
using namespace std;



/* The purpose of this program is to construct a smoother (consisting of a
  forward and a backward filter) for the problem
  of localisation of a mobile robot equipped with an ultrasonic sensor.
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

#define KALMAN      1
#define BOOTSTRAP   2

int main(int argc, char** argv)
{
  cerr << "==================================================" << endl
       << "Test of different smooters" << endl
       << "Mobile robot localisation example" << endl
       << "==================================================" << endl;

  int smoother_name;
  if (!(argc== 2 ))
  {
    cout << "Please provide one argument. Possible arguments are:" << endl
	 << "  kalman_smoother" << endl
	 << "  bootstrap_smoother" << endl;
    return 0;
  }
  else {
    string argument = argv[1];
    if (argument == "kalman_filter")             smoother_name = KALMAN;
    else if (argument == "bootstrap_filter")     smoother_name = BOOTSTRAP;
    else
      {
	cout << "Please provide another argument. Possible arguments are:" << endl
	     << "  kalman_filter" << endl
	     << "  bootstrap_filter" << endl;
	return 0 ;
      }
  }
  /***********************
   * PREPARE FILESTREAMS *
   **********************/
  ofstream fout_time, fout_E, fout_cov, fout_meas, fout_states, fout_particles, fout_numparticles, fout_E_smooth, fout_cov_smooth, fout_time_smooth, fout_particles_smooth, fout_numparticles_smooth;

  fout_time.open("time.out");
  fout_E.open("E.out");
  fout_cov.open("cov.out");
  fout_meas.open("meas.out");
  fout_states.open("states.out");
  fout_E_smooth.open("Esmooth.out");
  fout_cov_smooth.open("covsmooth.out");
  fout_time_smooth.open("timesmooth.out");

  if (smoother_name  == BOOTSTRAP )
    {
      fout_particles.open("particles.out");
      fout_numparticles.open("numparticles.out");
      fout_particles_smooth.open("particlessmooth.out");
      fout_numparticles_smooth.open("numparticlessmooth.out");
    }


  /****************************
   * Linear system model      *
   ***************************/

  // create gaussian
  ColumnVector sysNoise_Mu(3);
  sysNoise_Mu(1) = 0.0;
  sysNoise_Mu(2) = 0.0;
  sysNoise_Mu(3) = 0.0;

  SymmetricMatrix sysNoise_Cov(3);
  sysNoise_Cov(1,1) = pow(0.05,2);
  sysNoise_Cov(1,2) = 0.0;
  sysNoise_Cov(1,3) = 0.0;
  sysNoise_Cov(2,1) = 0.0;
  sysNoise_Cov(2,2) = pow(0.05,2);
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
  prior_Mu(1) = 1.0;//-1.0;
  prior_Mu(2) = 0.0;//1.0;
  prior_Mu(3) = 0.0;
  SymmetricMatrix prior_Cov(3);
  prior_Cov(1,1) = 0.1;
  prior_Cov(1,2) = 0.0;
  prior_Cov(1,3) = 0.0;
  prior_Cov(2,1) = 0.0;
  prior_Cov(2,2) = 1.0;
  prior_Cov(2,3) = 0.0;
  prior_Cov(3,1) = 0.0;
  prior_Cov(3,2) = 0.0;
  prior_Cov(3,3) = pow(0.8,2);
  Gaussian prior_cont(prior_Mu,prior_Cov);

  int NUM_SAMPLES =  1000;
  vector<Sample<ColumnVector> > prior_samples(NUM_SAMPLES);
  MCPdf<ColumnVector> prior_discr(NUM_SAMPLES,3);
  prior_cont.SampleFrom(prior_samples,NUM_SAMPLES,CHOLESKY,NULL);
  prior_discr.ListOfSamplesSet(prior_samples);

   cout<< "Prior initialised"<< "" << endl;
   cout << "Prior Mean = " << endl << prior_cont.ExpectedValueGet() << endl
        << "Covariance = " << endl << prior_cont.CovarianceGet() << endl;

  /******************************
   * Construction of the Filter *
   ******************************/
  //ExtendedKalmanFilter filter(&prior_cont);
  Filter<ColumnVector,ColumnVector> * filter = NULL;

  switch (smoother_name){
  case KALMAN:{
    cout << "Using the Extended Kalman Filter" << endl;
    filter = new ExtendedKalmanFilter(&prior_cont);
    break;}
  case BOOTSTRAP:{
    cout << "Using the bootstrapfilter, " << NUM_SAMPLES << " samples, dynamic resampling" << endl;
    filter = new BootstrapFilter<ColumnVector,ColumnVector> (&prior_discr, 0, NUM_SAMPLES/4.0);
    break;}
  default:{
    cout << "Type if filter not recognised on construction" <<endl;
    return 0 ;}
  }


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
  unsigned int num_timesteps = 3;

  /***************************
   * vector in which all posteriors will be stored*
   **************************/
  // in case of gaussian
  vector<Gaussian> posteriors_gauss(num_timesteps);
  vector<Gaussian>::iterator posteriors_it_gauss  = posteriors_gauss.begin();
  // in case of mcpdf
  vector<MCPdf<ColumnVector> > posteriors_mcpdf(num_timesteps);
  vector<MCPdf<ColumnVector> >::iterator posteriors_it_mcpdf  = posteriors_mcpdf.begin();

  /*******************
   * ESTIMATION LOOP *
   *******************/
  cout << "MAIN: Starting estimation" << endl;
  unsigned int time_step;
  for (time_step = 0; time_step < num_timesteps; time_step++)
    {
      // write date in files
      fout_time << time_step << ";" << endl;
      fout_meas << mobile_robot.Measure()(1) << ";" << endl;
      fout_states << mobile_robot.GetState()(1) << "," << mobile_robot.GetState()(2) << ","
                  << mobile_robot.GetState()(3) << ";" << endl;

      // write posterior to file
      Pdf<ColumnVector> * posterior = filter->PostGet();
      fout_E << posterior->ExpectedValueGet()(1) << "," << posterior->ExpectedValueGet()(2)<< ","
             << posterior->ExpectedValueGet()(3) << ";"  << endl;
      fout_cov << posterior->CovarianceGet()(1,1) << "," << posterior->CovarianceGet()(1,2) << ","
               << posterior->CovarianceGet()(1,3) << "," << posterior->CovarianceGet()(2,2) << ","
               << posterior->CovarianceGet()(2,3) << "," << posterior->CovarianceGet()(3,3) << ";" << endl;

      // write particles to file
      if (smoother_name == BOOTSTRAP)
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

      if ((time_step+1)%5 == 0){
        // DO ONE MEASUREMENT
        ColumnVector measurement = mobile_robot.Measure();
        // UPDATE FILTER
        filter->Update(&sys_model,input,&meas_model,measurement);
      }
      else{
        filter->Update(&sys_model,input);
      }

      // make copy of posterior
      //*posteriors_it = (Gaussian)(*posterior);
      switch (smoother_name){
      case KALMAN:{
        *posteriors_it_gauss = *(Gaussian*)(posterior);
        posteriors_it_gauss++;
        break;}
      case BOOTSTRAP:{
        *posteriors_it_mcpdf = *(MCPdf<ColumnVector>*)(posterior);
        posteriors_it_mcpdf++;
        break;}
      default:{
        cout << "Type if filter not recognised on construction" <<endl;
        return 0 ;}
        }
    } // estimation loop


    Pdf<ColumnVector> * posterior = filter->PostGet();
    // write data in files
    fout_time << time_step << ";" << endl;
    fout_meas << mobile_robot.Measure()(1) << ";" << endl;
    fout_states << mobile_robot.GetState()(1) << "," << mobile_robot.GetState()(2) << ","
                << mobile_robot.GetState()(3) << ";" << endl;

    // write posterior to file
    fout_E << posterior->ExpectedValueGet()(1) << "," << posterior->ExpectedValueGet()(2)<< ","
           << posterior->ExpectedValueGet()(3) << ";"  << endl;
    fout_cov << posterior->CovarianceGet()(1,1) << "," << posterior->CovarianceGet()(1,2) << ","
             << posterior->CovarianceGet()(1,3) << "," << posterior->CovarianceGet()(2,2) << ","
             << posterior->CovarianceGet()(2,3) << "," << posterior->CovarianceGet()(3,3) << ";" << endl;
    // write particles to file
    if (smoother_name == BOOTSTRAP)
	{
	  fout_numparticles << ((MCPdf<ColumnVector>*)posterior)->NumSamplesGet() << ";" << endl;
	  vector<WeightedSample<ColumnVector> > samples_list = ((MCPdf<ColumnVector>*)posterior)->ListOfSamplesGet() ;
	  vector<WeightedSample<ColumnVector> >::iterator sample;
          for ( sample = samples_list.begin() ; sample != samples_list.end() ; sample++ )
	      fout_particles << sample->ValueGet()(1) << "," << sample->ValueGet()(2) << ","
                             << sample->ValueGet()(3) << "," <<  sample->WeightGet() <<";"<<endl;
	}

  cout << "After " << time_step+1 << " timesteps " << endl;
  cout << " Posterior Mean = " << endl << posterior->ExpectedValueGet() << endl
       << " Covariance = " << endl << posterior->CovarianceGet() << "" << endl;


  cout << "======================================================" << endl
       << "End of the filter for mobile robot localisation" << endl
       << "======================================================"
       << endl;


  /***************************************
   * Construction of the Backward Filter *
   **************************************/
  BackwardFilter<ColumnVector> * backwardfilter = NULL;

  switch (smoother_name){
  case KALMAN:{
    cout << "Using the RauchTungStriebelSmoother" << endl;
    backwardfilter = new RauchTungStriebel((Gaussian*)posterior);
    break;}
  case BOOTSTRAP:{
    cout << "Using the particlesmoother, " << NUM_SAMPLES << " samples" << endl;
    backwardfilter = new ParticleSmoother<ColumnVector>((MCPdf<ColumnVector>*)posterior);
    break;}
  default:{
    cout << "Type if filter not recognised on construction" <<endl;
    return 0 ;}
  }
  fout_time_smooth << time_step << ";" << endl;
  // write posterior to file
  fout_E_smooth << posterior->ExpectedValueGet()(1) << "," << posterior->ExpectedValueGet()(2)<< ","
         << posterior->ExpectedValueGet()(3) << ";"  << endl;
  fout_cov_smooth << posterior->CovarianceGet()(1,1) << "," << posterior->CovarianceGet()(1,2) << ","
           << posterior->CovarianceGet()(1,3) << "," << posterior->CovarianceGet()(2,2) << ","
           << posterior->CovarianceGet()(2,3) << "," << posterior->CovarianceGet()(3,3) << ";" << endl;
  // write particles to file
  if (smoother_name == BOOTSTRAP)
  {
    fout_numparticles_smooth << ((MCPdf<ColumnVector>*)posterior)->NumSamplesGet() << ";" << endl;
    vector<WeightedSample<ColumnVector> > samples_list = ((MCPdf<ColumnVector>*)posterior)->ListOfSamplesGet() ;
    vector<WeightedSample<ColumnVector> >::iterator sample;
        for ( sample = samples_list.begin() ; sample != samples_list.end() ; sample++ )
        fout_particles_smooth << sample->ValueGet()(1) << "," << sample->ValueGet()(2) << ","
                           << sample->ValueGet()(3) << "," <<  sample->WeightGet() <<";"<<endl;
  }

  /*******************
   * ESTIMATION LOOP *
   *******************/
  cout << "======================================================" << endl
       << "Start of the smoother for mobile robot localisation" << endl
       << "======================================================"
       << endl;
  for (time_step = num_timesteps-1; time_step+1 > 0 ; time_step--)
    {
      //Pdf<ColumnVector> filtered;
      switch (smoother_name){
      case KALMAN:{
        posteriors_it_gauss--;
        Gaussian filtered = *posteriors_it_gauss;
        backwardfilter->Update(&sys_model,input, &filtered);
        break;}
      case BOOTSTRAP:{
        posteriors_it_mcpdf--;
        MCPdf<ColumnVector> filtered = *posteriors_it_mcpdf;
        cout << "before update" << endl;
        backwardfilter->Update(&sys_model,input, &filtered);
        cout << "after update" << endl;
        break;}
      default:{
        cout << "Type if filter not recognised on construction" <<endl;
        return 0 ;}
        }
      // UPDATE  BACKWARDFILTER
      Pdf<ColumnVector> * posterior = backwardfilter->PostGet();

      fout_time_smooth << time_step << ";" << endl;
      // write posterior to file
      fout_E_smooth << posterior->ExpectedValueGet()(1) << "," << posterior->ExpectedValueGet()(2)<< ","
             << posterior->ExpectedValueGet()(3) << ";"  << endl;
      fout_cov_smooth << posterior->CovarianceGet()(1,1) << "," << posterior->CovarianceGet()(1,2) << ","
               << posterior->CovarianceGet()(1,3) << "," << posterior->CovarianceGet()(2,2) << ","
               << posterior->CovarianceGet()(2,3) << "," << posterior->CovarianceGet()(3,3) << ";" << endl;
      // write particles to file
      if (smoother_name == BOOTSTRAP)
	  {
	    fout_numparticles_smooth << ((MCPdf<ColumnVector>*)posterior)->NumSamplesGet() << ";" << endl;
	    vector<WeightedSample<ColumnVector> > samples_list = ((MCPdf<ColumnVector>*)posterior)->ListOfSamplesGet() ;
	    vector<WeightedSample<ColumnVector> >::iterator sample;
            for ( sample = samples_list.begin() ; sample != samples_list.end() ; sample++ )
	        fout_particles_smooth << sample->ValueGet()(1) << "," << sample->ValueGet()(2) << ","
                               << sample->ValueGet()(3) << "," <<  sample->WeightGet() <<";"<<endl;
	  }

      // make copy of posterior
      //*posteriors_it = (Gaussian)(*posterior);
      switch (smoother_name){
      case KALMAN:{
        *posteriors_it_gauss = *(Gaussian*)(posterior);
        break;}
      case BOOTSTRAP:{
        *posteriors_it_mcpdf = *(MCPdf<ColumnVector>*)(posterior);
        break;}
      default:{
        cout << "Type if filter not recognised on construction" <<endl;
        return 0 ;}
        }

    } // estimation loop


  // delete the filter
  delete filter;
  delete backwardfilter;

  cout << "======================================================" << endl
       << "End of the backward filter for mobile robot localisation" << endl
       << "======================================================"
       << endl;

  /****************************
   * CLOSE FILESTREAMS
   ***************************/
  fout_time.close();
  fout_E.close();
  fout_cov.close();
  fout_meas.close();
  fout_states.close();
  fout_time_smooth.close();
  fout_E_smooth.close();
  fout_cov_smooth.close();
  if (smoother_name  == BOOTSTRAP )
    {
      fout_particles.close();
      fout_numparticles.close();
      fout_particles_smooth.close();
      fout_numparticles_smooth.close();
    }
  return 0;
}
