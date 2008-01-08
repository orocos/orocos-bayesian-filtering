// $Id: test_discrete_filter.cpp tdelaet $
// Copyright (C) 2007 Tinne De Laet <first dot last at mech dot kuleuven dot be>
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


#include <filter/histogramfilter.h>

#include <model/discretesystemmodel.h>

#include <pdf/discretepdf.h>
#include "conditionalUniformMeasPdf1d.h"
#include "../mobile_robot.h"

#include <iostream>
#include <fstream>

using namespace MatrixWrapper;
using namespace BFL;
using namespace std;



/* The purpose of this program is to construct a histogram filter for a very
   simple example, which is a 1d mobile robot localisation. 
   Furthermore the mobile robot is equipped with an ultrasonic sensor which
   directly measures the robot's position.
  
  The necessary SYSTEM MODEL is discrete.
  
  The MEASUREMENT MODEL is a home made one which takes into account the discrete nature
  of the state under consideration and some extra gaussian measurement noise
  
*/


int main(int argc, char** argv)
{
  cerr << "==================================================" << endl
       << "Test of histogram filter" << endl
       << "1D Mobile robot localisation example" << endl
       << "==================================================" << endl;


  int num_states = 20;
  int num_cond_args = 1;
  /****************************
  * Discrete system model      *
  ***************************/
  int cond_arg_dims[num_cond_args];
  cond_arg_dims[0] = num_states;
  DiscreteConditionalPdf sys_pdf(num_states,num_cond_args,cond_arg_dims);  // no  inputs
  std::vector<int> cond_args(num_cond_args);
  double prob_diag = 0.9;
  double prob_nondiag = (1-prob_diag)/(num_states-1);
  for (int state_kMinusOne = 0 ; state_kMinusOne < num_states ;  state_kMinusOne++)
    {
       cond_args[0] = state_kMinusOne;
       for (int state_k = 0 ; state_k < num_states ;  state_k++)
         {
           if (state_kMinusOne == state_k) sys_pdf.ProbabilitySet(prob_diag,state_k,cond_args);
           else sys_pdf.ProbabilitySet(prob_nondiag,state_k,cond_args);
         }
    }
  DiscreteSystemModel sys_model(&sys_pdf);
  std::cerr << "discrete system model created" << std::endl;

  /*********************************
   * Initialise measurement model *
   ********************************/

  // Construct the measurement noise (a scalar in this case)
  ColumnVector measNoise_Mu(1);
  measNoise_Mu(1) = 0.0;

  SymmetricMatrix measNoise_Cov(1);
  measNoise_Cov(1,1) = pow(1.0,2);
  Gaussian measurement_uncertainty(measNoise_Mu, measNoise_Cov);

  // create the model
  ConditionalUniformMeasPdf1d meas_pdf(measurement_uncertainty);
  MeasurementModel<MatrixWrapper::ColumnVector,int> meas_model(&meas_pdf);
  std::cout << "measurement model created" << std::endl;

  /****************************
  * Uniform prior DENSITY     *
  ***************************/
  DiscretePdf prior(num_states); //equal probability is set for all states
  std::cout << "uniform prior density created" << std::endl;

  /******************************
   * Construction of the Filter *
   ******************************/
  HistogramFilter filter(&prior);
  DiscretePdf * prior_test = filter.PostGet();
  std::cout << "filter created" << std::endl;
  
  /***************************
   * Initialise MOBILE ROBOT *
   **************************/
  // Model of mobile robot in world with one wall
  // The model is used to simultate the distance measurements.
  MobileRobot mobile_robot;
  ColumnVector input(2);
  input(1) = 0.1;
  input(2) = 0.0;


  /*******************
   * ESTIMATION LOOP *
   *******************/
  cout << "MAIN: Starting estimation" << endl;
  unsigned int time_step;
  for (time_step = 0; time_step < 100; time_step++)
    {
      // DO ONE STEP WITH MOBILE ROBOT
      mobile_robot.Move(input);

      // DO ONE MEASUREMENT
      ColumnVector measurement = mobile_robot.Measure();
     
      // UPDATE FILTER                                      
      filter.Update(&sys_model,&meas_model,measurement);
      //filter.Update(&sys_model);
      //filter.Update(&meas_model,measurement);

    } // estimation loop

  

  DiscretePdf *  posterior = filter.PostGet();
  cout << "After " << time_step+1 << " timesteps " << endl;
  cout << " Posterior probabilities = " << endl;
  for (int state = 0 ;  state< posterior->NumStatesGet() ; state++) cout << state << ": " << (double)(posterior->ProbabilityGet(state)) << endl;
  

  cout << "======================================================" << endl
       << "End of the Histogram filter for 1D mobile robot localisation" << endl
       << "======================================================"
       << endl;


  return 0;
}
