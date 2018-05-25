// Copyright (C) 2007 Wim Meeussen <wim DOT meeussen AT mech DOT kuleuven DOT be>
// Copyright (C) 2008 Tinne De Laet <tinne DOT delaet AT mech DOT kuleuven DOT be>
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
Complete_FilterTest::testComplete_FilterValue_Cont()
{
  double epsilon       = 0.015;
  double epsilon_large = 0.5;
  double epsilon_huge  = 2.0;

  /****************************
   * Initialise system model *
   ***************************/
  ColumnVector SysNoise_Mu(STATE_SIZE);
  SysNoise_Mu = 0.0;
  SysNoise_Mu(1) = MU_SYSTEM_NOISE_X;
  SysNoise_Mu(2) = MU_SYSTEM_NOISE_Y;
  SysNoise_Mu(3) = MU_SYSTEM_NOISE_THETA;

  SymmetricMatrix SysNoise_Cov(STATE_SIZE);
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
  prior_cont.SampleFrom(prior_samples,NUM_SAMPLES,SampleMthd::CHOLESKY,NULL);
  prior_discr.ListOfSamplesSet(prior_samples);

  // Mixture prior for the Mixture Boostrap filter
  ColumnVector prior_mu1(STATE_SIZE);
  SymmetricMatrix prior_sigma1(STATE_SIZE);
  prior_mu1(1) = PRIOR_MU_X1;
  prior_mu1(2) = PRIOR_MU_Y1;
  prior_mu1(STATE_SIZE) = PRIOR_MU_THETA1;
  prior_sigma1 = 0.0;
  prior_sigma1(1,1) = PRIOR_COV_X1;
  prior_sigma1(2,2) = PRIOR_COV_Y1;
  prior_sigma1(3,3) = PRIOR_COV_THETA1;
  Gaussian prior_cont1(prior_mu1,prior_sigma1);

  MCPdf<ColumnVector> mixcomp1(NUM_SAMPLES,STATE_SIZE);
  prior_cont1.SampleFrom(prior_samples,NUM_SAMPLES,SampleMthd::CHOLESKY,NULL);
  mixcomp1.ListOfSamplesSet(prior_samples);

  ColumnVector prior_mu2(STATE_SIZE);
  SymmetricMatrix prior_sigma2(STATE_SIZE);
  prior_mu2(1) = PRIOR_MU_X2;
  prior_mu2(2) = PRIOR_MU_Y2;
  prior_mu2(STATE_SIZE) = PRIOR_MU_THETA2;
  prior_sigma2 = 0.0;
  prior_sigma2(1,1) = PRIOR_COV_X2;
  prior_sigma2(2,2) = PRIOR_COV_Y2;
  prior_sigma2(3,3) = PRIOR_COV_THETA2;
  Gaussian prior_cont2(prior_mu2,prior_sigma2);

  MCPdf<ColumnVector> mixcomp2(NUM_SAMPLES,STATE_SIZE);
  prior_cont2.SampleFrom(prior_samples,NUM_SAMPLES,SampleMthd::CHOLESKY,NULL);
  mixcomp2.ListOfSamplesSet(prior_samples);

  ColumnVector prior_mu3(STATE_SIZE);
  SymmetricMatrix prior_sigma3(STATE_SIZE);
  prior_mu3(1) = PRIOR_MU_X3;
  prior_mu3(2) = PRIOR_MU_Y3;
  prior_mu3(STATE_SIZE) = PRIOR_MU_THETA3;
  prior_sigma3 = 0.0;
  prior_sigma3(1,1) = PRIOR_COV_X3;
  prior_sigma3(2,2) = PRIOR_COV_Y3;
  prior_sigma3(3,3) = PRIOR_COV_THETA3;
  Gaussian prior_cont3(prior_mu3,prior_sigma3);

  MCPdf<ColumnVector> mixcomp3(NUM_SAMPLES,STATE_SIZE);
  prior_cont3.SampleFrom(prior_samples,NUM_SAMPLES,SampleMthd::CHOLESKY,NULL);
  mixcomp3.ListOfSamplesSet(prior_samples);

  ColumnVector prior_mu4(STATE_SIZE);
  SymmetricMatrix prior_sigma4(STATE_SIZE);
  prior_mu4(1) = PRIOR_MU_X4;
  prior_mu4(2) = PRIOR_MU_Y4;
  prior_mu4(STATE_SIZE) = PRIOR_MU_THETA3;
  prior_sigma4 = 0.0;
  prior_sigma4(1,1) = PRIOR_COV_X4;
  prior_sigma4(2,2) = PRIOR_COV_Y4;
  prior_sigma4(3,3) = PRIOR_COV_THETA4;
  Gaussian prior_cont4(prior_mu4,prior_sigma4);

  MCPdf<ColumnVector> mixcomp4(NUM_SAMPLES,STATE_SIZE);
  prior_cont4.SampleFrom(prior_samples,NUM_SAMPLES,SampleMthd::CHOLESKY,NULL);
  mixcomp4.ListOfSamplesSet(prior_samples);

  vector<Pdf<ColumnVector>*> mixVec(3);
  mixVec[0] = &mixcomp1;
  mixVec[1] = &mixcomp2;
  mixVec[2] = &mixcomp3;
  //mixVec[3] = &mixcomp4;
  Mixture<ColumnVector> prior_mix(mixVec);

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
  // Model of mobile robot in world
  // The model is used to simultate the distance measurements.
  MobileRobot mobile_robot;
  ColumnVector input(INPUT_SIZE);
  input(1) = LIN_SPEED * DELTA_T;
  input(2)  = ROT_SPEED * DELTA_T;


  /******************************
   * Construction of the Filter *
   ******************************/
  Filter<ColumnVector,ColumnVector> *my_filter_extendedkalman, *my_filter_iteratedextendedkalman, *my_filter_bootstrap, *my_filter_ekparticle, *my_filter_mixtureBootstrap;
  my_filter_extendedkalman = new ExtendedKalmanFilter(&prior_cont);
  my_filter_iteratedextendedkalman = new IteratedExtendedKalmanFilter(&prior_cont,NUM_ITERATIONS);
  my_filter_bootstrap = new BootstrapFilter<ColumnVector,ColumnVector> (&prior_discr, RESAMPLE_PERIOD, RESAMPLE_THRESHOLD);
  my_filter_ekparticle = new EKParticleFilter(&prior_discr, 0, RESAMPLE_THRESHOLD);
  my_filter_mixtureBootstrap = new MixtureBootstrapFilter<ColumnVector,ColumnVector> (&prior_mix, RESAMPLE_PERIOD, RESAMPLE_THRESHOLD);

  /*******************
   * ESTIMATION LOOP *
   *******************/
  ColumnVector measurement ;
  ColumnVector mobile_robot_state ;
  Pdf<ColumnVector> * posterior_mixtureBootstrap;
  ofstream mixtureFile;
  if(OUTPUT_MIXTURE)
  {
    mixtureFile.open("mixtureOutput.txt");
  }

  cout << "Running 4 different filters. This may take a few minutes... " << endl;
  unsigned int time_step;
  for (time_step = 0; time_step < NUM_TIME_STEPS-1; time_step++)
    {
      // DO ONE STEP WITH MOBILE ROBOT
      mobile_robot.Move(input);

      // DO ONE MEASUREMENT
      measurement = mobile_robot.Measure();
      mobile_robot_state = mobile_robot.GetState();

      if(OUTPUT_MIXTURE)
      {
        posterior_mixtureBootstrap = my_filter_mixtureBootstrap->PostGet();
        vector<WeightedSample<ColumnVector> > los;
        vector<WeightedSample<ColumnVector> >::iterator los_it;
        int numComp = dynamic_cast<Mixture<ColumnVector> *>(posterior_mixtureBootstrap)->NumComponentsGet();
        mixtureFile << time_step << " " << numComp << " ";
        mixtureFile << mobile_robot_state(1) << " " << mobile_robot_state(2) << " " << mobile_robot_state(3) << " ";
        for(int i = 0 ; i<numComp ; i++ )
        {
            double componentWeight = ( dynamic_cast<Mixture<ColumnVector> *>(posterior_mixtureBootstrap)->WeightGet(i)) ;
            los = dynamic_cast<MCPdf<ColumnVector> *>( dynamic_cast<Mixture<ColumnVector> *>(posterior_mixtureBootstrap)->ComponentGet(i))->ListOfSamplesGet();
            mixtureFile << i << " " << componentWeight << " " << los.size()<< " " << STATE_SIZE << " ";
            for ( los_it=los.begin(); los_it != los.end() ; los_it++)
            {
                for (int j=0; j<STATE_SIZE ; j++)
                    mixtureFile << los_it->ValueGet()[j] << " ";
                mixtureFile<< los_it->WeightGet() << " ";
            }
        }
        mixtureFile<<endl;
       }
      // UPDATE FILTER
      my_filter_extendedkalman->Update(&sys_model,input,&meas_model, measurement);
      my_filter_iteratedextendedkalman->Update(&sys_model,input,&meas_model, measurement);
      my_filter_bootstrap->Update(&sys_model,input,&meas_model, measurement);
      //my_filter_ekparticle->Update(&sys_model,input,&meas_model, measurement);
      my_filter_mixtureBootstrap->Update(&sys_model,input,&meas_model, measurement);
    }

    if(OUTPUT_MIXTURE)
    {
      posterior_mixtureBootstrap = my_filter_mixtureBootstrap->PostGet();
      vector<WeightedSample<ColumnVector> > los;
      vector<WeightedSample<ColumnVector> >::iterator los_it;
      int numComp = dynamic_cast<Mixture<ColumnVector> *>(posterior_mixtureBootstrap)->NumComponentsGet();
      mixtureFile << time_step << " " << numComp << " ";
      mixtureFile << mobile_robot_state(1) << " " << mobile_robot_state(2) << " " << mobile_robot_state(3) << " ";
      for(int i = 0 ; i<numComp ; i++ )
      {
          double componentWeight = ( dynamic_cast<Mixture<ColumnVector> *>(posterior_mixtureBootstrap)->WeightGet(i)) ;
          los = dynamic_cast<MCPdf<ColumnVector> *>( dynamic_cast<Mixture<ColumnVector> *>(posterior_mixtureBootstrap)->ComponentGet(i))->ListOfSamplesGet();
          mixtureFile << i << " " << componentWeight << " " << los.size()<< " " << STATE_SIZE << " ";
          for ( los_it=los.begin(); los_it != los.end() ; los_it++)
          {
              for (int j=0; j<STATE_SIZE ; j++)
                  mixtureFile << los_it->ValueGet()[j] << " ";
              mixtureFile<< los_it->WeightGet() << " ";
          }
      }
      mixtureFile<<endl;
     }


  // ek_check
  Pdf<ColumnVector> * posterior_extendedkalman = my_filter_extendedkalman->PostGet();
  ColumnVector mean_ek_check(STATE_SIZE);
  mean_ek_check=mobile_robot.GetState();
  //mean_ek_check(1) = mobile_robot_state(1); mean_ek_check(2) = mobile_robot_state(2); mean_ek_check(3) = mobile_robot_state(3);
  SymmetricMatrix cov_ek_check(STATE_SIZE);
  cov_ek_check(1,1) = 0.0599729;   cov_ek_check(1,2) = 0.000291386; cov_ek_check(1,3) = 0.00223255;
  cov_ek_check(2,1) = 0.000291386; cov_ek_check(2,2) = 0.000277528; cov_ek_check(2,3) = 0.000644136;
  cov_ek_check(3,1) = 0.00223255;  cov_ek_check(3,2) = 0.000644136; cov_ek_check(3,3) = 0.00766009;
  CPPUNIT_ASSERT_EQUAL(approxEqual(posterior_extendedkalman->ExpectedValueGet(), mean_ek_check, epsilon_large),true);
  CPPUNIT_ASSERT_EQUAL(approxEqual(posterior_extendedkalman->CovarianceGet(), cov_ek_check, epsilon),true);

  // it_check
  Pdf<ColumnVector> * posterior_iteratedextendedkalman = my_filter_iteratedextendedkalman->PostGet();
  ColumnVector mean_it_check(STATE_SIZE);
  mean_it_check=mobile_robot.GetState();
  //mean_it_check(1) = mobile_robot_state(1); mean_it_check(2) = mobile_robot_state(2); mean_it_check(3) = mobile_robot_state(3);
  SymmetricMatrix cov_it_check(STATE_SIZE);
  cov_it_check = 0.0;
  cov_it_check(1,1) = 0.0611143;   cov_it_check(1,2) = 0.000315923; cov_it_check(1,3) = 0.00238938;
  cov_it_check(2,1) = 0.000315923; cov_it_check(2,2) = 0.000280736; cov_it_check(2,3) = 0.000665735;
  cov_it_check(3,1) = 0.00238938;  cov_it_check(3,2) = 0.000665735; cov_it_check(3,3) = 0.00775776;
  CPPUNIT_ASSERT_EQUAL(approxEqual(posterior_iteratedextendedkalman->ExpectedValueGet(), mean_it_check, epsilon_large),true);
  CPPUNIT_ASSERT_EQUAL(approxEqual(posterior_iteratedextendedkalman->CovarianceGet(), cov_it_check, epsilon),true);

  // bs_check
  Pdf<ColumnVector> * posterior_bootstrap = my_filter_bootstrap->PostGet();
  ColumnVector mean_bs_check(STATE_SIZE);
  mean_bs_check=mobile_robot.GetState();
  //mean_bs_check(1) = mobile_robot_state(1); mean_bs_check(2) = mobile_robot_state(2); mean_bs_check(3) = mobile_robot_state(3);
  SymmetricMatrix cov_bs_check(STATE_SIZE);
  cov_bs_check = 0.0;
  cov_bs_check(1,1) = PRIOR_COV_X;
  CPPUNIT_ASSERT_EQUAL(approxEqual(posterior_bootstrap->ExpectedValueGet(), mean_bs_check, epsilon_large),true);
  CPPUNIT_ASSERT_EQUAL(approxEqual(posterior_bootstrap->CovarianceGet(), cov_bs_check, epsilon),true);

  // ep_check
  /*
  Pdf<ColumnVector> * posterior_ekparticle = my_filter_ekparticle->PostGet();
  cout << " Posterior Mean = " << endl << posterior_ekparticle->ExpectedValueGet() << endl
       << " Covariance = " << endl << posterior_ekparticle->CovarianceGet() << "" << endl;
  ColumnVector mean_ep_check(STATE_SIZE);
  //mean_ep_check(1) = 6.64581; mean_ep_check(2) = -7.05499; mean_ep_check(3) = -0.76974;
  mean_ep_check=mobile_robot.GetState();
  SymmetricMatrix cov_ep_check(STATE_SIZE);
  cov_ep_check(1,1) = 0.0160492;   cov_ep_check(1,2) = 0.000193798; cov_ep_check(1,3) = 0.0013101;
  cov_ep_check(2,1) = 0.000193798; cov_ep_check(2,2) = 0.000289425; cov_ep_check(2,3) = 0.000701263;
  cov_ep_check(3,1) = 0.0013101;   cov_ep_check(3,2) = 0.000701263; cov_ep_check(3,3) = 0.00682061;
  cout << "mean_ep_check " << mean_ep_check << endl;
  cout << "cov_ep_check " << cov_ep_check << endl;
  CPPUNIT_ASSERT_EQUAL(approxEqual(posterior_ekparticle->ExpectedValueGet(), mean_ep_check, epsilon_huge),true);
  CPPUNIT_ASSERT_EQUAL(approxEqual(posterior_ekparticle->CovarianceGet(), cov_ep_check, epsilon_large),true);
  */
  // mixtureBoostrapFilter check
  posterior_mixtureBootstrap = my_filter_mixtureBootstrap->PostGet();
  ColumnVector mean_mbs_check(STATE_SIZE);
  //mean_mbs_check(1) = 6.64581; mean_mbs_check(2) = -7.05499; mean_mbs_check(3) = -0.76974;
  mean_mbs_check(1) = mobile_robot_state(1); mean_mbs_check(2) = mobile_robot_state(2); mean_mbs_check(3) = mobile_robot_state(3);
  //cout << "mixture weights:" << endl;
  vector<Probability> weights= dynamic_cast<Mixture<ColumnVector> *>(posterior_mixtureBootstrap)->WeightsGet();
  ColumnVector exp;
  for(int i = 0 ; i< dynamic_cast<Mixture<ColumnVector> *>(posterior_mixtureBootstrap)->NumComponentsGet(); i++ )
  {
        //cout << "weight component " << i << ": " << weights[i] << endl;
        exp= dynamic_cast<Mixture<ColumnVector> *>(posterior_mixtureBootstrap)->ComponentGet(i)->ExpectedValueGet();
        //cout << "expected value component " << i << ": " << exp << endl;
  }
  //cout << "expected value total: " << posterior_mixtureBootstrap->ExpectedValueGet() << endl;
  //cout << "should be : " << mean_mbs_check << endl;
  CPPUNIT_ASSERT_EQUAL(approxEqual(posterior_mixtureBootstrap->ExpectedValueGet(), mean_mbs_check, epsilon_huge),true);

  // closing file stream
  if(OUTPUT_MIXTURE)
  {
    mixtureFile.close();
  }

  // delete the filters
  delete my_filter_extendedkalman;
  delete my_filter_iteratedextendedkalman;
  delete my_filter_bootstrap;
  delete my_filter_ekparticle;
  delete my_filter_mixtureBootstrap;
}

void
Complete_FilterTest::testComplete_FilterValue_Discr()
{
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

  /*********************************
   * Initialise measurement model *
   ********************************/

  // Construct the measurement noise (a scalar in this case)
  ColumnVector measNoise_Mu(MEAS_SIZE);
  measNoise_Mu(1) = MU_MEAS_NOISE;

  SymmetricMatrix measNoise_Cov(MEAS_SIZE);
  measNoise_Cov(1,1) = SIGMA_MEAS_NOISE;
  Gaussian measurement_uncertainty(measNoise_Mu, measNoise_Cov);

  // create the model
  ConditionalUniformMeasPdf1d meas_pdf(measurement_uncertainty);
  MeasurementModel<MatrixWrapper::ColumnVector,int> meas_model(&meas_pdf);

  /****************************
  * Uniform prior DENSITY     *
  ***************************/
  DiscretePdf prior(num_states); //equal probability is set for all classed

  /******************************
   * Construction of the Filter *
   ******************************/
  HistogramFilter<ColumnVector> filter(&prior);
  DiscretePdf * prior_test = filter.PostGet();

  /***************************
   * initialise MOBILE ROBOT *
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
  unsigned int time_step;
  for (time_step = 0; time_step < NUM_TIME_STEPS-1; time_step++)
    {
      // DO ONE STEP WITH MOBILE ROBOT
      mobile_robot.Move(input);
      // DO ONE MEASUREMENT
      ColumnVector measurement = mobile_robot.Measure();
      // change sign of measurement (measurement model returns negative value)
      measurement(1) = 0-measurement(1);
      // UPDATE FILTER
      filter.Update(&sys_model,&meas_model,measurement);
    } // estimation loop

  // FIXME: This test needs more explanation...
  DiscretePdf *  posterior = filter.PostGet();
  for (int state=0; state< num_states; state++)
  {
    // std::cout << "state = " << state << " : " << "posterior->ProbabilityGet(state) = " << posterior->ProbabilityGet(state) << std::endl;
    if (state == (int)(round(mobile_robot.GetState()(2))) ){ // Y position  What does this comparison mean???
      CPPUNIT_ASSERT(posterior->ProbabilityGet(state) >0.9);
    }
    else {
      CPPUNIT_ASSERT(posterior->ProbabilityGet(state) <0.1);
    }
  }
}
