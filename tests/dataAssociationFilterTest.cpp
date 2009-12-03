// Copyright (C) 2008 Tinne De Laet <first DOT last AT mech DOT kuleuven DOT be>
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
 

#include "dataAssociationFilterTest.hpp"
#include "approxEqual.hpp"




// Registers the fixture into the 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION( DataAssociationFilterTest );

using namespace MatrixWrapper;
using namespace BFL;


#define NUM_FILTERS 3
#define STATE_SIZE 2 // Default Number of Samples
#define NUM_SAMPLES 100 // Default Number of Samples
#define RESAMPLE_PERIOD 0 // Default Resample Period
#define RESAMPLE_THRESHOLD (NUM_SAMPLES/4.0) // Threshold for Dynamic Resampling
#define PRIOR_MU_X1 1
#define PRIOR_MU_Y1 1
#define PRIOR_MU_X2 2
#define PRIOR_MU_Y2 2
#define PRIOR_MU_X3 3
#define PRIOR_MU_Y3 3
#define PRIOR_MU_X4 2
#define PRIOR_MU_Y4 2
#define PRIOR_MU_X5 5
#define PRIOR_MU_Y5 5
#define PRIOR_MU_X6 6
#define PRIOR_MU_Y6 6
#define PRIOR_COV_X 0.02
#define PRIOR_COV_Y 0.02


void 
DataAssociationFilterTest::setUp()
{

}

void 
DataAssociationFilterTest::testDataAssociationFilterMCPdf()
{
  cout<< "TESTDATAASSOCIATIONFILTERMCPDF" << endl;
  BootstrapFilter<ColumnVector,ColumnVector> *my_filter_bootstrap1, *my_filter_bootstrap2, *my_filter_bootstrap3, *my_filter_bootstrap4;

  // Continuous Gaussian prior (for Kalman filters)
  ColumnVector prior_mu(STATE_SIZE);  
  SymmetricMatrix prior_sigma(STATE_SIZE);
  prior_mu(1) = PRIOR_MU_X1;
  prior_mu(2) = PRIOR_MU_Y1;
  prior_sigma = 0.0;
  prior_sigma(1,1) = PRIOR_COV_X;
  prior_sigma(2,2) = PRIOR_COV_Y;
  Gaussian prior_cont(prior_mu,prior_sigma); 

  MCPdf<ColumnVector> prior_discr(NUM_SAMPLES,STATE_SIZE);
  vector<Sample<ColumnVector> > prior_samples(NUM_SAMPLES);
  prior_cont.SampleFrom(prior_samples,NUM_SAMPLES,CHOLESKY,NULL);
  prior_discr.ListOfSamplesSet(prior_samples);

  // first filter    
  my_filter_bootstrap1 = new BootstrapFilter<ColumnVector,ColumnVector> (&prior_discr, RESAMPLE_PERIOD, RESAMPLE_THRESHOLD);

  prior_mu(1) = PRIOR_MU_X2;
  prior_mu(2) = PRIOR_MU_Y2;
  prior_cont.ExpectedValueSet(prior_mu);
  prior_cont.SampleFrom(prior_samples,NUM_SAMPLES,CHOLESKY,NULL);
  prior_discr.ListOfSamplesSet(prior_samples);

  // second filter    
  my_filter_bootstrap2 = new BootstrapFilter<ColumnVector,ColumnVector> (&prior_discr, RESAMPLE_PERIOD, RESAMPLE_THRESHOLD);

  prior_mu(1) = PRIOR_MU_X3;
  prior_mu(2) = PRIOR_MU_Y3;
  prior_cont.ExpectedValueSet(prior_mu);
  prior_cont.SampleFrom(prior_samples,NUM_SAMPLES,CHOLESKY,NULL);
  prior_discr.ListOfSamplesSet(prior_samples);

  // second filter    
  my_filter_bootstrap3 = new BootstrapFilter<ColumnVector,ColumnVector> (&prior_discr, RESAMPLE_PERIOD, RESAMPLE_THRESHOLD);

  vector<Filter<ColumnVector,ColumnVector>* > filters(NUM_FILTERS); 
  filters[0] = my_filter_bootstrap1;
  filters[1] = my_filter_bootstrap2;
  filters[2] = my_filter_bootstrap3;

  // constructor
  double gamma = 0.1;
  double treshold = 0.2;
  DataAssociationFilterMCPdf<ColumnVector,ColumnVector> my_dataAssociationFilter(filters,gamma,treshold);
  vector<Pdf<ColumnVector>* > filters_post = my_dataAssociationFilter.PostGet();
  vector<int > filterID = my_dataAssociationFilter.FilterIDGet();

  for(int i = 0; i< NUM_FILTERS ; i++)
  {
     CPPUNIT_ASSERT_EQUAL( filters[i]->PostGet() , filters_post[i]);
     CPPUNIT_ASSERT_EQUAL( i+1 , filterID[i]);
     cout << "filterID[" << i << "]" << filterID[i] << endl;
  }
    
  CPPUNIT_ASSERT_EQUAL( NUM_FILTERS , my_dataAssociationFilter.NumFiltersGet());
  CPPUNIT_ASSERT_EQUAL( gamma , my_dataAssociationFilter.GammaGet());
  CPPUNIT_ASSERT_EQUAL( treshold , my_dataAssociationFilter.TresholdGet());

  // REMOVE FILTER
  vector<Filter<ColumnVector,ColumnVector>* > filters_bis(NUM_FILTERS); 
  filters_bis[0] = my_filter_bootstrap1;
  filters_bis[1] = my_filter_bootstrap3;

  CPPUNIT_ASSERT_EQUAL( false , my_dataAssociationFilter.RemoveFilter(-1));
  CPPUNIT_ASSERT_EQUAL( false , my_dataAssociationFilter.RemoveFilter(3));
  CPPUNIT_ASSERT_EQUAL( true , my_dataAssociationFilter.RemoveFilter(1));
  CPPUNIT_ASSERT_EQUAL( NUM_FILTERS-1 , my_dataAssociationFilter.NumFiltersGet());
  filters_post = my_dataAssociationFilter.PostGet();
  filterID = my_dataAssociationFilter.FilterIDGet();
  
  // check if other filters are still there
  CPPUNIT_ASSERT_EQUAL( filters_bis[0]->PostGet() , filters_post[0]);
  CPPUNIT_ASSERT_EQUAL( filters_bis[1]->PostGet() , filters_post[1]);
  CPPUNIT_ASSERT_EQUAL( 1 , filterID[0]);
  CPPUNIT_ASSERT_EQUAL( 3 , filterID[1]);
  

  // ADD FILTER
  prior_mu(1) = PRIOR_MU_X4;
  prior_mu(2) = PRIOR_MU_Y4;
  prior_cont.ExpectedValueSet(prior_mu);
  prior_cont.SampleFrom(prior_samples,NUM_SAMPLES,CHOLESKY,NULL);
  prior_discr.ListOfSamplesSet(prior_samples);
  my_filter_bootstrap4 = new BootstrapFilter<ColumnVector,ColumnVector> (&prior_discr, RESAMPLE_PERIOD, RESAMPLE_THRESHOLD);
  vector<Filter<ColumnVector,ColumnVector>* > filters_bis2(NUM_FILTERS); 
  filters_bis2[0] = my_filter_bootstrap1;
  filters_bis2[1] = my_filter_bootstrap3;
  filters_bis2[2] = my_filter_bootstrap4;

  my_dataAssociationFilter.AddFilter(my_filter_bootstrap4);
  CPPUNIT_ASSERT_EQUAL( NUM_FILTERS , my_dataAssociationFilter.NumFiltersGet());

  filters_post = my_dataAssociationFilter.PostGet();
  filterID = my_dataAssociationFilter.FilterIDGet();
  
  // check if other filters are still there
  CPPUNIT_ASSERT_EQUAL( filters_bis2[0]->PostGet() , filters_post[0]);
  CPPUNIT_ASSERT_EQUAL( filters_bis2[1]->PostGet() , filters_post[1]);
  CPPUNIT_ASSERT_EQUAL( filters_bis2[2]->PostGet() , filters_post[2]);
  CPPUNIT_ASSERT_EQUAL( 1 , filterID[0]);
  CPPUNIT_ASSERT_EQUAL( 3 , filterID[1]);
  CPPUNIT_ASSERT_EQUAL( 4 , filterID[2]);

  for (int teller_filters = 0;  teller_filters<my_dataAssociationFilter.NumFiltersGet() ; teller_filters++)
  {
    cout << "filters_post["<< teller_filters << "] expected value " << filters_post[teller_filters]->ExpectedValueGet() << endl;
    //cout << "filters_post[0] covariance " << ( (Matrix)((MCPdf<ColumnVector>*)(filters_post[teller_filters]))->CovarianceGet()) << endl;
    cout << "filterID[" << teller_filters << "]" << filterID[teller_filters] << endl;
  }

  Matrix probs(3,3);
  probs = 0.0;
  probs(1,1) = 0.4;
  //probs(1,2) = 0.4;
  //probs(1,3) = 0.4;
  //probs(2,1) = 0.4;
  probs(2,2) = 0.4;
  //probs(2,3) = 0.4;
  //probs(2,1) = 0.4;
  //probs(2,2) = 0.4;
  //probs(2,3) = 0.4;
  //probs(3,1) = 0.4;
  //probs(3,2) = 0.4;
  probs(3,3) = 0.4;

/*
  vector<vector<int > > associations = my_dataAssociationFilter.GetAssociations(probs,3);

  std::cout << "ASSOCIATIONS  " << std::endl;
  for (int i =0 ; i < 35 ; i++)
  {
       std::cout << "Association " << i+1 << std::endl;
       for (int j=0 ; j < 3 ; j++)
            std::cout << associations[i][j] << std::endl;
  } 

*/
 // create measurement model 
 Matrix H(2,2);
 H = 0.0;
 H(1,1) = 1.0;
 H(2,2) = 1.0;

 ColumnVector measNoise_Mu(2);
 measNoise_Mu(1) = 0.0;
 measNoise_Mu(2) = 0.0;

 SymmetricMatrix measNoise_Cov(2);
 measNoise_Cov = 0.0;
 measNoise_Cov(1,1) = pow(0.1,2);
 measNoise_Cov(2,2) = pow(0.1,2);

 Gaussian measurement_Uncertainty(measNoise_Mu, measNoise_Cov);

 LinearAnalyticConditionalGaussian meas_pdf(H,measurement_Uncertainty);

 LinearAnalyticMeasurementModelGaussianUncertainty meas_model(&meas_pdf);

 vector<ColumnVector> measurements(3);
 ColumnVector measurement(2);
 measurement(1) = 0.9;
 measurement(2) = 0.9;
 measurements[0] = measurement; 
 measurement(1) = 3.2;
 measurement(2) = 3.2;
 measurements[1] = measurement; 
 measurement(1) = 2.2;
 measurement(2) = 2.2;
 measurements[2] = measurement; 
 ColumnVector s;

 std::cout << "call GetMeasProbs " << endl;
 Matrix meas_probs = my_dataAssociationFilter.GetMeasProbs(&meas_model,measurements,s);
 std::cout << "MeasProbs " << meas_probs << endl;

 std::cout << "call getassociationprobs " << endl;
 vector<vector<double > > association_probs = my_dataAssociationFilter.GetAssociationProbs(&meas_model,measurements, s);
 std::cout << "ASSOCIATION_PROBS  " << std::endl;
 for (int i =0 ; i < my_dataAssociationFilter.NumFiltersGet() ; i++)
 {
      for (int j=0 ; j < measurements.size() ; j++)
      {
           std::cout << "Feature " << i+1 << " with  Object " << j+1 << std::endl;
           std::cout << association_probs[i][j] << std::endl;
      }
 } 

 Matrix A(2,2);
 A = 0.0;
 A(1,1) = 1.0;
 A(2,2) = 1.0;

 vector<Matrix> AB(1);
 AB[0] = A;

 ColumnVector sysNoise_Mu(2);
 sysNoise_Mu(1) = 0.0;
 sysNoise_Mu(2) = 0.0;

 SymmetricMatrix sysNoise_Cov(2);
 sysNoise_Cov = 0.0;
 sysNoise_Cov(1,1) = pow(0.2,2);
 sysNoise_Cov(2,2) = pow(0.2,2);

 Gaussian system_Uncertainty(sysNoise_Mu, sysNoise_Cov);
 LinearAnalyticConditionalGaussian sys_pdf(AB, system_Uncertainty);
 LinearAnalyticSystemModelGaussianUncertainty sys_model(&sys_pdf);

  // test system update
  cout<< "call system update!" << endl;
  my_dataAssociationFilter.Update(&sys_model);
  filters_post = my_dataAssociationFilter.PostGet();
  cout<< "SYSTEM UPDATE DONE!" << endl;
  for (int teller_filters = 0;  teller_filters<my_dataAssociationFilter.NumFiltersGet() ; teller_filters++)
  {
    cout << "filters_post["<< teller_filters << "] expected value " << filters_post[teller_filters]->ExpectedValueGet() << endl;
    //cout << "filters_post[0] covariance " << ( (Matrix)((MCPdf<ColumnVector>*)(filters_post[teller_filters]))->CovarianceGet()) << endl;
  }

  // test measurement update
  my_dataAssociationFilter.Update(&meas_model,measurements);
  my_dataAssociationFilter.Update(&meas_model,measurements);
  filters_post = my_dataAssociationFilter.PostGet();
  cout<< "MEASUREMENT UPDATE DONE!" << endl;
  for (int teller_filters = 0;  teller_filters<my_dataAssociationFilter.NumFiltersGet() ; teller_filters++)
  {
    cout << "filters_post["<< teller_filters << "] expected value " << filters_post[teller_filters]->ExpectedValueGet() << endl;
    cout << "filters_post["<< teller_filters << "] covariance " << (filters_post[teller_filters]->CovarianceGet()) << endl;
  }

  association_probs = my_dataAssociationFilter.GetAssociationProbs(&meas_model,measurements, s);
  my_dataAssociationFilter.Update(&meas_model,measurements);
}

void 
DataAssociationFilterTest::testDataAssociationFilterGaussian()
{
  cout<< "TESTDATAASSOCIATIONFILTERGAUSSIAN" << endl;
  ExtendedKalmanFilter *my_filter_kalman1, *my_filter_kalman2, *my_filter_kalman3, *my_filter_kalman4, *my_filter_kalman5, *my_filter_kalman6;

  // Continuous Gaussian prior (for Kalman filters)
  ColumnVector prior_mu(STATE_SIZE);  
  SymmetricMatrix prior_sigma(STATE_SIZE);
  prior_mu(1) = PRIOR_MU_X1;
  prior_mu(2) = PRIOR_MU_Y1;
  prior_sigma = 0.0;
  prior_sigma(1,1) = PRIOR_COV_X;
  prior_sigma(2,2) = PRIOR_COV_Y;
  Gaussian prior_cont(prior_mu,prior_sigma); 

  // first filter    
  my_filter_kalman1 = new ExtendedKalmanFilter(&prior_cont);

  prior_mu(1) = PRIOR_MU_X2;
  prior_mu(2) = PRIOR_MU_Y2;
  prior_cont.ExpectedValueSet(prior_mu);

  // second filter    
  my_filter_kalman2 = new ExtendedKalmanFilter(&prior_cont);

  prior_mu(1) = PRIOR_MU_X3;
  prior_mu(2) = PRIOR_MU_Y3;
  prior_cont.ExpectedValueSet(prior_mu);

  // third filter    
  my_filter_kalman3 = new ExtendedKalmanFilter(&prior_cont);

  vector<Filter<ColumnVector,ColumnVector>* > filters(NUM_FILTERS); 
  filters[0] = my_filter_kalman1;
  filters[1] = my_filter_kalman2;
  filters[2] = my_filter_kalman3;

  // constructor
  double gamma = 0.1;
  double treshold = 0.2;
  DataAssociationFilterGaussian my_dataAssociationFilter(filters,gamma,treshold);
  vector<Pdf<ColumnVector>* > filters_post = my_dataAssociationFilter.PostGet();
  vector<int > filterID = my_dataAssociationFilter.FilterIDGet();

  for(int i = 0; i< NUM_FILTERS ; i++)
  {
     CPPUNIT_ASSERT_EQUAL( filters[i]->PostGet() , filters_post[i]);
     CPPUNIT_ASSERT_EQUAL( i+1 , filterID[i]);
     cout << "filterID[" << i << "]" << filterID[i] << endl;
  }
    
  CPPUNIT_ASSERT_EQUAL( NUM_FILTERS , my_dataAssociationFilter.NumFiltersGet());
  CPPUNIT_ASSERT_EQUAL( gamma , my_dataAssociationFilter.GammaGet());
  CPPUNIT_ASSERT_EQUAL( treshold , my_dataAssociationFilter.TresholdGet());

  // REMOVE FILTER
  cout << "remove filter test started" << endl;
  vector<Filter<ColumnVector,ColumnVector>* > filters_bis(NUM_FILTERS); 
  filters_bis[0] = my_filter_kalman1;
  filters_bis[1] = my_filter_kalman3;

  CPPUNIT_ASSERT_EQUAL( false , my_dataAssociationFilter.RemoveFilter(-1));
  CPPUNIT_ASSERT_EQUAL( false , my_dataAssociationFilter.RemoveFilter(3));
  CPPUNIT_ASSERT_EQUAL( true , my_dataAssociationFilter.RemoveFilter(1));
  CPPUNIT_ASSERT_EQUAL( NUM_FILTERS-1 , my_dataAssociationFilter.NumFiltersGet());
  filters_post = my_dataAssociationFilter.PostGet();
  filterID = my_dataAssociationFilter.FilterIDGet();
  
  // check if other filters are still there
  CPPUNIT_ASSERT_EQUAL( filters_bis[0]->PostGet() , filters_post[0]);
  CPPUNIT_ASSERT_EQUAL( filters_bis[1]->PostGet() , filters_post[1]);
  CPPUNIT_ASSERT_EQUAL( 1 , filterID[0]);
  CPPUNIT_ASSERT_EQUAL( 3 , filterID[1]);

  cout << "remove filter test ended" << endl;
  
  // ADD FILTER
  cout << "add filter test started" << endl;
  prior_mu(1) = PRIOR_MU_X4;
  prior_mu(2) = PRIOR_MU_Y4;
  prior_cont.ExpectedValueSet(prior_mu);
  my_filter_kalman4 = new ExtendedKalmanFilter(&prior_cont);
  vector<Filter<ColumnVector,ColumnVector>* > filters_bis2(NUM_FILTERS); 
  filters_bis2[0] = my_filter_kalman1;
  filters_bis2[1] = my_filter_kalman3;
  filters_bis2[2] = my_filter_kalman4;

  my_dataAssociationFilter.AddFilter(my_filter_kalman4);
  CPPUNIT_ASSERT_EQUAL( NUM_FILTERS , my_dataAssociationFilter.NumFiltersGet());

  filters_post = my_dataAssociationFilter.PostGet();
  filterID = my_dataAssociationFilter.FilterIDGet();
  
  // check if other filters are still there
  CPPUNIT_ASSERT_EQUAL( filters_bis2[0]->PostGet() , filters_post[0]);
  CPPUNIT_ASSERT_EQUAL( filters_bis2[1]->PostGet() , filters_post[1]);
  CPPUNIT_ASSERT_EQUAL( filters_bis2[2]->PostGet() , filters_post[2]);
  CPPUNIT_ASSERT_EQUAL( 1 , filterID[0]);
  CPPUNIT_ASSERT_EQUAL( 3 , filterID[1]);
  CPPUNIT_ASSERT_EQUAL( 4 , filterID[2]);

  for (int teller_filters = 0;  teller_filters<my_dataAssociationFilter.NumFiltersGet() ; teller_filters++)
  {
    cout << "filters_post["<< teller_filters << "] expected value " << filters_post[teller_filters]->ExpectedValueGet() << endl;
    //cout << "filters_post[0] covariance " << ( (Matrix)((MCPdf<ColumnVector>*)(filters_post[teller_filters]))->CovarianceGet()) << endl;
    cout << "filterID[" << teller_filters << "]" << filterID[teller_filters] << endl;
  }

  prior_mu(1) = PRIOR_MU_X5;
  prior_mu(2) = PRIOR_MU_Y5;
  prior_cont.ExpectedValueSet(prior_mu);
  my_filter_kalman5 = new ExtendedKalmanFilter(&prior_cont);
  my_dataAssociationFilter.AddFilter(my_filter_kalman5);

  prior_mu(1) = PRIOR_MU_X6;
  prior_mu(2) = PRIOR_MU_Y6;
  prior_cont.ExpectedValueSet(prior_mu);
  my_filter_kalman6 = new ExtendedKalmanFilter(&prior_cont);
  my_dataAssociationFilter.AddFilter(my_filter_kalman6);

  cout << "make probs " << endl;
  Matrix probs(5,5);
  probs = 0.0;
  probs(1,1) = 0.4;
  //probs(1,2) = 0.4;
  //probs(1,3) = 0.4;
  //probs(2,1) = 0.4;
  probs(2,2) = 0.4;
  //probs(2,3) = 0.4;
  //probs(2,1) = 0.4;
  //probs(2,2) = 0.4;
  //probs(2,3) = 0.4;
  //probs(3,1) = 0.4;
  //probs(3,2) = 0.4;
  probs(3,3) = 0.4;
  probs(4,4) = 0.4;
  probs(5,5) = 0.4;
/*
  vector<vector<int > > associations = my_dataAssociationFilter.GetAssociations(probs,3);

  std::cout << "ASSOCIATIONS  " << std::endl;
  for (int i =0 ; i < 35 ; i++)
  {
       std::cout << "Association " << i+1 << std::endl;
       for (int j=0 ; j < 3 ; j++)
            std::cout << associations[i][j] << std::endl;
  } 

*/

 // create measurement model 
 Matrix H(2,2);
 H = 0.0;
 H(1,1) = 1.0;
 H(2,2) = 1.0;

 ColumnVector measNoise_Mu(2);
 measNoise_Mu(1) = 0.0;
 measNoise_Mu(2) = 0.0;

 SymmetricMatrix measNoise_Cov(2);
 measNoise_Cov(1,1) = pow(0.1,2);
 measNoise_Cov(2,2) = pow(0.1,2);

 Gaussian measurement_Uncertainty(measNoise_Mu, measNoise_Cov);

 LinearAnalyticConditionalGaussian meas_pdf(H,measurement_Uncertainty);

 LinearAnalyticMeasurementModelGaussianUncertainty meas_model(&meas_pdf);

 vector<ColumnVector> measurements(5);
 ColumnVector measurement(2);
 measurement(1) = 0.9;
 measurement(2) = 0.9;
 measurements[0] = measurement; 
 measurement(1) = 3.2;
 measurement(2) = 3.2;
 measurements[1] = measurement; 
 measurement(1) = 6.2;
 measurement(2) = 6.2;
 measurements[2] = measurement; 
 measurement(1) = 5.1;
 measurement(2) = 5.1;
 measurements[3] = measurement; 
 measurement(1) = 2.1;
 measurement(2) = 2.1;
 measurements[4] = measurement; 
 ColumnVector s;

 std::cout << "call GetMeasProbs " << endl;
 Matrix meas_probs = my_dataAssociationFilter.GetMeasProbs(&meas_model,measurements,s);
 std::cout << "MeasProbs " << meas_probs << endl;
 std::cout << "call getassociationprobs " << endl;
 vector<vector<double > > association_probs = my_dataAssociationFilter.GetAssociationProbs(&meas_model,measurements, s);
 std::cout << "ASSOCIATION_PROBS  " << std::endl;
 for (int i =0 ; i < my_dataAssociationFilter.NumFiltersGet() ; i++)
 {
      for (int j=0 ; j < measurements.size() ; j++)
      {
           std::cout << "Feature " << i+1 << " with  Object " << j+1 << std::endl;
           std::cout << association_probs[i][j] << std::endl;
      }
 } 

 Matrix A(2,2);
 A = 0.0;
 A(1,1) = 1.0;
 A(2,2) = 1.0;

 vector<Matrix> AB(1);
 AB[0] = A;

 ColumnVector sysNoise_Mu(2);
 sysNoise_Mu(1) = 0.0;
 sysNoise_Mu(2) = 0.0;

 SymmetricMatrix sysNoise_Cov(2);
 sysNoise_Cov = 0.0;
 sysNoise_Cov(1,1) = pow(0.2,2);
 sysNoise_Cov(2,2) = pow(0.2,2);

 Gaussian system_Uncertainty(sysNoise_Mu, sysNoise_Cov);
 LinearAnalyticConditionalGaussian sys_pdf(AB, system_Uncertainty);
 LinearAnalyticSystemModelGaussianUncertainty sys_model(&sys_pdf);

  // test system update
  cout<< "call system update!" << endl;
  my_dataAssociationFilter.Update(&sys_model);
  filters_post = my_dataAssociationFilter.PostGet();
  cout<< "SYSTEM UPDATE DONE!" << endl;
  for (int teller_filters = 0;  teller_filters<my_dataAssociationFilter.NumFiltersGet() ; teller_filters++)
  {
    cout << "filters_post["<< teller_filters << "] expected value " << filters_post[teller_filters]->ExpectedValueGet() << endl;
    //cout << "filters_post[0] covariance " << ( (Matrix)((MCPdf<ColumnVector>*)(filters_post[teller_filters]))->CovarianceGet()) << endl;
  }

  // test measurement update
  my_dataAssociationFilter.Update(&meas_model,measurements);
  my_dataAssociationFilter.Update(&meas_model,measurements);
  filters_post = my_dataAssociationFilter.PostGet();
  cout<< "MEASUREMENT UPDATE DONE!" << endl;
  for (int teller_filters = 0;  teller_filters<my_dataAssociationFilter.NumFiltersGet() ; teller_filters++)
  {
    cout << "filters_post["<< teller_filters << "] expected value " << filters_post[teller_filters]->ExpectedValueGet() << endl;
    cout << "filters_post["<< teller_filters << "] covariance " << (filters_post[teller_filters]->CovarianceGet()) << endl;
  }

  //association_probs = my_dataAssociationFilter.GetAssociationProbs(&meas_model,measurements, s);
  //my_dataAssociationFilter.Update(&meas_model,measurements);
}
void 
DataAssociationFilterTest::tearDown()
{
}
