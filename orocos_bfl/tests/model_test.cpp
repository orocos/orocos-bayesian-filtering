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
 
#include "model_test.hpp"
#include <cmath> // For sinus

#define NUM_DS 5
#define NUM_COND_ARGS 1 

//testLinearAnalyticSystemModelGaussianUncertainty
#define STATE_SIZE 3
#define INPUT_SIZE 2

#define SIGMA_NOISE 0.0001 // Noise variance (constant for every input here)
#define CORRELATION_NOISE 0.0 // Correlation between different noise (idem)

#define LINEAR_SPEED 1.0
#define ROT_SPEED 0.0

//testLinearAnalyticMeasurementModelGaussianUncertainty
#define MEASUREMENT_SIZE 1
// Coordinates of wall
#define RICO_WALL 0.5
#define OFFSET_WALL 30
#define DELTA_T 1

// Measurement noise
#define WALL_CT 1/(sqrt(pow(RICO_WALL,2) + 1))
#define MU_MEAS_NOISE OFFSET_WALL*WALL_CT
#define SIGMA_MEAS_NOISE 0.5
// Registers the fixture into the 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION( ModelTest );

using namespace BFL;

void 
ModelTest::setUp()
{
}

void 
ModelTest::tearDown()
{
}

void 
ModelTest::testDiscreteSystemModel()
{
  int cond_arg_dims[NUM_COND_ARGS] = { NUM_DS };
  int state_k;int state_kMinusOne;
  
  DiscreteConditionalPdf a_discretecondpdf(NUM_DS,NUM_COND_ARGS,cond_arg_dims);
  std::vector<int> cond_args(NUM_COND_ARGS);
  /* Set and Get all Probabilities*/
  double prob_diag = 0.98;
  double prob_nondiag = (1-prob_diag)/(NUM_DS-1);
  for (state_kMinusOne = 0 ; state_kMinusOne < NUM_DS ;  state_kMinusOne++)
    {
       cond_args[0] = state_kMinusOne;
       for (state_k = 0 ; state_k < NUM_DS ;  state_k++)
         {
           if (state_kMinusOne == state_k) a_discretecondpdf.ProbabilitySet(prob_diag,state_k,cond_args);
           else a_discretecondpdf.ProbabilitySet(prob_nondiag,state_k,cond_args);
         }
    }

  /* Construction */
  DiscreteSystemModel a_discreteSysModel(&a_discretecondpdf);
 
  /* Number of discrete states*/
  CPPUNIT_ASSERT_EQUAL( NUM_DS , (int)a_discreteSysModel.NumStatesGet());

  /* State size*/
  CPPUNIT_ASSERT_EQUAL( 1 , a_discreteSysModel.StateSizeGet());
  
  /* SystemWithoutInputs*/ 
  CPPUNIT_ASSERT_EQUAL( true , a_discreteSysModel.SystemWithoutInputs());
 
  /* ProbabilityGet without inputs*/
  for (state_kMinusOne = 0 ; state_kMinusOne < NUM_DS ;  state_kMinusOne++)
    {
       cond_args[0] = state_kMinusOne;
       for (state_k = 0 ; state_k < NUM_DS ;  state_k++)
         {
            CPPUNIT_ASSERT_EQUAL( (double)a_discretecondpdf.ProbabilityGet(state_k),(double)a_discreteSysModel.ProbabilityGet(state_k,state_kMinusOne));
         }
    }
 
  /* Simulate*/// ???
  //for (state_kMinusOne = 0 ; state_kMinusOne < NUM_DS ;  state_kMinusOne++)
  //   {
  //      a_discretecondpdf.ConditionalArgumentSet(0,state_kMinusOne);
  //      CPPUNIT_ASSERT_EQUAL( a_discretecondpdf.ExpectedValueGet() , a_discreteSysModel.Simulate(state_kMinusOne));
  //   }
  
  /* Copy constructor */
  DiscreteSystemModel b_discreteSysModel(a_discreteSysModel);
  for (state_kMinusOne = 0 ; state_kMinusOne < NUM_DS ;  state_kMinusOne++)
    {
       cond_args[0] = state_kMinusOne;
       for (state_k = 0 ; state_k < NUM_DS ;  state_k++)
         {
            CPPUNIT_ASSERT_EQUAL( (double)a_discreteSysModel.ProbabilityGet(state_k,state_kMinusOne),(double)b_discreteSysModel.ProbabilityGet(state_k,state_kMinusOne));
         }
    }

  // creating system model with input
  int num_cond_args_new = 2;
  int size_input = 2;
  int cond_args_dims_new[num_cond_args_new];
  cond_args_dims_new[0] =  NUM_DS; cond_args_dims_new[1] =  size_input;
  int input;
  
  DiscreteConditionalPdf c_discretecondpdf(NUM_DS,num_cond_args_new,cond_args_dims_new);
  std::vector<int> cond_args_new(num_cond_args_new);
  /* Set and Get all Probabilities*/
  for (state_kMinusOne = 0 ; state_kMinusOne < NUM_DS ;  state_kMinusOne++)
    {
       cond_args_new[0] = state_kMinusOne;
       for (input = 0 ; input < size_input ;  input++)
         {
            cond_args_new[1] = input;
            for (state_k = 0 ; state_k < NUM_DS ;  state_k++)
              {
                if (state_kMinusOne == state_k) c_discretecondpdf.ProbabilitySet(prob_diag,state_k,cond_args_new);
                else c_discretecondpdf.ProbabilitySet(prob_nondiag,state_k,cond_args_new);
              }
         }
    }

  DiscreteSystemModel c_discreteSysModel;
  c_discreteSysModel.SystemPdfSet(&c_discretecondpdf);

  /* SystemWithoutInputs*/ 
  CPPUNIT_ASSERT_EQUAL( false , c_discreteSysModel.SystemWithoutInputs());

  /* ProbabilityGet without inputs*/
  for (state_kMinusOne = 0 ; state_kMinusOne < NUM_DS ;  state_kMinusOne++)
    {
       cond_args[0] = state_kMinusOne;
       for (input = 0 ; input < size_input ;  input++)
         {
            cond_args_new[1] = input;
            for (state_k = 0 ; state_k < NUM_DS ;  state_k++)
              {
                 CPPUNIT_ASSERT_EQUAL( (double)c_discretecondpdf.ProbabilityGet(state_k),(double)c_discreteSysModel.ProbabilityGet(state_k,state_kMinusOne,input));
              }
         }
    }
}


void
ModelTest::testLinearAnalyticMeasurementModelGaussianUncertainty()
{
  Matrix A(STATE_SIZE,STATE_SIZE);
  Matrix B(STATE_SIZE,INPUT_SIZE);

  ColumnVector state(STATE_SIZE);
  ColumnVector initial_state(3);
  initial_state(1) = 1.0;
  initial_state(2) = 1.0;
  initial_state(3) = 0.5;
  ColumnVector input(INPUT_SIZE);
  ColumnVector initial_input(2);
  initial_input(1) = 1.0;
  initial_input(2) = 1.0;

  for (int i=1; i < STATE_SIZE+1; i++){state(i) = initial_state(i);}    
  for (int i=1; i < INPUT_SIZE+1; i++){input(i) = initial_input(i);}    

  // Uncertainty or Noice (Additive)
  ColumnVector Noise_Mu(STATE_SIZE);
  for (int row=0; row < STATE_SIZE; row++){Noise_Mu(row+1) = 0;}    
  SymmetricMatrix Noise_Cov(STATE_SIZE);
  for (int row=0; row < STATE_SIZE; row++)
    {
      for (int column=0; column < STATE_SIZE; column++)
	{
	  if (row == column) {Noise_Cov(row+1,column+1) = SIGMA_NOISE;}
	  else {Noise_Cov(row+1,column+1) = CORRELATION_NOISE;}
	}
    }
  Gaussian System_Uncertainty(Noise_Mu,Noise_Cov);
 
  // MATRIX A: unit matrix 
  for (int row=0; row < STATE_SIZE; row++)
    {
      for (int column=0; column < STATE_SIZE; column++)
	{
	  if (row == column) A(row+1,column+1)=1;
	  else A(row+1,column+1)=0;
	}
    }
  // MATRIX B
  B(STATE_SIZE,1) = 0.0; B(1,INPUT_SIZE) = 0.0; B(2,INPUT_SIZE) = 0.0;
  B(STATE_SIZE,INPUT_SIZE) = 1.0;
  B(1,1) = cos(state(STATE_SIZE)) * DELTA_T;
  B(2,1) = sin(state(STATE_SIZE)) * DELTA_T;

  vector<Matrix> v(2);
  v[0] = A;
  v[1] = B;
  LinearAnalyticConditionalGaussian pdf(v,System_Uncertainty);
  LinearAnalyticSystemModelGaussianUncertainty a_linAnSysModel(&pdf);

  /* State size Get */
  CPPUNIT_ASSERT_EQUAL( STATE_SIZE , a_linAnSysModel.StateSizeGet());

  /* System without inputs */
  CPPUNIT_ASSERT_EQUAL( false , a_linAnSysModel.SystemWithoutInputs());

  /* Get the system model matrices */
  CPPUNIT_ASSERT_EQUAL( A , a_linAnSysModel.AGet());
  CPPUNIT_ASSERT_EQUAL( B , a_linAnSysModel.BGet());

  /* df_dxGet */
  CPPUNIT_ASSERT_EQUAL( A , a_linAnSysModel.df_dxGet(input,state));
  CPPUNIT_ASSERT_EQUAL( B , a_linAnSysModel.BGet());

  /* PredictionGet */
  pdf.ConditionalArgumentSet(0,state);
  pdf.ConditionalArgumentSet(1,input);
  CPPUNIT_ASSERT_EQUAL( pdf.ExpectedValueGet() , a_linAnSysModel.PredictionGet(input,state) );
  
  /* CovarianceGet */
  CPPUNIT_ASSERT_EQUAL( pdf.CovarianceGet() , a_linAnSysModel.CovarianceGet(input,state) );
 
  /* Simulate */
  ColumnVector state_new(STATE_SIZE);
  state_new = a_linAnSysModel.Simulate(state, input);
  // TODO: can we check this in any way?
  
  /* ProbabilityGet */
  CPPUNIT_ASSERT_EQUAL((double)pdf.ProbabilityGet(state_new) , (double)a_linAnSysModel.ProbabilityGet(state_new,state,input) );

  // build new A and B matrices for new state
  // MATRIX B
  B(1,1) = cos(state_new(STATE_SIZE)) * DELTA_T;
  B(2,1) = sin(state_new(STATE_SIZE)) * DELTA_T;

  /* A and B set and get */
  a_linAnSysModel.ASet(A);
  a_linAnSysModel.BSet(B);
  CPPUNIT_ASSERT_EQUAL(A , a_linAnSysModel.AGet());
  CPPUNIT_ASSERT_EQUAL(B , a_linAnSysModel.BGet());
  
  //TODO: SystemPdfSet
}

void
ModelTest::testLinearAnalyticSystemModelGaussianUncertainty()
{
  Matrix H(MEASUREMENT_SIZE,STATE_SIZE);
  
  ColumnVector state(STATE_SIZE);
  ColumnVector initial_state(3);
  initial_state(1) = 1.0;
  initial_state(2) = 1.0;
  initial_state(3) = 0.5;
  for (int i=1; i < STATE_SIZE+1; i++){state(i) = initial_state(i);}    
  ColumnVector measurement(MEASUREMENT_SIZE);
  ColumnVector MeasNoise_Mu(MEASUREMENT_SIZE);
  SymmetricMatrix MeasNoise_Cov(MEASUREMENT_SIZE);

  // Fill up H
  H(1,1) = WALL_CT * RICO_WALL;
  H(1,2) = 0 - WALL_CT;
  H(1,STATE_SIZE) = 0;
  
  // Construct the measurement noise (a scalar in this case)
  MeasNoise_Mu(1) = MU_MEAS_NOISE;
  MeasNoise_Cov(1,1) = SIGMA_MEAS_NOISE;

  Gaussian Measurement_Uncertainty(MeasNoise_Mu,MeasNoise_Cov);
  LinearAnalyticConditionalGaussian pdf(H,Measurement_Uncertainty);
  LinearAnalyticMeasurementModelGaussianUncertainty a_linAnMeasModel(&pdf);

  /* Measurement size Get */
  CPPUNIT_ASSERT_EQUAL( MEASUREMENT_SIZE , a_linAnMeasModel.MeasurementSizeGet());

  /* System without sensor parameters */
  CPPUNIT_ASSERT_EQUAL( true , a_linAnMeasModel.SystemWithoutSensorParams());

  /* Get the measurement model matrices */
  CPPUNIT_ASSERT_EQUAL( H , a_linAnMeasModel.HGet());
  //CPPUNIT_ASSERT_EQUAL( J , a_linAnMeasModel.JGet());

  /* df_dxGet */
  CPPUNIT_ASSERT_EQUAL( H , a_linAnMeasModel.df_dxGet(0,state));

  /* PredictionGet */
  pdf.ConditionalArgumentSet(0,state);
  CPPUNIT_ASSERT_EQUAL( pdf.ExpectedValueGet() , a_linAnMeasModel.PredictionGet(0,state) );
  
  /* CovarianceGet */
  CPPUNIT_ASSERT_EQUAL( pdf.CovarianceGet() , a_linAnMeasModel.CovarianceGet(0,state) );
 
  /* Simulate */
  ColumnVector meas_new(MEASUREMENT_SIZE);
  meas_new = a_linAnMeasModel.Simulate(state);
  // TODO: can we check this in any way?
  
  /* ProbabilityGet */
  CPPUNIT_ASSERT_EQUAL((double)pdf.ProbabilityGet(meas_new) , (double)a_linAnMeasModel.ProbabilityGet(meas_new,state) );

  // build new H 
  H = 2;

  /* H  set and get */
  a_linAnMeasModel.HSet(H);
  CPPUNIT_ASSERT_EQUAL(H , a_linAnMeasModel.HGet());
  
  //TODO: SystemPdfSet


}
