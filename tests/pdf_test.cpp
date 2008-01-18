// Copyright (C) 2007 Klaas Gadeyne <first dot last at gmail dot com>
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
 
#include "pdf_test.hpp"
#include "approxEqual.hpp"

// Registers the fixture into the 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION( PdfTest );
using namespace BFL;

#define DIMENSION_INPUT 2
#define MU_0 2.8
#define MU_1 3.4
#define MU_2 3.4
#define MU_3 3.4
#define SIGMA 0.01
#define SIGMA2 0.00
#define WIDTH_0 1
#define WIDTH_1 2
const unsigned int NUM_COND_ARGS = 2;
const unsigned int DIMENSION = 2;
const unsigned int NUM_SAMPLES = 500;

//discrete pdf
#define NUM_DS 5
//discrete conditional pdf
#define NUM_COND_ARGS_DIS 1

void 
PdfTest::setUp()
{
  epsilon = 0.0001;

  _mu.resize(DIMENSION);
  _mu(1) = MU_0; _mu(2) = MU_1;
  _sigma.resize(DIMENSION);
  _sigma = 0.0;
  for (unsigned int rows=1; rows < DIMENSION + 1; rows++){ _sigma(rows,rows)=SIGMA; }
  _width.resize(DIMENSION);
  _width(1) = WIDTH_0 ;
  _width(2) = WIDTH_1 ;
}

void 
PdfTest::tearDown()
{
}

void 
PdfTest::testGaussian()
{
  Gaussian a_gaussian(_mu,_sigma);

  /* Sampling */
  // one sample
  Sample<ColumnVector> a_sample ;
    CPPUNIT_ASSERT_EQUAL( true, a_gaussian.SampleFrom(a_sample,DEFAULT,NULL));
  // 4 sigma test : REMARK: this test WILL occasionaly fail with probability
  // erf(4/sqrt(2)) for each sample
  for(int j = 1 ; j <= DIMENSION; j++)  
  {
    CPPUNIT_ASSERT( (a_sample.ValueGet())(j) > _mu(j) - 4.0 * sqrt(_sigma(j,j) ) );
    CPPUNIT_ASSERT( (a_sample.ValueGet())(j) < _mu(j) + 4.0 * sqrt(_sigma(j,j) ) );
  }
  // Box-Muller is not implemented yet
  CPPUNIT_ASSERT_EQUAL( false, a_gaussian.SampleFrom(a_sample,BOXMULLER,NULL));
  CPPUNIT_ASSERT_EQUAL( true, a_gaussian.SampleFrom(a_sample,CHOLESKY,NULL));
  // 4 sigma test : REMARK: this test WILL occasionaly fail with probability
  // erf(4/sqrt(2)) for each sample
  for(int j = 1 ; j <= DIMENSION; j++)  
  {
    CPPUNIT_ASSERT( (a_sample.ValueGet())(j) > _mu(j) - 4.0 * sqrt(_sigma(j,j) ) );
    CPPUNIT_ASSERT( (a_sample.ValueGet())(j) < _mu(j) + 4.0 * sqrt(_sigma(j,j) ) );
  }

  // list of samples
  vector<Sample<ColumnVector> > los(NUM_SAMPLES);
  CPPUNIT_ASSERT_EQUAL( true, a_gaussian.SampleFrom(los,NUM_SAMPLES,DEFAULT,NULL));
  // 4 sigma test : REMARK: this test WILL occasionaly fail with probability
  // erf(4/sqrt(2)) for each sample
  ColumnVector sample_mean(DIMENSION);
  sample_mean = 0.0;
  for(int i = 0 ; i < NUM_SAMPLES; i++)  
  {
    a_sample = los[i];
    sample_mean += a_sample.ValueGet();
    for(int j = 1 ; j <= DIMENSION; j++)  
    {
        CPPUNIT_ASSERT( (a_sample.ValueGet())(j) > _mu(j) - 4.0 * sqrt(_sigma(j,j) ) );
        CPPUNIT_ASSERT( (a_sample.ValueGet())(j) < _mu(j) + 4.0 * sqrt(_sigma(j,j) ) );
    }
  }
  // check whether mean and of samples is in withing the expected 3 sigma border
  // 3 sigma test : REMARK: this test WILL occasionaly fail with probability
  // erf(3/sqrt(2)) for each sample
  sample_mean = sample_mean / (double)NUM_SAMPLES;
  for(int j = 1 ; j <= DIMENSION; j++)  
  {
    CPPUNIT_ASSERT( sample_mean(j) <= _mu(j) + 3.0 * 1/sqrt(NUM_SAMPLES) * sqrt(_sigma(j,j)) );
    CPPUNIT_ASSERT( sample_mean(j) >= _mu(j) - 3.0 * 1/sqrt(NUM_SAMPLES) * sqrt(_sigma(j,j)) );
  }
  // Box-Muller is not implemented yet
  CPPUNIT_ASSERT_EQUAL( false, a_gaussian.SampleFrom(los,NUM_SAMPLES,BOXMULLER,NULL));
  CPPUNIT_ASSERT_EQUAL( true, a_gaussian.SampleFrom(los,NUM_SAMPLES,CHOLESKY,NULL));
  // 4 sigma test : REMARK: this test WILL occasionaly fail with probability
  // erf(4/sqrt(2)) for each sample
  sample_mean = 0.0;
  for(int i = 0 ; i < NUM_SAMPLES; i++)  
  {
    a_sample = los[i];
    sample_mean += a_sample.ValueGet();
    for(int j = 1 ; j <= DIMENSION; j++)  
    {
        CPPUNIT_ASSERT( (a_sample.ValueGet())(j) > _mu(j) - 4.0 * sqrt(_sigma(j,j) ) );
        CPPUNIT_ASSERT( (a_sample.ValueGet())(j) < _mu(j) + 4.0 * sqrt(_sigma(j,j) ) );
    }
  }
  // check whether mean and of samples is in withing the expected 3 sigma border
  // 3 sigma test : REMARK: this test WILL occasionaly fail with probability
  // erf(3/sqrt(2)) for each sample
  sample_mean = sample_mean / (double)NUM_SAMPLES;
  for(int j = 1 ; j <= DIMENSION; j++)  
  {
    CPPUNIT_ASSERT( sample_mean(j) <= _mu(j) + 3.0 * 1/sqrt(NUM_SAMPLES) * sqrt(_sigma(j,j)) );
    CPPUNIT_ASSERT( sample_mean(j) >= _mu(j) - 3.0 * 1/sqrt(NUM_SAMPLES) * sqrt(_sigma(j,j)) );
  }
  
  /* Setting and Getting mean/covariance */
  CPPUNIT_ASSERT_EQUAL( _mu, a_gaussian.ExpectedValueGet());
  CPPUNIT_ASSERT_EQUAL( _sigma, a_gaussian.CovarianceGet());

  /* Copy Constructor etc */
  Gaussian b_gaussian(a_gaussian);
  CPPUNIT_ASSERT_EQUAL( _mu, b_gaussian.ExpectedValueGet());
  CPPUNIT_ASSERT_EQUAL( _sigma, b_gaussian.CovarianceGet());

  /* Create Gaussian, allocate memory afterwards */
  Gaussian c_gaussian;
  c_gaussian.ExpectedValueSet(_mu);
  CPPUNIT_ASSERT_EQUAL( _mu, c_gaussian.ExpectedValueGet());
  c_gaussian.CovarianceSet(_sigma);
  CPPUNIT_ASSERT_EQUAL( _sigma, c_gaussian.CovarianceGet());
}  

void 
PdfTest::testUniform()
{
  Uniform a_uniform(_mu,_width);

  /* Setting and Getting center and widths */
  CPPUNIT_ASSERT_EQUAL(approxEqual( _mu, a_uniform.CenterGet(),epsilon),true);
  CPPUNIT_ASSERT_EQUAL(approxEqual( _width, a_uniform.WidthGet(),epsilon),true);

  /* Sampling */
  vector<Sample<ColumnVector> > los(NUM_SAMPLES);
  CPPUNIT_ASSERT_EQUAL( true, a_uniform.SampleFrom(los,NUM_SAMPLES,DEFAULT,NULL));
  // test if samples are located in the area
  vector<Sample<ColumnVector > >::iterator los_it;
  for (los_it = los.begin(); los_it!=los.end(); los_it++)
  { 
       ColumnVector current_sample = los_it->ValueGet();
       for (int i = 1; i < _mu.rows()+1 ; i++)
       {
            CPPUNIT_ASSERT(current_sample(i) > (_mu(i)-_width(i)/2) ) ;
            CPPUNIT_ASSERT(current_sample(i) < (_mu(i)+_width(i)/2) ) ;
       }
   }
  // Box-Muller is not implemented yet
  CPPUNIT_ASSERT_EQUAL( false, a_uniform.SampleFrom(los,NUM_SAMPLES,BOXMULLER,NULL));
  Sample<ColumnVector> a_sample ;
  CPPUNIT_ASSERT_EQUAL( true, a_uniform.SampleFrom(a_sample,DEFAULT,NULL));
  // Box-Muller is not implemented yet
  CPPUNIT_ASSERT_EQUAL( false, a_uniform.SampleFrom(a_sample,BOXMULLER,NULL));

  /* Getting the probability */
  double area = 1;
  for (int i =1 ; i < DIMENSION + 1 ; i++) area = area * _width(i); 
  CPPUNIT_ASSERT_EQUAL(1/area, (double)a_uniform.ProbabilityGet(_mu));
  CPPUNIT_ASSERT_EQUAL(0.0, (double)a_uniform.ProbabilityGet(_mu + _width));
  CPPUNIT_ASSERT_EQUAL(0.0, (double)a_uniform.ProbabilityGet(_mu - _width));
  ColumnVector test_prob(DIMENSION);
  test_prob = _mu;
  test_prob(DIMENSION) = _mu(DIMENSION) + _width(DIMENSION)*2/3;
  CPPUNIT_ASSERT_EQUAL(0.0, (double)a_uniform.ProbabilityGet(test_prob));
  
  /* Copy Constructor etc */
  Uniform b_uniform(a_uniform);
  CPPUNIT_ASSERT_EQUAL(approxEqual( _mu, b_uniform.CenterGet(),epsilon),true);
  CPPUNIT_ASSERT_EQUAL(approxEqual( _width, b_uniform.WidthGet(),epsilon),true);

  /* Create Uniform, allocate memory afterwards */
  Uniform c_uniform;
  c_uniform.UniformSet(_mu,_width);
  CPPUNIT_ASSERT_EQUAL(approxEqual( _mu, c_uniform.CenterGet(),epsilon),true);
  CPPUNIT_ASSERT_EQUAL(approxEqual( _width, c_uniform.WidthGet(),epsilon),true);
  /* Getting the probability */
  area = 1;
  for (int i =1 ; i < DIMENSION + 1 ; i++) area = area * _width(i); 
  CPPUNIT_ASSERT_EQUAL(1/area, (double)c_uniform.ProbabilityGet(_mu));
  CPPUNIT_ASSERT_EQUAL(0.0, (double)c_uniform.ProbabilityGet(_mu + _width));
  CPPUNIT_ASSERT_EQUAL(0.0, (double)c_uniform.ProbabilityGet(_mu - _width));
  test_prob = _mu;
  test_prob(DIMENSION) = _mu(DIMENSION) + _width(DIMENSION)*2/3;
  CPPUNIT_ASSERT_EQUAL(0.0, (double)c_uniform.ProbabilityGet(test_prob));
}  
void
PdfTest::testDiscretePdf()
{
  DiscretePdf a_discretepdf(NUM_DS);
  vector<Probability> uniform_vector(NUM_DS);
  for (int state = 0; state < NUM_DS ; state++) uniform_vector[state] = (Probability)(1.0/NUM_DS );

  /* Check initial uniform distribution */
  vector<Probability> result_proba = a_discretepdf.ProbabilitiesGet();
  for (int state = 0; state < NUM_DS ; state++) CPPUNIT_ASSERT_EQUAL( (double) (uniform_vector[state]), double(result_proba[state] )) ;

  /* Check dimension */
  CPPUNIT_ASSERT_EQUAL(1,(int)a_discretepdf.DimensionGet());

  /* Check the number of states */
  CPPUNIT_ASSERT_EQUAL(NUM_DS,(int)a_discretepdf.NumStatesGet());

  /* Copy constructor */
  DiscretePdf b_discretepdf(a_discretepdf);
  vector<Probability> result_probb = b_discretepdf.ProbabilitiesGet();
  for (int state = 0; state < NUM_DS ; state++) CPPUNIT_ASSERT_EQUAL( (double) (result_proba[state]), double(result_probb[state] )) ;
  CPPUNIT_ASSERT_EQUAL(a_discretepdf.DimensionGet(), b_discretepdf.DimensionGet());
  CPPUNIT_ASSERT_EQUAL(a_discretepdf.NumStatesGet(), b_discretepdf.NumStatesGet());

  /* Set and Get probabilities */
  // set one probability 
  double prob_new = 0.57;
  int new_el = NUM_DS-1; 
  CPPUNIT_ASSERT_EQUAL(true, a_discretepdf.ProbabilitySet(new_el,prob_new));
  CPPUNIT_ASSERT_EQUAL(prob_new, (double) a_discretepdf.ProbabilityGet(new_el));
  // check if the sum of the probabilities is still one!
  double sumProb = 0.0;
  for (int state = 0; state < NUM_DS ; state++)
  { 
    sumProb = sumProb + a_discretepdf.ProbabilityGet(state);
  }
  CPPUNIT_ASSERT_EQUAL(1.0, sumProb);

  // set all probabilities one by one
  for (int new_el_c = 0; new_el_c < NUM_DS; new_el_c++)
  {
     CPPUNIT_ASSERT_EQUAL(true, a_discretepdf.ProbabilitySet(new_el_c,prob_new));
     CPPUNIT_ASSERT_EQUAL(prob_new, (double) a_discretepdf.ProbabilityGet(new_el_c));
  }
  // check if the sum of the probabilities is still one!
  sumProb = 0.0;
  for (int state = 0; state < NUM_DS ; state++)
  { 
    sumProb = sumProb + a_discretepdf.ProbabilityGet(state);
  }
  CPPUNIT_ASSERT_EQUAL(1.0, sumProb);

  // all at the same time
  vector<Probability> prob_vec(NUM_DS);
  prob_vec[0] = 0.9;
  for (int state = 1; state < NUM_DS  ; state++)
  {
     prob_vec[state] = (Probability)( (1-0.9)/(NUM_DS-1) );
  }
  a_discretepdf.ProbabilitiesSet(prob_vec);
  result_proba = a_discretepdf.ProbabilitiesGet();
  for (int state = 0; state < NUM_DS ; state++) CPPUNIT_ASSERT_EQUAL( (double) (prob_vec[state]), double(result_proba[state] )) ;

  // check if the sum of the probabilities is still one!
  sumProb = 0.0;
  for (int state = 0; state < NUM_DS ; state++)
  { 
    sumProb = sumProb + a_discretepdf.ProbabilityGet(state);
  }
  CPPUNIT_ASSERT_EQUAL(1.0, sumProb);

  // special case for setting probability : set a new value for a probability
  // which was == 1
  DiscretePdf c_discretepdf(NUM_DS);
  vector<Probability> prob_vecc(NUM_DS);
  prob_vecc[0] = 1.0; 
  for (int state = 1; state < NUM_DS  ; state++)
  {
    prob_vecc[state] = 0.0; 
  }
  c_discretepdf.ProbabilitiesSet(prob_vecc);
  Probability new_prob0 = 0.5;
  double new_prob_other = (1-new_prob0)/(NUM_DS-1);
  c_discretepdf.ProbabilitySet(0,new_prob0);
  CPPUNIT_ASSERT_EQUAL((double)new_prob0, (double)(c_discretepdf.ProbabilityGet(0)));
  for (int state = 1; state < NUM_DS  ; state++)
  {
    CPPUNIT_ASSERT_EQUAL( new_prob_other, (double)(c_discretepdf.ProbabilityGet(state)));
  }
  
  /* Sampling */
  DiscretePdf d_discretepdf(NUM_DS);
  vector<Probability> prob_vecd(NUM_DS);
  prob_vecd[0] = 0.1; 
  prob_vecd[1] = 0.2; 
  prob_vecd[2] = 0.3; 
  prob_vecd[3] = 0.05; 
  prob_vecd[4] = 0.35; 
  d_discretepdf.ProbabilitiesSet(prob_vecd);
  vector<Sample<int> > los(NUM_SAMPLES);
  vector<Sample<int> >::iterator it;
  CPPUNIT_ASSERT_EQUAL( true, d_discretepdf.SampleFrom(los,NUM_SAMPLES,DEFAULT,NULL));
  // test if samples are distributed according to probabilities
  // remark that this test occasionally fails  
  vector<unsigned int> num_samples(NUM_DS);
  for (it = los.begin(); it!=los.end();it++)
  {
    num_samples[it->ValueGet()] +=1;
  }
  for (int i = 0 ; i< NUM_DS ; i++)
  { 
    //TODO: find some theoretical limits depending on the number of samples
    CPPUNIT_ASSERT( approxEqual( prob_vecd[i] , (double)num_samples[i] / NUM_SAMPLES, 0.05 ) );
  }
  
  // Ripley
  CPPUNIT_ASSERT_EQUAL( true, d_discretepdf.SampleFrom(los,NUM_SAMPLES,RIPLEY,NULL));
  // test if samples are distributed according to probabilities
  // remark that this test occasionally fails  
  vector<unsigned int> num_samples2(NUM_DS);
  for (it = los.begin(); it!=los.end();it++)
  {
    num_samples2[it->ValueGet()] +=1;
  }
  for (int i = 0 ; i< NUM_DS ; i++)
  { 
    //TODO: find some theoretical limits depending on the number of samples
   CPPUNIT_ASSERT( approxEqual( prob_vecd[i] , (double)num_samples2[i] / NUM_SAMPLES, 0.05 ) );
  }

  Sample<int> a_sample ;
  CPPUNIT_ASSERT_EQUAL( true, d_discretepdf.SampleFrom(a_sample,DEFAULT,NULL));
  // Ripley not implemented for one sample
  CPPUNIT_ASSERT_EQUAL( false, d_discretepdf.SampleFrom(a_sample,RIPLEY,NULL));
}

void
PdfTest::testLinearAnalyticConditionalGaussian()
{
  //creating additive gaussian noise
  Gaussian noise(_mu,_sigma);

  // Conditional gaussian with NUM_CONDITIONAL_ARGS conditional args and state vector of
  // dimension DIMENSION
  Matrix a(DIMENSION,DIMENSION); a = 1.0;
  Matrix b(DIMENSION,DIMENSION); b = 1.0;
  vector<Matrix> v(2); v[0] = a; v[1] = b;
  
  LinearAnalyticConditionalGaussian a_condgaussian(v, noise);

  /* Dimension check*/
  CPPUNIT_ASSERT_EQUAL( NUM_COND_ARGS, a_condgaussian.NumConditionalArgumentsGet());

  /* Matrix Check */
  for (unsigned int i=0; i < NUM_COND_ARGS; i++)
  {
    CPPUNIT_ASSERT_EQUAL( v[i], a_condgaussian.MatrixGet(i));
  }

  /* Getting mean/covariance  of the additive noise*/
  CPPUNIT_ASSERT_EQUAL( _mu, a_condgaussian.AdditiveNoiseMuGet());
  CPPUNIT_ASSERT_EQUAL( _sigma, a_condgaussian.AdditiveNoiseSigmaGet());

  /* Setting and getting the conditional args; one at a time */
  ColumnVector x2(DIMENSION); x2(1) = 1.0; x2(DIMENSION) = 0.5;
  ColumnVector u2(DIMENSION_INPUT); u2(1) = 1.0; u2(DIMENSION_INPUT) = 0.5;

  a_condgaussian.ConditionalArgumentSet(0,x2);
  a_condgaussian.ConditionalArgumentSet(1,u2);

  CPPUNIT_ASSERT_EQUAL( DIMENSION, a_condgaussian.DimensionGet());
  CPPUNIT_ASSERT_EQUAL( x2, a_condgaussian.ConditionalArgumentGet(0));
  CPPUNIT_ASSERT_EQUAL( u2, a_condgaussian.ConditionalArgumentGet(1));
  

  /* Setting and getting the conditional args; all together */
  ColumnVector x(DIMENSION); x(1) = 1.0; x(DIMENSION) = 0.5;
  ColumnVector u(DIMENSION_INPUT); u(1) = 1.0; u(DIMENSION_INPUT) = 0.5;

  std::vector<ColumnVector> cond_args(NUM_COND_ARGS);
  cond_args[0] = x; 
  cond_args[NUM_COND_ARGS-1] = u;

  a_condgaussian.ConditionalArgumentsSet(cond_args);

  CPPUNIT_ASSERT_EQUAL( DIMENSION, a_condgaussian.DimensionGet());
  for (unsigned int i=0; i < NUM_COND_ARGS; i++)
  {
    CPPUNIT_ASSERT_EQUAL( cond_args[i], a_condgaussian.ConditionalArgumentsGet()[i]);
  }

  /* Sampling */
  vector<Sample<ColumnVector> > los(NUM_SAMPLES);
  CPPUNIT_ASSERT_EQUAL( true, a_condgaussian.SampleFrom(los,NUM_SAMPLES,DEFAULT,NULL));
  // Box-Muller is not implemented yet
  CPPUNIT_ASSERT_EQUAL( false, a_condgaussian.SampleFrom(los,NUM_SAMPLES,BOXMULLER,NULL));
  CPPUNIT_ASSERT_EQUAL( true, a_condgaussian.SampleFrom(los,NUM_SAMPLES,CHOLESKY,NULL));
  Sample<ColumnVector> a_sample ;
  CPPUNIT_ASSERT_EQUAL( true, a_condgaussian.SampleFrom(a_sample,DEFAULT,NULL));
  // Box-Muller is not implemented yet
  CPPUNIT_ASSERT_EQUAL( false, a_condgaussian.SampleFrom(a_sample,BOXMULLER,NULL));
  CPPUNIT_ASSERT_EQUAL( true, a_condgaussian.SampleFrom(a_sample,CHOLESKY,NULL));

  /* Test dfGet */
  for (unsigned int i=0; i < NUM_COND_ARGS; i++)
  {
    CPPUNIT_ASSERT_EQUAL( v[i], a_condgaussian.dfGet(i));
  }
  
  /* Setting and Getting mean/covariance */
  // calculate expected value
  ColumnVector exp(DIMENSION); exp = 0.0;
  for (unsigned int i=0; i < NUM_COND_ARGS; i++)
  {
     exp = exp + v[i]*cond_args[i];
  }
  exp += _mu;
  CPPUNIT_ASSERT_EQUAL( exp, a_condgaussian.ExpectedValueGet());
  CPPUNIT_ASSERT_EQUAL( _sigma, a_condgaussian.CovarianceGet());

  /* Copy Constructor etc */
  LinearAnalyticConditionalGaussian b_condgaussian(a_condgaussian);
  CPPUNIT_ASSERT_EQUAL( _mu, b_condgaussian.AdditiveNoiseMuGet());
  CPPUNIT_ASSERT_EQUAL( _sigma, b_condgaussian.AdditiveNoiseSigmaGet());

  /* Setting and getting matrices */
  Matrix a2(DIMENSION,DIMENSION); a2 = 2.1;
  Matrix b2(DIMENSION,DIMENSION); b2 = 1.5;
  
  a_condgaussian.MatrixSet(0,a2);
  a_condgaussian.MatrixSet(1,b2);

  CPPUNIT_ASSERT_EQUAL( a2, a_condgaussian.MatrixGet(0));
  CPPUNIT_ASSERT_EQUAL( b2, a_condgaussian.MatrixGet(1));

  /* Setting and Getting mean/covariance  of the additive noise*/
  ColumnVector mu2;
  SymmetricMatrix sigma2;
  mu2.resize(DIMENSION);
  mu2(1) = MU_1; mu2(2) = MU_0;
  sigma2.resize(DIMENSION);
  sigma2 = 0.0;
  double sig2 = 0.05; 
  for (unsigned int rows=1; rows < DIMENSION + 1; rows++){ sigma2(rows,rows)=sig2; }

  a_condgaussian.AdditiveNoiseMuSet(mu2);
  a_condgaussian.AdditiveNoiseSigmaSet(sigma2);

  CPPUNIT_ASSERT_EQUAL( mu2, a_condgaussian.AdditiveNoiseMuGet());
  CPPUNIT_ASSERT_EQUAL( sigma2, a_condgaussian.AdditiveNoiseSigmaGet());

  /* What would be the best way to test ProbabilityGet? */

}

void
PdfTest::testDiscreteConditionalPdf()
{
  int cond_arg_dims[NUM_COND_ARGS_DIS] = { NUM_DS };
  int state_k;int state_kMinusOne;
  
  DiscreteConditionalPdf a_discretecondpdf(NUM_DS,NUM_COND_ARGS_DIS,cond_arg_dims);

  /* Get the dimension  */ 
  CPPUNIT_ASSERT_EQUAL( 1, (int)a_discretecondpdf.DimensionGet());
  
  /* Get the dimension  */ 
  CPPUNIT_ASSERT_EQUAL( NUM_DS, (int)a_discretecondpdf.NumStatesGet());

  /* Get number of conditional arguments  */ 
  CPPUNIT_ASSERT_EQUAL( NUM_COND_ARGS_DIS, (int)a_discretecondpdf.NumConditionalArgumentsGet());

  std::vector<int> cond_args(NUM_COND_ARGS_DIS);
  /* Set and Get all Probabilities*/
  double prob_diag = 0.9;
  double prob_nondiag = (1-0.9)/(NUM_DS-1);
  for (state_kMinusOne = 0 ; state_kMinusOne < NUM_DS ;  state_kMinusOne++)
    {
       cond_args[0] = state_kMinusOne;
       for (state_k = 0 ; state_k < NUM_DS ;  state_k++)
	     {
           if (state_kMinusOne == state_k) a_discretecondpdf.ProbabilitySet(prob_diag,state_k,cond_args);
           else a_discretecondpdf.ProbabilitySet(prob_nondiag,state_k,cond_args);
 	     }
    }

  /* Set and Get Conditional Arguments */
  int cond_arg = NUM_DS - 1;
  a_discretecondpdf.ConditionalArgumentSet(0, cond_arg);
  CPPUNIT_ASSERT_EQUAL( cond_arg, a_discretecondpdf.ConditionalArgumentGet(0));
  
  /* Get the probability for the states given the conditional argument set */
  for (state_k = 0 ; state_k < NUM_DS ;  state_k++)
    {
      if( state_k == cond_arg) CPPUNIT_ASSERT_EQUAL( prob_diag, (double)a_discretecondpdf.ProbabilityGet(state_k));
      else CPPUNIT_ASSERT_EQUAL( prob_nondiag, (double)a_discretecondpdf.ProbabilityGet(state_k));
    }

  /* Sampling */
  vector<Sample<int> > los(NUM_SAMPLES);
  CPPUNIT_ASSERT_EQUAL( true, a_discretecondpdf.SampleFrom(los,NUM_SAMPLES,DEFAULT,NULL));
  Sample<int> a_sample ;
  CPPUNIT_ASSERT_EQUAL( true, a_discretecondpdf.SampleFrom(a_sample,DEFAULT,NULL));

}

void
PdfTest::testMcpdf()
{
  /* Set and Get Dimension and number of samples*/
  MCPdf<ColumnVector> a_mcpdf(NUM_SAMPLES,DIMENSION);
  CPPUNIT_ASSERT_EQUAL( DIMENSION, a_mcpdf.DimensionGet());
  CPPUNIT_ASSERT_EQUAL( NUM_SAMPLES , a_mcpdf.NumSamplesGet());

  // Generating (exact) Samples from a Gaussian density with 
  // cholesky sampling
  Gaussian gaussian(_mu,_sigma);

  vector<Sample<ColumnVector> > exact_samples(NUM_SAMPLES);
  vector<Sample<ColumnVector> >::iterator it;
  gaussian.SampleFrom(exact_samples, NUM_SAMPLES,CHOLESKY,NULL);

  /* Getting and setting the list of samples (non-weighted)*/
  a_mcpdf.ListOfSamplesSet(exact_samples);
  const vector<WeightedSample<ColumnVector> > mcpdf_samples = a_mcpdf.ListOfSamplesGet();
  for (unsigned int i = 0; i < NUM_SAMPLES ; i++)
  {
     CPPUNIT_ASSERT_EQUAL( exact_samples[i].ValueGet(), mcpdf_samples[i].ValueGet());
     CPPUNIT_ASSERT_EQUAL( exact_samples[i].ValueGet(), a_mcpdf.SampleGet(i).ValueGet());
     CPPUNIT_ASSERT_EQUAL( 1.0/NUM_SAMPLES, mcpdf_samples[i].WeightGet());
  }

  /* List of samples update + getting and setting the list of samples (weighted)*/
  vector<WeightedSample<ColumnVector> > samples_weighted = mcpdf_samples;
  for (unsigned int i = 0 ; i < NUM_SAMPLES ; i++)  
  {
     //set a weight
     samples_weighted[i].WeightSet(i+1);
  }
  double tot_weight = (double)(NUM_SAMPLES+1)*((double)NUM_SAMPLES)/2.0;
  CPPUNIT_ASSERT_EQUAL( true , a_mcpdf.ListOfSamplesUpdate(samples_weighted) );
  for (unsigned int i = 0; i < NUM_SAMPLES ; i++)
  {
     CPPUNIT_ASSERT_EQUAL( samples_weighted[i].ValueGet(), a_mcpdf.SampleGet(i).ValueGet());
     CPPUNIT_ASSERT_EQUAL( (double)(samples_weighted[i].WeightGet())/tot_weight, a_mcpdf.SampleGet(i).WeightGet());
  }

  /* Copy Constructor etc */
  MCPdf<ColumnVector> b_mcpdf(a_mcpdf);
  for (unsigned int i = 0; i < NUM_SAMPLES ; i++)
  {
     CPPUNIT_ASSERT_EQUAL( a_mcpdf.SampleGet(i).ValueGet(), b_mcpdf.SampleGet(i).ValueGet());
     CPPUNIT_ASSERT_EQUAL( a_mcpdf.SampleGet(i).WeightGet(), b_mcpdf.SampleGet(i).WeightGet());
  }

  /* Sampling */
  vector<Sample<ColumnVector> > samples_test(NUM_SAMPLES);
  CPPUNIT_ASSERT_EQUAL( true, a_mcpdf.SampleFrom(samples_test,NUM_SAMPLES,RIPLEY,NULL));
  CPPUNIT_ASSERT_EQUAL( true, a_mcpdf.SampleFrom(samples_test,NUM_SAMPLES,DEFAULT,NULL));

  /* Expected Value*/
  vector<WeightedSample<ColumnVector> > los = a_mcpdf.ListOfSamplesGet();
  vector<WeightedSample<ColumnVector> >::iterator it2;
  ColumnVector cumSum(DIMENSION);
  cumSum=0.0;
  double sumWeights = 0.0;
  for ( it2 = los.begin() ; it2!= los.end() ; it2++ )
  {
    cumSum += ( it2->ValueGet() * it2->WeightGet() );
    sumWeights += it2->WeightGet();
  }
  CPPUNIT_ASSERT_EQUAL( approxEqual(cumSum/sumWeights, a_mcpdf.ExpectedValueGet(), epsilon),true);
 
  /* Covariance  + cumsumPDF*/
  ColumnVector mean(a_mcpdf.ExpectedValueGet());
  ColumnVector diff(DIMENSION); // Temporary storage
  Matrix diffsum(DIMENSION, DIMENSION);
  diffsum = 0.0;
  vector<double> cumPDF = a_mcpdf.CumulativePDFGet();
  vector<double>::iterator cumPDFit;
  cumPDFit = cumPDF.begin(); *cumPDFit = 0.0;
  double cumSumW = 0.0;
  for (it2 = los.begin(); it2 != los.end(); it2++)
  {
    diff = (it2->ValueGet() - mean);
    diffsum += diff * (diff.transpose() * it2->WeightGet());
    cumPDFit++;
    cumSumW += ( it2->WeightGet() / sumWeights);
    // test cumulative sum
    CPPUNIT_ASSERT_EQUAL(approxEqual(cumSumW, *cumPDFit, epsilon), true); 
  }
  CPPUNIT_ASSERT_EQUAL(approxEqual(diffsum/sumWeights, (Matrix)a_mcpdf.CovarianceGet(),epsilon),true);

  /* ProbabilityGet */ 

  /**************************
  // MCPDF with unsigned int
  *************************/

  /* Set and Get Dimension and number of samples*/
  MCPdf<unsigned int> a_mcpdf_uint(NUM_SAMPLES,1);
  unsigned int one = 1;
  CPPUNIT_ASSERT_EQUAL(one, a_mcpdf_uint.DimensionGet());
  CPPUNIT_ASSERT_EQUAL( NUM_SAMPLES , a_mcpdf_uint.NumSamplesGet());


  // Generating (exact) Samples from a discrete pdf 
  unsigned int num_states = 10;
  DiscretePdf discrete(num_states);

  vector<Sample<int> > samples_discrete_int(NUM_SAMPLES);
  vector<Sample<int> >::iterator it_discrete_int;
  discrete.SampleFrom(samples_discrete_int, NUM_SAMPLES);

  vector<Sample<unsigned int> > samples_discrete(NUM_SAMPLES);
  vector<Sample<unsigned int> >::iterator it_discrete;
  it_discrete = samples_discrete.begin(); 

  Sample<unsigned int> temp_sample;
  for(it_discrete_int = samples_discrete_int.begin(); it_discrete_int != samples_discrete_int.end(); it_discrete_int++)
  {
    temp_sample.ValueSet((*it_discrete_int).ValueGet());
    (*it_discrete)= temp_sample;
    it_discrete++;
}

  /* Getting and setting the list of samples (non-weighted)*/
  a_mcpdf_uint.ListOfSamplesSet(samples_discrete);
  const vector<WeightedSample<unsigned int> > mcpdf_samples_uint = a_mcpdf_uint.ListOfSamplesGet();
  for (unsigned int i = 0; i < NUM_SAMPLES ; i++)
  {
     CPPUNIT_ASSERT_EQUAL( samples_discrete[i].ValueGet(), mcpdf_samples_uint[i].ValueGet());
     CPPUNIT_ASSERT_EQUAL( samples_discrete[i].ValueGet(), a_mcpdf_uint.SampleGet(i).ValueGet());
     CPPUNIT_ASSERT_EQUAL( 1.0/NUM_SAMPLES, mcpdf_samples_uint[i].WeightGet());
  }

  /* List of samples update + getting and setting the list of samples (weighted)*/
  vector<WeightedSample<unsigned int > > samples_weighted_uint = mcpdf_samples_uint;
  for (unsigned int i = 0 ; i < NUM_SAMPLES ; i++)  
  {
     //set a weight
     samples_weighted_uint[i].WeightSet(i+1);
  }
  double tot_weight_uint = (double)(NUM_SAMPLES+1)*((double)NUM_SAMPLES)/2.0;
  CPPUNIT_ASSERT_EQUAL( true , a_mcpdf_uint.ListOfSamplesUpdate(samples_weighted_uint) );
  for (unsigned int i = 0; i < NUM_SAMPLES ; i++)
  {
     CPPUNIT_ASSERT_EQUAL( samples_weighted_uint[i].ValueGet(), a_mcpdf_uint.SampleGet(i).ValueGet());
     CPPUNIT_ASSERT_EQUAL( (double)(samples_weighted_uint[i].WeightGet())/tot_weight_uint, a_mcpdf_uint.SampleGet(i).WeightGet());
  }

  ///* Copy Constructor etc */
  MCPdf<unsigned int> b_mcpdf_uint(a_mcpdf_uint);
  for (unsigned int i = 0; i < NUM_SAMPLES ; i++)
  {
     CPPUNIT_ASSERT_EQUAL( a_mcpdf_uint.SampleGet(i).ValueGet(), b_mcpdf_uint.SampleGet(i).ValueGet());
     CPPUNIT_ASSERT_EQUAL( a_mcpdf_uint.SampleGet(i).WeightGet(), b_mcpdf_uint.SampleGet(i).WeightGet());
  }

  /* Sampling */
  vector<Sample<unsigned int> > samples_test_uint(NUM_SAMPLES);
  CPPUNIT_ASSERT_EQUAL( true, a_mcpdf_uint.SampleFrom(samples_test_uint,NUM_SAMPLES,DEFAULT,NULL));

  /* Expected Value*/
  vector<WeightedSample<unsigned int> > los_uint = a_mcpdf_uint.ListOfSamplesGet();
  vector<WeightedSample<unsigned int> >::iterator it2_uint;
  double cumSum_double;
  cumSum_double=0.0;
  sumWeights = 0.0;
  for ( it2_uint = los_uint.begin() ; it2_uint!= los_uint.end() ; it2_uint++ )
  {
    cumSum_double += ( (double)it2_uint->ValueGet() * it2_uint->WeightGet() );
    sumWeights += it2_uint->WeightGet();
  }
  CPPUNIT_ASSERT_EQUAL( approxEqual( (unsigned int ) (cumSum_double/sumWeights + 0.5) , a_mcpdf_uint.ExpectedValueGet(), epsilon),true);
 
  /* Covariance  + cumsumPDF*/
  unsigned int  mean_uint = a_mcpdf_uint.ExpectedValueGet();
  unsigned int  diff_uint;
  double diffsum_uint;
  diffsum_uint = 0.0;
  vector<double> cumPDF_uint = a_mcpdf_uint.CumulativePDFGet();
  vector<double>::iterator cumPDFit_uint;
  cumPDFit_uint = cumPDF_uint.begin(); *cumPDFit_uint = 0.0;
  double cumSumW_uint = 0.0;
  for (it2_uint = los_uint.begin(); it2_uint != los_uint.end(); it2_uint++)
  {
    diff_uint =  (it2_uint->ValueGet() - mean_uint);
    diffsum_uint += (double)(diff_uint * diff_uint)  * it2_uint->WeightGet();
    cumPDFit_uint++;
    cumSumW_uint += ( it2_uint->WeightGet() / sumWeights);
    // test cumulative sum
    CPPUNIT_ASSERT_EQUAL(approxEqual(cumSumW_uint, *cumPDFit_uint, epsilon), true); 
  }
  Matrix test_diff(1,1);
  test_diff(1,1) = diffsum_uint/sumWeights;
  CPPUNIT_ASSERT_EQUAL(approxEqual(test_diff, (Matrix)a_mcpdf_uint.CovarianceGet(),epsilon),true);
  /* ProbabilityGet */ 
}
