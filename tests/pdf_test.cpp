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
const unsigned int NUM_COND_ARGS = 2;
const unsigned int DIMENSION = 2;
const unsigned int NUM_SAMPLES = 100;

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
  vector<Sample<ColumnVector> > los(NUM_SAMPLES);
  CPPUNIT_ASSERT_EQUAL( true, a_gaussian.SampleFrom(los,NUM_SAMPLES,DEFAULT,NULL));
  // Box-Muller is not implemented yet
  CPPUNIT_ASSERT_EQUAL( false, a_gaussian.SampleFrom(los,NUM_SAMPLES,BOXMULLER,NULL));
  CPPUNIT_ASSERT_EQUAL( true, a_gaussian.SampleFrom(los,NUM_SAMPLES,CHOLESKY,NULL));
  Sample<ColumnVector> a_sample ;
  CPPUNIT_ASSERT_EQUAL( true, a_gaussian.SampleFrom(a_sample,DEFAULT,NULL));
  // Box-Muller is not implemented yet
  CPPUNIT_ASSERT_EQUAL( false, a_gaussian.SampleFrom(a_sample,BOXMULLER,NULL));
  CPPUNIT_ASSERT_EQUAL( true, a_gaussian.SampleFrom(a_sample,CHOLESKY,NULL));
  
  /* Setting and Getting mean/covariance */
  CPPUNIT_ASSERT_EQUAL( _mu, a_gaussian.ExpectedValueGet());
  CPPUNIT_ASSERT_EQUAL( _sigma, a_gaussian.CovarianceGet());

  /* Copy Constructor etc */
  Gaussian b_gaussian(a_gaussian);
  CPPUNIT_ASSERT_EQUAL( _mu, b_gaussian.ExpectedValueGet());
  CPPUNIT_ASSERT_EQUAL( _sigma, b_gaussian.CovarianceGet());
}  

void
PdfTest::testDiscretePdf()
{
  DiscretePdf a_discretepdf(NUM_DS);
  ColumnVector uniform_vector(NUM_DS);
  uniform_vector = a_discretepdf.ProbabilitiesGet();
 
  /* Check initial uniform distribution */
  for (int i = 0; i < NUM_DS ; i++)
    {
       CPPUNIT_ASSERT_EQUAL( 1.0/NUM_DS, (double)a_discretepdf.ProbabilityGet(i));
    }

  /* Set and Get probabilities */
  // one at a time
  double prob_new = 0.57;
  for (int new_el = 0; new_el < NUM_DS; new_el++)
    {
       CPPUNIT_ASSERT_EQUAL(true, a_discretepdf.ProbabilitySet(new_el,prob_new));
       CPPUNIT_ASSERT_EQUAL(prob_new, (double) a_discretepdf.ProbabilityGet(new_el));
    }
  ColumnVector prob_vec(NUM_DS);
  prob_vec(1) = 0.9;
  for (int i = 2; i < NUM_DS + 1 ; i++)
    {
       prob_vec(i) = (1-0.9)/(NUM_DS-1);
    }
  a_discretepdf.ProbabilitiesSet(prob_vec);
  CPPUNIT_ASSERT_EQUAL( prob_vec, a_discretepdf.ProbabilitiesGet());
 
  /* Sampling */
  vector<Sample<int> > los(NUM_SAMPLES);
  vector<Sample<int> >::iterator it;
  CPPUNIT_ASSERT_EQUAL( true, a_discretepdf.SampleFrom(los,NUM_SAMPLES,DEFAULT,NULL));
  // Ripley
  CPPUNIT_ASSERT_EQUAL( true, a_discretepdf.SampleFrom(los,NUM_SAMPLES,RIPLEY,NULL));

  Sample<int> a_sample ;
  CPPUNIT_ASSERT_EQUAL( true, a_discretepdf.SampleFrom(a_sample,DEFAULT,NULL));
  // Ripley not implemented for one sample
  CPPUNIT_ASSERT_EQUAL( false, a_discretepdf.SampleFrom(a_sample,RIPLEY,NULL));
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
  CPPUNIT_ASSERT_EQUAL( NUM_DS, (int)a_discretecondpdf.DimensionGet());
  
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
  //create mcpdf
  unsigned int num_samples = 100;
  unsigned int dimension = 10;
  MCPdf<ColumnVector> a_mcpdf(num_samples,dimension);
  
  /* Get Dimension and number of samples*/
  CPPUNIT_ASSERT_EQUAL( dimension, a_mcpdf.DimensionGet());
  CPPUNIT_ASSERT_EQUAL( num_samples , a_mcpdf.NumSamplesGet());

  /* Set and Get Dimension and number of samples*/
  a_mcpdf.DimensionSet(DIMENSION);
  a_mcpdf.NumSamplesSet(NUM_SAMPLES);
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
}
