// $Id$
// Copyright (C) 2003 Klaas Gadeyne <first dot last at gmail dot com>
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
//  

#include "conditionalgaussian.h"
#include <cmath>
#include "../wrappers/rng/rng.h"

namespace BFL
{
  using namespace MatrixWrapper;

  ConditionalGaussian::ConditionalGaussian(int dim, 
					   int num_conditional_arguments)
    : ConditionalPdf<ColumnVector,ColumnVector>(dim, num_conditional_arguments)
  {}

  /// Destructor
  ConditionalGaussian::~ConditionalGaussian(){}

  Probability 
  ConditionalGaussian::ProbabilityGet(const ColumnVector& input) const
  {
    // Update Mu
    ColumnVector Mu = ExpectedValueGet();
    ColumnVector diff = input - Mu;
    
    Probability temp = diff.transpose() * (ColumnVector)(CovarianceGet().inverse() * diff);
    Probability result = exp(-0.5 * temp) / sqrt(pow(M_PI*2,(double)DimensionGet()) * CovarianceGet().determinant());
    return result;
  }

  bool
  ConditionalGaussian::SampleFrom (vector<Sample<ColumnVector> >& samples, const int num_samples, int method, void * args) const
  {
    return Pdf<ColumnVector>::SampleFrom(samples, num_samples, method, args);
  }

  bool
  ConditionalGaussian::SampleFrom (Sample<ColumnVector>& sample, int method, void * args) const
  {
    // Sampling from a Gaussian is simple if DIMENSION = 1 or 2 (and the
    // 2 variables are independant!)
    // Then we can use inversion sampling (Box-Muller method)
    // So for 1D, we use Box-Muller, else we use the cholesky method
    // These are both methods that don't require any arguments
  
    // Update mu
    ColumnVector Mu = ExpectedValueGet();
  
    switch(method)
      {
      case DEFAULT: // Cholesky, see althere (bad implementation)
	{
	  Matrix Low_triangle;
	  bool result = CovarianceGet().cholesky_semidefinite(Low_triangle);
	  ColumnVector samples(DimensionGet()); ColumnVector SampleValue(DimensionGet());
	  for (unsigned int j=1; j < DimensionGet()+1; j++){samples(j) = rnorm(0,1);}
	  SampleValue = (Low_triangle * samples) + Mu;
	  sample.ValueSet(SampleValue);
	  return result;
	}
      case BOXMULLER: /// @todo Implement box-muller here
	{
	  cerr << "Box-Muller not implemented yet!" << endl;
	  return false;
	}
      case CHOLESKY: // Cholesky Sampling
	{
	  Matrix Low_triangle;
	  bool result = CovarianceGet().cholesky_semidefinite(Low_triangle);
	  /* For now we keep it simple, and use the scythe library
	     (although wrapped) with the uRNG that it uses itself only */
	  ColumnVector samples(DimensionGet()); ColumnVector SampleValue(DimensionGet());
	  /* Sample Gaussian._dimension samples from univariate
	     gaussian This could be done using several available
	     libraries, combined with different uniform RNG.  Both the
	     library to be used and the uRNG could be implemented as
	     #ifdef conditions, although I'm sure there must exist a
	     cleaner way to implement this!
	  */
	  for (unsigned int j=1; j < DimensionGet()+1; j++) samples(j) = rnorm(0,1);
	  SampleValue = (Low_triangle * samples) + Mu;
	  sample.ValueSet(SampleValue);
	  return result;
	}
      default:
	cerr << "Conditional Gaussian: Sampling method " << method
	     << "not implemented yet!" << endl;
	return false;
      }
  }

} // End namespace 
