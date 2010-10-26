// $Id$
// Copyright (C) 2003 Klaas Gadeyne <first dot last at gmail dot com>
// Copyright (C) 2008 Tinne De Laet <first dot last at mech dot kuleuven dot be>
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
    , _diff(dim)
    , _Mu(dim)
    , _Low_triangle(dim,dim)
    , _samples(dim)
    , _SampleValue(dim)
  {}

  /// Destructor
  ConditionalGaussian::~ConditionalGaussian(){}

  //Clone function
  ConditionalGaussian* ConditionalGaussian::Clone() const
  {
      return new ConditionalGaussian(*this);
  }

  Probability
  ConditionalGaussian::ProbabilityGet(const ColumnVector& input) const
  {
    // Update Mu
    _Mu = ExpectedValueGet();
    _diff = input - _Mu;

    Probability temp = _diff.transpose() * (ColumnVector)(CovarianceGet().inverse() * _diff);
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
    _Mu = ExpectedValueGet();

    switch(method)
      {
      case DEFAULT: // Cholesky, see althere (bad implementation)
	{
	  bool result = CovarianceGet().cholesky_semidefinite(_Low_triangle);
	  for (unsigned int j=1; j < DimensionGet()+1; j++){_samples(j) = rnorm(0,1);}
	  _SampleValue = (_Low_triangle * _samples) + _Mu;
	  sample.ValueSet(_SampleValue);
	  return result;
	}
      case BOXMULLER: /// @todo Implement box-muller here
	{
	  cerr << "Box-Muller not implemented yet!" << endl;
	  return false;
	}
      case CHOLESKY: // Cholesky Sampling
	{
	  bool result = CovarianceGet().cholesky_semidefinite(_Low_triangle);
	  /* For now we keep it simple, and use the scythe library
	     (although wrapped) with the uRNG that it uses itself only */
	  /* Sample Gaussian._dimension samples from univariate
	     gaussian This could be done using several available
	     libraries, combined with different uniform RNG.  Both the
	     library to be used and the uRNG could be implemented as
	     #ifdef conditions, although I'm sure there must exist a
	     cleaner way to implement this!
	  */
	  for (unsigned int j=1; j < DimensionGet()+1; j++) _samples(j) = rnorm(0,1);
	  _SampleValue = (_Low_triangle * _samples) + _Mu;
	  sample.ValueSet(_SampleValue);
	  return result;
	}
      default:
	cerr << "Conditional Gaussian: Sampling method " << method
	     << "not implemented yet!" << endl;
	return false;
      }
  }

} // End namespace
