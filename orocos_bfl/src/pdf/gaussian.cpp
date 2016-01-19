// $Id$
// Copyright (C) 2002 Klaas Gadeyne <first dot last at gmail dot com>
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
#include "gaussian.h"

#include "../wrappers/rng/rng.h" // Wrapper around several rng libraries

#include <cmath> // For sqrt and exp
#include <cassert>

namespace BFL
{
  using namespace MatrixWrapper;

  Gaussian::Gaussian (const ColumnVector& m, const SymmetricMatrix& s)
    : Pdf<ColumnVector> ( m.rows() )
    , _diff(DimensionGet())
    , _tempColumn(DimensionGet())
    , _samples(DimensionGet())
    , _sampleValue(DimensionGet())
    , _Low_triangle(DimensionGet(),DimensionGet())
  {
    // check if inputs are consistent
    assert (m.rows() == s.columns());
    _Mu = m;
    _Sigma = s;
    _Sigma_inverse.resize(DimensionGet());
    _Sigma_changed = true;
  }

  Gaussian::Gaussian (int dimension)
    : Pdf<ColumnVector>(dimension)
    , _diff(dimension)
    , _tempColumn(DimensionGet())
    , _samples(dimension)
    , _sampleValue(dimension)
    , _Low_triangle(dimension,dimension)
  {
    _Mu.resize(dimension);
    _Sigma.resize(dimension);
    _Sigma_inverse.resize(dimension);
    _Sigma_changed = true;
  }

  Gaussian::~Gaussian(){}

  std::ostream& operator<< (std::ostream& os, const Gaussian& g)
  {
    os << "\nMu:\n"    << g.ExpectedValueGet()
       << "\nSigma:\n" << g.CovarianceGet() << endl;
    return os;
  }

  //Clone function
  Gaussian* Gaussian::Clone() const
  {
      return new Gaussian(*this);
  }

  Probability Gaussian::ProbabilityGet(const ColumnVector& input) const
  {
    // only calculate these variables if sigma has changed
    if (_Sigma_changed){
      _Sigma_changed = false;
      _Sigma_inverse = _Sigma.inverse();
      _sqrt_pow = 1 / sqrt(pow(M_PI*2,(double)DimensionGet()) * _Sigma.determinant());
    }

    _diff = input;
    _diff -= _Mu;
    _Sigma_inverse.multiply(_diff,_tempColumn);
    //_tempColumn = _Sigma_inverse * _diff; 
    Probability temp = _diff.transpose() * _tempColumn;
    //Probability temp = _diff.transpose() * (_Sigma_inverse * _diff);
    Probability result = exp(-0.5 * temp) * _sqrt_pow;
    return result;
  }

  // Redefined for optimal performance.  Eg. do Cholesky decomposition
  // only once when drawing multiple samples at once!
  // See method below for more info regarding the algorithms
  bool
  Gaussian::SampleFrom (vector<Sample<ColumnVector> >& list_samples, const unsigned int num_samples, int method, void * args) const
  {
    list_samples.resize(num_samples); // will break real-timeness if list_samples.size()!=num_samples
    vector<Sample<ColumnVector> >::iterator rit = list_samples.begin();
    switch(method)
      {
      case DEFAULT: // Cholesky Sampling
      case CHOLESKY:
	{
	  bool result = _Sigma.cholesky_semidefinite(_Low_triangle);
	  while (rit != list_samples.end())
	    {
	      for (unsigned int j=1; j < DimensionGet()+1; j++) _samples(j) = rnorm(0,1);
	      _sampleValue = _Low_triangle * _samples ;
	      _sampleValue +=  this->_Mu;
	      rit->ValueSet(_sampleValue);
	      rit++;
	    }
	  return result;
	}
      case BOXMULLER: // Implement box-muller here
	// Only for univariate distributions.
	return false;
      default:
	return false;
      }
  }


  bool
  Gaussian::SampleFrom (Sample<ColumnVector>& one_sample, int method, void * args) const
  {
    /*  Exact i.i.d. samples from a Gaussian can be drawn in several
	ways:
	- if the DimensionGet() = 1 or 2 (and the 2 variables are
	independant), we can use inversion sampling (Box-Muller
	method)
        - For larger dimensions, we use can use the Cholesky method or
	an approached based on conditional distributions.
	(see ripley87, P.98 (bibtex below)).  The Cholesky method is
	generally preferred and the only one implemented for now.
    */
    switch(method)
      {
      case DEFAULT: // Cholesky Sampling, see eg.
      case CHOLESKY: // Cholesky Sampling, see eg.
	/*
	  @Book{		  ripley87,
	  author	= {Ripley, Brian D.},
	  title		= {Stochastic Simulation},
	  publisher	= {John Wiley and Sons},
	  year		= {1987},
	  annote	= {ISBN 0271-6356, WBIB 1 519.245}
	  }
	  p.98
	*/
	{
      bool result = _Sigma.cholesky_semidefinite(_Low_triangle);

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
	  _sampleValue = (_Low_triangle * _samples) + this->_Mu;
	  one_sample.ValueSet(_sampleValue);
	  return result;
	}
      case BOXMULLER: // Implement box-muller here
	// Only for univariate distributions.
	return false;
      default:
	return false;
      }
  }


  ColumnVector
  Gaussian::ExpectedValueGet (  ) const
  {
    return _Mu;
  }

  SymmetricMatrix
  Gaussian::CovarianceGet () const
  {
    return _Sigma;
  }

  void
  Gaussian::ExpectedValueSet (const ColumnVector& mu)
  {
    _Mu = mu;
    if (this->DimensionGet() == 0)
      {
	this->DimensionSet(mu.rows());
      }
    assert(this->DimensionGet() == mu.rows());
  }

  void
  Gaussian::CovarianceSet (const SymmetricMatrix& cov)
  {
    _Sigma = cov;
    _Sigma_changed = true;
    if (this->DimensionGet() == 0)
      {
	this->DimensionSet(cov.rows());
      }
    assert(this->DimensionGet() == cov.rows());
  }

  void
  Gaussian::DimensionSet ( unsigned int dim )
  {
    Pdf<ColumnVector>::DimensionSet(dim);
    _diff.resize(DimensionGet());
    _tempColumn.resize(DimensionGet());
    _samples.resize(DimensionGet());
    _sampleValue.resize(DimensionGet());
    _Low_triangle.resize(DimensionGet(),DimensionGet());
  }

} // End namespace BFL
