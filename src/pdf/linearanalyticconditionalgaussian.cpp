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

#include "linearanalyticconditionalgaussian.h"
#include "../wrappers/rng/rng.h" // Wrapper around several rng
                                 // libraries

namespace BFL
{
  using namespace MatrixWrapper;


  LinearAnalyticConditionalGaussian::LinearAnalyticConditionalGaussian(const vector<Matrix> & ratio,
						       const Gaussian& additiveNoise)
    : AnalyticConditionalGaussianAdditiveNoise(additiveNoise,ratio.size())
    , _ratio(ratio)
    , _mean_temp(DimensionGet())
    , _arg(DimensionGet())
  {
    // Initialise ConditionalArguments to 0
    ColumnVector arg;
    for (unsigned int i=0; i < NumConditionalArgumentsGet() ; i++)
      {
	arg.resize(_ratio[i].columns());
	arg = 0.0;
	ConditionalArgumentSet(i,arg);
      }
  }

  // Only one conditional argument
  LinearAnalyticConditionalGaussian::LinearAnalyticConditionalGaussian(const Matrix& a,
						       const Gaussian& additiveNoise)
    : AnalyticConditionalGaussianAdditiveNoise(additiveNoise,1)
    , _mean_temp(DimensionGet())
    , _arg(DimensionGet())
  {
    _ratio.resize(1);
    _ratio[0] = a;
    // Initialise ConditionalArguments to 0
    ColumnVector x(a.columns()); x = 0.0;
    ConditionalArgumentSet(0,x);
  }

  LinearAnalyticConditionalGaussian::~LinearAnalyticConditionalGaussian(){}

  //Clone function
  LinearAnalyticConditionalGaussian* LinearAnalyticConditionalGaussian::Clone() const
  {     
      return new LinearAnalyticConditionalGaussian(*this);
  }

  ColumnVector
  LinearAnalyticConditionalGaussian::ExpectedValueGet() const
  {
    _mean_temp = 0.0;
    for (unsigned int i=0; i < NumConditionalArgumentsGet() ; i++)
      {
	_arg = ConditionalArgumentGet(i);
	_mean_temp += (ColumnVector) (MatrixGet(i) * _arg);
      }
    _mean_temp += AdditiveNoiseMuGet();
    return _mean_temp;
  }

  Matrix
  LinearAnalyticConditionalGaussian::dfGet(unsigned int i) const
  {
    assert(i < NumConditionalArgumentsGet());
    return _ratio[i];
  }

  void
  LinearAnalyticConditionalGaussian::NumConditionalArgumentsSet(unsigned int numconditionalarguments)
  {
    ConditionalPdf<ColumnVector,ColumnVector>::NumConditionalArgumentsSet(numconditionalarguments);
    _ratio.resize(numconditionalarguments);
  }

  void
  LinearAnalyticConditionalGaussian::MatrixSet(unsigned int i, const Matrix & m)
  {
    assert(i < NumConditionalArgumentsGet());
    _ratio[i] = m;
  }

  const Matrix&
  LinearAnalyticConditionalGaussian::MatrixGet(unsigned int i) const
  {
    assert(i < NumConditionalArgumentsGet());
    return _ratio[i];
  }

} // End namespace BFL
