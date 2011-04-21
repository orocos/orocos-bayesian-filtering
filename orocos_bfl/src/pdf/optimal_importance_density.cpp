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

#include "optimal_importance_density.h"
#include "../wrappers/rng/rng.h" // Wrapper around several rng
                                 // libraries

namespace BFL
{
#define OptImpDensity OptimalImportanceDensity

  OptImpDensity::OptImpDensity(AnalyticConditionalGaussian * SystemPdf,
			       LinearAnalyticConditionalGaussian * MeasPdf)
    : AnalyticConditionalGaussian(SystemPdf->DimensionGet(),
			  SystemPdf->NumConditionalArgumentsGet()
			  + MeasPdf->NumConditionalArgumentsGet()),
      _SystemPdf(SystemPdf),
      _MeasPdf(MeasPdf)
  {
    Matrix tmp((SystemPdf->DimensionGet()),(SystemPdf->DimensionGet()));
    tmp = (SystemPdf->AdditiveNoiseSigmaGet().inverse()
	   + ( MeasPdf->MatrixGet(0).transpose() * (Matrix) (MeasPdf->AdditiveNoiseSigmaGet().inverse()) * MeasPdf->MatrixGet(0) )).inverse();
    tmp.convertToSymmetricMatrix(this->_additiveNoise_Sigma);
    this->_additiveNoise_Mu = 0.0;
  }

  OptImpDensity::~OptImpDensity(){};

  ColumnVector
  OptImpDensity::ExpectedValueGet() const
  {
    ColumnVector mean(DimensionGet()); mean = 0.0;
    ColumnVector arg;
    for (int i=0; i < NumConditionalArgumentsGet() ; i++)
      {
	arg = ConditionalArgumentGet(i);
	mean += (ColumnVector) (this->MatrixGet(i) * arg);
      }
    mean += AdditiveNoiseMuGet();
    return mean;
  }

  SymmetricMatrix
  OptImpDensity::CovarianceGet() const
  {
    return AdditiveNoiseSigmaGet();
  }

  Matrix
  OptImpDensity::dfGet(int i) const
  {
    assert(i < NumConditionalArgumentsGet() && i >= 0);
    return _ratio[i];
  }

  void
  OptImpDensity::NumConditionalArgumentsSet(int numconditionalarguments)
  {
    cerr << "You probably don't want to use this one, do you?" << endl;
  }

  void
  OptImpDensity::MatrixSet(int i, const Matrix & m)
  {
    assert(i < NumConditionalArgumentsGet());
    _ratio[i] = m;
  }

  const Matrix&
  OptImpDensity::MatrixGet(int i) const
  {
    assert(i < NumConditionalArgumentsGet());
    return _ratio[i];
  }

} // End namespace BFL
