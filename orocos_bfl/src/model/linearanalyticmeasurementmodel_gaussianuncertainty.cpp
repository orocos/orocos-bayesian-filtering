// $Id$
// Copyright (C) 2002 Klaas Gadeyne <first dot last at gmail dot com>
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
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
//

#include "linearanalyticmeasurementmodel_gaussianuncertainty.h"
#include "../sample/sample.h"
#include <cassert>

namespace BFL
{
  using namespace MatrixWrapper;



  LinearAnalyticMeasurementModelGaussianUncertainty::LinearAnalyticMeasurementModelGaussianUncertainty
  ( LinearAnalyticConditionalGaussian* pdf)
    : AnalyticMeasurementModelGaussianUncertainty( pdf )
  {}

  /*
  LinearAnalyticMeasurementModelGaussianUncertainty::LinearAnalyticMeasurementModelGaussianUncertainty
  ( const LinearAnalyticMeasurementModelGaussianUncertainty& l )
    : AnalyticMeasurementModelGaussianUncertainty()
  {
    this->_MeasurementPdf = new LinearAnalyticConditionalGaussian(*l.MeasurementPdfGet());
  }
  */

  LinearAnalyticMeasurementModelGaussianUncertainty::~LinearAnalyticMeasurementModelGaussianUncertainty()
  {}


  ColumnVector
  LinearAnalyticMeasurementModelGaussianUncertainty::PredictionGet(const ColumnVector& u,
								   const ColumnVector& x)
  {
    MeasurementPdfGet()->ConditionalArgumentSet(0,x);
    if (MeasurementPdfGet()->NumConditionalArgumentsGet() == 2) MeasurementPdfGet()->ConditionalArgumentSet(1,u);
    return MeasurementPdfGet()->ExpectedValueGet();
  }


  SymmetricMatrix
  LinearAnalyticMeasurementModelGaussianUncertainty::CovarianceGet(const ColumnVector& u,
								   const ColumnVector& x)
  {
    MeasurementPdfGet()->ConditionalArgumentSet(0,x);
    if (MeasurementPdfGet()->NumConditionalArgumentsGet() == 2) MeasurementPdfGet()->ConditionalArgumentSet(1,u);
    return  MeasurementPdfGet()->CovarianceGet();
  }


  Matrix
  LinearAnalyticMeasurementModelGaussianUncertainty::df_dxGet(const ColumnVector& u,
							      const ColumnVector& x)
  {
    MeasurementPdfGet()->ConditionalArgumentSet(0,x);
    if (MeasurementPdfGet()->NumConditionalArgumentsGet() == 2) MeasurementPdfGet()->ConditionalArgumentSet(1,u);
    return dynamic_cast<AnalyticConditionalGaussian *>(MeasurementPdfGet())->dfGet(0);
  }


  void
  LinearAnalyticMeasurementModelGaussianUncertainty::HSet(const Matrix& h)
  {
    dynamic_cast<LinearAnalyticConditionalGaussian *>(MeasurementPdfGet())->MatrixSet(0,h);
  }


  void
  LinearAnalyticMeasurementModelGaussianUncertainty::JSet(const Matrix& j)
  {
    dynamic_cast<LinearAnalyticConditionalGaussian *>(MeasurementPdfGet())->MatrixSet(1,j);
  }


  const Matrix&
  LinearAnalyticMeasurementModelGaussianUncertainty::HGet() const
  {
    return dynamic_cast<LinearAnalyticConditionalGaussian *>(_MeasurementPdf)->MatrixGet(0);
  }


  const Matrix&
  LinearAnalyticMeasurementModelGaussianUncertainty::JGet() const
  {
    return dynamic_cast<LinearAnalyticConditionalGaussian *>(_MeasurementPdf)->MatrixGet(1);
  }


}
