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
// Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
//

#include "analyticsystemmodel_gaussianuncertainty.h"

namespace BFL
{

  using namespace MatrixWrapper;


  // Constructor
  AnalyticSystemModelGaussianUncertainty::AnalyticSystemModelGaussianUncertainty
  (AnalyticConditionalGaussian* Systempdf)
    : SystemModel<ColumnVector>(Systempdf)
  {}

  // Copy constructor
  /*
  AnalyticSystemModelGaussianUncertainty::AnalyticSystemModelGaussianUncertainty
  (const AnalyticSystemModelGaussianUncertainty& model)
    : SystemModel<ColumnVector>(&(model.SystemPdfGet()))
  {}
  */

  // Destructor
  AnalyticSystemModelGaussianUncertainty::~AnalyticSystemModelGaussianUncertainty()
  {}


  Matrix
  AnalyticSystemModelGaussianUncertainty::df_dxGet(const ColumnVector& u, const ColumnVector& x)
  {
    SystemPdfGet()->ConditionalArgumentSet(0,x);
    if (SystemPdfGet()->NumConditionalArgumentsGet() == 2) SystemPdfGet()->ConditionalArgumentSet(1,u);
    return dynamic_cast<AnalyticConditionalGaussian *>(SystemPdfGet())->dfGet(0);
  }


  ColumnVector
  AnalyticSystemModelGaussianUncertainty::PredictionGet(const ColumnVector& u, const ColumnVector& x)
  {
    SystemPdfGet()->ConditionalArgumentSet(0,x);
    if (SystemPdfGet()->NumConditionalArgumentsGet() == 2) SystemPdfGet()->ConditionalArgumentSet(1,u);
    return SystemPdfGet()->ExpectedValueGet();
  }


  SymmetricMatrix
  AnalyticSystemModelGaussianUncertainty::CovarianceGet(const ColumnVector& u, const ColumnVector& x)
  {
    SystemPdfGet()->ConditionalArgumentSet(0,x);
    if (SystemPdfGet()->NumConditionalArgumentsGet() == 2) SystemPdfGet()->ConditionalArgumentSet(1,u);
    return dynamic_cast<AnalyticConditionalGaussian *>(SystemPdfGet())->CovarianceGet();
  }


} // End namespace BFL
