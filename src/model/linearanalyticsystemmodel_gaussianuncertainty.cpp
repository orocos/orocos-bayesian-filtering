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

#include "../sample/sample.h"
#include "linearanalyticsystemmodel_gaussianuncertainty.h"
#include "../pdf/linearanalyticconditionalgaussian.h"
#include <cassert>

namespace BFL
{
  using namespace MatrixWrapper;


  LinearAnalyticSystemModelGaussianUncertainty::LinearAnalyticSystemModelGaussianUncertainty
  ( LinearAnalyticConditionalGaussian* pdf)
    : AnalyticSystemModelGaussianUncertainty( pdf )
  {
  }

  LinearAnalyticSystemModelGaussianUncertainty::~LinearAnalyticSystemModelGaussianUncertainty()
  {}


  void LinearAnalyticSystemModelGaussianUncertainty::ASet(const Matrix& a)
  {
    dynamic_cast<LinearAnalyticConditionalGaussian *>(SystemPdfGet())->MatrixSet(0,a);
  }

  void LinearAnalyticSystemModelGaussianUncertainty::BSet(const Matrix& b)
  {
    dynamic_cast<LinearAnalyticConditionalGaussian *>(SystemPdfGet())->MatrixSet(1,b);
  }

  const Matrix&
  LinearAnalyticSystemModelGaussianUncertainty::AGet() const
  {
    return dynamic_cast<LinearAnalyticConditionalGaussian *>(_SystemPdf)->MatrixGet(0);
  }

  const Matrix&
  LinearAnalyticSystemModelGaussianUncertainty::BGet() const
  {
    return dynamic_cast<LinearAnalyticConditionalGaussian *>(_SystemPdf)->MatrixGet(1);
  }


} // End namespace BFL
