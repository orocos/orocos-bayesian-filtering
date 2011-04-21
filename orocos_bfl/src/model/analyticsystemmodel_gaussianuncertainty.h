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
#ifndef __SYSTEM_MODEL_GAUSSIANUNCERTAINTY__
#define __SYSTEM_MODEL_GAUSSIANUNCERTAINTY__

#include "../pdf/analyticconditionalgaussian.h"
#include "systemmodel.h"

namespace BFL
{

  /// Class for analytic system models with additive Gauss. uncertainty
  /** Class representing all analytic system Models with Additive
      Gaussian noise
  */
  class AnalyticSystemModelGaussianUncertainty: public SystemModel<MatrixWrapper::ColumnVector>
    {
    public:
      /// Constructor
      /** @param Systempdf AnalyticConditionalGaussian representing \f$
	  P(X_k | X_{k-1}, U_{k}) \f$
      */
      AnalyticSystemModelGaussianUncertainty(AnalyticConditionalGaussian* Systempdf);

      /// Default copy Constructor, interface class
      /* @param model The Analytic System Model with additive Gaussian
	  uncertainty to be copied
      */
      // AnalyticSystemModelGaussianUncertainty(const AnalyticSystemModelGaussianUncertainty& model);

      /// Destructor
      virtual ~AnalyticSystemModelGaussianUncertainty();

      /// Returns F-matrix
      /** \f[ F = \frac{df}{dx} \mid_{u,x} \f] used by kalman filter variants
	  @param u The value of the input in which the derivate is evaluated
	  @param x The value in the state in which the derivate is
	  evaluated
	  @bug Should actually be defined for _any_ continuous system
	  model!  There should be a class between this one and system
	  model tout court, not assuming gaussian uncertainty!
      */
      MatrixWrapper::Matrix df_dxGet(const MatrixWrapper::ColumnVector& u, const MatrixWrapper::ColumnVector& x);

      /// Returns prediction of state
      MatrixWrapper::ColumnVector PredictionGet(const MatrixWrapper::ColumnVector& u, const MatrixWrapper::ColumnVector& x);

      /// Covariance of system noise
      MatrixWrapper::SymmetricMatrix CovarianceGet(const MatrixWrapper::ColumnVector& u, const MatrixWrapper::ColumnVector& x);
    };

} // End namespace BFL

#endif // __SYSTEM_MODEL_GAUSSIANUNCERTAINTY__


