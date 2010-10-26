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
#ifndef __MEASUREMENT_MODEL_GAUSSIANUNCERTAINTY__
#define __MEASUREMENT_MODEL_GAUSSIANUNCERTAINTY__

#include "../pdf/analyticconditionalgaussian.h"
#include "measurementmodel.h"

namespace BFL
{

  /** Class representing all continuous analytic Measurement Models with
      additive Gaussian Uncertainty
  */
  class AnalyticMeasurementModelGaussianUncertainty: public MeasurementModel<MatrixWrapper::ColumnVector,MatrixWrapper::ColumnVector>
    {
    public:
      /// Constructor
      /** @param Measurementpdf ConditionalPdf<S,T> representing \f$ P(Z_k |
	  X_{k} (, U_{k})) \f$
	  @see MEASUREMENT_SIZE, STATE_SIZE, INPUT_SIZE, _MeasurementPdf
      */
      AnalyticMeasurementModelGaussianUncertainty(AnalyticConditionalGaussian* Measurementpdf = NULL);

      // default Copy Constructor, interface class

      /// Destructor
      virtual ~AnalyticMeasurementModelGaussianUncertainty();

      /// Returns H-matrix
      /** \f[ H = \frac{df}{dx} \mid_{u,x} \f] used by extended kalman filter
	  @param u The value of the input in which the derivate is evaluated
	  @param x The value in the state in which the derivate is evaluated
      */
      virtual MatrixWrapper::Matrix df_dxGet(const MatrixWrapper::ColumnVector& u, const MatrixWrapper::ColumnVector& x);

      /// Returns estimation of measurement
      virtual MatrixWrapper::ColumnVector PredictionGet(const MatrixWrapper::ColumnVector& u, const MatrixWrapper::ColumnVector& x);

      /// Returns covariance on the measurement
      virtual MatrixWrapper::SymmetricMatrix CovarianceGet(const MatrixWrapper::ColumnVector& u, const MatrixWrapper::ColumnVector& x);


    };

} // End namespace BFL

#endif // __MEASUREMENT_MODEL_GAUSSIANUNCERTAINTY__
