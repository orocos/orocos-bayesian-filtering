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
#ifndef __LINEAR_MEASUREMENT_MODEL_GAUSSIAN_UNCERTAINTY__
#define __LINEAR_MEASUREMENT_MODEL_GAUSSIAN_UNCERTAINTY__

#include "analyticmeasurementmodel_gaussianuncertainty.h"
#include "../pdf/gaussian.h"
#include "../pdf/linearanalyticconditionalgaussian.h"

namespace BFL
{

  /// Class for linear analytic measurementmodels with additive gaussian noise
  /** This class represents all measurementmodels of the form
      \f[ z_k = H \times x_k + J \times s_{k} + N(\mu,\Sigma) \f]
  */
  class LinearAnalyticMeasurementModelGaussianUncertainty :
    public AnalyticMeasurementModelGaussianUncertainty
    {
    public:
      /// Constructor
      /** @param pdf Conditional pdf, with Gaussian uncertainty
       */
      LinearAnalyticMeasurementModelGaussianUncertainty( LinearAnalyticConditionalGaussian* pdf = NULL);

      // Default Copy constructor will do
      // LinearAnalyticMeasurementModelGaussianUncertainty(const LinearAnalyticMeasurementModelGaussianUncertainty& l);

      // Destructor
      virtual ~LinearAnalyticMeasurementModelGaussianUncertainty();

      // redefinition of virtual functions
      virtual MatrixWrapper::Matrix          df_dxGet     (const MatrixWrapper::ColumnVector& u, const MatrixWrapper::ColumnVector& x);
      virtual MatrixWrapper::ColumnVector    PredictionGet(const MatrixWrapper::ColumnVector& u, const MatrixWrapper::ColumnVector& x);
      virtual MatrixWrapper::SymmetricMatrix CovarianceGet(const MatrixWrapper::ColumnVector& u, const MatrixWrapper::ColumnVector& x);

      /// Set Matrix H
      /** This can be particularly useful for time-varying systems
	  @param h Matrix H
      */
      void HSet(const MatrixWrapper::Matrix& h);

      /// Set Matrix J
      /** This can be particularly useful for time-varying systems
	  @param j Matrix J
      */
      void JSet(const MatrixWrapper::Matrix& j);

      /// Get Matrix H
      const MatrixWrapper::Matrix& HGet() const;

      /// Get Matrix J
      const MatrixWrapper::Matrix& JGet() const;


    protected:

    };

} // End namespace BFL

#endif // __LINEAR_MEASUREMENT_MODEL_GAUSSIAN_UNCERTAINTY__
