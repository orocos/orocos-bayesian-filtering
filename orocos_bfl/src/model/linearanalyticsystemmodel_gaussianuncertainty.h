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
#ifndef __LINEAR_ANALYTIC_SYSTEM_MODEL_GAUSSIAN_UNCERTAINTY__
#define __LINEAR_ANALYTIC_SYSTEM_MODEL_GAUSSIAN_UNCERTAINTY__

#include "systemmodel.h"
#include "../pdf/gaussian.h"
#include "analyticsystemmodel_gaussianuncertainty.h"
#include "../pdf/linearanalyticconditionalgaussian.h"

namespace BFL
{

  /// Class for linear analytic systemmodels with additive gaussian noise
  /** This class represents all systemmodels of the form
      \f[ x_k = A \times x_{k-1} + B \times u_{k} + N(\mu,\Sigma) \f]
  */
  class LinearAnalyticSystemModelGaussianUncertainty : public AnalyticSystemModelGaussianUncertainty
    {
    public:
      /// Constructor
      /**@pre LinearAnalyticConditionalGaussian should have 1/2 conditional
	 Arguments (checked) and the first conditional argument should
	 be x!
	 @param pdf Conditional pdf with Gaussian uncertainty
       */
      LinearAnalyticSystemModelGaussianUncertainty( LinearAnalyticConditionalGaussian* pdf);

      // Default Copy Constructor will do

      /// Destructor
      virtual ~LinearAnalyticSystemModelGaussianUncertainty();

      /// Set Matrix A
      /** This can be particularly useful for time-varying systems
	  @param a Matrix a
      */
      void ASet(const MatrixWrapper::Matrix & a);
      /// Set Matrix B
      /** This can be particularly useful for time-varying systems
	  @param b Matrix b
      */
      void BSet(const MatrixWrapper::Matrix & b);

      /// Get Matrix A
      const MatrixWrapper::Matrix& AGet() const;

      /// Get Matrix B
      const MatrixWrapper::Matrix& BGet() const;

    };

} // End namespace BFL

#endif // __LINEAR_ANALYTIC_SYSTEM_MODEL_GAUSSIAN_UNCERTAINTY__
