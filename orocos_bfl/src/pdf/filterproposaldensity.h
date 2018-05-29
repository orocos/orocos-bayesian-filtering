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
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
//

#ifndef __FILTER_PROP_DENSITY__
#define __FILTER_PROP_DENSITY__

#include "analyticconditionalgaussian.h"
#include "gaussian.h"
#include "../filter/filter.h"
#include "../model/analyticmeasurementmodel_gaussianuncertainty.h"
#include "../model/analyticsystemmodel_gaussianuncertainty.h"

namespace BFL
{
  /// Proposal Density for non-linear systems with additive Gaussian Noise (using a (analytic) Filter)
  /** Calculates an importance density for all systems of the
      form
      \f[ x_k = f(x_{k-1}[,u_k]) + v_k, \quad v_k \sim N(0, \Sigma_v) \f]
      \f[ z_k = h(x_k[,s_k]) + w_k, \quad w_k \sim N(0, \Sigma_w) \f]

     This means all systems with a system equation and measurement equation that use a
     AnalyticConditionalGaussian Class.
     It uses a Filter to generate a proposal
  */
  class FilterProposalDensity : public AnalyticConditionalGaussian
    {
    public:
      /// Constructor
      /** @param SysModel
	  @param MeasModel
      */
      FilterProposalDensity(AnalyticSystemModelGaussianUncertainty * SysModel,
			    AnalyticMeasurementModelGaussianUncertainty * MeasModel);

      /// Copy constructor
      /** @param fpd
	  @bug Not implemented yet
       */
      FilterProposalDensity(const FilterProposalDensity & fpd);

      /// Destructor
      virtual ~FilterProposalDensity();

      // redefine pure virtual functions
      virtual MatrixWrapper::ColumnVector    ExpectedValueGet() const;
      virtual MatrixWrapper::SymmetricMatrix CovarianceGet()    const;
      virtual MatrixWrapper::Matrix          dfGet(unsigned int i)       const;

      /// Set SystemModel
      /** @param SysModel
       */
      void SystemModelSet(AnalyticSystemModelGaussianUncertainty * SysModel);

      /// Set Measurementmodel
      /** @param MeasModel
       */
      void MeasurementModelSet(AnalyticMeasurementModelGaussianUncertainty * MeasModel);

      /// Set SampleCov
      /**
	 @param cov
      */
      void SampleCovSet(MatrixWrapper::SymmetricMatrix & cov);

    protected:
      mutable Gaussian * _TmpPrior;
      mutable Filter<MatrixWrapper::ColumnVector,MatrixWrapper::ColumnVector> * _filter;

      AnalyticSystemModelGaussianUncertainty * _sysmodel;
      AnalyticMeasurementModelGaussianUncertainty * _measmodel;

      MatrixWrapper::SymmetricMatrix _sample_cov;

      /// internal method
      virtual void FilterStep() const;

    };

} // End namespace BFL

#endif //  __FILTER_PROP_DENSITY__
