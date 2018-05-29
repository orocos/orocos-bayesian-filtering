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

#ifndef __EKF_PROP_DENSITY__
#define __EKF_PROP_DENSITY__

#include "filterproposaldensity.h"
#include "../filter/extendedkalmanfilter.h"

namespace BFL
{
  /// Proposal Density for non-linear systems with additive Gaussian Noise (using a EKF Filter)
  /** Calculates an importance density for all systems of the
      form
      \f[ x_k = f(x_{k-1}[,u_k]) + v_k, \quad v_k \sim N(0, \Sigma_v) \f]
      \f[ z_k = h(x_k[,s_k]) + w_k, \quad w_k \sim N(0, \Sigma_w) \f]

     This means all systems with a system equation and measurement equation that use a
     AnalyticConditionalGaussian Class.
  */
  class EKFProposalDensity : public FilterProposalDensity
    {
    public:
      /// Constructor
      /** @param SysModel
	  @param MeasModel
      */
      EKFProposalDensity(AnalyticSystemModelGaussianUncertainty * SysModel,
			 AnalyticMeasurementModelGaussianUncertainty * MeasModel);

      /// Destructor
      virtual ~EKFProposalDensity();

    };

} // End namespace BFL

#endif //  __FILTER_PROP_DENSITY__
