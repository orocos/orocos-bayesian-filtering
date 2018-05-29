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

#ifndef __OPTIMAL_IMPORTANCE_DENSITY__
#define __OPTIMAL_IMPORTANCE_DENSITY__

#include "analyticconditionalgaussian.h"

namespace BFL
{
  /// Optimal importance density for Nonlinear Gaussian SS Models
  /**
     Describes the optimal importance density for all systems of the
     form
     \f[ x_k = f(x_{k-1}) + v_k, \quad v_k \sim N(0, \Sigma_v) \f]
     \f[ z_k = H x_k + w_k, \quad w_k \sim N(0, \Sigma_w) \f]

     This means all systems with a system equation that uses a
     AnalyticConditionalGaussian Class and a measurement equation that uses a
     LinearAnalyticConditionalGaussian class
  */
  class OptimalImportanceDensity : public AnalyticConditionalGaussian
    {
    public:
      /// Constructor
      /** @param SystemPdf
	  @param MeasPdf
      */
      OptimalImportanceDensity(AnalyticConditionalGaussian * SystemPdf,
			       LinearAnalyticConditionalGaussian * MeasPdf);

      // Default copy constructor

      /// Destructor
      virtual ~OptimalImportanceDensity();

      // redefine pure virtual functions
      virtual ColumnVector    ExpectedValueGet() const;
      virtual SymmetricMatrix CovarianceGet()    const;
      virtual Matrix          dfGet(int i)       const;

    private:
      AnalyticConditionalGaussian * _SystemPdf;
      LinearAnalyticConditionalGaussian * _MeasPdf;

    };

} // End namespace BFL

#include "optimal_importance_density.cpp"

#endif //  __OPTIMAL_IMPORTANCE_DENSITY__
