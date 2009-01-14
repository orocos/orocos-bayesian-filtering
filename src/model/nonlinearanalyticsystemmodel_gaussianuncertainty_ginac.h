// $Id$
// Copyright (C) 2002 Klaas Gadeyne <first dot last at gmail dot com>
//                    Wim Meeussen  <wim dot meeussen at mech dot kuleuven dot ac dot be>
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

#ifndef __NON_LINEAR_SYSTEM_MODEL_GAUSSIAN_UNCERTAINTY_GINAC__
#define __NON_LINEAR_SYSTEM_MODEL_GAUSSIAN_UNCERTAINTY_GINAC__

#include "analyticsystemmodel_gaussianuncertainty.h"
#include "../pdf/gaussian.h"
#include "../pdf/nonlinearanalyticconditionalgaussian_ginac.h"
#include <ginac/ginac.h>
#include <vector>
#include <iostream>

namespace BFL
{

  using namespace std;

  /// Class for nonlinear analytic systemmodels with additive gaussian noise
  /** This class represents all measurementmodels of the form
      \f[ x_k = f(x_{k-}) + N(\mu, \Sigma) \f]

  */
  class NonLinearAnalyticSystemModelGaussianUncertainty_Ginac
    : public AnalyticSystemModelGaussianUncertainty
    {
    public:
      /// Constructor
      /** @param pdf conditional pdf, gaussian uncertainty
       */
      NonLinearAnalyticSystemModelGaussianUncertainty_Ginac(NonLinearAnalyticConditionalGaussian_Ginac* const pdf);

      /// Destructor
      virtual ~NonLinearAnalyticSystemModelGaussianUncertainty_Ginac();

      /// output stream for system model
      // Not yet implemented
      /*
      friend std::ostream& operator<< (std::ostream& os,
				       NonLinearAnalyticSystemModelGaussianUncertainty_Ginac& m);
      */

      /// Get function
      GiNaC::matrix FunctionGet();

      /// Get State symbols
      vector<GiNaC::symbol> StateGet();

      /// Get input symbols
      vector<GiNaC::symbol> InputGet();

    };

} // End namespace BFL

#endif // __NON_LINEAR_SYSTEM_MODEL_GAUSSIAN_UNCERTAINTY_GINAC__
