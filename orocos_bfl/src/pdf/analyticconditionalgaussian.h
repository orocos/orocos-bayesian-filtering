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

#ifndef __ANALYTIC_CONDITIONAL_GAUSSIAN__
#define __ANALYTIC_CONDITIONAL_GAUSSIAN__

#include "conditionalgaussian.h"

namespace BFL
{

  /// Abstract Class representing all _FULL_ Analytical Conditional gaussians
  /** So this class represents all Pdf's of the type
      \f[ P ( A | B, C, D, ... ) \f] where
      \f[ \mu_A = f(B,C,D, ...) \f] and
      \f[ \Sigma_A = g(B,C,D, ...) \f] and

      \f[ A = N(\mu_A, \Sigma_A) \f]

  */
  class AnalyticConditionalGaussian : public ConditionalGaussian
    {
    public:
      /// Constructor
      /**
	 @param dim Dimension of state
	 @param num_conditional_arguments The number of conditional
	 arguments.
      */
      AnalyticConditionalGaussian(int dim = 0, int num_conditional_arguments=0);

      // Default Copy constructor will do

      /// Destructor
      virtual ~AnalyticConditionalGaussian();

      /// returns derivative from function to n-th conditional variable
      /** @param i Number of the conditional variable to use for
	  partial derivation
	  @return Partial derivative with respect to conditional
	  variable i
      */
      virtual MatrixWrapper::Matrix dfGet(unsigned int i) const;

    };

} // End namespace BFL

#endif // __ANALYTIC_CONDITIONAL_GAUSSIAN__

