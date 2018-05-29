// $Id$
// Copyright (C) 2003 Klaas Gadeyne <first dot last at gmail dot com>
//                    Wim Meeussen  <wim dot meeussen at mech dot kuleuven dot be>
//                    Tinne De Laet  <tinne dot delaet at mech dot kuleuven dot be>
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

#ifndef __LINEAR_SYSTEM_CONDITIONAL_GAUSSIAN__
#define __LINEAR_SYSTEM_CONDITIONAL_GAUSSIAN__

#include "analyticconditionalgaussian_additivenoise.h"

namespace BFL
{
  /// Linear Conditional Gaussian
  /**
     - \f$ \mu = Matrix[1] . ConditionalArguments[0] +
     Matrix[2]. ConditionalArguments[1]  + ... + Noise.\mu \f$
     - Covariance is independent of the ConditionalArguments, and is
     the covariance of the Noise pdf
  */
  class LinearAnalyticConditionalGaussian : public AnalyticConditionalGaussianAdditiveNoise
    {
    public:
      /// Constructor
      /** @pre:  Every Matrix should have the same amount of rows!
	  This is currently not checked.  The same goes for the number
	  of columns, which should be equal to the number of rows of
	  the corresponding conditional argument!
	  @param ratio: vector containing the different matrices of
	  the linear relationship between the conditional arguments
	  and \f$\mu\f$
	  @param additiveNoise Pdf representing the additive Gaussian uncertainty
      */
      LinearAnalyticConditionalGaussian(const vector<MatrixWrapper::Matrix> & ratio,
				const Gaussian& additiveNoise);

      /// Constructor (overloaded)
      /** @pre There is only 1 conditional argument.
	  @param a Matrix for calculation of \f$\mu\f$:
	  \f$ \mu = a . ConditionalArguments[0] + Noise.\mu \f$
	  @param additiveNoise Pdf representing the additive Gaussian uncertainty
      */
      LinearAnalyticConditionalGaussian(const MatrixWrapper::Matrix& a, const Gaussian& additiveNoise);

      // Default copy constructor will do

      /// Destructor
      virtual ~LinearAnalyticConditionalGaussian();

      ///Clone function
      virtual LinearAnalyticConditionalGaussian* Clone() const;

      // implement virtual functions
      virtual MatrixWrapper::ColumnVector    ExpectedValueGet() const;
      virtual MatrixWrapper::Matrix          dfGet(unsigned int i)       const;

      /// Be careful: you don't want to use this one: Redefined.
      /** @bug This method is not implemented, we can ReSize the
	       std::vector<BFL::Matrix>, but we don't know the
	       dimensions of the matrices self.  So this will most
	       certainly result in a segfault.  Anyway, why would you
	       need this?
      */
      virtual void NumConditionalArgumentsSet(unsigned int numconditionalarguments);

      /// Set the i-th Matrix for calculation of \f$ \mu \f$
      /** Set the i-th Matrix of the \f$ \mu \f$ calculation in the
	  conditonal gaussian class
	  @pre i < Numconditionalarg
	  @param i index determining which conditional Arg. will be
	  multiplied with the given matrix
	  @param m Matrix for calculation of \f$ \mu \f$:
	  \f$ \mu = ... m . ConditionalArguments[i] + ... \f$
      */
      void MatrixSet(unsigned int i, const MatrixWrapper::Matrix& m);

      /// Get the i-th matrix of the system
      /**
	 @param i index determining which conditional Arg. multiplier
	 matrix will returned
	 @return the n-th Matrix of the system-equation
      */
      const MatrixWrapper::Matrix& MatrixGet(unsigned int i) const;

    private:

      vector<MatrixWrapper::Matrix> _ratio;
      // variables to avoid allocation during expectedValueGet call
      mutable MatrixWrapper::ColumnVector _mean_temp;
      mutable MatrixWrapper::ColumnVector _arg;

    };

} // End namespace BFL

#endif //  __LINEAR_SYSTEM_CONDITIONAL_GAUSSIAN__
