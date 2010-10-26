// $Id$
// Copyright (C) 2003 Klaas Gadeyne <first dot last at gmail dot com>
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

#ifndef __NONLINEAR_SYSTEM_CONDITIONAL_GAUSSIAN_GINAC__
#define __NONLINEAR_SYSTEM_CONDITIONAL_GAUSSIAN_GINAC__

#include "analyticconditionalgaussian_additivenoise.h"
#include <ginac/ginac.h>
#include <iostream>

namespace BFL
{
  /// Conditional Gaussian for an analytic nonlinear system using Ginac:
  /**
     Describes classes of the type \f[ P(z | subs) \f]
     with
     \f[ z=f(subs) + N(\mu,\Sigma) \f]
     or
     \f[ z=f(subs,c+N(\mu, \Sigma)) \f]

     Constructor for the first type:
     \f[ NonLinearAnalyticConditionalGaussian_Ginac(f(subs), subs, N(\mu, \Sigma) ) \f]

     Constructor for the second type:
     \f[ NonLinearAnalyticConditionalGaussian_Ginac(f(subs,z), subs, N(\mu, \Sigma) ,c) \f]

     When the second type is used, the additive noise on c will be
     converted to additive noise on f, by locally linearising the
     function.
     @bug:  This class is higly biased towards filtering applications.
  */
  class NonLinearAnalyticConditionalGaussian_Ginac : public AnalyticConditionalGaussianAdditiveNoise
    {
    public:
      /// constructor
      /**
	 @param func function to be evaluated for expected value
	 @param u symbols to be substituted (by numeric values) for
	 evaluation.  These can be system inputs or sensor parameters
	 @param x symbols representing state
	 @param additiveNoise Gaussian representing additive noise
	 @param cond parameters where additive noise applies to
      */
      NonLinearAnalyticConditionalGaussian_Ginac( const GiNaC::matrix& func,
						  const vector<GiNaC::symbol>& u,
						  const vector<GiNaC::symbol>& x,
						  const Gaussian& additiveNoise,
						  const vector<GiNaC::symbol>& cond );

      /// constructor
      /**
	 @param func function to be evaluated for expected value
	 @param u symbols to be substituted (by numeric values) for
	 evaluation. These can be system inputs or sensor parameters
	 @param x symbols representing state
	 @param additiveNoise Gaussian representing additive noise
	 on function output
      */
      NonLinearAnalyticConditionalGaussian_Ginac( const GiNaC::matrix& func,
						  const vector<GiNaC::symbol>& u,
						  const vector<GiNaC::symbol>& x,
						  const Gaussian& additiveNoise );
      /// copy constructor
      NonLinearAnalyticConditionalGaussian_Ginac( const NonLinearAnalyticConditionalGaussian_Ginac& g);

      /// Destructor
      virtual ~NonLinearAnalyticConditionalGaussian_Ginac();

      /// output stream for measurement model
      friend std::ostream& operator<< (std::ostream& os, NonLinearAnalyticConditionalGaussian_Ginac& p);

      /// return function
      GiNaC::matrix FunctionGet();

      /// return substitution symbols
      vector<GiNaC::symbol> InputGet();

      /// return state symbols
      vector<GiNaC::symbol> StateGet();

      /// Get conditional arguments
      vector<GiNaC::symbol> ConditionalGet();

      // redefinition of virtual functions
      virtual MatrixWrapper::ColumnVector    ExpectedValueGet() const;
      virtual MatrixWrapper::SymmetricMatrix CovarianceGet()    const;

      // Redefinition of dfGet
      /**
       * @bug only implemented for i = 0 for now (so in a filter
       * context, only the derivative with respect to x is implemented
       **/
      virtual MatrixWrapper::Matrix dfGet(unsigned int i)       const;


    private:
      GiNaC::matrix 		  func_sym;
      vector<GiNaC::symbol> 	  cond_sym, u_sym, x_sym;
      unsigned int			  cond_size, u_size, x_size, func_size;
      vector<GiNaC::ex> 		  dfunc_dcond, dfunc_dx;


    };

} // End namespace

#endif //  __NONLINEAR_SYSTEM_CONDITIONAL_GAUSSIAN_GINAC__
