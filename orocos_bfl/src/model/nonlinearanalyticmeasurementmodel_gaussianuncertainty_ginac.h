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
#ifndef __NON_LINEAR_MEASUREMENT_MODEL_GAUSSIAN_UNCERTAINTY_GINAC__
#define __NON_LINEAR_MEASUREMENT_MODEL_GAUSSIAN_UNCERTAINTY_GINAC__

#include "analyticmeasurementmodel_gaussianuncertainty.h"
#include "../pdf/gaussian.h"
#include "../pdf/nonlinearanalyticconditionalgaussian_ginac.h"
#include <ginac/ginac.h>
#include <vector>
#include <iostream>

namespace BFL
{

  using namespace std;

  /// Class for nonlinear analytic measurementmodels with additive gaussian noise
  /** This class represents all measurementmodels of the form
      \f[ h(x)=z \ or \ h(x,z)=0 \f]

  */
  class NonLinearAnalyticMeasurementModelGaussianUncertainty_Ginac :
    public AnalyticMeasurementModelGaussianUncertainty
    {

    public:
      /// Constructor
      /** @param pdf conditional pdf, gaussian uncertainty
       */
      NonLinearAnalyticMeasurementModelGaussianUncertainty_Ginac(NonLinearAnalyticConditionalGaussian_Ginac* const pdf);

      /// copy constructor
      //  NonLinearAnalyticMeasurementModelGaussianUncertainty_Ginac(const NonLinearAnalyticMeasurementModelGaussianUncertainty_Ginac& model);

      /// Destructor
      virtual ~NonLinearAnalyticMeasurementModelGaussianUncertainty_Ginac();

      /// output stream for measurement model
      // Not yet implemented
      /*
	 friend std::ostream& operator<< (std::ostream& os,
	 NonLinearAnalyticMeasurementModelGaussianUncertainty_Ginac&
      //m);
      */

      // redefinition of virtual functions
      virtual MatrixWrapper::Matrix df_dxGet (const MatrixWrapper::ColumnVector& u,
					      const MatrixWrapper::ColumnVector& x);
      virtual MatrixWrapper::ColumnVector PredictionGet(const MatrixWrapper::ColumnVector& u,
							const MatrixWrapper::ColumnVector& x);
      virtual MatrixWrapper::SymmetricMatrix CovarianceGet(const MatrixWrapper::ColumnVector& u,
					    const MatrixWrapper::ColumnVector& x);

      /// Get function
      GiNaC::matrix FunctionGet();

      /// Get State symbols
      vector<GiNaC::symbol> StateGet();

      /// Get input symbols
      vector<GiNaC::symbol> InputGet();

      /// Get conditional arguments
      vector<GiNaC::symbol> ConditionalGet();

    };

} // End namespace BFL

#endif // __NON_LINEAR_MEASUREMENT_MODEL_GAUSSIAN_UNCERTAINTY_GINAC__
