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

#include "../sample/sample.h"
#include "nonlinearanalyticmeasurementmodel_gaussianuncertainty_ginac.h"

namespace BFL
{

  NonLinearAnalyticMeasurementModelGaussianUncertainty_Ginac::NonLinearAnalyticMeasurementModelGaussianUncertainty_Ginac
  (NonLinearAnalyticConditionalGaussian_Ginac* const pdf)
    : AnalyticMeasurementModelGaussianUncertainty( new NonLinearAnalyticConditionalGaussian_Ginac( *pdf ) )
  {}


  // copy constructor: Not implemented yet
  /*
  NonLinearAnalyticMeasurementModelGaussianUncertainty_Ginac::NonLinearAnalyticMeasurementModelGaussianUncertainty_Ginac
  (const NonLinearAnalyticMeasurementModelGaussianUncertainty_Ginac& m)
    : AnalyticMeasurementModelGaussianUncertainty( new NonLinearAnalyticConditionalGaussian_Ginac( (NonLinearAnalyticConditionalGaussian_Ginac*) m.MeasurementPdfGet()) )
  {}
  */

  NonLinearAnalyticMeasurementModelGaussianUncertainty_Ginac::~NonLinearAnalyticMeasurementModelGaussianUncertainty_Ginac() {}

  /*
  std::ostream& operator<< (std::ostream& os, NonLinearAnalyticMeasurementModelGaussianUncertainty_Ginac& m)
  {
    os << "\nMeasurementModel:"  << endl;
    os << *(m.MeasurementPdfGet());
    return os;
  }
  */


  MatrixWrapper::ColumnVector
  NonLinearAnalyticMeasurementModelGaussianUncertainty_Ginac::PredictionGet(const MatrixWrapper::ColumnVector& u,
								      const MatrixWrapper::ColumnVector& x)
  {
    MeasurementPdfGet()->ConditionalArgumentSet(1,u);
    MeasurementPdfGet()->ConditionalArgumentSet(0,x);
    return MeasurementPdfGet()->ExpectedValueGet();
  }


  MatrixWrapper::SymmetricMatrix
  NonLinearAnalyticMeasurementModelGaussianUncertainty_Ginac::CovarianceGet(const MatrixWrapper::ColumnVector& u,
								      const MatrixWrapper::ColumnVector& x)
  {
    MeasurementPdfGet()->ConditionalArgumentSet(1,u);
    MeasurementPdfGet()->ConditionalArgumentSet(0,x);
    return MeasurementPdfGet()->CovarianceGet();
  }


  MatrixWrapper::Matrix
  NonLinearAnalyticMeasurementModelGaussianUncertainty_Ginac::df_dxGet(const MatrixWrapper::ColumnVector& u,
								       const MatrixWrapper::ColumnVector& x)
  {
    MeasurementPdfGet()->ConditionalArgumentSet(1,u);
    MeasurementPdfGet()->ConditionalArgumentSet(0,x);
    return ((AnalyticConditionalGaussian *) MeasurementPdfGet())->dfGet(0);
  }

  GiNaC::matrix
  NonLinearAnalyticMeasurementModelGaussianUncertainty_Ginac::FunctionGet()
  {
    return ((NonLinearAnalyticMeasurementModelGaussianUncertainty_Ginac *) MeasurementPdfGet())->FunctionGet();
  }

  vector<GiNaC::symbol>
  NonLinearAnalyticMeasurementModelGaussianUncertainty_Ginac::StateGet()
  {
    return ((NonLinearAnalyticMeasurementModelGaussianUncertainty_Ginac *)MeasurementPdfGet())->StateGet();
  }

  vector<GiNaC::symbol>
  NonLinearAnalyticMeasurementModelGaussianUncertainty_Ginac::InputGet()
  {
    return ((NonLinearAnalyticMeasurementModelGaussianUncertainty_Ginac *) MeasurementPdfGet())->InputGet();
  }

  vector<GiNaC::symbol>
  NonLinearAnalyticMeasurementModelGaussianUncertainty_Ginac::ConditionalGet()
  {
    return ((NonLinearAnalyticMeasurementModelGaussianUncertainty_Ginac *) MeasurementPdfGet())->ConditionalGet();
  }

} // End namespace BFL
