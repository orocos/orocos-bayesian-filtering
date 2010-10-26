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
#include "nonlinearanalyticsystemmodel_gaussianuncertainty_ginac.h"
#include <cassert>

namespace BFL
{

  using namespace std;

  NonLinearAnalyticSystemModelGaussianUncertainty_Ginac::NonLinearAnalyticSystemModelGaussianUncertainty_Ginac
  (NonLinearAnalyticConditionalGaussian_Ginac* const pdf)
    : AnalyticSystemModelGaussianUncertainty( new NonLinearAnalyticConditionalGaussian_Ginac( *pdf ) )
  {}


  NonLinearAnalyticSystemModelGaussianUncertainty_Ginac::~NonLinearAnalyticSystemModelGaussianUncertainty_Ginac()
  {}

  /*
  std::ostream& operator<< (std::ostream& os, NonLinearAnalyticSystemModelGaussianUncertainty_Ginac& m)
  {
    os << "\nSystemModel:"  << endl;
    os << *(SystemPdfGet());
    return os;
  }
  */

  GiNaC::matrix
  NonLinearAnalyticSystemModelGaussianUncertainty_Ginac::FunctionGet()
  {
    return ((NonLinearAnalyticSystemModelGaussianUncertainty_Ginac *) SystemPdfGet())->FunctionGet();
  }

  vector<GiNaC::symbol>
  NonLinearAnalyticSystemModelGaussianUncertainty_Ginac::StateGet()
  {
    return ((NonLinearAnalyticSystemModelGaussianUncertainty_Ginac *)SystemPdfGet())->StateGet();
  }

  vector<GiNaC::symbol>
  NonLinearAnalyticSystemModelGaussianUncertainty_Ginac::InputGet()
  {
    return ((NonLinearAnalyticSystemModelGaussianUncertainty_Ginac *) SystemPdfGet())->InputGet();
  }


} // End namespace BFL
