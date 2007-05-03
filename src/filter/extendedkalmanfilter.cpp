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
#include "extendedkalmanfilter.h"

namespace BFL
{
  using namespace MatrixWrapper;
  


#define AnalyticSys    AnalyticSystemModelGaussianUncertainty
#define AnalyticMeas   AnalyticMeasurementModelGaussianUncertainty


ExtendedKalmanFilter::ExtendedKalmanFilter(Gaussian* prior)
: KalmanFilter(prior)
{}

ExtendedKalmanFilter::~ExtendedKalmanFilter()
{}

void
ExtendedKalmanFilter::SysUpdate(SystemModel<ColumnVector>* const sysmodel,
                                const ColumnVector& u)
{
  ColumnVector    x = _post->ExpectedValueGet();
  ColumnVector    J = ((AnalyticSys*)sysmodel)->PredictionGet(u,x); 
  Matrix          F = ((AnalyticSys*)sysmodel)->df_dxGet(u,x);
  SymmetricMatrix Q = ((AnalyticSys*)sysmodel)->CovarianceGet(u,x);
  
  CalculateSysUpdate(J, F, Q);
}

void
ExtendedKalmanFilter::MeasUpdate(MeasurementModel<ColumnVector,ColumnVector>* const measmodel,
                                 const ColumnVector& z, 
			         const ColumnVector& s)
{
  ColumnVector    x = _post->ExpectedValueGet();
  ColumnVector    Z = ((AnalyticMeas*)measmodel)->PredictionGet(s,x); 
  Matrix          H = ((AnalyticMeas*)measmodel)->df_dxGet(s,x);
  SymmetricMatrix R = ((AnalyticMeas*)measmodel)->CovarianceGet(s,x);
  
  CalculateMeasUpdate(z, Z, H, R);
}

} // end namespace BFL
