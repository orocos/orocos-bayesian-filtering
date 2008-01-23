// $Id$
// Copyright (C) 2003 Klaas Gadeyne <first dot last at gmail dot com>
//                    Wim Meeussen  <wim dot meeussen at mech dot kuleuven dot be>
//                    Tinne De Laet <tinne dot delaet at mech dot kuleuven dot be>
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
, _x(prior->DimensionGet())
, _J(prior->DimensionGet())
, _F(prior->DimensionGet(),prior->DimensionGet())
, _Q(prior->DimensionGet())
{}

ExtendedKalmanFilter::~ExtendedKalmanFilter()
{}

void
ExtendedKalmanFilter::AllocateMeasModelExt(const vector<unsigned int>& meas_dimensions)
{
  unsigned int meas_dimension;
  for(int i = 0 ; i< meas_dimensions.size(); i++)
  {
      // find if variables with size meas_sizes[i] are already allocated
      meas_dimension = meas_dimensions[i];
      _mapMeasUpdateVariablesExt_it =  _mapMeasUpdateVariablesExt.find(meas_dimension);
      if( _mapMeasUpdateVariablesExt_it == _mapMeasUpdateVariablesExt.end())
      {
          //variables with size z.rows() not allocated yet
          _mapMeasUpdateVariablesExt_it = (_mapMeasUpdateVariablesExt.insert
              (std::pair<unsigned int, MeasUpdateVariablesExt>( meas_dimension,MeasUpdateVariablesExt(meas_dimension,_x.rows()) ))).first;
       }
   }
}

void
ExtendedKalmanFilter::AllocateMeasModelExt(const unsigned int& meas_dimension)
{
   // find if variables with size meas_sizes[i] are already allocated
   _mapMeasUpdateVariablesExt_it =  _mapMeasUpdateVariablesExt.find(meas_dimension);
   if( _mapMeasUpdateVariablesExt_it == _mapMeasUpdateVariablesExt.end())
   {
       //variables with size z.rows() not allocated yet
       _mapMeasUpdateVariablesExt_it = (_mapMeasUpdateVariablesExt.insert
           (std::pair<unsigned int, MeasUpdateVariablesExt>( meas_dimension,MeasUpdateVariablesExt(meas_dimension,_x.rows()) ))).first;
    }
}

void
ExtendedKalmanFilter::SysUpdate(SystemModel<ColumnVector>* const sysmodel,
                                const ColumnVector& u)
{
  _x = _post->ExpectedValueGet();
  _J = ((AnalyticSys*)sysmodel)->PredictionGet(u,_x); 
  _F = ((AnalyticSys*)sysmodel)->df_dxGet(u,_x);
  _Q = ((AnalyticSys*)sysmodel)->CovarianceGet(u,_x);
  
  CalculateSysUpdate(_J, _F, _Q);
}

void
ExtendedKalmanFilter::MeasUpdate(MeasurementModel<ColumnVector,ColumnVector>* const measmodel,
                                 const ColumnVector& z, 
			         const ColumnVector& s)
{
  // allocate measurement for z.rows() if needed
  AllocateMeasModelExt(z.rows());
  
  (_mapMeasUpdateVariablesExt_it->second)._Z = ((AnalyticMeas*)measmodel)->PredictionGet(s,_x); 
  (_mapMeasUpdateVariablesExt_it->second)._H = ((AnalyticMeas*)measmodel)->df_dxGet(s,_x);
  (_mapMeasUpdateVariablesExt_it->second)._R = ((AnalyticMeas*)measmodel)->CovarianceGet(s,_x);
  
  CalculateMeasUpdate(z, (_mapMeasUpdateVariablesExt_it->second)._Z, (_mapMeasUpdateVariablesExt_it->second)._H, (_mapMeasUpdateVariablesExt_it->second)._R);
}

} // end namespace BFL
