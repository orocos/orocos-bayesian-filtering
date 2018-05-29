// $Id$
// Copyright (C) 2003 Klaas Gadeyne <first dot last at gmail dot com>
//                    Wim Meeussen  <wim dot meeussen at mech dot kuleuven dot ac dot be>
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
#include "iteratedextendedkalmanfilter.h"

namespace BFL
{
  using namespace MatrixWrapper;


#define AnalyticSys    AnalyticSystemModelGaussianUncertainty
#define AnalyticMeas   AnalyticMeasurementModelGaussianUncertainty


  IteratedExtendedKalmanFilter::IteratedExtendedKalmanFilter(Gaussian* prior, unsigned int nr_it, InnovationCheck* I)
    : KalmanFilter(prior)
    , _nr_iterations(nr_it)
    , _innovationChecker(I)
    , _x(prior->DimensionGet())
    , _x_i(prior->DimensionGet())
    , _x_i_prev(prior->DimensionGet())
    , _J(prior->DimensionGet())
    , _innovation(prior->DimensionGet())
    , _F(prior->DimensionGet(),prior->DimensionGet())
    , _Q(prior->DimensionGet())
    , _P_Matrix(prior->DimensionGet())
  {}

  IteratedExtendedKalmanFilter::~IteratedExtendedKalmanFilter()
  {}

  void
  IteratedExtendedKalmanFilter::AllocateMeasModelIExt(const vector<unsigned int>& meas_dimensions)
  {
    unsigned int meas_dimension;
    for(int i = 0 ; i< meas_dimensions.size(); i++)
    {
        // find if variables with size meas_sizes[i] are already allocated
        meas_dimension = meas_dimensions[i];
        _mapMeasUpdateVariablesIExt_it =  _mapMeasUpdateVariablesIExt.find(meas_dimension);
        if( _mapMeasUpdateVariablesIExt_it == _mapMeasUpdateVariablesIExt.end())
        {
            //variables with size z.rows() not allocated yet
            _mapMeasUpdateVariablesIExt_it = (_mapMeasUpdateVariablesIExt.insert
                (std::pair<unsigned int, MeasUpdateVariablesIExt>( meas_dimension,MeasUpdateVariablesIExt(meas_dimension,_x.rows()) ))).first;
         }
     }
  }

  void
  IteratedExtendedKalmanFilter::AllocateMeasModelIExt(const unsigned int& meas_dimension)
  {
     // find if variables with size meas_sizes[i] are already allocated
     _mapMeasUpdateVariablesIExt_it =  _mapMeasUpdateVariablesIExt.find(meas_dimension);
     if( _mapMeasUpdateVariablesIExt_it == _mapMeasUpdateVariablesIExt.end())
     {
         //variables with size z.rows() not allocated yet
         _mapMeasUpdateVariablesIExt_it = (_mapMeasUpdateVariablesIExt.insert
             (std::pair<unsigned int, MeasUpdateVariablesIExt>( meas_dimension,MeasUpdateVariablesIExt(meas_dimension,_x.rows()) ))).first;
      }
  }

  void
  IteratedExtendedKalmanFilter::SysUpdate(SystemModel<ColumnVector>* const sysmodel,
					  const ColumnVector& u)
  {
    _x = _post->ExpectedValueGet();
    _J = ((AnalyticSys*)sysmodel)->PredictionGet(u,_x);
    _F = ((AnalyticSys*)sysmodel)->df_dxGet(u,_x);
    _Q = ((AnalyticSys*)sysmodel)->CovarianceGet(u,_x);

    /*
      cout << "x :\n" << x << endl;
      cout << "J:\n " << J << endl;
      cout << "F:\n " << F << endl;
      cout << "Q:\n " << Q << endl;
    */

    CalculateSysUpdate(_J, _F, _Q);
  }

  void
  IteratedExtendedKalmanFilter::MeasUpdate(MeasurementModel<ColumnVector,ColumnVector>* const measmodel,
					   const ColumnVector& z,
					   const ColumnVector& s)
  {
    // allocate measurement for z.rows() if needed
    AllocateMeasModelIExt(z.rows());

    _x = _post->ExpectedValueGet();
    _P_Matrix = _post->CovarianceGet();
    _x_i = _post->ExpectedValueGet();

    bool            test_innovation = true;

    for (unsigned int i=0; i<_nr_iterations && test_innovation; i++)
      {
    _x_i_prev = _x_i;
	(_mapMeasUpdateVariablesIExt_it->second)._H_i  = ((AnalyticMeas*)measmodel)->df_dxGet(s,_x_i);
	(_mapMeasUpdateVariablesIExt_it->second)._R_i  = ((AnalyticMeas*)measmodel)->CovarianceGet(s,_x_i);
	_S_i  = ( (_mapMeasUpdateVariablesIExt_it->second)._H_i * (Matrix)_P_Matrix * ((_mapMeasUpdateVariablesIExt_it->second)._H_i.transpose()) ) + (Matrix)((_mapMeasUpdateVariablesIExt_it->second)._R_i);
	(_mapMeasUpdateVariablesIExt_it->second)._K_i  = (Matrix)_P_Matrix * ((_mapMeasUpdateVariablesIExt_it->second)._H_i.transpose()) * (_S_i.inverse());
	(_mapMeasUpdateVariablesIExt_it->second)._Z_i  = ((AnalyticMeas*)measmodel)->PredictionGet(s,_x_i) + ( (_mapMeasUpdateVariablesIExt_it->second)._H_i * (_x - _x_i) );
	_x_i  = _x + (_mapMeasUpdateVariablesIExt_it->second)._K_i * (z - (_mapMeasUpdateVariablesIExt_it->second)._Z_i);
    _innovation = (_x_i - _x_i_prev);
    if (_innovationChecker != NULL)
        test_innovation = _innovationChecker->check(_innovation) ; //test if the innovation is not too small

	/*
	  cout << "H_i :\n" << H_i << endl;
	  cout << "R_i:\n " << R_i << endl;
	  cout << "S_i:\n " << S_i << endl;
	  cout << "S_i inverse:\n " << S_i.inverse() << endl;
	  cout << "innovation:\n" << z- ( ((NLinMeas*)measmodel)->PredictionGet(s,x_i) )<< endl;
	  cout << "x_i:\n" << x_i << endl;
	*/
      }

    CalculateMeasUpdate(z, (_mapMeasUpdateVariablesIExt_it->second)._Z_i, (_mapMeasUpdateVariablesIExt_it->second)._H_i, (_mapMeasUpdateVariablesIExt_it->second)._R_i);
  }

}
