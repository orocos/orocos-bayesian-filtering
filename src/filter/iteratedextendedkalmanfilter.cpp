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
// Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
//  
#include "iteratedextendedkalmanfilter.h"

namespace BFL
{
  using namespace MatrixWrapper;
  

#define AnalyticSys    AnalyticSystemModelGaussianUncertainty
#define AnalyticMeas   AnalyticMeasurementModelGaussianUncertainty


  IteratedExtendedKalmanFilter::IteratedExtendedKalmanFilter(Gaussian* prior, unsigned int nr_it, InnovationCheck* I)
    : KalmanFilter(prior),
      nr_iterations(nr_it),
      innovationChecker(I)
  {}

  IteratedExtendedKalmanFilter::~IteratedExtendedKalmanFilter()
  {}

  void
  IteratedExtendedKalmanFilter::SysUpdate(SystemModel<ColumnVector>* const sysmodel,
					  const ColumnVector& u)
  {
    ColumnVector    x = _post->ExpectedValueGet();
    ColumnVector    J = ((AnalyticSys*)sysmodel)->PredictionGet(u,x); 
    Matrix          F = ((AnalyticSys*)sysmodel)->df_dxGet(u,x);
    SymmetricMatrix Q = ((AnalyticSys*)sysmodel)->CovarianceGet(u,x);
  
    /*
      cout << "x :\n" << x << endl;
      cout << "J:\n " << J << endl;
      cout << "F:\n " << F << endl;
      cout << "Q:\n " << Q << endl;
    */
  
    CalculateSysUpdate(J, F, Q);
  }

  void
  IteratedExtendedKalmanFilter::MeasUpdate(MeasurementModel<ColumnVector,ColumnVector>* const measmodel,
					   const ColumnVector& z, 
					   const ColumnVector& s)
  {

    ColumnVector    x_k = _post->ExpectedValueGet();
    SymmetricMatrix P_k = _post->CovarianceGet();
    ColumnVector    x_i = _post->ExpectedValueGet();

    ColumnVector    x_i_prev;
    Matrix          H_i;
    SymmetricMatrix R_i;
    Matrix          S_i;
    Matrix          K_i;
    ColumnVector    Z_i;
    ColumnVector    innovation;
    bool            test_innovation = true;
  
    for (unsigned int i=0; i<nr_iterations && test_innovation; i++)
      {
    x_i_prev = x_i;
	H_i  = ((AnalyticMeas*)measmodel)->df_dxGet(s,x_i);
	R_i  = ((AnalyticMeas*)measmodel)->CovarianceGet(s,x_i);
	S_i  = ( H_i * (Matrix)P_k * (H_i.transpose()) ) + (Matrix)R_i;
	K_i  = (Matrix)P_k * (H_i.transpose()) * (S_i.inverse());
	Z_i  = ((AnalyticMeas*)measmodel)->PredictionGet(s,x_i) + ( H_i * (x_k - x_i) );  
	x_i  = x_k + K_i * (z - Z_i);
    innovation = (x_i - x_i_prev);
    if (innovationChecker != NULL)
        test_innovation = innovationChecker->check(innovation) ; //test if the innovation is not too small

    
	/*
	  cout << "H_i :\n" << H_i << endl;
	  cout << "R_i:\n " << R_i << endl;
	  cout << "S_i:\n " << S_i << endl;
	  cout << "S_i inverse:\n " << S_i.inverse() << endl;
	  cout << "innovation:\n" << z- ( ((NLinMeas*)measmodel)->PredictionGet(s,x_i) )<< endl;
	  cout << "x_i:\n" << x_i << endl;
	*/
      }
  
    CalculateMeasUpdate(z, Z_i, H_i, R_i);
  }

}
