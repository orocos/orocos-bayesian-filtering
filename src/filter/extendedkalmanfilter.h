// $Id$
// Copyright (C) 2002 Klaas Gadeyne <first dot last at gmail dot com>
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

#ifndef __EXTENDED_KALMAN_FILTER__
#define __EXTENDED_KALMAN_FILTER__

#include "kalmanfilter.h"
#include "../pdf/conditionalpdf.h"
#include "../pdf/gaussian.h"
# include <map>

namespace BFL
{

/** This is a class implementing the Kalman Filter (KF) class for 
    Extended Kalman Filters.
    
    The System- and MeasurementUpdate equasions are not linear, and 
    will be approximated by local linearisations.

    @see KalmanFilter
    @note that if the system/measurement model that you pass to the
    update calls is not analytical with additive gaussian noise, you
    will get rubbish results
*/
class ExtendedKalmanFilter : public KalmanFilter
{
public:
  /** Constructor
      @pre you created the prior
      @param prior pointer to the Gaussian prior density
  */
  ExtendedKalmanFilter(Gaussian* prior);
   
  /// Destructor
  virtual ~ExtendedKalmanFilter();

  /// Function to allocate memory needed during the measurement update,
  //  For realtime use, this function should be called before calling measUpdate
  /*  @param vector containing the dimension of the measurement models which are
      going to be used
  */
  void AllocateMeasModelExt( const vector<unsigned int>& meas_dimensions);

  /// Function to allocate memory needed during the measurement update
  //  For realtime use, this function should be called before calling measUpdate
  /*  @param dimension of the measurement models which is
      going to be used
  */
  void AllocateMeasModelExt( const unsigned int& meas_dimensions);

private:
  struct MeasUpdateVariablesExt
  {
    SymmetricMatrix _R;
    Matrix _H;
    ColumnVector _Z;
    MeasUpdateVariablesExt() {};
    MeasUpdateVariablesExt(unsigned int meas_dimension, unsigned int state_dimension):
      _R(meas_dimension) 
    , _H(meas_dimension,state_dimension)
    , _Z(meas_dimension)
{};
  }; //struct

 protected:
  virtual void SysUpdate(SystemModel<MatrixWrapper::ColumnVector>* const sysmodel,
                         const MatrixWrapper::ColumnVector& u);
  virtual void MeasUpdate(MeasurementModel<MatrixWrapper::ColumnVector,MatrixWrapper::ColumnVector>* const measmodel,
                          const MatrixWrapper::ColumnVector& z, 
			  const MatrixWrapper::ColumnVector& s);
  // variables to avoid allocation on the heap
  ColumnVector _x;
  ColumnVector _J;
  Matrix    _F;
  SymmetricMatrix _Q;
  std::map<unsigned int, MeasUpdateVariablesExt> _mapMeasUpdateVariablesExt;
  std::map<unsigned int, MeasUpdateVariablesExt>::iterator _mapMeasUpdateVariablesExt_it;


};  // class

} // End namespace BFL

#endif // __EXTENDED_KALMAN_FILTER__
