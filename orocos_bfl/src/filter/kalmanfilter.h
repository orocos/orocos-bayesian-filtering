// $Id$
// Copyright (C) 2002 Klaas Gadeyne <first dot last at gmail dot com>
//                    Wim Meeussen  <wim dot meeussen at mech dot kuleuven dot ac dot be>
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
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
//

#ifndef __KALMAN_FILTER__
#define __KALMAN_FILTER__

#include "filter.h"
#include "../pdf/gaussian.h"
#include "../model/analyticmeasurementmodel_gaussianuncertainty.h"
#include "../model/analyticsystemmodel_gaussianuncertainty.h"
# include <map>

namespace BFL
{

/// Class representing the family of all Kalman Filters (EKF, IEKF, ...)
/** This is a class representing the family of all Kalman
    Filter (KF).  Kalman filters are filters in which the Posterior
    density is represented by a Gaussian density.  Kalman filters are
    only applicable to continuous systems.

    The system of updating the Posterior density is implemented in this
    base class. However, the parameters used for this update differ for
    different KFs (Simple KF,EKF,IEKF): that's why the xUpdate members
    are still pure virtual functions.

    This class is the base class for all sorts of KFs.

    @see Gaussian
    @see LinearAnalyticSystemModelGaussianUncertainty
*/
class KalmanFilter : public Filter<MatrixWrapper::ColumnVector,MatrixWrapper::ColumnVector>
{
public:
  /// Constructor
  /** @pre you created the prior
      @param prior pointer to the Gaussian Pdf prior density
  */
  KalmanFilter(Gaussian* prior);

  /// Destructor
  virtual ~KalmanFilter();

  // implement virtual function
  virtual Gaussian* PostGet();

  /// Function to allocate memory needed during the measurement update,
  //  For realtime use, this function should be called before calling measUpdate
  /*  @param vector containing the dimension of the measurement models which are
      going to be used
  */
  void AllocateMeasModel( const vector<unsigned int>& meas_dimensions);

  /// Function to allocate memory needed during the measurement update
  //  For realtime use, this function should be called before calling measUpdate
  /*  @param dimension of the measurement models which is
      going to be used
  */
  void AllocateMeasModel( const unsigned int& meas_dimensions);

private:
  struct MeasUpdateVariables
  {
    Matrix _S_Matrix;
    Matrix _K;
    ColumnVector _innov;
    Matrix _postHT;
    MeasUpdateVariables() {};
    MeasUpdateVariables(unsigned int meas_dimension, unsigned int state_dimension):
      _S_Matrix(meas_dimension,meas_dimension)
    , _K(state_dimension,meas_dimension)
    , _innov(meas_dimension)
    , _postHT(state_dimension,meas_dimension)
{};
  }; //struct

protected:
  // variables to avoid allocation during update calls
  ColumnVector  _Mu_new;
  SymmetricMatrix _Sigma_new;
  Matrix _Sigma_temp;
  Matrix _Sigma_temp_par;
  Matrix _SMatrix;
  Matrix _K;
  double _Nis;
  std::map<unsigned int, MeasUpdateVariables> _mapMeasUpdateVariables;
  std::map<unsigned int, MeasUpdateVariables>::iterator _mapMeasUpdateVariables_it;


  /** Very dirty hack to avoid ugly methods PostSigmaSet
      and PostMuSet to be public!
      NonMinimalKalmanFilter should be redesigned though!
  */
  friend class NonminimalKalmanFilter;

  /// Set covariance of posterior estimate
  void PostSigmaSet( const MatrixWrapper::SymmetricMatrix& s);

  /// Set expected value of posterior estimate
  void PostMuSet( const MatrixWrapper::ColumnVector& c);

  /** Calculate Kalman filter System Update
      \f[ x_k = J \f]
      \f[ P_k = F.P_{k-}.F' + Q \f]
  */
  void CalculateSysUpdate(const MatrixWrapper::ColumnVector& J, const MatrixWrapper::Matrix& F, const MatrixWrapper::SymmetricMatrix& Q);

  /** Calculate Kalman filter Measurement Update
      \f[ x_k = x_{k-} + K.(z - Z) \f]
      \f[ P_k = (I-K.H).P_{k-} \f]
      with
      \f[ K = P_{k-}.H'.(H.P_{k-}.H'+R)^{-1} \f]
  */
  void CalculateMeasUpdate(const MatrixWrapper::ColumnVector& z, const MatrixWrapper::ColumnVector& Z, const MatrixWrapper::Matrix& H, const MatrixWrapper::SymmetricMatrix& R);

  /// System Update
  /** Update the filter's Posterior density using the deterministic
      inputs to the system and the system model
      @param sysmodel pointer to the system model the filter should use
      @param u input to the system
  */
  virtual void SysUpdate(SystemModel<MatrixWrapper::ColumnVector>* const sysmodel,
			 const MatrixWrapper::ColumnVector& u) = 0;

  /// Measurement Update (overloaded)
  /** Update the filter's Posterior density using the sensor
      measurements, an input and the measurement model.  This method is
      used when the measurements depend on the inputs too (doesn't
      happen very often, does it?)
      BEWARE: the first time the measurment update is called with a new size of measurement, new allocations are done
      @param measmodel pointer to the measurement model the filter
      should use
      @param z sensor measurement
      @param s input to the system (must be of the same type as u
      for now, since this was not yet implemented in ConditionalPdf
  */
  virtual void MeasUpdate(MeasurementModel<MatrixWrapper::ColumnVector,MatrixWrapper::ColumnVector>* const measmodel,
			  const MatrixWrapper::ColumnVector& z,
			  const MatrixWrapper::ColumnVector& s) = 0;

  virtual bool UpdateInternal(SystemModel<MatrixWrapper::ColumnVector>* const sysmodel,
			      const MatrixWrapper::ColumnVector& u,
			      MeasurementModel<MatrixWrapper::ColumnVector,MatrixWrapper::ColumnVector>* const measmodel,
			      const MatrixWrapper::ColumnVector& z,
			      const MatrixWrapper::ColumnVector& s);

  // Calculate Normalised Innovation Squared (NIS) value
  void CalculateNis(const MatrixWrapper::ColumnVector& z, const MatrixWrapper::ColumnVector& Z, const MatrixWrapper::Matrix& H, const MatrixWrapper::SymmetricMatrix& R);

}; // class



} // End namespace BFL

#endif // __KALMAN_FILTER__
