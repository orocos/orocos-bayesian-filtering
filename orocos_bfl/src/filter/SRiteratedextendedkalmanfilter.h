// $Id$
// Copyright (C) 2002 Klaas Gadeyne <first dot last at gmail dot com>
//                    Peter Slaets  <peter dot slaets at mech dot kuleuven dot ac dot be>
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

#ifndef __SR_ITERATED_EXTENDED_KALMAN_FILTER__
#define __SR_ITERATED_EXTENDED_KALMAN_FILTER__

#include "kalmanfilter.h"
#include "../pdf/conditionalpdf.h"
#include "../pdf/gaussian.h"

namespace BFL
{

  /** This is a class implementing the Kalman Filter (KF) class for
      Square Root Iterated Extended Kalman Filters.
      this is a possible implementation of a Kalman filter, which in
       will yield better numerical stable results, since the covariance matrix
       is defined as the Square root of the Covariance matrix of the state estimation.
       See
     @verbatim
     @Book{   Anderson_auxiliary,
      author    = {Anderson, B.D.O. and Moore, J.B.},
      title     = {Optimal filtering},
      publisher = {Prentice-Hall, Englewood Cliffs, NJ },
      year      = {1979}
    }
     @endverbatim
      for more details.
     Note that this particular implementation:
	        - Currently only works for implicit measurement model
		(LinearAnalyticMeasurementModelGaussianUncertainty_Implicit)
    		the implicit measurement model is described as :\f[  0 = f (x,z) \f]
    		and internally defined as a linear measurement model with virtual measurement z_k^{virtual}
    		\f[ z_k^{virtual} = H(x_k,z_k) \times x_k  + N(\mu(x_{k},z_k) ,\Sigma(x_k,z_k)) \f]
		@see KalmanFilter
  */
  class SRIteratedExtendedKalmanFilter : public KalmanFilter
    {
    public:
      /** Constructor
	  @pre you created the prior
	  @param prior pointer to the Monte Carlo Pdf prior density
	  @param nr_it the number of iterations in one update
      */
      SRIteratedExtendedKalmanFilter(Gaussian* prior, unsigned int nr_it=1);

      /// Destructor
      virtual ~SRIteratedExtendedKalmanFilter();
      /// Perform a system update with the current measurement model ans system model using an input u
      virtual void SysUpdate(SystemModel<MatrixWrapper::ColumnVector>* const sysmodel,
			     const MatrixWrapper::ColumnVector& u);
      /// Perform a system update with the current measurement model and system model
      virtual void SysUpdate(SystemModel<MatrixWrapper::ColumnVector>* const sysmodel);
      /// Perform a measurement update use a measurement model, measurements z and virutal measurements s
      virtual void MeasUpdate(MeasurementModel<MatrixWrapper::ColumnVector,MatrixWrapper::ColumnVector>* const measmodel,
			      const MatrixWrapper::ColumnVector& z,
			      const MatrixWrapper::ColumnVector& s);
       /// Perform a measurement update use a measurement model, measurements z
      virtual void MeasUpdate(MeasurementModel<MatrixWrapper::ColumnVector,MatrixWrapper::ColumnVector>* const measmodel,const MatrixWrapper::ColumnVector& z);
      /// Returns a square root of the covariance of the measurement input u
      MatrixWrapper::Matrix SRCovarianceGet() const;
      /// Set the square root covariance to a specific value
      void  	SRCovarianceSet( MatrixWrapper::Matrix JP_new);
      /// Set mean and covariance of the state estimation to a specific value
      void 	PriorSet(MatrixWrapper::ColumnVector& X, MatrixWrapper::SymmetricMatrix& P);

      /** Calculate Kalman filter Measurement Update
      \f[ x_k = x_{k-} + K.(z - Z) \f]
      \f[ S=(H_k*J_k*J_k^{'}*H_k+D*R*D^{'})^{'} \f]
      \f[ Sr=cholesky(S)^{'}\f]
      \f[ J_k = J_k-J_k*J_k^{'}*H^{'}*inv(S)^{'}*inv(R+Sr)*H_k*J_k  \f]
      with
      \f[ K = J_{k}*J_{k}^{'}*H'*inv(S)^{'}*inv(S) \f]
      */
      virtual void CalculateMeasUpdate(MatrixWrapper::ColumnVector z, MatrixWrapper::ColumnVector Z, MatrixWrapper::Matrix H, MatrixWrapper::SymmetricMatrix R);

      /// Calculate K_i , invS and  Sr
      /**
       \f[ Sr= cholesky(H_i*JP*JP^{'}*H_i^{'}+R_i)^{'}  \f]
       \f[ invS=inv(Sr)\f]
       \f[ K_i=JP*JP^{'}*H_i^{'}*invS^{'}*invS \f]
       */
      virtual void  CalculateMatrix(MatrixWrapper::Matrix& H_i  , MatrixWrapper::SymmetricMatrix& R_i , MatrixWrapper::Matrix& invS  , MatrixWrapper::Matrix& Sr  , MatrixWrapper::Matrix& K_i );

      /// Calculate the new state estimate
      /**
	 \f[ x_k = x_{k-} + K.(z - Z) \f]
      */
      virtual void CalculateMean(MatrixWrapper::ColumnVector& x_k, const MatrixWrapper::ColumnVector& z, MatrixWrapper::ColumnVector& Z_i ,MatrixWrapper::Matrix& K_i);

      /// Calculate the covariance of the new state estimate (P_k)
      /**
        \f[ JP = JP-JP*JP^{'}*H_k^{'}*inv(S)^{'}*inv(R+Sr)*H_k*JP  \f]
  	\f[ P_k=JP*JP^{'} \f]
       */
      virtual void CalculateCovariance( MatrixWrapper::Matrix& R_vf, MatrixWrapper::Matrix& H_i, MatrixWrapper::Matrix& invS, MatrixWrapper::Matrix& SR );

    private:
      /// Variable indicating the number of iterations of the filter
      unsigned int nr_iterations;
      /// the upper triangular matrix cholesky decompostion of the state covariance (\f$ JP*JP^{'}=P \f$)
      MatrixWrapper::Matrix JP;
      /// inv(S) matrix
      //BFL::Matrix invS;
      /// define \f$ S=H*JP*JP^{'}*H^{'}+D*R*D^{'} \f$
      //BFL::Matrix Sr;
      /// the Kalman gain factors
      //BFL::Matrix K_i;
    };  // class

} // End namespace BFL

#endif // __ITERATED_EXTENDED_KALMAN_FILTER__

