// $Id$
// Copyright (C) 2002 Klaas Gadeyne <first dot last at gmail dot com>
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
#ifndef __LINEAR_MEASUREMENT_MODEL_GAUSSIAN_UNCERTAINTY_IMPLICIT__
#define __LINEAR_MEASUREMENT_MODEL_GAUSSIAN_UNCERTAINTY_IMPLICIT__


#include "../pdf/gaussian.h"
#include "../pdf/linearanalyticconditionalgaussian.h"
#include "linearanalyticmeasurementmodel_gaussianuncertainty.h"

namespace BFL
{

  /// Class for linear analytic measurementmodels with additive gaussian noise
  // derived from an implicit measurement equation $h(x,z)=0$
  //
  /** This class represents all measurement models of the form
    \f[  0 = f (x,z) \f]
    as a linear measurement model with virtual measurement z_k^{virtual}
    \f[ z_k^{virtual} = H(x_k,z_k) \times x_k  + N(\mu(x_{k},z_k) ,\Sigma(x_k,z_k)) \f]
  */
  class LinearAnalyticMeasurementModelGaussianUncertainty_Implicit : public LinearAnalyticMeasurementModelGaussianUncertainty
    {
    public:
      /// Constructor
      /** @param pdf Conditional pdf, with Gaussian uncertainty
       */
      LinearAnalyticMeasurementModelGaussianUncertainty_Implicit( LinearAnalyticConditionalGaussian* pdf);
      /// Constructor
      LinearAnalyticMeasurementModelGaussianUncertainty_Implicit();

      // Default Copy constructor will do
      // LinearAnalyticMeasurementModelGaussianUncertainty_Implicit(const LinearAnalyticMeasurementModelGaussianUncertainty_Implicit& l);

      /// Destructor
      virtual ~LinearAnalyticMeasurementModelGaussianUncertainty_Implicit();

      // redefinition of virtual functions
      //
      virtual const  	MatrixWrapper::ColumnVector&    fGet     () const =0;
      virtual const 	int        TypeGet     ()  const=0;
      virtual 		MatrixWrapper::Matrix&          dfGet     (int number)  =0 ;
      /// Returns H-matrix calculated with measurement z and state x
       /** \f[ H = \frac{df}{dx} \mid_{ z, x} \f] used to determine the covariance of noise on the linear measurement equation
	  @param  z The value of the input in which the derivate is evaluated
	  @param  x The value in the state in which the derivate is evaluated
      */
      virtual MatrixWrapper::Matrix          df_dxGet     (const MatrixWrapper::ColumnVector& z, const MatrixWrapper::ColumnVector& x) =0;
      /// Returns D-matrix calculated with measurement z and state x
      /** \f[ D = \frac{df}{dz} \mid_{ z, x} \f] used to determine the covariance of noise on the linear measurement equation
	  @param z  The value of the input in which the derivate is evaluated
	  @param x The value in the state in which the derivate is evaluated
      */
      virtual MatrixWrapper::Matrix&          df_dzGet     (const MatrixWrapper::ColumnVector& z, const MatrixWrapper::ColumnVector& x)=0;
      /// Return a prediction for the mean of the noise on the linear measurement equation, calculated with measurements z and state x
      virtual MatrixWrapper::ColumnVector    PredictionGet(const MatrixWrapper::ColumnVector& z, const MatrixWrapper::ColumnVector& x)=0;
      /// Return a prediction for the mean of the noise on the linear measurement equation, using the current x and z
      virtual MatrixWrapper::ColumnVector    ExpectedValueGet()=0;
      /// Returns  covariance of the noise on the linearised  measurement model evaluated using measurements z and states x
      /** The linearised measurement equation look like:
 	  \f[ z_k^{virtual} = H(x_{k},z_k) \times x_k  + N(\mu(x_{k},z_k) ,\Sigma(x_k,z_k)) \f]
	  with noise
	  \f[ =N(\mu(x_{k},z_k), \Sigma(x_k,z_k))\f]
	  and covariance
	  \f[  \Sigma(x_k,z_k)= D(x_k,z_k)*R*D(x_k,z_k)'  \f]
	  and R the noise on the measurements z .
       */
      virtual MatrixWrapper::SymmetricMatrix& CovarianceGet()=0;
      /// Returns  covariance of the noise on the linearised  measurement model evaluated using current z and states x
      /** The linearised measurement equation look like:
 	  \f[ z_k^{virtual} = H(x_{k},z_k) \times x_k  + N(\mu(x_{k},z_k) ,\Sigma(x_k,z_k)) \f]
	  with noise
	  \f[ =N(\mu(x_{k},z_k), \Sigma(x_k,z_k))\f]
	  and covariance
	  \f[  \Sigma(x_k,z_k)= D(x_k,z_k)*R*D(x_k,z_k)'  \f]
	  and R the noise on the measurements z .
       */
      virtual MatrixWrapper::SymmetricMatrix CovarianceGet(const MatrixWrapper::ColumnVector& z, const MatrixWrapper::ColumnVector& x)=0;
      virtual void            Calculate(const MatrixWrapper::ColumnVector& x ,const  MatrixWrapper::ColumnVector& z,const MatrixWrapper::Matrix& R)=0;

      /// Returns square root of the covariance of the measurements z
      virtual const 	     MatrixWrapper::Matrix&          SRCovariance() const=0;
      /// Returns 1 if D-matrix equals the identity matrix else 0
      virtual const   	     int&     	    		     Is_Identity()  const=0;
    };
} // End namespace BFL

#endif // __LINEAR_MEASUREMENT_MODEL_GAUSSIAN_UNCERTAINTY__
