// $Id$
// Copyright (C) 2003 Klaas Gadeyne <first dot last at gmail dot com>
//                    Wim Meeussen  <wim dot meeussen at mech dot kuleuven dot be>
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
#include "kalmanfilter.h"
#include <cmath>

namespace BFL
{
  using namespace MatrixWrapper;


  KalmanFilter::KalmanFilter(Gaussian * prior)
    : Filter<ColumnVector,ColumnVector>(prior)
    , _Mu_new(prior->DimensionGet())
    , _Sigma_new(prior->DimensionGet())
    , _Sigma_temp(prior->DimensionGet(),prior->DimensionGet())
    , _Sigma_temp_par(prior->DimensionGet(),prior->DimensionGet())
  {
    // create posterior dencity
    _post = new Gaussian(*prior);
  }

  KalmanFilter::~KalmanFilter()
  {
    delete _post;
  }

  void
  KalmanFilter::AllocateMeasModel(const vector<unsigned int>& meas_dimensions)
  {
    unsigned int meas_dimension;
    for(int i = 0 ; i< meas_dimensions.size(); i++)
    {
        // find if variables with size meas_sizes[i] are already allocated
        meas_dimension = meas_dimensions[i];
        _mapMeasUpdateVariables_it =  _mapMeasUpdateVariables.find(meas_dimension);
        if( _mapMeasUpdateVariables_it == _mapMeasUpdateVariables.end())
        {
            //variables with size z.rows() not allocated yet
            _mapMeasUpdateVariables_it = (_mapMeasUpdateVariables.insert
                (std::pair<unsigned int, MeasUpdateVariables>( meas_dimension,MeasUpdateVariables(meas_dimension,_Mu_new.rows()) ))).first;
         }
     }
  }

  void
  KalmanFilter::AllocateMeasModel(const unsigned int& meas_dimension)
  {
     // find if variables with size meas_sizes[i] are already allocated
     _mapMeasUpdateVariables_it =  _mapMeasUpdateVariables.find(meas_dimension);
     if( _mapMeasUpdateVariables_it == _mapMeasUpdateVariables.end())
     {
         //variables with size z.rows() not allocated yet
         _mapMeasUpdateVariables_it = (_mapMeasUpdateVariables.insert
             (std::pair<unsigned int, MeasUpdateVariables>( meas_dimension,MeasUpdateVariables(meas_dimension,_Mu_new.rows()) ))).first;
      }
  }

  void
  KalmanFilter::CalculateSysUpdate(const ColumnVector& J, const Matrix& F, const SymmetricMatrix& Q)
  {
    _Sigma_temp = F * ( (Matrix)_post->CovarianceGet() * F.transpose());
    _Sigma_temp += (Matrix)Q;
    _Sigma_temp.convertToSymmetricMatrix(_Sigma_new);

    // set new state gaussian
    PostMuSet   ( J );
    PostSigmaSet( _Sigma_new );
  }

  void
  KalmanFilter::CalculateMeasUpdate(const ColumnVector& z, const ColumnVector& Z, const Matrix& H, const SymmetricMatrix& R)
  {
    // allocate measurement for z.rows() if needed
    AllocateMeasModel(z.rows());

    (_mapMeasUpdateVariables_it->second)._postHT =   (Matrix)(_post->CovarianceGet()) * H.transpose() ;
    (_mapMeasUpdateVariables_it->second)._S_Matrix =  H * (_mapMeasUpdateVariables_it->second)._postHT;
    (_mapMeasUpdateVariables_it->second)._S_Matrix += (Matrix)R;

    // _K = covariance * H' * S(-1)
    (_mapMeasUpdateVariables_it->second)._K =  (_mapMeasUpdateVariables_it->second)._postHT * ( (_mapMeasUpdateVariables_it->second)._S_Matrix.inverse());

    // calcutate new state gaussian
    // Mu = expectedValue + K*(z-Z)
    (_mapMeasUpdateVariables_it->second)._innov = z-Z;
     _Mu_new  =  (_mapMeasUpdateVariables_it->second)._K * (_mapMeasUpdateVariables_it->second)._innov  ;
     _Mu_new  +=  _post->ExpectedValueGet() ;
    // Sigma = post - K*H*post
    _Sigma_temp = (_post->CovarianceGet());
    _Sigma_temp_par = (_mapMeasUpdateVariables_it->second)._K * H ;
    _Sigma_temp -=  _Sigma_temp_par * (Matrix)(_post->CovarianceGet());
    // convert to symmetric matrix
    _Sigma_temp.convertToSymmetricMatrix(_Sigma_new);

    // set new state gaussian
    PostMuSet   ( _Mu_new );
    PostSigmaSet( _Sigma_new );

    /*
      cout << "H:\n" << H << endl;
      cout << "R:\n" << R << endl;
      cout << "Z:\n" << Z << endl;
      cout << "inov:\n" << z-Z << endl;
      cout << "S:\n" << S << endl;
      cout << "S.inverse:\n" << S.inverse() << endl;
      cout << "K:\n" << K << endl;
      cout << "Mu_new:\n" << Mu_new << endl;
      cout << "sigma_new\n" << Sigma_new << endl;
    */
  }

  void KalmanFilter::CalculateNis(const ColumnVector& z, const ColumnVector& Z, const Matrix& H, const SymmetricMatrix& R){
    /* The Normalised Innovation Squared (NIS) value
    * 
    *  innov^T * S^-1 * innov
    */
    // allocate measurement for z.rows() if needed
    AllocateMeasModel(z.rows());

    // Calculate Innovation Covariance
    (_mapMeasUpdateVariables_it->second)._postHT =   (Matrix)(_post->CovarianceGet()) * H.transpose() ;
    (_mapMeasUpdateVariables_it->second)._S_Matrix =  H * (_mapMeasUpdateVariables_it->second)._postHT;
    (_mapMeasUpdateVariables_it->second)._S_Matrix += (Matrix)R;
    
    // Calculate Innovation
    (_mapMeasUpdateVariables_it->second)._innov = z-Z;

    _Nis = ((_mapMeasUpdateVariables_it->second)._innov.transpose()) * (((_mapMeasUpdateVariables_it->second)._S_Matrix.inverse()) * ((_mapMeasUpdateVariables_it->second)._innov));
  }

  bool
  KalmanFilter::UpdateInternal(SystemModel<ColumnVector>* const sysmodel,
			       const ColumnVector& u,
			       MeasurementModel<ColumnVector,ColumnVector>* const measmodel,
			       const ColumnVector& z, const ColumnVector& s)
  {
    if (sysmodel != NULL)
      {
	SysUpdate(sysmodel,u);
      }
    if (measmodel != NULL)
      {
	MeasUpdate(measmodel,z,s);
      }
    return true;
  }

  void
  KalmanFilter::PostSigmaSet( const SymmetricMatrix& s)
  {
    dynamic_cast<Gaussian *>(_post)->CovarianceSet(s);
  }

  void
  KalmanFilter::PostMuSet( const ColumnVector& c)
  {
    dynamic_cast<Gaussian *>(_post)->ExpectedValueSet(c);
  }


  Gaussian*
  KalmanFilter::PostGet()
  {
    return (Gaussian*)Filter<ColumnVector,ColumnVector>::PostGet();
  }

}
