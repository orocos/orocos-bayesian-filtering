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
#include "kalmanfilter.h"
#include <cmath>

namespace BFL
{
  using namespace MatrixWrapper;


  KalmanFilter::KalmanFilter(Gaussian * prior)
    : Filter<ColumnVector,ColumnVector>(prior)
  {
    // create posterior dencity
    _post = new Gaussian(*prior);
  }

  KalmanFilter::~KalmanFilter()
  {
    delete _post;
  }

  void 
  KalmanFilter::CalculateSysUpdate(ColumnVector J, Matrix F, SymmetricMatrix Q)
  {
    Matrix temp = F * (Matrix)_post->CovarianceGet() * F.transpose() + (Matrix)Q;
    SymmetricMatrix Sigma_new(_post->DimensionGet());
    temp.convertToSymmetricMatrix(Sigma_new);
 
    // set new state gaussian
    PostMuSet   ( J );    
    PostSigmaSet( Sigma_new );
  }

  void 
  KalmanFilter::CalculateMeasUpdate(ColumnVector z, ColumnVector Z, Matrix H, SymmetricMatrix R)
  {
    // build K matrix
    Matrix S = ( H * (Matrix)(_post->CovarianceGet()) * (H.transpose()) ) + (Matrix)R;
    Matrix K = (Matrix)(_post->CovarianceGet()) * (H.transpose()) * (S.inverse());
  
    // calcutate new state gaussian
    ColumnVector Mu_new  = ( _post->ExpectedValueGet() + K * (z - Z)  );  
    Matrix Sigma_new_matrix = (Matrix)(_post->CovarianceGet()) - K * H * (Matrix)(_post->CovarianceGet());
    // convert to symmetric matrix
    SymmetricMatrix Sigma_new(_post->DimensionGet());
    Sigma_new_matrix.convertToSymmetricMatrix(Sigma_new);
  
    // set new state gaussian
    PostMuSet   ( Mu_new );
    PostSigmaSet( Sigma_new );
  
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
