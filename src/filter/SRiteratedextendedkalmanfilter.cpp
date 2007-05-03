// $Id$
// Copyright (C) 2003 Klaas Gadeyne <first dot last at gmail dot com>
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
// Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
//  
#include "SRiteratedextendedkalmanfilter.h"
#include "../model/linearanalyticmeasurementmodel_gaussianuncertainty_implicit.h"

namespace BFL
{
using namespace MatrixWrapper;
#define AnalyticSys                   AnalyticSystemModelGaussianUncertainty
#define LinearAnalyticMeas_Implicit   LinearAnalyticMeasurementModelGaussianUncertainty_Implicit
#define Numerical_Limitation          100*100
  
SRIteratedExtendedKalmanFilter::SRIteratedExtendedKalmanFilter(Gaussian* prior, unsigned int nr_it)
    : KalmanFilter(prior),
      nr_iterations(nr_it), JP(prior->CovarianceGet().rows(),prior->CovarianceGet().rows())
  {
	  (prior->CovarianceGet()).cholesky(JP);
  }
	
  SRIteratedExtendedKalmanFilter::~SRIteratedExtendedKalmanFilter(){}

  void
  SRIteratedExtendedKalmanFilter::SysUpdate(SystemModel<MatrixWrapper::ColumnVector>* const sysmodel,const MatrixWrapper::ColumnVector& u)
  {

    MatrixWrapper::ColumnVector    x = _post->ExpectedValueGet();
    MatrixWrapper::ColumnVector    J = ((AnalyticSys*)sysmodel)->PredictionGet(u,x); 
    MatrixWrapper::Matrix          F = ((AnalyticSys*)sysmodel)->df_dxGet(u,x);
    MatrixWrapper::SymmetricMatrix Q = ((AnalyticSys*)sysmodel)->CovarianceGet(u,x);
  //  cout<<"JP1\n"<<JP<<endl;

    CalculateSysUpdate(J, F, Q);
  //  cout<<"JP2\n"<<JP<<endl;
   // cout<<"post_covar\n"<<_post->CovarianceGet()<<endl;
    
    ((_post->CovarianceGet()).cholesky(JP));
    JP = JP.transpose();
    // cout<<"JP3\n"<<JP<<endl;

  }

  void
  SRIteratedExtendedKalmanFilter::SysUpdate(SystemModel<MatrixWrapper::ColumnVector>* const sysmodel)
  {
    MatrixWrapper::ColumnVector u(0);
    SysUpdate(sysmodel, u);
  }

  void
  SRIteratedExtendedKalmanFilter::MeasUpdate(MeasurementModel<MatrixWrapper::ColumnVector,MatrixWrapper::ColumnVector>* const measmodel,const MatrixWrapper::ColumnVector& z,const MatrixWrapper::ColumnVector& s)
  {

    MatrixWrapper::Matrix invS(z.rows(),z.rows());
    MatrixWrapper::Matrix Sr(z.rows(),z.rows());
    MatrixWrapper::Matrix K_i(_post->CovarianceGet().rows(),z.rows());
    
    MatrixWrapper::ColumnVector    x_k = _post->ExpectedValueGet();
    MatrixWrapper::SymmetricMatrix P_k = _post->CovarianceGet();
    MatrixWrapper::ColumnVector    x_i = _post->ExpectedValueGet();
  
    MatrixWrapper::Matrix       	H_i;    	MatrixWrapper::SymmetricMatrix  	 R_i;
    MatrixWrapper::Matrix  		R_vf;		MatrixWrapper::Matrix       		 SR_vf;
    MatrixWrapper::ColumnVector     	Z_i;
    MatrixWrapper::Matrix 	   	U; 		MatrixWrapper::ColumnVector              V;   	MatrixWrapper::Matrix             W;   
    MatrixWrapper::Matrix         	JP1;      	int change;				
    

    Matrix diag(JP.rows(),JP.columns());
    Matrix invdiag(JP.rows(),JP.columns());
    diag=0;invdiag=0;change=0;
    V=0;U=0;W=0;

    // matrix determining the numerical limitations of covariance matrix:
    for(unsigned int j=1;j<JP.rows()+1;j++){diag(j,j)=100; invdiag(j,j)=0.01;}
    
    
    for (unsigned int i=1; i<nr_iterations+1; i++)
      {	
	x_i = _post->ExpectedValueGet();

	H_i  = ((LinearAnalyticMeas_Implicit*)measmodel)->df_dxGet(s,x_i);
	Z_i  = ((LinearAnalyticMeas_Implicit*)measmodel)->ExpectedValueGet() + ( H_i * (x_k - x_i) );  
	
	R_i  = ((LinearAnalyticMeas_Implicit*)measmodel)->CovarianceGet();	
	SR_vf  = ((LinearAnalyticMeas_Implicit*)measmodel)->SRCovariance();
	
	// check two different types of Kalman filters:
	 if(((LinearAnalyticMeas_Implicit*)measmodel)->Is_Identity()==1)
         {
                     R_vf = SR_vf.transpose();
         }
         else
         {
                     R_i.cholesky(R_vf);
                     R_vf = R_vf.transpose();
         }
	
	// numerical limitations
	// The singular values of the Square root covariance matrix are limited the the value of 10e-4 
	// because of numerical stabilisation of the Kalman filter algorithm.  
	 JP.SVD(V,U,W);	
	 MatrixWrapper::Matrix V_matrix(U.columns(),W.columns());
	for(unsigned int k=1;k<JP.rows()+1;k++)
	{
   	        V_matrix(k,k) = V(k);
		V(k)=max(V(k),1.0/(Numerical_Limitation));
		if(V(k)==1/(Numerical_Limitation)){change=1;}
	}
	if(change==1)
	{	
		JP   = U*V_matrix*(W.transpose());
	}

	// end limitations	

	CalculateMatrix(H_i,  R_i , invS , K_i , Sr );

	CalculateMean(x_k, z, Z_i , K_i);

	if (i==nr_iterations)
	{
		CalculateCovariance( R_vf, H_i, invS, Sr );
	}

    }  
  }


  void
  SRIteratedExtendedKalmanFilter::MeasUpdate(MeasurementModel<ColumnVector,ColumnVector>* const measmodel,
					   const ColumnVector& z)
  {
    ColumnVector s(0);
    MeasUpdate(measmodel, z, s);
  }
  
  Matrix  SRIteratedExtendedKalmanFilter::SRCovarianceGet() const
 {
		 return (Matrix) JP;
 }
 
 void SRIteratedExtendedKalmanFilter::SRCovarianceSet(Matrix JP_new)
 {
		 JP=JP_new;
 }
 
 void SRIteratedExtendedKalmanFilter::PriorSet(ColumnVector& X_prior,SymmetricMatrix& P_prior)
 {
 	PostMuSet( X_prior );
	PostSigmaSet( P_prior );
 }
void 
  SRIteratedExtendedKalmanFilter::CalculateMeasUpdate(ColumnVector z, ColumnVector Z, Matrix H, SymmetricMatrix R)
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
    
      }
void
  SRIteratedExtendedKalmanFilter::CalculateMatrix(Matrix& H_i, SymmetricMatrix& R_i, Matrix& invS, Matrix& K_i, Matrix& Sr)
  {
	MatrixWrapper::Matrix S_i1,S_i2,S_temp1;
	MatrixWrapper::SymmetricMatrix S_temp2,S_temp;
	S_i1 = ( H_i * (Matrix)JP * (Matrix) (JP.transpose())* (H_i.transpose()) );
	S_i2 = (Matrix) R_i;
	S_temp1 = (S_i1 + S_i2).transpose();
	S_temp1.convertToSymmetricMatrix(S_temp);
	
	S_temp.cholesky(Sr);
	Sr  = Sr.transpose();
	
	invS = Sr.inverse();
	K_i  = JP*(JP.transpose())*(H_i.transpose())*(invS.transpose())*invS;
	  
/*	  cout<<"H_i\n"<<H_i<<endl;
	  cout<<"JP\n"<<JP<<endl;
	  cout<<"S_i1\n"<<S_i1<<endl;
	  cout<<"S_i1\n"<<S_i1<<endl;
	  cout<<"S_i2\n"<<S_i2<<endl;
	  cout<<"K_i\n"<<K_i<<endl;
	  */
    }  

void 
  SRIteratedExtendedKalmanFilter::CalculateMean(ColumnVector& x_k, const  ColumnVector& z, ColumnVector& Z_i ,Matrix& K_i)
  {
	MatrixWrapper::ColumnVector x_i;
	x_i  = x_k + K_i * (z - Z_i);
	PostMuSet( x_i );
  }

void 
  SRIteratedExtendedKalmanFilter::CalculateCovariance(Matrix& R_vf, Matrix& H_i, Matrix& invS ,Matrix& Sr)
  {
	MatrixWrapper::Matrix temp;
	temp = (Matrix)R_vf+(Matrix)Sr;	
	JP   = (Matrix)JP -(Matrix)JP*(Matrix)(JP.transpose()) * (H_i.transpose()) * (Matrix)(invS.transpose())*(Matrix)(temp.inverse())*H_i*(Matrix)JP;
	MatrixWrapper::SymmetricMatrix Sigma;
	MatrixWrapper::Matrix Sigma1;
	Sigma1=(JP*(JP.transpose())).transpose();
	Sigma1.convertToSymmetricMatrix(Sigma);
	PostSigmaSet(Sigma);
  }
}
