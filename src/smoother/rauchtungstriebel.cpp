// $Id: rauchtungstriebel.cpp 6736 2006-12-22 11:24:42Z tdelaet $
// Copyright (C) 2006  Tinne De Laet <first dot last at mech dot kuleuven dot be>
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
#include "rauchtungstriebel.h"
#include <cmath>

namespace BFL
{
  using namespace MatrixWrapper;

#define AnalyticSys    AnalyticSystemModelGaussianUncertainty

  RauchTungStriebel::RauchTungStriebel(Gaussian * prior)
    : BackwardFilter<ColumnVector>(prior)
  {
    // create posterior dencity
    _post = new Gaussian(*prior);
  }

  RauchTungStriebel::~RauchTungStriebel()
  {
    delete _post;
  }

  void
  RauchTungStriebel::SysUpdate(SystemModel<ColumnVector>* const sysmodel, const ColumnVector& u,  Pdf<ColumnVector>* const filtered_post)
  {
    ColumnVector x = _post->ExpectedValueGet();
    ColumnVector xf = filtered_post->ExpectedValueGet();
    Matrix F = ((AnalyticSys*)sysmodel)->df_dxGet(u,x);
    SymmetricMatrix Q = ((AnalyticSys*)sysmodel)->CovarianceGet(u,x);
  
    Matrix Ppred = F * (Matrix)filtered_post->CovarianceGet() * F.transpose() + (Matrix)Q;
    Matrix Pxx = (Matrix)filtered_post->CovarianceGet() * F.transpose();
    Matrix K = Pxx * Ppred.inverse();
    ColumnVector xpred = ((AnalyticSys*)sysmodel)->PredictionGet(u,xf);
    ColumnVector xsmooth =  xf + K * (x - xpred);

    Matrix Psmooth = (Matrix)filtered_post->CovarianceGet() - K * ( Ppred - (Matrix)_post->CovarianceGet() ) * K.transpose();
    SymmetricMatrix Sigma_new(_post->DimensionGet());
    Psmooth.convertToSymmetricMatrix(Sigma_new);

    // set new state gaussian
    PostMuSet   ( xsmooth );    
    PostSigmaSet( Sigma_new );
  }

  bool 
  RauchTungStriebel::UpdateInternal(SystemModel<ColumnVector>* const sysmodel, const ColumnVector& u,Pdf<ColumnVector>* const  filtered_post)
  {
	SysUpdate(sysmodel,u,(Gaussian*)filtered_post);
    return true;
  }

  void
  RauchTungStriebel::PostSigmaSet( const SymmetricMatrix& s)
  {
    dynamic_cast<Gaussian *>(_post)->CovarianceSet(s);
  }

  void
  RauchTungStriebel::PostMuSet( const ColumnVector& c)
  {
    dynamic_cast<Gaussian *>(_post)->ExpectedValueSet(c);
  }

}
