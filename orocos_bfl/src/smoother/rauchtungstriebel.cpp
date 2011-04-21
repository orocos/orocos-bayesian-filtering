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
    , _x(prior->DimensionGet())
    , _xf(prior->DimensionGet())
    , _xpred(prior->DimensionGet())
    , _xsmooth(prior->DimensionGet())
    , _F(prior->DimensionGet(),prior->DimensionGet())
    , _Ppred(prior->DimensionGet(),prior->DimensionGet())
    , _Pxx(prior->DimensionGet(),prior->DimensionGet())
    , _K(prior->DimensionGet(),prior->DimensionGet())
    , _Psmooth(prior->DimensionGet(),prior->DimensionGet())
    , _Q(prior->DimensionGet())
    , _Sigma_new(prior->DimensionGet())
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
    _x = _post->ExpectedValueGet();
    _xf = filtered_post->ExpectedValueGet();
    _F = ((AnalyticSys*)sysmodel)->df_dxGet(u,_x);
    _Q = ((AnalyticSys*)sysmodel)->CovarianceGet(u,_x);
  
    _Ppred = _F * (Matrix)filtered_post->CovarianceGet() * _F.transpose() + (Matrix)_Q;
    _Pxx = (Matrix)filtered_post->CovarianceGet() * _F.transpose();
    _K = _Pxx * _Ppred.inverse();
    _xpred = ((AnalyticSys*)sysmodel)->PredictionGet(u,_xf);
    _xsmooth =  _xf + _K * (_x - _xpred);

    _Psmooth = (Matrix)filtered_post->CovarianceGet() - _K * ( _Ppred - (Matrix)_post->CovarianceGet() ) * _K.transpose();
    _Psmooth.convertToSymmetricMatrix(_Sigma_new);

    // set new state gaussian
    PostMuSet   ( _xsmooth );    
    PostSigmaSet( _Sigma_new );
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
