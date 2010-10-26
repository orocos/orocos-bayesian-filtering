// $Id: rauchtungstriebel.h 6736 2006-12-22 11:24:42Z tdelaet $
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

#ifndef __RAUCHTUNGSTRIEBEL__
#define __RAUCHTUNGSTRIEBEL__

#include "backwardfilter.h"
#include "../pdf/gaussian.h"
#include "../pdf/conditionalpdf.h"
#include "../model/analyticsystemmodel_gaussianuncertainty.h"

namespace BFL
{

/// Class representing  all Rauch-Tung-Striebel backward filters
/** This is a class representing the Rauch-Tung-Striebel backward filter.
 *  It is a backward filter in which the Posterior
    density is represented by a Gaussian density.  Rauch-Tung-Striebel backward filter  are
    only applicable to continuous systems.

    The system of updating the Posterior density is implemented in this
    base class.

    @see Gaussian
    @see LinearAnalyticSystemModelGaussianUncertainty
*/
class RauchTungStriebel : public BackwardFilter<MatrixWrapper::ColumnVector>
{
public:
  /// Constructor
  /** @pre you created the prior
      @param prior pointer to the Gaussian Pdf prior density
  */
  RauchTungStriebel(Gaussian* prior);

  /// Destructor
  virtual ~RauchTungStriebel();

protected:

  /// Set covariance of posterior estimate
  void PostSigmaSet( const MatrixWrapper::SymmetricMatrix& s);

  /// Set expected value of posterior estimate
  void PostMuSet( const MatrixWrapper::ColumnVector& c);

  /// System Update
  /** Update the filter's Posterior density using the deterministic
      inputs to the system and the system model
      @param sysmodel pointer to the system model the filter should use
      @param u input to the system
      @param filtered_post posterior from forward Bayesian filter
  */
  virtual void SysUpdate(SystemModel<MatrixWrapper::ColumnVector>* const sysmodel, const MatrixWrapper::ColumnVector& u , Pdf<ColumnVector>* const filtered_post);

  virtual bool UpdateInternal(SystemModel<ColumnVector>* const sysmodel, const ColumnVector& u,  Pdf<ColumnVector>* const filtered_post);

private:
    // Variables to avoid allocation during sysupdate call
    ColumnVector _x, _xf, _xpred, _xsmooth;
    Matrix _F, _Ppred, _Pxx, _K, _Psmooth;
    SymmetricMatrix _Q, _Sigma_new;
}; // class

} // End namespace BFL

#endif //__RAUCHTUNGSTRIEBEL__
