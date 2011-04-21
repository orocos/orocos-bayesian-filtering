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
#include "nonminimalkalmanfilter.h"
#include "nonminimal_state/linearise.h"
#include "../pdf/nonlinearanalyticconditionalgaussian.h"

namespace BFL
{

#define NLinSys	   NonLinearAnalyticSystemModelGaussianUncertainty
#define NLinMeas   NonLinearAnalyticMeasurementModelGaussianUncertainty


  NonminimalKalmanFilter::NonminimalKalmanFilter(Gaussian* prior, 
						 unsigned int NrIterations,
						 vector<NLSysModel*>   minimalsysmodels,
						 vector<NLMeasModel*>  minimalmeasmodels,
						 vector<GiNaC::symbol> makelinear)
    : KalmanFilter(prior)
  {
    // create linearise
    Linear = new Linearise(minimalsysmodels, minimalmeasmodels, makelinear);
  
    // create STATE
    MinimalState    = Linear->NonlinearStateGet();
    NonminimalState = Linear->LinearStateGet();
    
    // create PRIOR
    ColumnVector mu_prior(NonminimalState.size());  mu_prior = 0;
    SymmetricMatrix sigma_prior(NonminimalState.size());  sigma_prior = 0;
    for (unsigned int i=0; i< NonminimalState.size(); i++)
      sigma_prior(i+1,i+1) = 333*333;
    NonminimalPrior = new Gaussian(mu_prior,sigma_prior);
    MinimalPrior    = prior;
    cout << "nonminimal prior " << *NonminimalPrior << endl;
  
    // create FILTER
    NonminimalFilter = new IteratedExtendedKalmanFilter(NonminimalPrior, NrIterations);
    MinimalFilter    = new IteratedExtendedKalmanFilter(MinimalPrior,    NrIterations);
  

    // create MODEL
    vector<GiNaC::symbol> empty_sym(0);
    ColumnVector mu_add(NonminimalState.size()); mu_add = 0;
    Gaussian additiveNoise(mu_add,NonminimalFilter->PostGet().CovarianceGet());
    NonLinearConditionalGaussian pdf(Linear->SubstitutionGet(), empty_sym, MinimalState, additiveNoise);
    MinimalMeasModel = new NLinMeas( &pdf );
  
    // show substitutions
    for (unsigned int i=0; i<NonminimalState.size() ; i++)
      cout << NonminimalState[i] << " -> " << Linear->SubstitutionGet()[i] << endl;
  }


  NonminimalKalmanFilter::~NonminimalKalmanFilter()
  {
    delete Linear;
    delete NonminimalPrior;
    delete NonminimalFilter;
    delete MinimalFilter;
    delete MinimalMeasModel;
  }

  void
  NonminimalKalmanFilter::SysUpdate(SystemModel<ColumnVector>* const sysmodel,
				    const ColumnVector& u)
  {
    NonminimalFilter->SysUpdate(Linear->LinearSysModelGet((NLinSys*)sysmodel),u);
  }

  void
  NonminimalKalmanFilter::MeasUpdate(MeasurementModel<ColumnVector,ColumnVector>* const measmodel,
				     const ColumnVector& z, 
				     const ColumnVector& s)
  {
    NonminimalFilter->MeasUpdate(Linear->LinearMeasModelGet((NLinMeas*)measmodel),z,s);
  }

}// End namespace
