// $Id$
// Copyright (C) 2002 Klaas Gadeyne <first dot last at gmail dot com>
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

#ifndef __NONMINIMAL_KALMAN_FILTER__
#define __NONMINIMAL_KALMAN_FILTER__

#include "kalmanfilter.h"
#include "../pdf/conditionalpdf.h"
#include "../pdf/gaussian.h"
#include "../model/nonlinearanalyticmeasurementmodel_gaussianuncertainty.h"
#include "../model/nonlinearanalyticsystemmodel_gaussianuncertainty.h"
#include "nonminimal_state/linearise.h"
#include "../filter/iteratedextendedkalmanfilter.h"

namespace BFL
{

#define NLSysModel	NonLinearAnalyticSystemModelGaussianUncertainty
#define NLMeasModel	NonLinearAnalyticMeasurementModelGaussianUncertainty


  /** This is a class implementing the Kalman Filter (KF) class for
      Non Minimal State Kalman Filters.

      The System- and MeasurementUpdate equasions are not linear. Substituting the
      state by a non-minimal state will make the System- and MeasurementUpdate
      linear in the non-minimal state
      @see KalmanFilter
      @todo Seriously reimplement this class!
  */
  class NonminimalKalmanFilter : public KalmanFilter
    {
    public:
      /** Constructor
	  @pre you created the prior
	  @param prior pointer to the Gaussian prior density
	  @param NrIterations Number of iterations that will be used (this class uses 2 IEKF:
	  one for the minimal state and one for the nonminimal state)
	  @param minimalsysmodels vector of measurement models that can be used for measurement updates
	  @param minimalmeasmodels vector of system models that can be used for system updates
	  @param nonlinearstate symbols of state in wich the filter should be linearised
      */
      NonminimalKalmanFilter(Gaussian* prior,
			     unsigned int NrIterations,
			     vector<NLSysModel*>   minimalsysmodels,
			     vector<NLMeasModel*>  minimalmeasmodels,
			     vector<GiNaC::symbol> nonlinearstate = *(new vector<GiNaC::symbol>));

      /// Destructor
      virtual ~NonminimalKalmanFilter();

      // virtual functions
      virtual void SysUpdate(SystemModel<ColumnVector>* const sysmodel,
			     const ColumnVector& u);
      virtual void MeasUpdate(MeasurementModel<ColumnVector,ColumnVector>* const measmodel,
			      const ColumnVector& z,
			      const ColumnVector& s);

    private:
      vector<GiNaC::symbol>        MinimalState, NonminimalState;
      Linearise                    *Linear;
      IteratedExtendedKalmanFilter *NonminimalFilter, *MinimalFilter;
      Gaussian                     *NonminimalPrior, *MinimalPrior;
      NLMeasModel                  *MinimalMeasModel;

    };  // class

} // End namespace

#endif // __NONMINIMAL_KALMAN_FILTER__
