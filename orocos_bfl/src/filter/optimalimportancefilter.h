// Copyright (C) 2003 Klaas Gadeyne <first dot last at gmail dot com>

 /***************************************************************************
 *   This library is free software; you can redistribute it and/or         *
 *   modify it under the terms of the GNU General Public                   *
 *   License as published by the Free Software Foundation;                 *
 *   version 2 of the License.                                             *
 *                                                                         *
 *   As a special exception, you may use this file as part of a free       *
 *   software library without restriction.  Specifically, if other files   *
 *   instantiate templates or use macros or inline functions from this     *
 *   file, or you compile this file and link it with other files to        *
 *   produce an executable, this file does not by itself cause the         *
 *   resulting executable to be covered by the GNU General Public          *
 *   License.  This exception does not however invalidate any other        *
 *   reasons why the executable file might be covered by the GNU General   *
 *   Public License.                                                       *
 *                                                                         *
 *   This library is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU     *
 *   Lesser General Public License for more details.                       *
 *                                                                         *
 *   You should have received a copy of the GNU General Public             *
 *   License along with this library; if not, write to the Free Software   *
 *   Foundation, Inc., 51 Franklin Street, Fifth Floor,                                    *
 *   Boston, MA  02110-1301  USA                                *
 *                                                                         *
 ***************************************************************************/

// $Id$

#ifndef __OPTIMALIMPORTANCE_FILTER__
#define __OPTIMALIMPORTANCE_FILTER__

#include "particlefilter.h"

namespace BFL
{

/// Particular particle filter: Proposal PDF = Optimal Importance function
/** This is one (simple) particular implementation of a particle
    filter, in which the proposal density is equal to the pdf
    \f$ P(x_k | x_{k-1}, z_k) \f$.  Note that this pdf can only be
    easily (analytically) determined for a limited class of models!
    The current implementation focusses on systems with a linear
    measurement model.  @see

    @verbatim
    @Article{         doucet98bis,
    author =       {Doucet, Arnaud and Godsill, Simon and Andrieu, Christophe},
    title =        {On Sequential Monte Carlo Sampling Methods for
                    Bayesian Filtering},
    journal =      {Statistics and Computing},
    year =         {2000},
    volume =       {10},
    number =       {3},
    pages =        {197--208},
    }
    @endverbatim

    for a more thorough discussion about all these issues and the
    possible suboptimal alternatives in case one is not able to sample
    from the optimal importance function.

 */
template <typename StateVar, typename MeasVar> class Optimalimportancefilter
  : public ParticleFilter<StateVar,MeasVar>
{
 protected:
  /// Construct Optimal importance density from a sys and meas. model
  /**
     @param sysmodel system model to use
     @param measmodel measurement model to use for proposal construction
   */
  virtual void ConstructProposal(SystemModel<StateVar>* const sysmodel,
				 MeasurementModel<MeasVar,StateVar>* const measmodel);

 public:
  /// Constructor
  /** @pre you created the necessary models and the prior
      @param prior pointer to the Monte Carlo Pdf prior density
      @param resampleperiod fixed resampling period (if desired)
      @param resamplethreshold threshold used when dynamic resampling
      @param resamplescheme resampling scheme, see header file for
      different defines and their meaning
  */
  OptimalImportanceFilter(MCPdf<StateVar> * prior,
			  int resampleperiod = 0,
			  double resamplethreshold = 0,
			  int resamplescheme = DEFAULT_RS);

  /// Destructor
  virtual ~OptimalImportanceFilter();
  /// Copy constructor
  OptimalImportanceFilter(const OptimalImportanceFilter<StateVar,MeasVar> & filt);

  virtual void Update(SystemModel<StateVar>* const sysmodel,
		      const StateVar& u,
		      MeasurementModel<MeasVar,StateVar>* const measmodel,
		      const MeasVar& z,
		      const StateVar& s);
  virtual void Update(SystemModel<StateVar>* const sysmodel,
		      MeasurementModel<MeasVar,StateVar>* const measmodel,
		      const MeasVar& z,
		      const StateVar& s);
  virtual void Update(SystemModel<StateVar>* const sysmodel,
		      MeasurementModel<MeasVar,StateVar>* const measmodel,
		      const MeasVar& z);
  virtual void Update(SystemModel<StateVar>* const sysmodel,
		      const StateVar& u,
		      MeasurementModel<MeasVar,StateVar>* const measmodel,
		      const MeasVar& z);

  /// Only sysupdate
  virtual void Update(SystemModel<StateVar>* const sysmodel,
		      const StateVar& u);
  virtual void Update(SystemModel<StateVar>* const sysmodel);

  /// Only measupdate
  virtual void Update(MeasurementModel<MeasVar,StateVar>* const measmodel,
		      const MeasVar& z,
		      const StateVar& s);
  virtual void Update(MeasurementModel<MeasVar,StateVar>* const measmodel,
		      const MeasVar& z);
};

#include "optimalimportancefilter.cpp"

}

#endif // __OPTIMALIMPORTANCE_FILTER__
