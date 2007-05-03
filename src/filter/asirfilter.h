// Copyright (C) 2003 Klaas Gadeyne <first dot last at gmail dot com>
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
// $Id$

#ifndef __ASIR_FILTER__
#define __ASIR_FILTER__

#include "particlefilter.h"

namespace BFL
{

  /// ASIR: Auxiliary Particle Filter 
  /** This is a possible implementation of a particle filter, which in
      some/many  cases will yield better results, since measures are
      taken to make the proposal density more similar to the posterior.  
      See 

      @verbatim
      @Article{	  pitt_auxiliary,
      author	= {Pitt, M. and Shephard, N.},
      title	= {Filtering via simulation: auxiliary particle filter},
      journal	= {Journal of the American Statistical Association},
      year	= {1999},
      note	= {forthcoming}
      }
      @endverbatim

      for more details.

      Note that this particular implementation:
      - Uses currently a particular version of important sampling
      - Still uses the the system pdf as "actual" proposal density,
      which amongs others means the results will be suboptimal in case
      of a log concave measurement density.

      Particular issue with of proposalstep in case of ASIR
      filter... 
      The current implementation uses the approximation (see p. 11 of
      the Pitt and Shephard paper, we use their notation here ---
      \f[ \alpha \f] denoting the state and \f[ y \f] denoting the
      measurement)
      \f[ f(y_{t+1} | \alpha_{t+1}) \approx f(y_{t+1} | \mu_{t+1}) \f]
      where \f[ \mu \f] is the \em mean value of the system pdf.

      Note that the ASIR needs the measurementmodel for its
      proposalstep, to obtain a better proposal.  Note also that
      normally, the ASIR performs better that the standard SIR filter
      in case of outliers, but worse for "normal data" (due to the
      extra resampling stage).
  */
  template <typename StateVar, typename MeasVar> class ASIRFilter
    : public ParticleFilter<StateVar,MeasVar>
    {
    protected:
      /// Actual implementation of updateinternal
      virtual void UpdateInternal(SystemModel<StateVar>* const sysmodel,
				  const StateVar& u,
				  MeasurementModel<MeasVar,StateVar>* const measmodel,
				  const MeasVar& z,
				  const StateVar& s);

    public:
      /// Constructor
      /** @pre you created the necessary models and the prior
	  @param prior pointer to the Monte Carlo Pdf prior density
	  @param resampleperiod fixed resampling period (if desired)
	  @param resamplethreshold threshold used when dynamic resampling
	  @param resamplescheme resampling scheme, see header file for
	  different defines and their meaning
      */
      ASIRFilter(MCPdf<StateVar> * prior,
		 int resampleperiod = 0,
		 double resamplethreshold = 0,
		 int resamplescheme = DEFAULT_RS);

      /// Destructor
      virtual ~ASIRFilter();

      // Default Copy constructor will do
    };

#include "asirfilter.cpp"

}

#endif // __ASIR_FILTER__
