// Copyright (C) 2009 Tinne De Laet <first dot last at gmail dot com>
// $Id: mixtureBoostrapFilter.h 2009-02-03 tdelaet $
//
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
 *   Foundation, Inc., 59 Temple Place,                                    *
 *   Suite 330, Boston, MA  02111-1307  USA                                *
 *                                                                         *
 ***************************************************************************/

#ifndef __MIXTURE_BOOTSTRAP_FILTER__
#define __MIXTURE_BOOTSTRAP_FILTER__

#include "mixtureParticleFilter.h"

namespace BFL
{

/// Particular mixture particle filter : Proposal PDF = SystemPDF
/** This is one (simple) particular implementation of a particle
    filter, in which the proposal density is equal to the pdf
    describing the system model (aka as SystemPdf), and involving a
    resampling step

    The reason why I chose the name bootstrap filter is the fact that
    this is the name used in the book by Doucet et al.

    @verbatim
    @Book{		  doucet_book,
    editor	= {Doucet, Arnaud and de Freytas, Nando and Gordon, Neil},
    title		= {{S}equential {M}onte {C}arlo {M}ethods in {P}ractice},
    publisher	= {Springer--Verlag},
    year		= {2001},
    series	= {Statistics for engineering and information science},
    month		= {january},
    annote	= {see http://www-sigproc.eng.cam.ac.uk/~ad2/book.html}
    }
    @endverbatim

    (and I presume this will become a/the standard book about particle
    filtering).  Typical for the bootstrap filter is the fact that the
    proposal density is chosen to be the SystemPdf of the SystemModel.
    So there is no proposal density in the constructor here

    @todo The implementation is very slow for the moment.  It would
    probably be much faster to add a vector<WeightedSample> to the
    private members of this class.
    @see Pdf
 */
template <typename StateVar, typename MeasVar> class MixtureBootstrapFilter
  : public MixtureParticleFilter<StateVar,MeasVar>
{
 protected:
  /// Actual implementation of updateinternal
  virtual bool UpdateInternal(SystemModel<StateVar>* const sysmodel,
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
	  @param maintainMixturePeriod fixed mixture maintainance period
  */
  MixtureBootstrapFilter(Mixture<StateVar> * prior,
		 int resampleperiod = 0,
		 double resamplethreshold = 0,
		 int resamplescheme = DEFAULT_RS,
         int maintainMixturePeriod = 1 );

  /// Constructor
  /** @pre you created the necessary models and the prior
      @param prior pointer to the Monte Carlo Pdf prior density
      @param post pointer to the Monte Carlo Pdf post density
      @param resampleperiod fixed resampling period (if desired)
      @param resamplethreshold threshold used when dynamic resampling
      @param resamplescheme resampling scheme, see header file for
      different defines and their meaning
	  @param maintainMixturePeriod fixed mixture maintainance period
  */
  MixtureBootstrapFilter(Mixture<StateVar> * prior,
		  Mixture<StateVar> * post,
		 int resampleperiod = 0,
		 double resamplethreshold = 0,
		 int resamplescheme = DEFAULT_RS,
         int maintainMixturePeriod = 1 );

  /// Destructor
  virtual ~MixtureBootstrapFilter();

  // Default Copy constructor will do

};

#include "mixtureBootstrapFilter.cpp"

}

#endif // __MIXTURE_BOOTSTRAP_FILTER__
