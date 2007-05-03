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

#ifndef __BOOTSTRAP_FILTER__
#define __BOOTSTRAP_FILTER__

#include "particlefilter.h"

namespace BFL
{

/// Particular particle filter : Proposal PDF = SystemPDF
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
template <typename StateVar, typename MeasVar> class BootstrapFilter
  : public ParticleFilter<StateVar,MeasVar>
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
  */
  BootstrapFilter(MCPdf<StateVar> * prior,
		 int resampleperiod = 0,
		 double resamplethreshold = 0,
		 int resamplescheme = DEFAULT_RS);

  /// Constructor
  /** @pre you created the necessary models and the prior
      @param prior pointer to the Monte Carlo Pdf prior density
      @param post pointer to the Monte Carlo Pdf post density
      @param resampleperiod fixed resampling period (if desired)
      @param resamplethreshold threshold used when dynamic resampling
      @param resamplescheme resampling scheme, see header file for
      different defines and their meaning
  */
  BootstrapFilter(MCPdf<StateVar> * prior,
		  MCPdf<StateVar> * post,
		 int resampleperiod = 0,
		 double resamplethreshold = 0,
		 int resamplescheme = DEFAULT_RS);

  /// Destructor
  virtual ~BootstrapFilter();

  // Default Copy constructor will do

};

#include "bootstrapfilter.cpp"

}

#endif // __BOOTSTRAP_FILTER__
