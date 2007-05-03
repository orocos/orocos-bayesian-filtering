// $Id: backwardfilter.h 6736 2006-12-21 11:24:42Z tdelaet $
// Copyright (C) 2006 Tinne De Laet <first dot last at mech dot kuleuven dot be>
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

#ifndef __BACKWARDFILTER__
#define __BACKWARDFILTER__

#include "../model/systemmodel.h"
#include "../pdf/pdf.h"

namespace BFL
{
  using namespace std;

  /// Virtual Baseclass representing all bayesian backward filters
  /** This is the virtual baseclass representing all simularities between
       Bayesian backward filters.
      One class of smoothers corresponds to a forward pass (with a
      classical Bayesian filter) and a backward pass (with a backward
      Bayesian filter)
      The backward filters related with these smoothers are related to a
      System Model (they don't need a measurement model!(for as far I see now))

      This class is the base class for a rauch tung striebel backward filter, a 
      backward particle filters, ...

      @see Pdf SystemModel ConditionalPdf
    
      StateVar represents the form of the states and inputs
  */

  template <typename StateVar> class BackwardFilter
    {
    protected:
 
      /// prior Pdf
      Pdf<StateVar> * _prior;

      /// Pointer to the Posterior Pdf.
      /** The Posterior Pdf represents the subjective belief of the person
	  applying the filter AFTER processing a backwards step.**/
      Pdf<StateVar> * _post;
  
      /// Represents the current timestep of the filter
      int _timestep;

      /// Actual implementation of Update, varies along filters
      /** @param sysmodel pointer to the used system model
	  @param u input param for proposal density
      @param filtered_post is the posterior obtained by filtering of the timestep you want to smooth 
      */
      virtual bool UpdateInternal(SystemModel<StateVar>* const sysmodel,
				  const StateVar& u,
				  Pdf<StateVar>* const filtered_post)=0;

    public:
      /// Constructor
      /** @pre you created the prior
	  @param prior pointer to the prior Pdf
      */
      BackwardFilter(Pdf<StateVar> * prior);

      /// copy constructor
      BackwardFilter(const BackwardFilter<StateVar>& filt);

      /// destructor
      virtual ~BackwardFilter();

      /// Reset Filter 
      virtual void Reset(Pdf<StateVar> * prior);

      /// Full Update (system with inputs)
      /** @param sysmodel pointer to the system model to use for update
	  @param u input to the system
      @param filtered_post filtered posterior
       */
      virtual bool Update(SystemModel<StateVar>* const sysmodel,
              const StateVar& u,
              Pdf<StateVar>* const filtered_post);

      /// Full Update (system without inputs)
      /** @param sysmodel pointer to the system model to use for
	  update
      @param filtered_post filtered posterior
       */
      virtual bool Update(SystemModel<StateVar>* const sysmodel,
              Pdf<StateVar>* const filtered_post);

    /// Get Posterior density
      /** Get the current Posterior density
	  @return a pointer to the current posterior
      */
      virtual Pdf<StateVar> * PostGet();  

      /// Get current time
      /** Get the current time of the filter
	  @return the current timestep
      */
      int TimeStepGet() const;
    };    

  // For template instantiation
#include "backwardfilter.cpp"

} // End namespace BFL

#endif // __BACKWARDFILTER__
