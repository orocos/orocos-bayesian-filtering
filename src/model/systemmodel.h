// $Id$
// Copyright (C) 2002 Klaas Gadeyne <first dot last at gmail dot com>
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
#ifndef __SYSTEM_MODEL__
#define __SYSTEM_MODEL__

#include "../pdf/conditionalpdf.h"

namespace BFL
{

  /** Template class representing all possible (continu and discrete)
      System Models
      @todo Check if there should be a "model" base class...
      @bug Currently supports only systemmodels of the form P(x | x, u), 
      where both u and x are continu or discrete.  So it
      lacks support for mixed systems () and systems with extra
      parameters.  You are welcome to provide an API and implementation
      for this :-)
  */
  template<typename T> class SystemModel
    {
    protected:

      /// ConditionalPdf representing \f$ P(X_k | X_{k-1}, U_{k}) \f$
      /* @bug Since, for now, the library only supports only conditional
	 arguments of the same type, both X and U have to be of the same
	 type (ie both continu or both discrete!).  I imagine there must
	 be systems for which this approach is not general enough @see
	 ConditionalPdf
      */
      ConditionalPdf<T,T>* _SystemPdf;

      /// System with no inputs?
      bool _systemWithoutInputs;

    public:
      /// Constructor
      /** @param systempdf ConditionalPdf<T,T> representing \f$ P(X_k |
	  X_{k-1}, U_{k}) \f$
	  @see STATE_SIZE, INPUT_SIZE, _SystemPdf
      */
      SystemModel(ConditionalPdf<T,T>* systempdf=NULL);
  
      /// Destructor
      virtual ~SystemModel();
  
      /// Copy constructor
      /// SystemModel(const SystemModel<T>& model);

      /// Get State Size
      /** @return the statesize of the system
       */
      int StateSizeGet() const;

      /// Number of Conditional Arguments..
      bool SystemWithoutInputs() const;

      // NO LONGER RELEVANT
      // void StateSizeSet(int); // necessary??
      // Get Input Size
      /* @return the statesize of the system
       */
      // int InputSizeGet() const;
      // void InputSizeSet(int); // necessary??

      /// Get the SystemPDF
      /** @return a reference to the ConditionalPdf describing the system
       */
      ConditionalPdf<T,T>* SystemPdfGet();
  
      /// Set the SystemPDF
      /** @param pdf a reference to the ConditionalPdf describing the system
       */
      void SystemPdfSet(ConditionalPdf<T,T>* pdf);

      /// Simulate the system
      /** @param x current state of the system
	  @param u input to the system
	  @return State where we arrive by simulating the system model for
	  1 step
	  @param sampling_method the sampling method to be used while
	  sampling from the Conditional Pdf describing the system (if not
	  specified = DEFAULT)
	  @param sampling_args Sometimes a sampling method can have some
	  extra parameters (eg mcmc sampling)
	  @note Maybe the return value would better be a Sample<T> instead
	  of a T 
      */
      T Simulate (const T& x, const T& u, int sampling_method = DEFAULT, void * sampling_args = NULL);
      /// Simulate the system (no input system)
      /** @param x current state of the system
	  @return State where we arrive by simulating the system model for
	  1 step
	  @note Maybe the return value would better be a Sample<T> instead
	  of a T 
	  @param sampling_method the sampling method to be used while
	  sampling from the Conditional Pdf describing the system (if not
	  specified = DEFAULT)
	  @param sampling_args Sometimes a sampling method can have some
	  extra parameters (eg mcmc sampling)
      */
  
      T Simulate (const T& x, int sampling_method = DEFAULT, void * sampling_args = NULL);

      /// Get the probability of arriving in a next state
      /** @param x_k the next state (at time k)
	  @param x_kminusone the current state (at time k-1)
	  @param u  the input
	  @return the probability value
      */
  
      Probability ProbabilityGet(const T& x_k, const T& x_kminusone, const T& u );

      /// Get the probability of arriving in a next state
      /** (no-input-system)
	  @param x_k the next state (at time k)
	  @param x_kminusone the current state (at time k-1)
	  @return the probability value
      */
      Probability ProbabilityGet(const T& x_k, const T& x_kminusone );
  

    };

#include "systemmodel.cpp"

} // End namespace BFL

#endif // __SYSTEM_MODEL__
