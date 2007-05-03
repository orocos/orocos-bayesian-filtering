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

#ifndef __FILTER__
#define __FILTER__

#include "../model/systemmodel.h"
#include "../model/measurementmodel.h"
#include "../pdf/pdf.h"

namespace BFL
{
  using namespace std;

  /// Virtual Baseclass representing all bayesian filters
  /** This is the a virtual baseclass representing all similarities
      between Bayesian filters.  All Bayesian filter are related to i) a
      System Model, ii) a Measurement Model and iii) a Prior density
      reflecting the subjective belief of the person applying the filter
      BEFORE getting sensor or any other form of information about the
      system.  

      This class is the base class for particle filters, kalman filters,
      ...

      This class is also a template class with 2 templates, in this way
      it allows for "semi-discrete" systems, eg.  a system with a fixed
      number of states but with continuous observations, as eg needed in
      ASR.
    
      @see Pdf SystemModel MeasurementModel ConditionalPdf
      @bug For now, due to a "bug" (= non-existence of a feature :-) in
      the ConditionalPdf class, STATES AND INPUTS MUST BE OF THE SAME
      TYPE (both discrete, or both continuous!  This means that you can
      this class for the following type of systems:
      - States, inputs and measurements continuous (most frequently
      used?)
      - States and inputs continous, Measurements discrete
      - States and inputs discrete, Measurements continously
      - States, inputs and measurements discrete
    
      StateVar represents the form of the states and inputs
      MeasVar represents the form of the measurements
      BEWARE: THIS IS CONTRARY TO THE NOTATION USED IN "measurementmodel.h"
  */
  template <typename StateVar, typename MeasVar> class Filter
    {
    protected:
 
      /// prior Pdf
      Pdf<StateVar> * _prior;

      /// Pointer to the Posterior Pdf.
      /** The Posterior Pdf represents the subjective belief of the person
	  applying the filter AFTER processing inputs and measurements.
	  A filter does not maintain the beliefs at all timesteps t, since
	  this leads to non-constant (or ever growing if you prefer)
	  memory requirements.
	  However, it is possible, to copy the Posterior density at all
	  timesteps in your application by means of the PostGet() member
	  function
	  @see PostGet()
      */
      Pdf<StateVar> * _post;
  
      /// Represents the current timestep of the filter
      /** @todo Check wether this really belongs here
       */
      int _timestep;

      /// Actual implementation of Update, varies along filters
      /** @param sysmodel pointer to the used system model
	  @param u input param for proposal density
	  @param measmodel pointer to the used measurementmodel
	  @param z measurement param for proposal density
	  @param s sensor param for proposal density
      */
      virtual bool UpdateInternal(SystemModel<StateVar>* const sysmodel,
				  const StateVar& u,
				  MeasurementModel<MeasVar,StateVar>* const measmodel,
				  const MeasVar& z,
				  const StateVar& s)=0;

    public:
      /// Constructor
      /** @pre you created the prior
	  @param prior pointer to the prior Pdf
      */
      Filter(Pdf<StateVar> * prior);

      /// copy constructor
      /** @bug we should make a copy of the prior
       */
      Filter(const Filter<StateVar,MeasVar>& filt);

      /// destructor
      virtual ~Filter();

      /// Reset Filter 
      virtual void Reset(Pdf<StateVar> * prior);

      /// Full Update (system with inputs/sensing params)
      /** @param sysmodel pointer to the system model to use for update
	  @param u input to the system
	  @param measmodel pointer to the measurement model to use for update
	  @param z measurement
	  @param s "sensing parameter"
       */
      virtual bool Update(SystemModel<StateVar>* const sysmodel,
			  const StateVar& u,
			  MeasurementModel<MeasVar,StateVar>* const measmodel,
			  const MeasVar& z,
			  const StateVar& s);

      /// Full Update (system without inputs, with sensing params)
      /** @param sysmodel pointer to the system model to use for
	  update
	  @param measmodel pointer to the measurement model to use for
	  update
	  @param z measurement
	  @param s "sensing parameter"
       */
      virtual bool Update(SystemModel<StateVar>* const sysmodel,
			  MeasurementModel<MeasVar,StateVar>* const measmodel,
			  const MeasVar& z,
			  const StateVar& s);
      /// Full Update (system without inputs/sensing params)
      /** @param sysmodel pointer to the system model to use for
	  update
	  @param measmodel pointer to the measurement model to use for
	  update
	  @param z measurement
       */
      virtual bool Update(SystemModel<StateVar>* const sysmodel,
			  MeasurementModel<MeasVar,StateVar>* const measmodel,
			  const MeasVar& z);
      /// Full Update (system with inputs, without sensing params)
      /** @param sysmodel pointer to the system model to use for update
	  @param u input to the system
	  @param measmodel pointer to the measurement model to use for
	  update
	  @param z measurement
       */
      virtual bool Update(SystemModel<StateVar>* const sysmodel,
			  const StateVar& u,
			  MeasurementModel<MeasVar,StateVar>* const measmodel,
			  const MeasVar& z);

      /// System Update (system with inputs)
      /** @param sysmodel pointer to the system model to use for update
	  @param u input to the system
       */
      virtual bool Update(SystemModel<StateVar>* const sysmodel,
			  const StateVar& u);
      /// System Update (system without inputs)
      /** @param sysmodel pointer to the system model to use for update
       */
      virtual bool Update(SystemModel<StateVar>* const sysmodel);

      /// Measurement Update (system with "sensing params")
      /** @param measmodel pointer to the measurement model to use for
	  update 
	  @param z measurement
	  @param s "sensing parameter"
       */
      virtual bool Update(MeasurementModel<MeasVar,StateVar>* const measmodel,
			  const MeasVar& z,
			  const StateVar& s);
      /// Measurement Update (system without "sensing params")
      /** @param measmodel pointer to the measurement model to use for
	  update
	  @param z measurement
       */
      virtual bool Update(MeasurementModel<MeasVar,StateVar>* const measmodel,
			  const MeasVar& z);

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
#include "filter.cpp"

} // End namespace BFL

#endif // __FILTER__
