// $Id$
// Copyright (C) 2002 Klaas Gadeyne <first dot last at gmail dot com>
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
 *   Foundation, Inc., 51 Franklin Street, Fifth Floor,                                    *
 *   Boston, MA  02110-1301  USA                                *
 *                                                                         *
 ***************************************************************************/

#ifndef __FILTER__
#define __FILTER__

#include "../model/systemmodel.h"
#include "../model/measurementmodel.h"
#include "../pdf/pdf.h"

namespace BFL
{
  using namespace std;

  /// Abstract class representing an interface for Bayesian Filters
  /** This is the Abstract interface class that defines the interface
      of Bayesian filters.  These filters are all related to i) a
      System Model, ii) a Measurement Model and iii) a Prior density
      reflecting the subjective belief of the person applying the filter
      BEFORE getting sensor or any other form of information about the
      modeled system.

      This class is the base class for particle filters, kalman filters,
      ...

      This class is a template class with 2 templates.  In this way
      it allows filtering for "semi-discrete" models, eg. models
      with a fixed number of states (discrete states) but with
      continuous observations, as needed in Automatic Speech
      Recognition.

      @see Pdf SystemModel MeasurementModel ConditionalPdf
      @bug For now, due to a "bug" (= non-existence of a feature :-) in
      the ConditionalPdf class, STATES AND INPUTS MUST BE OF THE SAME
      TYPE (both discrete, or both continuous!  This means that you can
      use this class for the following model types:
      - States, inputs and measurements continuous (most frequently
      used?)
      - States and inputs continous, Measurements discrete
      - States and inputs discrete, Measurements continous
      - States, inputs and measurements discrete

      StateVar represents the nature of the states and inputs
      MeasVar represents the nature of the measurements

      BEWARE: The order of the template arguments is reversed with
      respect to the notation used in "measurementmodel.h"
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
