// $Id: histogramfilter.h 14935 2007-12-17  $
// Copyright (C) 2007 Tinne De Laet  <tinne dot delaet at mech dot kuleuven dot be>
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

#ifndef __HISTOGRAM_FILTER__
#define __HISTOGRAM_FILTER__

#include "filter.h"
#include "../pdf/discretepdf.h"
#include "../model/measurementmodel.h"
#include "../model/discretesystemmodel.h"

namespace BFL
{

/// Class representing the histogram filter
/** This is a class representing the histogram filter.
    A histogram filter is the basic discrete state filter for histogram
    representations of the state.
    The implementation is based upon Probabilistic Robotics book of Thrun, Burgard, Fox

    @Book{              ThrunBurgardFox2005,
      author          = {Thrun, S. and Burgard, W. and Fox, D.},
      title           = {Probabilistic Robotics},
      publisher       = {MIT Press},
      year            = {2005},
      issn_isbn       = {0-262-20162-3},
      annote          = {\url{http://www.probabilistic-robotics.org}},
      keywords        = {Bayes theory, estimation}
    }
    The system of updating the Posterior density is implemented in this
    class.
*/

template <typename MeasVar> class HistogramFilter : public Filter<int,MeasVar>
{
public:
  /// Constructor
  /** @pre you created the prior
      @param prior pointer to the Discrete Pdf prior density
  */
  HistogramFilter(DiscretePdf* prior);

  /// Destructor
  virtual ~HistogramFilter();

  // implement virtual function
  virtual DiscretePdf* PostGet();

protected:
  /// While updating store list of old probabilities
  vector<Probability > _old_prob;
  /// While updating store list of new probabilities
  vector<Probability > _new_prob;

  /** Calculate Discrete filter System Update
    @param sysmodel pointer to the system model the filter should use
    @param u input to the system
  */
  void SysUpdate(SystemModel<int>* const sysmodel,
			 const int& u);

  /// Measurement Update
  /** Update the filter's Posterior density using the sensor
      measurements, an input and the measurement model.
      @param measmodel pointer to the measurement model the filter
      should use
      @param z sensor measurement
      @param s input to the system (must be of the same type as u
      for now, since this was not yet implemented in ConditionalPdf
  */
  void MeasUpdate(MeasurementModel<MeasVar,int>* const measmodel,
			  const MeasVar& z,
			  const int& s);

  bool UpdateInternal(SystemModel<int>* const sysmodel,
			      const int& u,
			      MeasurementModel<MeasVar,int>* const measmodel,
			      const MeasVar& z,
			      const int& s);
}; // class


#include "histogramfilter.cpp"

} // End namespace BFL

#endif // __HISTOGRAM_FILTER__
