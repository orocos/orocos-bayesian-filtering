// $Id$
// Copyright (C) 2002 Klaas Gadeyne <first dot last at gmail dot com>
//               2008 Tinne De Laet <first dot last at mech dot kuleuven dot be>
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
#ifndef DISCRETEPDF_H
#define DISCRETEPDF_H

#include "pdf.h"
#include "../wrappers/matrix/vector_wrapper.h"
#include "../wrappers/matrix/matrix_wrapper.h"
#include <vector>

namespace BFL
{
  /// Class representing a PDF on a discrete variable
  /** This class is an instantation from the template class Pdf, with
      added methods to get a set the probability of a certain discrete
      value (methods only relevant for discrete pdfs)
  */
  class DiscretePdf : public Pdf<int> // inherit abstract_template_class
    {
    protected:
      /// The number of discrete state
      unsigned int _num_states;

      /// Pointer to the discrete PDF-values, the sum of the elements = 1
      vector<Probability> *_Values_p;

      /// Normalize all the probabilities (eg. after setting a probability)
      bool NormalizeProbs();

      /// STL-vector containing the Cumulative PDF (for efficient sampling)
      vector<double> _CumPDF;

      /// Updates the cumPDF
      bool CumPDFUpdate();

    public:
      /// Constructor (dimension = number of classes) An equal probability is set for all classes
      /** @param num_states number of different classes or states
       */
      DiscretePdf(unsigned int num_states=0);

      /// Copy Constructor
      DiscretePdf(const DiscretePdf &);

      /// Destructor
      virtual ~DiscretePdf();

      /// Get the number of discrete States
      unsigned int NumStatesGet()const;

      /// Implementation of virtual base class method
      Probability ProbabilityGet(const unsigned int& state) const;

      /// Function to change/set the probability of a single state
      /** Changes the probabilities such that AFTER normalization the
          probability of the state "state" is equal to the probability a
           @param state number of state of which the probability will be set
           @param a probability value to which the probability of state "state"
            will be set (must be <= 1)
      */
      bool ProbabilitySet(unsigned int state, Probability a);

      bool SampleFrom (vector<Sample<int> >& list_samples,
		       const unsigned int num_samples,
		       int method = DEFAULT,
		       void * args = NULL) const;
      bool SampleFrom (Sample<int>& one_sample, int method = DEFAULT, void * args = NULL) const;

      /// Get all probabilities
      vector<Probability> ProbabilitiesGet() const;

      /// Set all probabilities
      /**  @param values vector<Probability> containing the new probability values.
           The sum of the probabilities of this list is not required to be one
           since the normalization is automatically carried out.
      */
      bool ProbabilitiesSet(vector<Probability> & values);

      /// Get the index of the most probable state
      int MostProbableStateGet();

    };

} // End namespace

#endif // DISCRETEPDF_H
