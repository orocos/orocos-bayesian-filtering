// $Id$
// Copyright (C) 2002 Klaas Gadeyne <first dot last at gmail dot com>
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
      /// Pointer to the discrete PDF-values
      MatrixWrapper::ColumnVector *_Values_p;
      /// Total sum of weighs associated with the variables
      double _SumWeights;
      /// Update the Sum of weights (eg. after setting a weight)
      bool SumWeightsUpdate();

      /// STL-vector containing the Cumulative PDF (for efficient sampling)
      vector<double> _CumPDF;
      /// After updating weights, we have to update the cumPDF
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
      Probability ProbabilityGet(const unsigned int& input) const;
      /// Only relevant for discrete Pdf's
      bool ProbabilitySet(unsigned int input, Probability a);

      bool SampleFrom (vector<Sample<int> >& list_samples,
		       const unsigned int num_samples,
		       int method = DEFAULT, 
		       void * args = NULL) const;
      bool SampleFrom (Sample<int>& one_sample, int method = DEFAULT, void * args = NULL) const;

      /// Get all probabilities
      MatrixWrapper::ColumnVector ProbabilitiesGet() const;
      /// Set all probabilities
      bool ProbabilitiesSet(MatrixWrapper::ColumnVector & values);

    };

} // End namespace

#endif // DISCRETEPDF_H
