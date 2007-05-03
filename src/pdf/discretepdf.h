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
      /// Pointer to the discrete PDF-values
      MatrixWrapper::ColumnVector *_Values_p;
      /// Total sum of weighs associated with the variables
      double _SumWeights;
      /// Update the Sum of weights (eg. after setting a weight)
      void SumWeightsUpdate();

      /// STL-vector containing the Cumulative PDF (for efficient sampling)
      vector<double> _CumPDF;
      /// After updating weights, we have to update the cumPDF
      void CumPDFUpdate();
    public:
      /// Constructor (dimension = number of classes) An equal probability is set for all classes 
      /** @param num_states number of different classes or states
       */
      DiscretePdf(int num_states=0);

      /// Copy Constructor
      DiscretePdf(const DiscretePdf &);

      /// Destructor
      virtual ~DiscretePdf();

      /// Implementation of virtual base class method
      Probability ProbabilityGet(const int& input) const;
      /// Only relevant for discrete Pdf's
      bool ProbabilitySet(int input, Probability a);

      bool SampleFrom (vector<Sample<int> >& list_samples,
		       const int num_samples,
		       int method = DEFAULT, 
		       void * args = NULL) const;
      bool SampleFrom (Sample<int>& one_sample, int method = DEFAULT, void * args = NULL) const;

      /// Get all probabilities
      MatrixWrapper::ColumnVector ProbabilitiesGet() const;
      /// Set all probabilities
      void ProbabilitiesSet(MatrixWrapper::ColumnVector & values);

    };

} // End namespace

#endif // DISCRETEPDF_H
