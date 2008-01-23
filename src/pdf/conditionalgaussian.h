// $Id$
// Copyright (C) 2002 Klaas Gadeyne <first dot last at gmail dot com>
// Copyright (C) 2008 Tinne De Laet <first dot last at mech dot kuleuven dot be>
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

#ifndef __CONDITIONAL_GAUSSIAN__
#define __CONDITIONAL_GAUSSIAN__

#include "conditionalpdf.h"

namespace BFL
{
  /// Abstract Class representing all Conditional gaussians
  /** This class inherits only from ConditionalPdf<ColumnVector,
      ColumnVector>.

      So this class represents all Pdf's of the type
      \f[ P ( A | B, C, D, ... ) \f] where
      \f[ \mu_A = f(B,C,D, ...) \f] and
      \f[ \Sigma_A = g(B,C,D, ...) \f] and

      \f[ A = N(\mu_A, \Sigma_A) \f]
      
      f and g are not necessarily analytical functions
  */
  class ConditionalGaussian : public ConditionalPdf<MatrixWrapper::ColumnVector, MatrixWrapper::ColumnVector>
    {
    public:
      /// Constructor
      /**
	 @param dim Dimension of state
	 @param num_conditional_arguments The number of conditional
	 arguments.
      */
      ConditionalGaussian(int dim = 0, int num_conditional_arguments=0);

      // Default Copy constructor will do

      /// Destructor
      virtual ~ConditionalGaussian();

      // implemented virtuals!
      virtual Probability ProbabilityGet(const MatrixWrapper::ColumnVector& input) const;
      virtual bool SampleFrom (Sample<MatrixWrapper::ColumnVector>& sample, int method=DEFAULT, void * args=NULL) const;
      virtual bool SampleFrom (std::vector<Sample<MatrixWrapper::ColumnVector> >& samples, const int num_samples,
			       int method=DEFAULT, void * args=NULL) const;

    protected:
      // variables to avoid allocation on the heap during sampling
      mutable ColumnVector _diff;
      mutable ColumnVector _Mu;
      mutable Matrix _Low_triangle;
      mutable ColumnVector _samples;
      mutable ColumnVector _SampleValue;

    };

} // End namespace BFL

#endif // __CONDITIONAL_GAUSSIAN__

