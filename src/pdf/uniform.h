// $Id: uniform.h tdelaet$
// Copyright (C) 2007 Tinne De Laet <first dot last at gmail dot com>
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
#ifndef UNIFORM_H
#define UNIFORM_H

#include "pdf.h"

namespace BFL
{
  /// Class representing uniform density
  class Uniform: public Pdf<MatrixWrapper::ColumnVector>
    {
    private:
      /// Lower border of the uniform distribution
      MatrixWrapper::ColumnVector _Lower;
      /// Upper border of the uniform distribution
      MatrixWrapper::ColumnVector _Higher;
      /// Height of the uniform distribution
      double _Height; //the height of the uniform distribution

      // variables to avoid allocation on the heap during resampling
      mutable ColumnVector _samples;

    public:
      /// Constructor
      /**
	 @param Center center of the uniform distribution
	 @param Width width of the uniform distribution
      */
      Uniform (const MatrixWrapper::ColumnVector& Center, const MatrixWrapper::ColumnVector& Width);

      /// constructor with only dimensions or nothing
      Uniform (int dimension = 0);

      /// Default Copy Constructor will do

      /// Destructor
      virtual ~Uniform();

      /// output stream for Uniform distribution
      friend std::ostream& operator<< (std::ostream& os, const Uniform& u);

      ///Clone function
      virtual Uniform* Clone() const;

      // Redefinition of pure virtuals
      virtual Probability ProbabilityGet(const MatrixWrapper::ColumnVector& input) const;
      bool SampleFrom (vector<Sample<MatrixWrapper::ColumnVector> >& list_samples,
		       const int num_samples,
		       int method=DEFAULT,
		       void * args=NULL) const;
      virtual bool SampleFrom (Sample<MatrixWrapper::ColumnVector>& one_sample, int method=DEFAULT, void * args=NULL) const;

      /// Get the center of the uniform
      /** Get the center of the uniform
      */
      virtual MatrixWrapper::ColumnVector CenterGet() const;

      /// Get the Width of the uniform distribution
      /** Get the Width of the uniform distribution
      */
      virtual MatrixWrapper::ColumnVector WidthGet() const;

      /// Set the center and width of the uniform
      /** Set the center and width of the uniform
	  @param center The new center of uniform distribution
	  @param width The new width of the uniform distribution
      */
      void UniformSet (const MatrixWrapper::ColumnVector& center, const MatrixWrapper::ColumnVector& width);

    };

} // end namespace
#endif
