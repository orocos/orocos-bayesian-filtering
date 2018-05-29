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
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
//
#ifndef GAUSSIAN_H
#define GAUSSIAN_H

#include "pdf.h"

namespace BFL
{
  /// Class representing Gaussian (or normal density)
  class Gaussian: public Pdf<MatrixWrapper::ColumnVector>
    {
    private:
      MatrixWrapper::ColumnVector _Mu;
      MatrixWrapper::SymmetricMatrix _Sigma;

      // variables to avoid recalculation of inverse
      mutable bool _Sigma_changed;
      mutable MatrixWrapper::SymmetricMatrix _Sigma_inverse;
      mutable double _sqrt_pow;
      mutable ColumnVector _diff; //needed in probabilityGet
      mutable ColumnVector _tempColumn; //needed in probabilityGet
      // variables to avoid allocation on the heap during resampling
      mutable ColumnVector _samples;
      mutable ColumnVector _sampleValue;
      mutable Matrix _Low_triangle;

    public:
      /// Constructor
      /**
	 @param Mu Mean Vector of the Gaussian
	 @param Sigma Covariance Matrix of the Gaussian
      */
      Gaussian (const MatrixWrapper::ColumnVector& Mu, const MatrixWrapper::SymmetricMatrix& Sigma);

      /// constructor with only dimensions or nothing
      Gaussian (int dimension = 0);

      /// Default Copy Constructor will do

      /// Destructor
      virtual ~Gaussian();

      /// output stream for Gaussian
      friend std::ostream& operator<< (std::ostream& os, const Gaussian& g);

      ///Clone function
      virtual Gaussian* Clone() const;

      // Redefinition of pure virtuals
      virtual Probability ProbabilityGet(const MatrixWrapper::ColumnVector& input) const;
      bool SampleFrom (vector<Sample<MatrixWrapper::ColumnVector> >& list_samples,
		       const unsigned int num_samples,
		       const SampleMthd method=SampleMthd::DEFAULT,
		       void * args=NULL) const;
      virtual bool SampleFrom (Sample<MatrixWrapper::ColumnVector>& one_sample, const SampleMthd method=SampleMthd::DEFAULT, void * args=NULL) const;

      virtual MatrixWrapper::ColumnVector ExpectedValueGet() const;
      virtual MatrixWrapper::SymmetricMatrix CovarianceGet() const;
      virtual void DimensionSet(unsigned int dim);

      // For a Gaussian this should be possible
      /// Set the Expected Value
      /** Set the Expected Value
	  @param mu The new Expected Value
      */
      void ExpectedValueSet (const MatrixWrapper::ColumnVector& mu);

      /// Set the Covariance Matrix
      /** Set the Covariance Matrix
	  @param cov The new Covariance matrix
      */
      void CovarianceSet (const MatrixWrapper::SymmetricMatrix& cov);
    };

} // end namespace
#endif
