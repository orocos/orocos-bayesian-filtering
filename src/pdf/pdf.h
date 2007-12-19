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
#ifndef PDF_H
#define PDF_H

#include <vector>
#include <iostream>
#include "../wrappers/matrix/vector_wrapper.h"
#include "../wrappers/matrix/matrix_wrapper.h"
#include "../sample/sample.h"
#include "../bfl_err.h"
#include "../bfl_constants.h"


namespace BFL
{
  using namespace std;

  // Defines for different sampling methods
#define DEFAULT 0 // Default sampling method, must be valid for every PDF!!
#define BOXMULLER 1
#define CHOLESKY 2
#define RIPLEY 3 // For efficient sampling from discrete/mcpdfs

  /// Class PDF: Virtual Base class representing Probability Density Functions
  template <typename T> class Pdf 
    {

    public:
      /// Constructor
      /** @param dimension int representing the number of rows of the state
	  vector for the continuous case, or the discrete number of states
	  for the discrete case
      */
      Pdf(unsigned int dimension=0);
  
      // Default Copy Constructor will do the job
  
      /// Destructor
      virtual ~Pdf();

      /// Draw multiple samples from the Pdf (overloaded)
      /** @param list_samples list of samples that will contain result of sampling
          @param num_samples Number of Samples to be drawn (iid)
	  @param method Sampling method to be used.  Each sampling method
	  is currently represented by a #define statement, eg. 
	  #define BOXMULLER 1
	  @param args Pointer to a struct representing extra sample
	  arguments.
	  "Sample Arguments" can be anything (the number of steps a
	  gibbs-iterator should take, the interval width in MCMC, ... (or
	  nothing), so it is hard to give a meaning to what exactly
	  Sample Arguments should represent...
	  @todo replace the C-call "void * args" by a more object-oriented
	  structure: Perhaps something like
	  virtual Sample * Sample (const int num_samples,class Sampler)
	  @bug Sometimes the compiler doesn't know which method to choose!
      */
      virtual bool SampleFrom (vector<Sample<T> >& list_samples,
			       const unsigned int num_samples,
			       int method = DEFAULT, 
			       void * args = NULL) const;

      /// Draw 1 sample from the Pdf:
      /** There's no need to create a list for only 1 sample!
	  @param one_sample sample that will contain result of sampling
	  @param method Sampling method to be used.  Each sampling method
	  is currently represented by a #define statement, eg. 
	  #define BOXMULLER 1
	  @param args Pointer to a struct representing extra sample
	  arguments
	  @see SampleFrom()
	  @bug Sometimes the compiler doesn't know which method to choose!
      */
      virtual bool  SampleFrom (Sample<T>& one_sample,
				int method = DEFAULT, 
				void * args = NULL) const;

      /// Get the probability of a certain argument
      /** @param input T argument of the Pdf
	  @return the probability value of the argument
      */
      virtual Probability ProbabilityGet(const T& input) const;

      /// Get the dimension of the argument
      /** @return the dimension of the argument
       */
      unsigned int DimensionGet() const;
  
      /// Set the dimension of the argument
      /** @param dim the dimension
       */
      virtual void DimensionSet(unsigned int dim);

      /// Get the expected value E[x] of the pdf 
      /** Get low order statistic (Expected Value) of this AnalyticPdf
	  @return The Expected Value of the Pdf (a ColumnVector with
	  DIMENSION rows)
	  @note No set functions here!  This can be useful for analytic
	  functions, but not for sample based representations!
	  @note For certain discrete Pdfs, this function has no
	  meaning, what is the average between yes and no?
      */
      virtual T ExpectedValueGet() const;

      /// Get the Covariance Matrix E[(x - E[x])^2] of the Analytic pdf 
      /** Get first order statistic (Covariance) of this AnalyticPdf
	  @return The Covariance of the Pdf (a SymmetricMatrix of dim
	  DIMENSION)
	  @todo extend this more general to n-th order statistic
	  @bug Discrete pdfs should not be able to use this!
      */
      virtual MatrixWrapper::SymmetricMatrix CovarianceGet() const;
  
    protected:
      /// Dimension of the argument x of P(x | ...).
      /** In case of a discrete pdf, this is the number of discrete
	  classes
      */
      unsigned int _dimension;

    };

template<typename T> 
Pdf<T>::Pdf(unsigned int dim)
{
  assert((int)dim >= 0);

  _dimension = dim;
#ifdef __CONSTRUCTOR__
  cout << "Pdf constructor" << endl;
#endif
}

template<typename T> 
Pdf<T>::~Pdf()
{
#ifdef __DESTRUCTOR__
  cout << "Pdf destructor" << endl;
#endif
}

template<typename T> inline unsigned int 
Pdf<T>::DimensionGet () const
{
  return _dimension;
}

template<typename T> void 
Pdf<T>::DimensionSet ( unsigned int dim )
{
  assert((int)dim >= 0);  
  _dimension = dim;
}

template<typename T> bool
Pdf<T>::SampleFrom (vector<Sample<T> >& list_samples, 
		    const unsigned int num_samples, 
		    int method, 
		    void * args) const
{
  list_samples.resize(num_samples);
  typename vector<Sample<T> >::iterator sample_it;
  for (sample_it = list_samples.begin(); sample_it != list_samples.end() ; sample_it++)
    if (!this->SampleFrom(*sample_it, method,args)) 
      return false;

  return true;
}

template<typename T> bool
Pdf<T>::SampleFrom (Sample<T>& one_sample,
		    int method, 
		    void * args) const
{
  cerr << "Error: The SampleFrom function was called, but you didn't implement it!\n";
  exit(-BFL_ERRMISUSE);
  return false;
}

template<typename T> Probability 
Pdf<T>::ProbabilityGet (const T& input) const
{
  cerr << "Error: The ProbabilityGet function was called, but you didn't implement it!\n";
  exit(-BFL_ERRMISUSE);
  return 1;
}

template<typename T>  T 
Pdf<T>::ExpectedValueGet () const
{
  cerr << "Error: The ExpectedValueGet function was called, but you didn't implement it!\n";
  exit(-BFL_ERRMISUSE);
  T t;
  return t;
}


template<typename T> MatrixWrapper::SymmetricMatrix 
Pdf<T>::CovarianceGet () const
{
  cerr << "Error: The CovarianceGet function was called, but you didn't implement it!\n";
  exit(-BFL_ERRMISUSE);
  MatrixWrapper::SymmetricMatrix m;
  return m;
}



} // End namespace
#endif
