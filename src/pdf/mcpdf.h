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
#ifndef MCPDF_H
#define MCPDF_H

#include "pdf.h"
#include "../sample/weightedsample.h"
#include "../wrappers/rng/rng.h"
#include "../bfl_err.h"
#include <list>
#include <vector>
#include <cassert>

namespace BFL
{

  /// Monte Carlo Pdf: Sample based implementation of Pdf
  /** Class Monte Carlo Pdf:  This is a sample based representation of a
      Pdf P(x), which can both be continu or discrete
      @todo This class can and should be made far more efficient!!!
  */
  template <typename T> class MCPdf: public Pdf<T> 
    {
    protected:

      /// Sum of all weights: used for normalising purposes
      double _SumWeights; 
      /// STL-list containing the list of samples
      vector<WeightedSample<T> > _listOfSamples;
      /// STL-iterator
      //  vector<WeightedSample<T> >::iterator _it;
      /// STL-list containing the Cumulative PDF (for efficient sampling)
      vector<double> _CumPDF;
      /// STL-iterator for cumulative PDF list
      //  typename vector<double>::iterator CumPDFit;

  
      /// After updating weights, we have to recalculate the sum of weights
      void SumWeightsUpdate();
      /// Normalizing the weights
      bool NormalizeWeights();
      /// After updating weights, we have to update the cumPDF
      void CumPDFUpdate();

    public:
      /// Constructor
      /** @param num_samples the number of samples this pdf has
	  @param dimension the dimension of these samples (necessary for
	  if you want to avoid too much memory allocation at runtime)
      */
      MCPdf(unsigned int num_samples = 0, unsigned int dimension=0);
      /// destructor
      virtual ~MCPdf();
      /// copy constructor
      MCPdf(const MCPdf<T> &);

      // implemented virtual functions
      bool SampleFrom (Sample<T>& one_sample, int method = DEFAULT, void * args = NULL) const;
      bool SampleFrom (vector<Sample<T> >& list_samples, unsigned int num_samples, int method = DEFAULT, 
		       void * args = NULL) const;
      T ExpectedValueGet() const;
      MatrixWrapper::SymmetricMatrix CovarianceGet() const;


      /// Set number of samples
      /** @param num_samples the number of samples offcourse :-)
	  @see sample, weightedsample
      */
      void NumSamplesSet(unsigned int num_samples);
  

      /// Get number of samples
      /** @return the number of samples
       */
      unsigned int NumSamplesGet() const;

      /// Get one sample
      /** @return sample
          @param i the ith sample
       */
      const WeightedSample<T>& SampleGet(unsigned int i) const;
  
      /// Set the list of weighted samples
      /** @param list_of_samples an STL-list containing the list of all
	  weighted samples
      */
      bool ListOfSamplesSet(const vector<WeightedSample<T> > & list_of_samples);
      /// Overloading: Set the list of Samples (uniform weights)
      /** @param list_of_samples an STL-list containing the list of all
	  samples
      */
      bool ListOfSamplesSet(const vector<Sample<T> > & list_of_samples);

      /// Get the list of weighted samples
      /** @return an STL-list with the list of weighted samples
       */
      const vector<WeightedSample<T> > & ListOfSamplesGet() const;

      /// Update the list of samples (overloaded)
      /** @param list_of_samples the list of weighted samples
	  @pre list_of_samples must contain exactly as many elements as
	  this->NumSamplesGet() returns
      */
      bool ListOfSamplesUpdate(const vector<WeightedSample<T> > & list_of_samples);

      /// Update the list of samples (overloaded)
      /** @param list_of_samples the list of samples
	  @pre list_of_samples must contain exactly as many elements as
	  this->NumSamplesGet() returns
      */
      bool ListOfSamplesUpdate(const vector<Sample<T> > & list_of_samples);

      /// Add a sample to the list
      /** @param sample the sample to be added
	  @todo what's the best way to remove some samples?
      */
      // void SampleAdd(WeightedSample<T>  sample);

      /// Get the Cumulative Pdf
      /** @return a vector of doubles representing the CumulativePDF
       */
      vector<double> & CumulativePDFGet();
      
    };

  /////////////////////////////////////////////////////////////////
  // Template Code here
  /////////////////////////////////////////////////////////////////
  
  // Constructor
  template <typename T> MCPdf<T>::MCPdf(unsigned int num_samples, unsigned int dimension) : 
    Pdf<T>(dimension)
    {
      _SumWeights = 0;
      WeightedSample<T> my_sample(dimension);
      _listOfSamples.insert(_listOfSamples.begin(),num_samples,my_sample);
      _CumPDF.insert(_CumPDF.begin(),num_samples+1,0.0);
#ifdef __CONSTRUCTOR__
      // if (num_samples > 0)
      cout << "MCPDF Constructor: NumSamples = " << _listOfSamples.size()
	   << ", CumPDF Samples = " << _CumPDF.size()
	   << ", _SumWeights = " << _SumWeights << endl; 
#endif // __CONSTRUCTOR__
    }

  // Destructor
  template <typename T> 
    MCPdf<T>::~MCPdf()
    {
#ifdef __DESTRUCTOR__
      cout << "MCPDF::Destructor" << endl;
#endif // __DESTRUCTOR__
    }

  // Copy constructor
  template <typename T> 
    MCPdf<T>::MCPdf(const MCPdf & pdf) : Pdf<T>(pdf)
    {
      this->_listOfSamples = pdf._listOfSamples;
      this->_CumPDF = pdf._CumPDF;
      _SumWeights = pdf._SumWeights;
#ifdef __CONSTRUCTOR__
      cout << "MCPDF Copy Constructor: NumSamples = " << _listOfSamples.size()
	   << ", CumPDF Samples = " << _CumPDF.size()
	   << ", SumWeights = " << _SumWeights << endl;
#endif // __CONSTRUCTOR__
    }


  template <typename T> bool
    MCPdf<T>::SampleFrom (vector<Sample<T> >& list_samples,
			  const unsigned int numsamples,
			  int method, 
			  void * args) const
    {
      list_samples.resize(numsamples);
      switch(method)
	{
	case DEFAULT: // O(N log(N) efficiency)
	  {
	    return Pdf<T>::SampleFrom(list_samples, numsamples,method,args);
	  }
	case RIPLEY: // Only possible here ( O(N) efficiency )
	  /* See 
	     @Book{		  ripley87,
	     author	= {Ripley, Brian D.},
	     title		= {Stochastic Simulation},
	     publisher	= {John Wiley and Sons},
	     year		= {1987},
	     annote	= {ISBN 0271-6356, WBIB 1 519.245}
	     }
	  */
	  // GENERATE N ORDERED IID UNIFORM SAMPLES
	  {
	    std::vector<double> unif_samples(numsamples);
	    for ( unsigned int i = 0; i < numsamples ; i++)
	      unif_samples[i] = runif();

	    /* take n-th racine of u_N */
	    unif_samples[numsamples-1] = pow(unif_samples[numsamples-1], double (1.0/numsamples));
	    /* rescale other samples */
	    // only resample if more than one sample
	    if (numsamples > 1)
	      for ( int i = numsamples-2; i >= 0 ; i--)
		unif_samples[i] = pow(unif_samples[i],double (1.0/(i+1))) * unif_samples[i+1];

	    // CHECK WHERE THESE SAMPLES ARE IN _CUMPDF
	    unsigned int index = 0;
	    unsigned int size;
	    size = _listOfSamples.size();
	    typename vector<WeightedSample<T> >::const_iterator it = _listOfSamples.begin();
	    typename vector<double>::const_iterator CumPDFit = _CumPDF.begin();
	    typename vector<Sample<T> >::iterator sit = list_samples.begin();

	    for ( unsigned int i = 0; i < numsamples ; i++)
	      {
		while ( unif_samples[i] > *CumPDFit )
		  {
		    assert(index <= size);
		    index++; it++; CumPDFit++;
		  }
		it--; 
		*sit = *it; 
		it++; 
		sit++;
	      }
	    return true;
	  }
	default:
	  {
	    cerr << "MCPdf::Samplefrom(int, void *): No such sampling method" << endl;
	    return false;
	  }
	}
    }

  template <typename T> bool
    MCPdf<T>::SampleFrom(Sample<T>& one_sample, int method, void * args) const
    {
      switch(method)
	{
	case DEFAULT:
	  {
	    // Sample from univariate uniform rng between 0 and 1;
	    double unif_sample; unif_sample = runif();
	    // Compare where we should be:
	    unsigned int index = 0; 
	    unsigned int size; size = _listOfSamples.size();
	    typename vector<WeightedSample<T> >::const_iterator it;
	    it = _listOfSamples.begin(); 
	    typename vector<double>::const_iterator CumPDFit;
	    CumPDFit = _CumPDF.begin();
	    
	    while ( unif_sample > *CumPDFit )
	      {
		// check for internal error
		assert(index <= size);
		index++; it++; CumPDFit++;
	      }
	    it--;
	    one_sample = *it;
	    return true;
	  }
	default:
	  {
	    cerr << "MCPdf::Samplefrom(int, void *): No such sampling method" << endl;
	    return false;
	  }
	}
    }


  template <typename T> unsigned int MCPdf<T>::NumSamplesGet() const 
    {
      return _listOfSamples.size();
    }

  template <typename T> const WeightedSample<T>& 
    MCPdf<T>::SampleGet(unsigned int i) const
    {
      assert(i < NumSamplesGet());
      return _listOfSamples[i];
    }

  // Get and set number of samples
  template <typename T> void 
    MCPdf<T>::NumSamplesSet(unsigned int num_samples)
    {
#ifdef __MCPDF_DEBUG__
      cout << "MCPDF::NumSamplesSet BEFORE:  NumSamples " << _listOfSamples.size() << endl;
      cout << "MCPDF::NumSamplesSet BEFORE:  CumPDF Samples " << _CumPDF.size() << endl;
#endif // __MCPDF_DEBUG__
      unsigned int ns = num_samples;
      unsigned int size = _listOfSamples.size();
      static typename vector<double>::iterator CumPDFit;
      static typename vector<WeightedSample<T> >::iterator it;
      if (size < ns) // Add samples
	{
	  WeightedSample<T> ws;
	  _listOfSamples.insert(_listOfSamples.end(),(ns - size),ws);
	  _CumPDF.insert(_CumPDF.end(),(ns - size),0.0);
	}
      else if (size > ns) // Delete some samples
	{
	  it = _listOfSamples.begin(); CumPDFit = _CumPDF.begin();
	  for ( unsigned int index = 0; index < (size-ns); index++ )
	    {
	      it = _listOfSamples.erase(it);
	      CumPDFit = _CumPDF.erase(CumPDFit);
	    }
#ifdef __MCPDF_DEBUG__
	  cout << "MCPDF::NumSamplesSet: WARNING DELETING SAMPLES!!" << endl;
#endif // __MCPDF_DEBUG__
	}
      else {;} // Do nothing (number of samples are equal)
#ifdef __MCPDF_DEBUG__
      cout << "MCPDF::NumSamplesSet: Setting NumSamples to " << _listOfSamples.size() << endl;
      cout << "MCPDF::NumSamplesSet: Setting CumPDF Samples to " << _CumPDF.size() << endl;
#endif // __MCPDF_DEBUG__
    }


  // Get and set the list of samples
  template <typename T> bool
    MCPdf<T>::ListOfSamplesSet(const vector<WeightedSample<T> > & los)
    {
      // Allocate necessary memory
      this->NumSamplesSet(los.size());
      _listOfSamples = los;
#ifdef __MCPDF_DEBUG__
      cout << "MCPDF::ListOfSamplesSet: NumSamples = " << ListOfSamples.size() << endl;
#endif // __MCPDF_DEBUG__
      return this->NormalizeWeights();
    }


  template <typename T> bool
    MCPdf<T>::ListOfSamplesSet(const vector<Sample<T> > & los)
    {
      unsigned int numsamples = los.size();
      typename vector<Sample<T> >::const_iterator lit; lit=los.begin();
      static typename vector<WeightedSample<T> >::iterator it;
      // Allocate necessary memory
      this->NumSamplesSet(numsamples);
      // Update the list of samples
      for ( it = _listOfSamples.begin() ; it != _listOfSamples.end() ; it++ )
	{
	  *it = *lit; ; 
	  it->WeightSet(1.0/numsamples);
	  lit++;
	}
      _SumWeights = 1.0;
      // Update Cum PDF
      this->CumPDFUpdate();

#ifdef __MCPDF_DEBUG__
      cout << "MCPDF ListOfSamplesSet: NumSamples = " << _listOfSamples.size()
	   << " SumWeights = " << _SumWeights << endl;
#endif // __MCPDF_DEBUG__

      return true;
    }

  template <typename T> const vector<WeightedSample<T> > & 
    MCPdf<T>::ListOfSamplesGet() const
    {
      return _listOfSamples;
    }


  template <typename T> bool
    MCPdf<T>::ListOfSamplesUpdate(const vector<WeightedSample<T> > & los)
    {
      assert (los.size() == _listOfSamples.size());
      if (los.size() != 0)
	{
	  _listOfSamples = los;
	  return this->NormalizeWeights();
	}
      return true;
    }

  template <typename T> bool
    MCPdf<T>::ListOfSamplesUpdate(const vector<Sample<T> > & los)
    {
      unsigned int numsamples = _listOfSamples.size();
      if ((numsamples = los.size()) == _listOfSamples.size())
	{
	  assert (numsamples != 0);
	  typename vector<Sample<T> >::const_iterator lit; lit=los.begin();
	  static typename vector<WeightedSample<T> >::iterator it;
	  // Allocate necessary memory
	  this->NumSamplesSet(numsamples);
	  // Update the sumweights
	  for ( it = _listOfSamples.begin() ; it != _listOfSamples.end() ; it++ )
	    {
	      *it = *lit; ; 
	      it->WeightSet(1.0/numsamples);
	      lit++;
	    }
	  _SumWeights = 1.0;
	  this->CumPDFUpdate();
	}
      return true;
    }


  template <typename T> void 
    MCPdf<T>::SumWeightsUpdate()
    {
      double SumOfWeights = 0.0; 
      double current_weight;
      static typename vector<WeightedSample<T> >::iterator it;
      for ( it = _listOfSamples.begin() ; it != _listOfSamples.end() ; it++ )
	{
	  current_weight = it->WeightGet();
	  SumOfWeights += current_weight;
	}
      assert(SumOfWeights > 0);
      this->_SumWeights = SumOfWeights;
      //  CumPDFUpdate();
#ifdef __MCPDF_DEBUG__
      cout << "MCPDF::SumWeightsUpdate: SumWeights = " << _SumWeights << endl;
#endif // __MCPDF_DEBUG__
    }

  template <typename T> bool
    MCPdf<T>::NormalizeWeights()
    {
      static typename vector<WeightedSample<T> >::iterator it;
      this->SumWeightsUpdate();
      if (this->_SumWeights != 0)
	{
	  for ( it = _listOfSamples.begin() ; it != _listOfSamples.end() ; it++ )
	    {
	      it->WeightSet(it->WeightGet() / _SumWeights);
	    }
	  this->_SumWeights = 1.0;
	  this->CumPDFUpdate();
	  return true;
	}
      // all samples have weight = 0
      else
	return false;
    }


  template <typename T> void 
    MCPdf<T>::CumPDFUpdate()
    {
      double CumSum=0.0; 
      static typename vector<double>::iterator CumPDFit;
      static typename vector<WeightedSample<T> >::iterator it;
      CumPDFit = _CumPDF.begin(); *CumPDFit = 0.0;

      // Calculate the Cumulative PDF
      for ( it = _listOfSamples.begin() ; it != _listOfSamples.end() ; it++ )
	{
	  CumPDFit++;
	  // Calculate the __normalised__ Cumulative sum!!!
	  CumSum += ( it->WeightGet() / _SumWeights);
	  *CumPDFit = CumSum;
	}
    }


  template <typename T>
    T MCPdf<T>::ExpectedValueGet (  ) const
    {
      cerr << "MCPDF ExpectedValueGet: not implemented for the template parameters you use." 
	   << endl << "Use template specialization as shown in mcpdf.cpp " << endl; 

      assert(0);
      T result;
      return result;
    }


  template <typename T>
    MatrixWrapper::SymmetricMatrix MCPdf<T>::CovarianceGet (  ) const
    {
      cerr << "MCPDF CovarianceGet: not implemented for the template parameters you use." 
	   << endl << "Use template specialization as shown in mcpdf.cpp " << endl; 

      assert(0);
      MatrixWrapper::SymmetricMatrix result;
      return result;
    }



  template <typename T>
    vector<double> & MCPdf<T>::CumulativePDFGet()
    {
      return _CumPDF;
    }


  
} // End namespace BFL

#include "mcpdf.cpp"

#endif
