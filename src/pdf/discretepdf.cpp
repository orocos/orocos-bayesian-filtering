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
#include "discretepdf.h"
#include "../bfl_err.h"
#include "../wrappers/rng/rng.h"
#include <vector>
#include <iostream>


namespace BFL
{
  using namespace std;
  using namespace MatrixWrapper;
  

  DiscretePdf::DiscretePdf(unsigned int num_states): Pdf<int>(1)
        ,_num_states(num_states)
  {
    //discrete pdf has dimension 1
    _Values_p = new vector<Probability>(num_states);
    for (int i=0; i<NumStatesGet();i++)
    {
        (*_Values_p)[i] = (Probability)(1.0/NumStatesGet());
    }
    _CumPDF.insert(_CumPDF.begin(),num_states+1,0.0);
    CumPDFUpdate();
#ifdef __CONSTRUCTOR__
    cout << "DiscretePdf constructor\n";
#endif // __CONSTRUCTOR__
  }

  DiscretePdf::DiscretePdf(const DiscretePdf & my_dpdf):Pdf<int>(my_dpdf)
        ,_num_states(my_dpdf.NumStatesGet())
  { 
    _Values_p = new vector<Probability>(this->NumStatesGet());
    (*_Values_p) = my_dpdf.ProbabilitiesGet();
    _CumPDF.insert(_CumPDF.begin(),NumStatesGet()+1,0.0);
    CumPDFUpdate();
#ifdef __CONSTRUCTOR__
    cout << "DiscretePdf copy constructor\n";
#endif // __CONSTRUCTOR__
  }

  DiscretePdf::~DiscretePdf()
  {
#ifdef __CONSTRUCTOR__
    cout << "DiscretePdf destructor\n";
#endif
    // Release memory!
    delete _Values_p;
  }
  
  unsigned int DiscretePdf::NumStatesGet()const
  {
    return _num_states;
  }


  Probability DiscretePdf::ProbabilityGet(const unsigned int& state) const
  {
    assert((int)state >= 0 && state < NumStatesGet());

    return (*_Values_p)[state];
  }

  bool DiscretePdf::ProbabilitySet(unsigned int state, Probability a) 
  {
    assert((int)state >= 0 && state < NumStatesGet());
    assert(a<=1);
    
    // renormalize other probabilities such that sum of probabilities will be
    // one after the probability of state is set to a
    // This should keep the probabilities normalized
    Probability old_prob_state = ProbabilityGet(state);
    if (old_prob_state!=1) {
        double normalization_factor = (1-a)/(1-old_prob_state);
        for (int i=0; i<NumStatesGet();i++)
        {
            (*_Values_p)[i] = (Probability)( (double)( (*_Values_p)[i] )* normalization_factor);
        }
    }
    else{
        for (int i=0; i<NumStatesGet();i++)
        {
            (*_Values_p)[i] = (Probability)( (1-a)/(NumStatesGet()-1));
        }
    }
    (*_Values_p)[state] = a;
    return CumPDFUpdate();
  }

  vector<Probability> DiscretePdf::ProbabilitiesGet() const
  {
    return *_Values_p;
  }

  bool DiscretePdf::ProbabilitiesSet(vector<Probability> & v)
  {
    assert(v.size() == NumStatesGet());

    *_Values_p = v;
    //normalize the probabilities and update the cumulative sum
    return (NormalizeProbs() && CumPDFUpdate());
  }

  // For optimal performance!
  bool
  DiscretePdf::SampleFrom (vector<Sample<int> >& list_samples,
			   const unsigned int num_samples,
			   int method, 
			   void * args) const
  {
    switch(method)
    {
      case DEFAULT: // O(N log(N) efficiency)
	  return Pdf<int>::SampleFrom(list_samples, num_samples,method,args);

      case RIPLEY: // See mcpdf.cpp for more explanation
	  {
	    list_samples.resize(num_samples);
	    // GENERATE N ORDERED IID UNIFORM SAMPLES
	    std::vector<double> unif_samples(num_samples);
	    for ( unsigned int i = 0; i < num_samples ; i++)
	      unif_samples[i] = runif();

	    /* take n-th racine of u_N */
	    unif_samples[num_samples-1] = pow(unif_samples[num_samples-1],
	  				   double (1.0/num_samples));
	    /* rescale samples */
	    for ( int i = num_samples-2; i >= 0 ; i--)
	        unif_samples[i] = pow(unif_samples[i], double (1.0/(i+1))) * unif_samples[i+1];

	    // CHECK WHERE THESE SAMPLES ARE IN _CUMPDF
	    unsigned int index = 0;
	    unsigned int num_states = NumStatesGet();
	    vector<double>::const_iterator CumPDFit = _CumPDF.begin();
	    vector<Sample<int> >::iterator sit = list_samples.begin();

	    for ( unsigned int i = 0; i < num_samples ; i++)
	      {
	        while ( unif_samples[i] > *CumPDFit )
	        {
	  	// check for internal error
	  	assert(index <= num_states);
	  	index++; CumPDFit++;
	        }
          int a = index - 1;
	      sit->ValueSet(a); 
	      sit++;
	      }
	    return true;
	  }
      default:
	    cerr << "DiscretePdf::Samplefrom(int, void *): No such sampling method" << endl;
	    return false;
    }
  }



  bool DiscretePdf::SampleFrom (Sample<int>& one_sample, int method, void * args) const
  {
    switch(method)
      {
      case DEFAULT:
	{
	  // Sample from univariate uniform rng between 0 and 1;
	  double unif_sample; unif_sample = runif();
	  // Compare where we should be
	  unsigned int index = 0;
	  while ( unif_sample > _CumPDF[index] )
	    {
	      assert(index <= NumStatesGet());
	      index++;
	    }
	  int a = index - 1;
	  one_sample.ValueSet(a);
	  return true;
	}
      default:
	cerr << "DiscretePdf::Samplefrom(int, void *): No such sampling method" 
	     << endl;
	return false;
      }
  }

  bool DiscretePdf::NormalizeProbs()
  {
    double SumOfProbs = 0.0; 
    for ( unsigned int i = 0; i < NumStatesGet() ; i++){
        SumOfProbs += (*_Values_p)[i];
    } 
    if (SumOfProbs > 0){
      for ( unsigned int i = 0; i < NumStatesGet() ; i++){
          (*_Values_p)[i] = (Probability)( (double) ( (*_Values_p)[i]) /SumOfProbs);
      }
      return true;
    }
    else{
      cerr << "DiscretePdf::NormalizeProbs(): SumOfProbs = " << SumOfProbs << endl;
      return false;
    }
  }

  bool DiscretePdf::CumPDFUpdate()
  {
    // precondition: sum of probabilities should be 1
    double CumSum=0.0; 
    static vector<double>::iterator CumPDFit;
    CumPDFit = _CumPDF.begin(); 
    *CumPDFit = 0.0;

    // Calculate the Cumulative PDF
    for ( unsigned int i = 0; i < NumStatesGet(); i++)
    {
	    CumPDFit++;
	    // Calculate the __normalised__ Cumulative sum!!!
	    CumSum += ( (*_Values_p)[i] );
	    *CumPDFit = CumSum;
    }
    // Check if last element of valuelist is +- 1
    assert( (_CumPDF[NumStatesGet()] >= 1.0 - NUMERIC_PRECISION) &&
	    (_CumPDF[NumStatesGet()] <= 1.0 + NUMERIC_PRECISION) );

    _CumPDF[NumStatesGet()]=1;

    return true;
  } 

  
} // End namespace
