// $Id$
// Copyright (C) 2003 Klaas Gadeyne <first dot last at gmail dot com>
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
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
//
#include "discreteconditionalpdf.h"
#include "../wrappers/rng/rng.h"
#include <vector>
#include <iostream>

namespace BFL
{
  using namespace std;

  DiscreteConditionalPdf::DiscreteConditionalPdf(int num_states,
						 int num_conditional_arguments,
						 int cond_arg_dimensions[])
    : ConditionalPdf<int,int>(1,num_conditional_arguments)
    , _num_states(num_states)
    , _probs(num_states)
    , _valuelist(num_states+1)
  {
    _cond_arg_dims_p = new int[num_conditional_arguments];
    int total_dim = 1;
    for ( int arg = 0 ; arg < num_conditional_arguments; arg++)
      {
	_cond_arg_dims_p[arg] = cond_arg_dimensions[arg];
	total_dim *= cond_arg_dimensions[arg];
      }
    _total_dimension = total_dim * num_states;
#ifdef __CONSTRUCTOR__
    cout << "DCPdf Constructor: total dimension allocated: "
	 << _total_dimension << endl;
#endif // __CONSTRUCTOR__
    _probability_p = new double[_total_dimension];
  }

  DiscreteConditionalPdf::~DiscreteConditionalPdf()
  {
#ifdef __DCPDFDEBUG__
    cout << "DCPdf Destructor:" << endl;
#endif // __DCPDFDEBUG__
    //delete _probability_p;
    //delete _cond_arg_dims_p;
  }


  DiscreteConditionalPdf::DiscreteConditionalPdf(const DiscreteConditionalPdf & pdf)
    : ConditionalPdf<int,int>(pdf)
    ,_num_states(pdf.NumStatesGet())
    , _probs(pdf.NumStatesGet())
    , _valuelist(pdf.NumStatesGet()+1)
  {
    _cond_arg_dims_p = new int[pdf.NumConditionalArgumentsGet()];
    int total_dim = 1;
    for ( unsigned int arg = 0 ; arg < NumConditionalArgumentsGet(); arg++)
      {
	_cond_arg_dims_p[arg] = pdf._cond_arg_dims_p[arg];
	total_dim *= _cond_arg_dims_p[arg];
      }
    total_dim *= _num_states;
    _total_dimension = total_dim;
    _probability_p = new double[total_dim];
    for (int prob = 0 ; prob < total_dim ; prob++)
      {
	_probability_p[prob] = pdf._probability_p[prob];
      }
  }

  //Clone function
  DiscreteConditionalPdf* DiscreteConditionalPdf::Clone() const
  {
      return new DiscreteConditionalPdf(*this);
  }

  // Get the number of discrete states
  unsigned int DiscreteConditionalPdf::NumStatesGet() const
  {
     return _num_states;
  }

  // Calculate index (used by ProbabilityGet and ProbabilitySet)
  int DiscreteConditionalPdf::IndexGet(const int& input,
				       const std::vector<int>& condargs) const
  {
    int index = 0;
    int blocksize = 1;

    // The first hyperdimension is that of input itself
    index += input * blocksize;
    blocksize *= NumStatesGet();
    // The other ons are those of the conditional args
    for (unsigned int arg = 0 ; arg < NumConditionalArgumentsGet() ; arg++ )
      {
	index += condargs[arg] * blocksize;
	blocksize *= (this->_cond_arg_dims_p)[arg];
      }
#ifdef __INDEXDEBUG__
    cout << "DCPdf::IndexGet -> Index = " << index << endl;
#endif // __INDEXDEBUG__
    return index;
  }



  Probability DiscreteConditionalPdf::ProbabilityGet(const int& input) const
  {
    unsigned int index = IndexGet(input, ConditionalArgumentsGet());
    double prob = (this->_probability_p)[index];
    return prob;
  }

  // Typical for discrete Pdf's
  void DiscreteConditionalPdf::ProbabilitySet(const double& prob,
					      const int& input,
					      const std::vector<int>& condargs) const
  {
    int index = this->IndexGet(input, condargs);
    _probability_p[index] = prob;
  }

  bool DiscreteConditionalPdf::SampleFrom(Sample<int>& one_sample, const SampleMthd method, void * args) const
  {
    // Get the elements of which to sample from
    int startindex = IndexGet(0,ConditionalArgumentsGet());
    double SumWeights = 0.0; double CumSum=0.0;
    unsigned int index;

    for ( index = 0; index < NumStatesGet() ; index++ )
      {
	_probs[index] = _probability_p[startindex+index];
	CumSum += _probs[index];
#ifdef __DCPDFDEBUG__
#define TABWIDTH 10
	cout << setw(TABWIDTH) << _probs[index];
#endif // __DCPDFDEBUG__
      }
    SumWeights = CumSum;
    _valuelist[0] = 0.0; CumSum = 0.0;
    for ( index = 1; index <= NumStatesGet() ; index++ )
      {
	CumSum += _probs[index-1]/SumWeights;
	_valuelist[index] = CumSum;
      }
    // Check if last element of valuelist is +- 1
    assert ( (_valuelist[NumStatesGet()] >= 1.0 - NUMERIC_PRECISION) &&
	     (_valuelist[NumStatesGet()] <= 1.0 + NUMERIC_PRECISION) );

    _valuelist[NumStatesGet()]=1;

    // Sample from univariate uniform rng between 0 and 1;
    double unif_sample; unif_sample = runif();
    // Compare where we should be: THIS CAN BE MADE FASTER: TODO
    index = 0;
    while ( unif_sample > _valuelist[index] )
      {
	assert(index <= NumStatesGet());
	index++;
      }
    one_sample.ValueSet(index-1);
    return true;
  }


  bool
  DiscreteConditionalPdf::SampleFrom(vector<Sample<int> >& list_samples,
				     unsigned int num_samples, const SampleMthd method,
				     void * args) const
  {
    list_samples.resize(num_samples); // will break real-timeness if list_samples.size()!=num_samples
    return Pdf<int>::SampleFrom(list_samples, num_samples,method,args);
  }


} // End namespace BFL
