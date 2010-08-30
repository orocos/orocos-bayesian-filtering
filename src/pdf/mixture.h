// $Id: mixture.h 2009-01-21 tdelaet $
// Copyright (C)  2009 Tinne De Laet <first dot last at mech dot kuleuven dot be>
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
#ifndef MIXTURE_H
#define MIXTURE_H

#include "pdf.h"
#include "discretepdf.h"
#include "../wrappers/matrix/vector_wrapper.h"
#include "../wrappers/matrix/matrix_wrapper.h"
#include "../wrappers/rng/rng.h"
#include <vector>

namespace BFL
{
  /// Class representing a mixture of PDFs, the mixture can contain different
  //kind of pdfs but they should all share the same state space i.e. they all
  //have to be of the same template type Pdf<T> 
  /** This class is an instantation from the template class Pdf, with
      added methods to get a set the individual components and mixture weights
  */
  template <typename T> class Mixture : public Pdf<T> // inherit abstract_template_class
    {
    protected:
      /// The number of components
      unsigned int _numComponents;

      /// Pointer to the vector of mixture weights, the sum of the elements = 1 
      vector<Probability> *_componentWeights;

      /// Pointer to the vector of component pdfs 
      vector<Pdf<T>* > *_componentPdfs;

      /// Normalize the component weigths (eg. after setting a component weight)
      bool NormalizeWeights();

      /// Vector containing the cumulative component weights (for efficient sampling)
      // BEWARE: this vector has length _numComponents +1 (first element is
      // always 0, last element is always 1)
      vector<double> _cumWeights;

      /// Updates the cumWeights
      bool CumWeightsUpdate();

      /// Called when a the number of components=0 and if method is called which
      //requires number of components>0
      void TestNotInit() const;

    public:
      /// Constructor: An equal weight is set for all components
      /** @param dimension the dimension of the state space
       */
      Mixture(const unsigned int dimension=0);

      /// Constructor: An equal weights is set for all components
      /** @param componentVector vector of components (Pdf<T>*)
       */
      template <typename U> Mixture(const U &componentVector);

      /// Copy Constructor

      Mixture(const Mixture & my_mixture);

      /// Destructor
      virtual ~Mixture();

      ///Clone function
      virtual Mixture* Clone() const;

      /// Get the number of components
      unsigned int NumComponentsGet()const;

      /// Implementation of virtual base class method
      Probability ProbabilityGet(const T& state) const;

      bool SampleFrom (vector<Sample<T> >& list_samples,
		       const unsigned int num_samples,
		       int method = DEFAULT, 
		       void * args = NULL) const;

      bool SampleFrom (Sample<T>& one_sample, int method = DEFAULT, void * args = NULL) const;

      T ExpectedValueGet() const;

      MatrixWrapper::SymmetricMatrix CovarianceGet (  ) const;

      /// Get all component weights
      vector<Probability> WeightsGet() const;

      /// Get the component weight of component "componentNumber"
      /**  @param componentNumber number of the component 
         (must be >=0 and <_numComponents)
      */
      Probability WeightGet(unsigned int componentNumber) const;

      /// Set all component weights
      /**  @weights values vector<Probability> containing the new component weights . 
           The sum of the probabilities of this list is not required to be one
           since the normalization is automatically carried out. 
           The size of weights should be equal to the number of components
      */
      bool WeightsSet(vector<Probability> & weights);

      /// Function to change/set the weigth of a single component
      /** Changes the component weights such that AFTER normalization the
          weight of the component "componentNumber" is equal to the weight w 
           @param componentNumber number of the component of which the weight will be set
           (must be >=0 and <_numComponents)
           @param w Probability to which the weight of component
           "componentNumber" 
            will be set (must be <= 1)
      */

      bool WeightSet(unsigned int componentNumber, Probability w);

      /// Get the index of the most probable component, if a few component are
      //equally probable, the index of the first most probable one is returned
      int MostProbableComponentGet() const;

      /// Add a component pdf: THIS IS A NON-REALTIME OPERATION
       // the weight of the new component is set to 0 (except when the number of
       // components is zero, then the weight is set to 1)!!
      /** @param pdf Component pdf which will be added
      */
      bool AddComponent(Pdf<T>& pdf );

      /// Add a component pdf with weight w: THIS IS A NON-REALTIME OPERATION
      /** @param pdf Component pdf which will be added
         @param weight the weight of the new component
      */
      bool AddComponent(Pdf<T>& pdf, Probability w);

      /// Delete a component pdf: THIS IS A NON_REALTIME OPERATION
      /**  @param componentNumber number of the component which will be deleted
         (must be >=0 and <_numComponents)
      */
      bool DeleteComponent(unsigned int componentNumber );

      /// Get the vector of pointers to the component pdfs
      vector<Pdf<T>* > ComponentsGet() const;

      /// Get the pointer to the component pdf of component "componentNumber"
      /**  @param componentNumber number of the component 
         (must be >=0 and <_numComponents)
      */
      Pdf<T>*  ComponentGet(unsigned int componentNumber) const;
    };

/////////////////////////////////////////////////////////////////
// Template Code here
/////////////////////////////////////////////////////////////////

// Constructor
//TODO: is this usefull because pointers to components point to nothing!
template<typename T>
Mixture<T>::Mixture(const unsigned int dimension):
  Pdf<T>(dimension)
  , _numComponents(0)
  , _cumWeights(_numComponents+1)
  {
    //create pointer to vector of component weights
    _componentWeights = new vector<Probability>(this->NumComponentsGet());

    //create pointer to vector of pointers to pdfs
    _componentPdfs = new vector< Pdf<T>* >(NumComponentsGet());
#ifdef __CONSTRUCTOR__
    cout << "Mixture constructor\n";
#endif // __CONSTRUCTOR__
  }

// Constructor
template<typename T> template <typename U>
Mixture<T>::Mixture(const U &componentVector): Pdf<T>(componentVector[0]->DimensionGet() )
    , _numComponents(componentVector.size())
{
    //create pointer to vector of component weights
    _componentWeights = new vector<Probability>(NumComponentsGet());
    for (int i=0; i<NumComponentsGet();i++)
    {
        (*_componentWeights)[i] = (Probability)(1.0/NumComponentsGet());
    }
    _cumWeights.insert(_cumWeights.begin(),NumComponentsGet()+1,0.0);
    CumWeightsUpdate();

    //create pointer to vector of pointers to weights
    _componentPdfs = new vector< Pdf<T>* >(NumComponentsGet());
    //TODO: take copy or point to same???
    for (int i=0; i<NumComponentsGet();i++)
    {
        //TODO: will this call the constructor of e.g. Gaussian or always the
        //general one?
        (*_componentPdfs)[i] = (componentVector[i])->Clone();
    }
#ifdef __CONSTRUCTOR__
    cout << "Mixture constructor\n";
#endif // __CONSTRUCTOR__
}

template<typename T > 
Mixture<T>::Mixture(const Mixture & my_mixture):Pdf<T>(my_mixture.DimensionGet() )
        ,_numComponents(my_mixture.NumComponentsGet())
  {
    //create pointer to vector of component weights
    _componentWeights = new vector<Probability>(this->NumComponentsGet());
    (*_componentWeights) = my_mixture.WeightsGet();
    _cumWeights.insert(_cumWeights.begin(),NumComponentsGet()+1,0.0);
    CumWeightsUpdate();

    //create pointer to vector of pointers to weights
    _componentPdfs = new vector< Pdf<T>* >(NumComponentsGet());
    for (int i=0; i<NumComponentsGet();i++)
    {
        (*_componentPdfs)[i] = (my_mixture.ComponentGet(i))->Clone();
    }
#ifdef __CONSTRUCTOR__
    cout << "Mixture copy constructor\n";
#endif // __CONSTRUCTOR__
  }

template<typename T>
  Mixture<T>::~Mixture()
  {
#ifdef __CONSTRUCTOR__
    cout << "Mixture destructor\n";
#endif
    // Release memory!
    delete _componentWeights;
    for (int i=0; i<NumComponentsGet();i++)
    {
        delete (*_componentPdfs)[i];
    }
    delete _componentPdfs;
  }

template<typename T>
  Mixture<T>* Mixture<T>::Clone() const
  {
      return new Mixture(*this);
  }

template<typename T>
  unsigned int Mixture<T>::NumComponentsGet() const
  {
    return _numComponents;
  }

template<typename T>
  Probability Mixture<T>::ProbabilityGet(const T& state) const
  {
    TestNotInit(); 
    Probability prob(0.0);
    for (int i=0; i<NumComponentsGet();i++)
    {
        prob= prob + (*_componentWeights)[i] * (*_componentPdfs)[i]->ProbabilityGet(state);
    }
    return prob;
  }

template<typename T>
  bool Mixture<T>::SampleFrom (vector<Sample<T> >& list_samples,
			   const unsigned int num_samples,
			   int method,
			   void * args) const
  {
    TestNotInit(); 
    switch(method)
    {
      case DEFAULT: // O(N log(N) efficiency)
	  return Pdf<T>::SampleFrom(list_samples, num_samples,method,args);

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

	    // CHECK WHERE THESE SAMPLES ARE IN _CUMWEIGHTS
	    unsigned int index = 0;
	    unsigned int num_states = NumComponentsGet();
	    vector<double>::const_iterator CumPDFit = _cumWeights.begin();
	    typename vector<Sample<T> >::iterator sit = list_samples.begin();

	    for ( unsigned int i = 0; i < num_samples ; i++)
	      {
	        while ( unif_samples[i] > *CumPDFit )
	        {
	  	// check for internal error
	  	assert(index <= num_states);
	  	index++; CumPDFit++;
	        }
          // index-1 is a sample of the discrete pdf of the mixture weights
          // get a sample from the index-1'th mixture component
	      (*_componentPdfs)[index-1]->SampleFrom(*sit,method,args);
	      sit++;
	      }
	    return true;
	  }
      default:
	    cerr << "Mixture::Samplefrom(T, void *): No such sampling method" << endl;
	    return false;
    }
  }
template<typename T>
  bool Mixture<T>::SampleFrom (Sample<T>& one_sample, int method, void * args) const
  {
    TestNotInit(); 
    switch(method)
      {
      case DEFAULT:
	{
	  // Sample from univariate uniform rng between 0 and 1;
	  double unif_sample; unif_sample = runif();
	  // Compare where we should be
	  unsigned int index = 0;
	  while ( unif_sample > _cumWeights[index] )
	    {
	      assert(index <= NumComponentsGet());
	      index++;
	    }
         // index-1 is a sample of the discrete pdf of the mixture weights
      // get a sample from the index-1'th mixture component
	  (*_componentPdfs)[index-1]->SampleFrom(one_sample,method,args);
	  return true;
	}
      default:
	cerr << "Mixture::Samplefrom(T, void *): No such sampling method"
	     << endl;
	return false;
      }
  }

template<typename T>
  T Mixture<T>::ExpectedValueGet() const
  {
    cerr << "Mixture ExpectedValueGet: not implemented for the template parameters you use." 
	 << endl << "Use template specialization as shown in mixture.cpp " << endl; 
    assert(0);
    T expectedValue;
    return expectedValue;
  }

template <typename T>
  MatrixWrapper::SymmetricMatrix Mixture<T>::CovarianceGet (  ) const
  {
    TestNotInit(); 
    cerr << "Mixture CovarianceGet: not implemented since so far I don't believe its usefull"
     << endl << "If you decide to implement is: Use template specialization as shown in mcpdf.cpp " << endl; 

    assert(0);
    MatrixWrapper::SymmetricMatrix result;
    return result;
  }

template<typename T>
  vector<Probability> Mixture<T>::WeightsGet() const
  {
    TestNotInit(); 
    return *_componentWeights;
   }

template<typename T>
  Probability Mixture<T>::WeightGet(unsigned int componentNumber) const
  {
    TestNotInit(); 
    assert((int)componentNumber >= 0 && componentNumber < NumComponentsGet());
    return (*_componentWeights)[componentNumber];
   }

template<typename T>
  bool Mixture<T>::WeightsSet(vector<Probability> & weights)
  {
    TestNotInit(); 
    assert(weights.size() == NumComponentsGet());
    *_componentWeights = weights;
    //normalize the probabilities and update the cumulative sum
    return (NormalizeWeights() && CumWeightsUpdate());
  }

template<typename T>
  bool Mixture<T>::WeightSet(unsigned int componentNumber, Probability weight)
  {
    TestNotInit(); 
    assert((int)componentNumber >= 0 && componentNumber < NumComponentsGet());
    assert((double)weight<=1.0);
    
    if (NumComponentsGet() == 1)
    {
        (*_componentWeights)[0] = weight;
    }
    else
    {
        // renormalize other weights such that sum of probabilities will be
        // one after the weight of the component is set to weight
        // This should keep the probabilities normalized
        Probability old_weight = WeightGet(componentNumber);
        if ((double)old_weight!=1.0) {
            double normalization_factor = (1-weight)/(1-old_weight);
            for (int i=0; i<NumComponentsGet();i++)
            {
                (*_componentWeights)[i] = (Probability)( (double)( (*_componentWeights)[i] )* normalization_factor);
            }
        }
        else{
            for (int i=0; i<NumComponentsGet();i++)
            {
                (*_componentWeights)[i] = (Probability)( (1-weight)/(NumComponentsGet()-1));
            }
        }
        (*_componentWeights)[componentNumber] = weight;
    }
    return CumWeightsUpdate();
  }

template<typename T>
  int Mixture<T>::MostProbableComponentGet() const
  {
    TestNotInit(); 
    int index_mostProbable= -1;
    Probability prob_mostProbable= 0.0;
    for (int component = 0 ; component < NumComponentsGet() ; component++)
    {
       if ( (*_componentWeights)[component] > prob_mostProbable)
       {
            index_mostProbable= component;
            prob_mostProbable= (*_componentWeights)[component];
       }
    }
    return index_mostProbable;
  }

template<typename T>
  bool Mixture<T>::AddComponent(Pdf<T>& pdf) 
  {
    if (NumComponentsGet()==0)
        return AddComponent(pdf, Probability(1.0));
    else
    {
        _numComponents++;
        (*_componentPdfs).push_back(pdf.Clone() );

        (*_componentWeights).push_back(Probability(0.0));
	    _cumWeights.push_back(0.0);
        //assert length of vectors
        assert(NumComponentsGet()==(*_componentPdfs).size());
        assert(NumComponentsGet()==(*_componentWeights).size());
        assert(NumComponentsGet()+1==_cumWeights.size());
        return (NormalizeWeights() && CumWeightsUpdate());
    }
  }

template<typename T>
  bool Mixture<T>::AddComponent(Pdf<T>& pdf, Probability w) 
  {
    if (NumComponentsGet()==0 && w!=1.0)
        return AddComponent(pdf, Probability(1.0));
    else
    {
        _numComponents++;
        (*_componentPdfs).push_back(pdf.Clone() );
        (*_componentWeights).push_back(Probability(0.0));
        _cumWeights.resize(NumComponentsGet()+1);
        //assert length of vectors
        assert(NumComponentsGet()==(*_componentPdfs).size());
        assert(NumComponentsGet()==(*_componentWeights).size());
        assert(NumComponentsGet()+1==_cumWeights.size());
        WeightSet(_numComponents-1,w);
        return (NormalizeWeights() && CumWeightsUpdate());
     }
  }

template<typename T>
  bool Mixture<T>::DeleteComponent(unsigned int componentNumber) 
  {
    //assert length of vectors
    assert(NumComponentsGet()==(*_componentPdfs).size());
    assert(NumComponentsGet()==(*_componentWeights).size());
    assert(NumComponentsGet()+1==_cumWeights.size());

    TestNotInit(); 
    assert((int)componentNumber >= 0 && componentNumber < NumComponentsGet());
    _numComponents--;
    Pdf<T>* pointer = (*_componentPdfs)[componentNumber];
    delete pointer;
    (*_componentPdfs).erase((*_componentPdfs).begin()+componentNumber);
    (*_componentWeights).erase((*_componentWeights).begin()+componentNumber);
    _cumWeights.resize(NumComponentsGet()+1);
    //assert length of vectors
    assert(NumComponentsGet()==(*_componentPdfs).size());
    assert(NumComponentsGet()==(*_componentWeights).size());
    assert(NumComponentsGet()+1==_cumWeights.size());
    if(_numComponents==0) //don't do normalization if numComponents == 0
        return true;
    else
        return (NormalizeWeights() && CumWeightsUpdate());
  }

template<typename T>
  vector<Pdf<T>*> Mixture<T>::ComponentsGet() const
  {
    TestNotInit();
    return _componentPdfs;
  }

template<typename T>
  Pdf<T>* Mixture<T>::ComponentGet(unsigned int componentNumber) const
  {
    TestNotInit();
    return (*_componentPdfs)[componentNumber];
  }

template<typename T>
  void Mixture<T>::TestNotInit() const
  {
    if (NumComponentsGet() == 0)
    {
    cerr << "Mixture method called which requires that the number of components is not zero"
	 << endl << "Current number of components: " << NumComponentsGet() << endl;
    assert(0);
    }
  }

template<typename T>
  bool Mixture<T>::NormalizeWeights()
  {
    double SumOfWeights = 0.0;
    for ( unsigned int i = 0; i < NumComponentsGet() ; i++){
        SumOfWeights += (*_componentWeights)[i];
    }
    if (SumOfWeights > 0){
      for ( unsigned int i = 0; i < NumComponentsGet() ; i++){
          (*_componentWeights)[i] = (Probability)( (double) ( (*_componentWeights)[i]) /SumOfWeights);
      }
      return true;
    }
    else{
      cerr << "Mixture::NormalizeProbs(): SumOfWeights = " << SumOfWeights << endl;
      return false;
    }
  }

template<typename T>
  bool Mixture<T>::CumWeightsUpdate()
  {
    // precondition: sum of probabilities should be 1
    double CumSum=0.0;
    static vector<double>::iterator CumWeightsit;
    CumWeightsit = _cumWeights.begin();
    *CumWeightsit = 0.0;

    // Calculate the Cumulative PDF
    for ( unsigned int i = 0; i < NumComponentsGet(); i++)
    {
	    CumWeightsit++;
	    // Calculate the __normalised__ Cumulative sum!!!
	    CumSum += ( (*_componentWeights)[i] );
	    *CumWeightsit = CumSum;
    }
    // Check if last element of valuelist is +- 1
    assert( (_cumWeights[NumComponentsGet()] >= 1.0 - NUMERIC_PRECISION) &&
	    (_cumWeights[NumComponentsGet()] <= 1.0 + NUMERIC_PRECISION) );

    _cumWeights[NumComponentsGet()]=1;
    return true;
  }

} // End namespace

#include "mixture.cpp"

#endif // MIXTURE_H
