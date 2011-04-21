// $Id: dataAssociationFilter.h 
// Copyright (C) 2008 Tinne De Laet <first dot last at mech dot kuleuven dot be>
//  
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

#ifndef __DATA_ASSOCIATION_FILTER__
#define __DATA_ASSOCIATION_FILTER__

#include "../model/systemmodel.h"
#include "../model/measurementmodel.h"
#include "../pdf/pdf.h"
#include "../filter/filter.h"
#include <vector>
#include "../wrappers/matrix/matrix_wrapper.h"

namespace BFL
{
  using namespace std;
  using namespace MatrixWrapper;

  /// Class for data association filters
  /** This is a clas that defines the data association filters
      These filters are all related to an underlying set of filters and are used
      data with the correct filters.
  */
  template <typename StateVar, typename MeasVar> class DataAssociationFilter
    {
    protected:
 
      /// Vector of pointers to filters
      std::vector< Filter<StateVar,MeasVar>* > _filters; 

      /// iterator for vector of pointers to filters
      typename std::vector< Filter<StateVar,MeasVar>* >::iterator _iter_filters; 

      /// Vector of filter numbers
      std::vector< int > _filterID; 

      /// iterator for vector of filterID
      std::vector< int >::iterator _iter_filterID; 
    
      int _timestep;

      /// Probability of false positive (feature measurement without object)
      double _gamma;

      /// Treshold of probability to take into account associations
      double _treshold;

      /// The effective number of filters
      int _numFilters;

      /// vector of pointers to posteriors of filters
      std::vector< Pdf<StateVar>* > _posts; 

      /// iterator for vector of pointers to posteriors of filters
      typename std::vector<Pdf<StateVar>* >::iterator _iter_posts; 

      /// vector containing prob(z|object) 
      Matrix _probs;
    
      /// helper variable for probabilities
      Matrix _probs_rest;

      /// maximum number of filters
      int _maxFilters ; 
      /// maximum number of features
      int _maxFeatures; 
      /// maximum number of associations
      int _maxAssociations; 
      /// maximum calls of GetAssociations to calculate associations
      int _maxCalls; 
      /// counter for the number of different filters used so far to obtain filterID
      int _filterCounter; 

      /// vector genest associations
      vector< vector<vector<int> >  > _vec_ass; 
      /// vector of genest number of associations
      vector< int > _vec_number_associations;
      /// Vector of associations
      vector<vector<int> > _associations; 
      /// Vector of possible objects associated with the considered feature
      vector<int> _objects_ass;  
      /// Vector of genest possible objects associated with the considered feature
      vector<vector<int> > _vec_objects_ass;   
      /// Vector of genest number associations
      vector< int > _vec_num_objects_ass; 
      /// Vector of genest probabilities
      vector< Matrix > _vec_probs; 
      /// number of associated objects
      int _num_objects_ass ; 
      /// number of assocations
      int _num_ass ; 
      /// The number of times GetAssociations is called
      int _numberGetAssociationsCalls; 
      /// The probabilites that the data is associated with each of the filters
      vector<vector<double> > _association_probs;  
    
      /// Implementation of Update
      // calls update on the underlying filters
      /** @param sysmodel pointer to the used system model
	  @param u input param for proposal density
	  @param measmodel pointer to the used measurementmodel
	  @param z measurement param for proposal density
	  @param s sensor param for proposal density
      */
      virtual bool UpdateInternal(SystemModel<StateVar>* const sysmodel,
				  const StateVar& u,
				  MeasurementModel<MeasVar,StateVar>* const measmodel,
				  const vector<MeasVar>& z,
				  const StateVar& s)=0;

      // Get associations 
      // Get all possible associations of the data with the filters
      /** @param probs matrix containg the probability of each of the  measurements given each of the filters 
	  @param number_features the number of features
      */
      void GetAssociations(const MatrixWrapper::Matrix& probs, int number_features);

    public:
      /// Constructor
      /** @pre you created the prior
	  @param filters pointer to the prior vector of filters
	  @param gamma probability of a false positive (measurement is not caused by any of the objects)
	  @param treshold threshold on the probability of a association from which the association is taken into account 
	  @param maxFilters the maximum number of filters (for allocation)
	  @param maxFeatures the maximum number of features (for allocation)
	  @param maxAssociations the maximum number of associations (for allocation)
	  @param maxCalls the maximum number of calls of getAssociations (for allocation)
      */
      DataAssociationFilter (vector<Filter<StateVar,MeasVar> *> filters, double gamma = 0.0, double treshold= 0.0, int maxFilters =20, int maxFeatures=20, int maxAssociations=2000,int maxCalls=2000);

      /// copy constructor
      /** @pre you have another DataAssociationFilter 
	  @param DataAssociationFilter you want to copy 
      */
      DataAssociationFilter (const DataAssociationFilter<StateVar,MeasVar> & filters);

      /// Get the maximum number of filters
      /** @return the maximum number of filters
      */
      int GetMaxFilters(){return _maxFilters;};

      /// Get the maximum number of features
      /** @return the maximum number of features
      */
      int GetMaxFeatures(){return _maxFeatures;};

      /// Get the maximum number of associations
      /** @return the maximum number of associations
      */
      int GetMaxAssociations(){return _maxAssociations;};

      /// Get the maximum number of getAssociation calls
      /** @return the maximum number of getAssociation calls
      */
      int GetMaxCalls(){return _maxCalls;};

      /// destructor
      virtual ~DataAssociationFilter();

      /// Get the probabilities of the measurements for each of the filters 
      /** @param sysmodel pointer to the used system model
	  @param measmodel pointer to the used measurementmodel
	  @param z vector of measurements 
	  @param s sensor param for proposal density
      @return Matrix containing the probability of each of the measurements given each of the filters
      */
      virtual MatrixWrapper::Matrix GetMeasProbs(MeasurementModel<MeasVar,StateVar>* const measmodel , const vector<MeasVar> &z, const StateVar& s)=0;

      /// Get the filters of the DataAssociationFilter
      /** @return vector of Filters underlying the dataAssociationFilter
      */
      vector< Filter<StateVar,MeasVar> > FiltersGet();

      /// Get Posterior density
      /** Get the current Posterior density
	  @return a vector of pointers to the current posteriors
      */
      vector< Pdf<StateVar> * > PostGet();  

      /// Get current time
      /** Get the current time of the filter
	  @return the current timestep
      */
      int TimeStepGet() const;

      /// Get the filtersIDs of the DataAssociationFilter
      /** @return vector of ints indicating the ID (identification number) of the filter 
      */
      vector< int > FilterIDGet();

      /// Get the filterCounter of the DataAssociationFilter
      /** @return the filter counter (maximum ID of the filters so far created)
      */
      int FilterCounterGet();

      /// Get the probability of false positive
      /** @return the probability of a false positive 
      */
      double GammaGet();

      /// Get the treshold on the probability from which associations are taken into account
      /** @return the treshold of the probability from which associations are  taken into account 
      */
      double TresholdGet();

      /// Reset Filters
      /** @param filters new vector of filters for the DataAssociationFilter
      */
      void Reset(vector< Filter<StateVar,MeasVar> *> filters);

      /// Get number of Filters
      /** @return the current number of filters 
      */
      int NumFiltersGet();

      /// Add a filter
      /** @param filter pointer to the filter that will be added
      */
      void AddFilter(Filter<StateVar,MeasVar>* filter);

      /// Remove a filter
      /** @param index index of the filter that will be removed. This index
         should be between 0 and NumFitlersGet()-1, else false is returned
          @return bool indicating of the removal was succesfull or not
      */
      bool RemoveFilter(int index);

      // Get the association probabilities of the measurements with each
      // probable association
      /** 
	  @param measmodel pointer to the used measurementmodel
	  @param z vector of measurements 
	  @param s sensor param for proposal density
      @return vector of vector of doubles containing the probability of each of the possible associations
      */
      vector<vector<double> > GetAssociationProbs(MeasurementModel<MeasVar,StateVar>* const measmodel , const vector<MeasVar> &z, const StateVar& s);

      /// Full Update (system with inputs/sensing params)
      /** @param sysmodel pointer to the system model to use for update
	  @param u input to the system
	  @param measmodel pointer to the measurement model to use for update
	  @param z measurement
	  @param s "sensing parameter"
       */
      bool Update(SystemModel<StateVar>* const sysmodel,
			  const StateVar& u,
			  MeasurementModel<MeasVar,StateVar>* const measmodel,
			  const vector<MeasVar>& z,
			  const StateVar& s);

      /// Full Update (system without inputs, with sensing params)
      /** @param sysmodel pointer to the system model to use for
	  update
	  @param measmodel pointer to the measurement model to use for
	  update
	  @param z measurement
	  @param s "sensing parameter"
       */
      bool Update(SystemModel<StateVar>* const sysmodel,
			  MeasurementModel<MeasVar,StateVar>* const measmodel,
			  const vector<MeasVar>& z,
			  const StateVar& s);

      /// Full Update (system without inputs/sensing params)
      /** @param sysmodel pointer to the system model to use for
	  update
	  @param measmodel pointer to the measurement model to use for
	  update
	  @param z measurement
       */
      bool Update(SystemModel<StateVar>* const sysmodel,
			  MeasurementModel<MeasVar,StateVar>* const measmodel,
			  const vector<MeasVar>& z);

      /// Full Update (system with inputs, without sensing params)
      /** @param sysmodel pointer to the system model to use for update
	  @param u input to the system
	  @param measmodel pointer to the measurement model to use for
	  update
	  @param z measurement
       */
      bool Update(SystemModel<StateVar>* const sysmodel,
			  const StateVar& u,
			  MeasurementModel<MeasVar,StateVar>* const measmodel,
			  const vector<MeasVar>& z);

      /// System Update (system with inputs)
      /** @param sysmodel pointer to the system model to use for update
	  @param u input to the system
       */
      bool Update(SystemModel<StateVar>* const sysmodel,
			  const StateVar& u);

      /// System Update (system without inputs)
      /** @param sysmodel pointer to the system model to use for update
       */
      bool Update(SystemModel<StateVar>* const sysmodel);

      /// Measurement Update (system with "sensing params")
      /** @param measmodel pointer to the measurement model to use for
	  update 
	  @param z measurement
	  @param s "sensing parameter"
       */
      bool Update(MeasurementModel<MeasVar,StateVar>* const measmodel,
			  const vector<MeasVar>& z,
			  const StateVar& s);

      /// Measurement Update (system without "sensing params")
      /** @param measmodel pointer to the measurement model to use for
	  update
	  @param z measurement
       */
      bool Update(MeasurementModel<MeasVar,StateVar>* const measmodel,
			  const vector<MeasVar>& z);

    };


////////////// IMPLEMENTATION/////////////////////////////////////

    #define StateVar SVar
    #define MeasVar MVar

    // Constructor
    template<typename SVar, typename MVar> 
    DataAssociationFilter<SVar,MVar>::DataAssociationFilter(vector<Filter<SVar,MVar> * > filters, double gamma, double treshold, int maxFilters, int maxFeatures, int maxAssociations, int maxCalls)
      : _timestep(0),
        _gamma(gamma),
        _treshold(treshold),
        _maxFilters(maxFilters),
        _maxFeatures(maxFeatures),
        _maxAssociations(maxAssociations),
        _maxCalls(maxCalls)
    {
      _numFilters = filters.size();

      _probs.resize(_maxFilters,_maxFeatures);
      _probs = 0.0;
      _posts.resize(_maxFilters);
      _filters.resize(_maxFilters);
      _filterID.resize(_maxFilters);

      typename std::vector< Filter<StateVar,MeasVar>* >::iterator iter_filters_local; 
      _iter_posts = _posts.begin();
      _iter_filters = _filters.begin();
      int counter = 0;
      for(iter_filters_local = filters.begin() ; iter_filters_local != filters.end() ; iter_filters_local++) 
      {
        *_iter_filters =(*iter_filters_local); 
        *_iter_posts = (*_iter_filters)->PostGet();
        _filterID[counter] = counter+1;
        _iter_posts++;
        _iter_filters++;
        counter ++;
      }
      _filterCounter = this->NumFiltersGet();
      _iter_posts = _posts.begin();
      _iter_filters = _filters.begin(); 

      _associations.resize(_maxAssociations);
      _associations.assign(_maxAssociations,vector<int>(_maxFeatures));
      _vec_ass.resize(_maxCalls);
      _vec_ass.assign(_maxCalls,_associations);
      _objects_ass.resize(_maxFilters+1);  // vector of possible objects associated with the considered feature
      _num_objects_ass = 0;
      _num_ass = 0;
      _numberGetAssociationsCalls = 0;
      _vec_number_associations.resize(_maxCalls);
      _vec_objects_ass.resize(_maxCalls);  // vector of genest possible objects associated with the considered feature
      _vec_objects_ass.assign(_maxCalls,_objects_ass);
      _vec_num_objects_ass.resize(_maxCalls); // vector of genest number associations

      _probs_rest.resize(_maxFilters,_maxFeatures);
      _probs_rest = 0.0;
      _vec_probs.resize(_maxCalls);
      _vec_probs.assign(_maxCalls,_probs_rest);

      _association_probs.resize(_maxFeatures);
      _association_probs.assign(_maxFeatures,vector<double>(_maxFilters));
    }
     
    template<typename SVar, typename MVar> 
    DataAssociationFilter<SVar,MVar>::~DataAssociationFilter()
    {
        // throw away filters
        for(_iter_filters = _filters.begin() ; _iter_filters != _filters.end() ; _iter_filters++) 
        {
          delete *_iter_filters;
        }
    }

    template<typename SVar, typename MVar> 
    DataAssociationFilter<SVar,MVar>::DataAssociationFilter(const DataAssociationFilter& filters)
        : _timestep(0)
    {   
        _maxFilters = filters.GetMaxFilters();
        _maxFeatures = filters.GetMaxFeatures();
        _maxAssociations = filters.GetMaxAssociations();
        _maxCalls = filters.GetMaxCalls();
        _probs.resize(_maxFilters,_maxFeatures);
        _probs = 0.0;
        _posts.resize(_maxFilters);
        _filters.resize(_maxFilters);
        _filterID.resize(_maxFilters);

        _numFilters = filters.NumFiltersGet();
        _gamma = filters.GammaGet();
        _treshold = filters.TresholdGet();
        _filters = filters.FiltersGet();
        _filterID = filters.FilterIDGet();
        _filterCounter = filters.FilterCounterGet();
        _iter_posts = _posts.begin();
        int counter = 0;
        for(_iter_filters = _filters.begin() ; _iter_filters != _filters.end() ; _iter_filters++) 
        {
          *_iter_posts = (*_iter_filters)->PostGet();
          _iter_posts++;
        }
        _iter_posts = _posts.begin();
        _iter_filters = _filters.begin(); 
        _probs.resize(_maxFilters,_maxFeatures);
        _probs = 0.0;

        _associations.resize(_maxAssociations);
        _associations.assign(_maxAssociations,vector<int>(_maxFeatures));
        _vec_ass.resize(_maxCalls);
        _vec_ass.assign(_maxCalls,_associations);
        _objects_ass.resize(_maxFilters+1);  // vector of possible objects associated with the considered feature
        _num_objects_ass = 0;
        _num_ass = 0;
        _numberGetAssociationsCalls = 0;
        _vec_number_associations.resize(_maxCalls);
        _vec_objects_ass.resize(_maxCalls);  // vector of genest possible objects associated with the considered feature
        _vec_objects_ass.assign(_maxCalls,_objects_ass);
        _vec_num_objects_ass.resize(_maxCalls); // vector of genest number associations

        _probs_rest.resize(_maxFilters,_maxFeatures);
        _probs_rest = 0.0;

        _association_probs.resize(_maxFeatures);
        _association_probs.assign(_maxFeatures,vector<double>(_maxFilters));
    }

    template<typename SVar, typename MVar> vector<Filter<SVar,MVar> >
    DataAssociationFilter<SVar,MVar>:: FiltersGet()
    {
        return _filters;
    }

    template<typename SVar, typename MVar> vector<int>
    DataAssociationFilter<SVar,MVar>:: FilterIDGet()
    {
        return _filterID;
    }

    template<typename SVar, typename MVar> int
    DataAssociationFilter<SVar,MVar>:: FilterCounterGet()
    {
        return _filterCounter;
    }
    
    template<typename SVar, typename MVar> double 
    DataAssociationFilter<SVar,MVar>:: GammaGet()
    {
        return _gamma;
    }

    template<typename SVar, typename MVar> double 
    DataAssociationFilter<SVar,MVar>:: TresholdGet()
    {
        return _treshold;
    }

    template<typename SVar, typename MVar>  void
    DataAssociationFilter<SVar,MVar>::Reset(vector< Filter<SVar,MeasVar> *> filters)
    {
      _numFilters = filters.size();
      typename std::vector< Filter<StateVar,MeasVar>* >::iterator iter_filters_local; 
      _iter_posts = _posts.begin();
      _iter_filters = _filters.begin();
      int counter = 0;
      for(iter_filters_local = filters.begin() ; iter_filters_local != filters.end() ; iter_filters_local++) 
      {
        *_iter_filters =(*iter_filters_local); 
        *_iter_posts = (*_iter_filters)->PostGet();
        _filterID[counter] = counter+1;
        _iter_posts++;
        _iter_filters++;
        counter ++ ;
      }
      _filterCounter = this->NumFiltersGet();
      _iter_posts = _posts.begin();
      _iter_filters = _filters.begin(); 
    }

    template<typename SVar, typename MVar>  int
    DataAssociationFilter<SVar,MVar>::NumFiltersGet()
    {
      return _numFilters;
    }

    template<typename SVar, typename MVar>  void
    DataAssociationFilter<SVar,MVar>::AddFilter(Filter<SVar,MVar>* filter)
    {
        _filters[_numFilters] = filter;
        _posts[_numFilters] = filter->PostGet();
        _filterID[_numFilters] = _filterCounter + 1;
        _filterCounter ++;
        _numFilters++;
    }
    
    template<typename SVar, typename MVar>  bool
    DataAssociationFilter<SVar,MVar>::RemoveFilter(int index)
    {
        if(index< 0 || index > _numFilters-1)
        {
            return false;
        }
        else
        {
            // get iterator to the index'th element
            _iter_filters = _filters.begin();
            _iter_posts = _posts.begin();
            _iter_filterID = _filterID.begin();
            for(int i = 0 ; i<index ; i++ )
            {
                _iter_posts++;
                _iter_filters++;
                _iter_filterID++;
            }
            // delete the filter itself!
            Filter<SVar,MVar>* temp_pointer = *_iter_filters;
            _posts.erase(_iter_posts);
            _filters.erase(_iter_filters);
            _filterID.erase(_iter_filterID);
            delete temp_pointer;
            // TODO: nothing more efficient to shift elements of vector forward?
            _posts.resize(_maxFilters);
            _filters.resize(_maxFilters);
            _numFilters--;
            return true;
        }
    }

    template<typename SVar, typename MVar> int 
    DataAssociationFilter<SVar,MVar>::TimeStepGet() const
    {
      return _timestep;
    }
    
    //template<typename SVar, typename MVar> bool
    //DataAssociationFilter<SVar,MVar>::UpdateInternal(SystemModel<SVar>* const sysmodel,
    //			  const SVar& u,
    //			  MeasurementModel<MVar,SVar>* const measmodel,
    //			  const MVar& z,
    //			  const SVar& s)
    //{
    //  bool return_bool = true;
    //  _iter_posts = _posts.begin();
    //  for(_iter_filters = _filters.begin() ; _iter_filters != _filters.end() ; _iter_filters++) 
    //  {
    //    return_bool = (return_bool && _iter_filters->Update(sysmodel,u,measmodel,z,s) ); 
    //    *_iter_posts = (*_iter_filters)->PostGet();
    //    _iter_posts++;
    //  }
    //  return return_bool;
    //}
    
    template<typename SVar, typename MVar> bool
    DataAssociationFilter<SVar,MVar>::Update(SystemModel<SVar>* const sysmodel,
    			  const SVar& u,
    			  MeasurementModel<MVar,SVar>* const measmodel,
    			  const vector<MVar>& z,
    			  const SVar& s)
    {
      return this->UpdateInternal(sysmodel,u,measmodel,z,s);
    }
    
    template<typename SVar, typename MVar> bool
    DataAssociationFilter<SVar,MVar>::Update(SystemModel<SVar>* const sysmodel,
    			  MeasurementModel<MVar,SVar>* const measmodel,
    			  const vector<MVar>& z,
    			  const SVar& s)
    {
      SVar u;
      return this->UpdateInternal(sysmodel,u,measmodel,z,s);
    }
    
    template<typename SVar, typename MVar> bool
    DataAssociationFilter<SVar,MVar>::Update(SystemModel<SVar>* const sysmodel,
    			  MeasurementModel<MVar,SVar>* const measmodel,
    			  const vector<MVar>& z)
    {
      SVar s; SVar u;
      return this->UpdateInternal(sysmodel,u,measmodel,z,s);
    }
    
    template<typename SVar, typename MVar> bool
    DataAssociationFilter<SVar,MVar>::Update(SystemModel<SVar>* const sysmodel,
    			  const SVar& u,
    			  MeasurementModel<MVar,SVar>* const measmodel,
    			  const vector<MVar>& z)
    {
      SVar s;
      return this->UpdateInternal(sysmodel,u,measmodel,z,s);
    }
    
    template<typename SVar, typename MVar> bool
    DataAssociationFilter<SVar,MVar>::Update(SystemModel<SVar>* const sysmodel,
    			  const SVar& u)
    {
      SVar s; vector<MVar> z(1);
      return this->UpdateInternal(sysmodel,u,NULL,z,s);
    }
    
    template<typename SVar, typename MVar> bool
    DataAssociationFilter<SVar,MVar>::Update(SystemModel<SVar>* const sysmodel)
    {
      SVar s;  SVar u; vector<MVar> z;
      return this->UpdateInternal(sysmodel,u,NULL,z,s);
    }
    
    template<typename SVar, typename MVar> bool
    DataAssociationFilter<SVar,MVar>::Update(MeasurementModel<MVar,SVar>* const measmodel,
    			  const vector<MVar>& z,
    			  const SVar& s)
    {
      SVar u; 
      return this->UpdateInternal(NULL,u,measmodel,z,s);
    }
    
    template<typename SVar, typename MVar> bool
    DataAssociationFilter<SVar,MVar>:: Update(MeasurementModel<MVar,SVar>* const measmodel,
    			  const vector<MVar>& z)
    {
      SVar u; SVar s;
      return this->UpdateInternal(NULL,u,measmodel,z,s);
    }
    
    template<typename SVar, typename MVar> vector< Pdf<SVar> * >
    DataAssociationFilter<SVar,MVar>::PostGet()
    {
      //for(_iter_filters = _filters.begin() ; _iter_filters != _filters.end() ; _iter_filters++) 
      //{
      //  *_iter_posts = (*_iter_filters)->PostGet();
      //  _iter_posts++;
      //}
      return _posts;
    }

    /*
    // update with weighting between all different possible associations
    template<typename SVar, typename MVar> vector<vector<double> > 
    DataAssociationFilter<SVar,MVar>::GetAssociationProbs(MeasurementModel<MeasVar,StateVar>* const measmodel , const vector<MeasVar> &z, const StateVar& s)
    {
      int number_features = z.size();
      int number_objects = this->NumFiltersGet();

#ifndef NDEBUG
      //cout << "Get meas probabilities " << endl;
#endif
      this->GetMeasProbs(measmodel , z , s);
#ifndef NDEBUG
      std::cout << "_probs" << _probs << std::endl;
#endif
      // stores result in _probs

      // get associations
      //cout << "Get associations " << endl;
      this->GetAssociations(_probs,z.size());
      // stores result in _associations
      // number of associations in _num_ass
      for (int k = 0 ; k< _num_ass; k++)
      {
        std::cout << "_association[" << k << "]" << std::endl;
        for(int l = 0 ; l<number_features; l++) 
            std::cout<< _associations[k][l] << std::endl;

      }

      //cout << "Get assocation probs " << endl;

      double prod_prob = 1.0;
    
      vector<vector<double> > association_probs(number_features);
      association_probs.assign(number_features,vector<double>(number_objects,0.0));

      for (int j=0 ; j<number_features ; j++)
      {
        for(int i=0; i<number_objects; i++)
        {
            // calculate association_probs[j][i]
            // check all associations for association of feature j with object i
            for (int k = 0 ; k< _num_ass; k++)
            {
                //check if feature j is assigned to object i in association k
                if (_associations[k][j] == i+1)
                {
                    // check this association
                    prod_prob = 1.0;
                    for (int l = 0 ; l< number_features  ; l++) // see with which object each feature is associated
                    {
                        // feature l is associated with object _associations[k][l]
                        if (_associations[k][l]  == 0 ) // false alarm
                            prod_prob = prod_prob * _gamma; 
                        else
                            prod_prob = prod_prob * _probs(_associations[k][l],l+1);
                    }
                    association_probs[j][i] = association_probs[j][i] + prod_prob;
                }
            }
            cout<< "association_probs[" << j << "][" << i << "] " << association_probs[j][i] << endl;
         }
       }
      return association_probs;
    }
   */ 

    // update with only taking into account most probable association
    template<typename SVar, typename MVar> vector<vector<double> > 
    DataAssociationFilter<SVar,MVar>::GetAssociationProbs(MeasurementModel<MeasVar,StateVar>* const measmodel , const vector<MeasVar> &z, const StateVar& s)
    {
#ifndef NDEBUG
      cout << "(DataAssociationFilter) entered GetAssociationProbs " << endl;
#endif
      int number_features = z.size();
      int number_objects = this->NumFiltersGet();

#ifndef NDEBUG
      cout << "Get meas probabilities " << endl;
#endif
      this->GetMeasProbs(measmodel , z , s);
#ifndef NDEBUG
      cout << "Got meas probabilities " << endl;
#endif
#ifndef NDEBUG
      //std::cout << "_probs" << _probs << std::endl;
#endif
#ifndef NDEBUG
      //std::cout << "_gamma" << _gamma << std::endl;
#endif
#ifndef NDEBUG
      //std::cout << "_treshold" << _treshold << std::endl;
#endif
      // stores result in _probs

      // get associations
#ifndef NDEBUG
      cout << "Get associations " << endl;
#endif
      this->GetAssociations(_probs,z.size());
#ifndef NDEBUG
      cout << "Got associations jeeej" << endl;
#endif
      // stores result in _associations
      // number of associations in _num_ass
      //for (int k = 0 ; k< _num_ass; k++)
      //{
#ifndef NDEBUG
      //  std::cout << "_association[" << k << "]" << std::endl;
#endif
      //  for(int l = 0 ; l<number_features; l++) 
#ifndef NDEBUG
      //      std::cout<< _associations[k][l] << std::endl;
#endif

      //}

#ifndef NDEBUG
      cout << "Get assocation probs " << endl;
#endif
      // find most probable association prob_ass[k] is probability of
      // association k
      vector<double> prob_ass(_num_ass);
      double temp;
      double max_prob_ass = -1.0;
      int best_association = -1;
      for (int k = 0 ; k< _num_ass; k++) //loop over all associations
      {
          temp = 1.0;
          for (int l = 0 ; l< number_features  ; l++) // see with which object each feature is associated
          {
              if (_associations[k][l]  == 0 ) // false alarm
                  temp = temp *  _gamma; 
              else
                  temp = temp * _probs(_associations[k][l],l+1);
          }
          prob_ass[k] = temp;
          if( prob_ass[k] > max_prob_ass )
          {
            max_prob_ass = prob_ass[k];
            best_association = k;
          }
      }
#ifndef NDEBUG
      //cout << "best_association " << best_association << endl;
#endif
      //for(int l = 0 ; l<number_features; l++) 
#ifndef NDEBUG
      //      std::cout<< _associations[best_association][l] << std::endl;
#endif

      int k = best_association;

      //vector<vector<double> > association_probs(number_features);
      //association_probs.assign(number_features,vector<double>(number_objects,0.0));
      _association_probs.assign(_maxFeatures,vector<double>(_maxFilters,0.0));

      //check if feature j is assigned to object i in best association
      int object = -1;
      for (int l = 0 ; l< number_features  ; l++) // see with which object each feature is associated
      {
          // feature l is associated with object _associations[k][l]
          if (_associations[k][l]  == 0 ) // false alarm
          {
              //prod_prob =  _gamma; 
          }
          else
          {
              object = _associations[k][l]-1;
              _association_probs[l][object] = _probs(_associations[k][l],l+1);
          }
      }
      return _association_probs;
    }
/*

    template<typename SVar, typename MVar> vector<vector <int> > 
    DataAssociationFilter<SVar,MVar>::GetAssociations(const MatrixWrapper::Matrix& probs,  int number_features)
    {
        // restricted associations: a feature can only belong to one object and each object only has one feature.
        // probs: matrix of probabilities p(zj|xi), probability of measurement zj given object xi
        // size probs: [number of objects x number of features]

#ifndef NDEBUG
        //std::cout << "*****************************"  <<std::endl;
#endif
#ifndef NDEBUG
        //std::cout << "number_features " << number_features <<std::endl;
#endif
        vector<vector<int> > associations(1);

        vector<int> objects_ass(1);  // vector of possible objects associated with the considered feature
        objects_ass[0] = 0;         // 0: feature is not caused by any of the objects

        associations[0] = objects_ass;

        int number_ass = objects_ass.size();

        // base case: number_features = 1;
        if (number_features == 1)
        {
#ifndef NDEBUG
            //std::cout << "Base case" << std::endl;
#endif
            vector<int> objects_vec(1);
            for (int objects = 1; objects <= number_objects; objects++)
            {
                if(probs(objects,1) > _treshold) // if feature is probably caused by object
                {
                    objects_vec[0] = objects;
                    objects_ass.push_back(objects); // add object to object association list
                    associations.push_back(objects_vec); // add object to object association list
                }
            }
            number_ass = objects_ass.size();
            //std::cout << "Calculated associations " << std::endl;
            //for (int i =0 ; i < associations.size() ; i++)
            //{
            //     std::cout << "Association " << i+1 << std::endl;
            //     for (int j=0 ; j < associations[i].size() ; j++)
            //          std::cout << associations[i][j] << std::endl;
            //}
            //std::cout << "END *****************************"  <<std::endl;
            return associations; 
        }
        else // not yet in base case
        {
            // get all the possibilities for the first feature
            for (int objects = 1; objects <= number_objects; objects++)
            {
                if(probs(objects,1) > _treshold) // if feature is probably caused by object
                    objects_ass.push_back(objects); // add object to object association list
            }
            number_ass = objects_ass.size();
            //std::cout << "Number associations for first feature " << number_ass << std::endl;
            int counter = 0;
            int object;
            for (int i=0 ; i < number_ass ; i++)
            {
                object = objects_ass[i];
                //std::cout << "Considering object " << object << std::endl;
                // delete object from other possibilities if object !=0
                // probs_rest = probabilities for features if object for feature
                // 1 is chosen
                MatrixWrapper::Matrix probs_rest(number_objects,number_features-1);
                //std::cout << "probs " << probs << std::endl;
                probs_rest = probs.sub(1,number_objects,2,number_features);
                //std::cout << "probs_rest " << probs_rest << std::endl;
                if (object != 0)
                {
                    for (int k = 1 ; k<=number_features - 1 ; k++)
                        probs_rest(object,k) = -1;
                    //probs_rest.sub(object,object,1,number_features-1) = -1;// if you set to -1 will never be selected if _treshold >= 0
                    //std::cout << "probs_rest AFTER DELETING DOUBLES " << probs_rest << std::endl;
                }
                //std::cout << "getasociations for  " << probs_rest << std::endl;
                std::vector<std::vector<int> > associations_other = GetAssociations(probs_rest,number_features-1);
                //std::cout << "gotassociations" << std::endl;
                int num_associations_other = associations_other.size();
                //std::cout << "num_associations_other " << num_associations_other << std::endl;
                associations.resize(counter+num_associations_other);
                //std::cout << "associations.size() " << associations.size() << std::endl;
                for (int teller_other = 0 ; teller_other < num_associations_other ; teller_other ++)
                {
                    associations[counter+teller_other].resize(number_features);
                    //std::cout << "associations[ " << counter+teller_other << "] size " << associations[counter+teller_other].size() << std::endl;
                    //cout << "teller_other " << teller_other << endl;
                    associations[counter + teller_other][0] = object;
                    //cout << "associations[" << counter+teller_other << "][ 0] " << associations[counter + teller_other][0] << endl;
                    for (int teller_features = 1 ; teller_features <= number_features-1 ; teller_features++ )
                    {
                        //cout << "teller_features " << teller_features << endl;
                        associations[counter + teller_other][teller_features] = associations_other[teller_other][teller_features-1];
                        //cout << "associations[" << counter+teller_other << " ][ " << teller_features << "] " << associations[counter + teller_other][teller_features] << endl;
                    }
                }
                //std::cout << "Updated associations " << std::endl;
                //for (int i =0 ; i < associations.size() ; i++)
                //  {
                //       std::cout << "Association " << i+1 << std::endl;
                //       for (int j=0 ; j < associations[i].size() ; j++)
                //            std::cout << associations[i][j] << std::endl;
                //  }
                counter = counter + num_associations_other;
            }
            std::cout << "Calculated " << associations.size() << " associations " << std::endl;
            //for (int i =0 ; i < associations.size() ; i++)
            //  {
            //       std::cout << "Association " << i+1 << std::endl;
            //       for (int j=0 ; j < associations[i].size() ; j++)
            //            std::cout << associations[i][j] << std::endl;
            //  }

            //std::cout << "END *****************************"  <<std::endl;
            return associations;
        }

    }
*/

    template<typename SVar, typename MVar> void 
    DataAssociationFilter<SVar,MVar>::GetAssociations(const MatrixWrapper::Matrix& probs, int number_features)
    {
        // restricted associations: a feature can only belong to one object and each object only has one feature.
        // probs: matrix of probabilities p(zj|xi), probability of measurement zj given object xi
        // size probs: [number of objects x number of features]
#ifndef NDEBUG
        cout << "(DataAssociationFilter) entered GetAssociations " << endl;
#endif

        int number_objects = this->NumFiltersGet();
#ifndef NDEBUG
        cout << "number of objects " << number_objects << endl;
#endif
#ifndef NDEBUG
        cout << "number of features " << number_features << endl;
#endif
        assert(probs.rows() > number_objects -1);
        // probs must have the correct size

        _numberGetAssociationsCalls ++ ;
        int numberCall = _numberGetAssociationsCalls;
        _vec_probs[numberCall-1] = probs;
        // in this call _vec_associations[numberCall- 1] will be
        // filled in

#ifndef NDEBUG
        //std::cout << "*****************************"  <<std::endl;
#endif
#ifndef NDEBUG
        //std::cout << "Called GetAssociations the  " <<  numberCall << "'th time"  <<std::endl;
#endif
#ifndef NDEBUG
        //std::cout << "number_features " << number_features <<std::endl;
#endif
#ifndef NDEBUG
        //std::cout << "probs " << probs <<std::endl;
#endif
        
        //_objects_ass
        _objects_ass[0] = 0;         // 0: feature is not caused by any of the objects
        _associations[0][0] = 0;
        _num_objects_ass = 1;
        _num_ass = 1;

        // base case: number_features = 1;
        if (number_features == 1)
        {
#ifndef NDEBUG
            std::cout << "Base case" << std::endl;
#endif
            for (int objects = 1; objects <= number_objects; objects++)
            {
#ifndef NDEBUG
                std::cout << "object "  << objects << std::endl;
#endif
#ifndef NDEBUG
                std::cout << "probs(objects,1) "  << probs(objects,1) << std::endl;
#endif
                if(probs(objects,1) > _treshold) // if feature is probably caused by object
                {
                    _objects_ass[_num_objects_ass] = objects; // add object to object association list
                    _vec_objects_ass[numberCall-1] = _objects_ass;
                    _vec_ass[numberCall-1][_num_ass][0] = objects;
                    if (numberCall == 1) // the loop ended!!!!
                    {
                        _associations[_num_ass][0] = objects;
                    }
                    _num_objects_ass ++;
                    _num_ass ++;
                }
            }
            _vec_objects_ass[numberCall-1] = _objects_ass;
            _vec_num_objects_ass[numberCall-1] = _num_objects_ass;

#ifndef NDEBUG
            std::cout << "Calculated associations " << std::endl;
#endif
#ifndef NDEBUG
            //std::cout << "numberCall " << numberCall << std::endl;
#endif
           // for (int i =0 ; i < _num_ass ; i++)
           // {
#ifndef NDEBUG
           //      std::cout << "Association " << i+1 << std::endl;
#endif
           //      for (int j=0 ; j < number_features ; j++)
#ifndef NDEBUG
           //           std::cout << _vec_ass[numberCall-1][i][j] << std::endl;
#endif
           // }
#ifndef NDEBUG
            //std::cout << "_vec_ass.size() " << _vec_ass.size() <<  std::endl;
#endif
#ifndef NDEBUG
            //std::cout << "numberCall " << numberCall << std::endl;
#endif
            _vec_number_associations[numberCall-1] = _num_ass;
#ifndef NDEBUG
            std::cout << "number associations " << _num_ass <<  std::endl;
#endif
            if(numberCall == 1)
            {
#ifndef NDEBUG
                cout << "RESETTING" << endl;
#endif
                _numberGetAssociationsCalls = 0;
            }
#ifndef NDEBUG
            std::cout << "end basecase"  <<std::endl;
#endif
            //return _associations;
        }
        else // not yet in base case
        {
            // get all the possibilities for the first feature
            for (int objects = 1; objects <= number_objects; objects++)
            {
                if(probs(objects,1) > _treshold) // if feature is probably caused by object
                {
                    _objects_ass[_num_objects_ass] = objects; // add object to object association list
                    _num_objects_ass ++;
#ifndef NDEBUG
                    //std::cout << "object " << objects << "added"  <<std::endl;
#endif
                }
            }
            _vec_objects_ass[numberCall-1] = _objects_ass;
            _vec_num_objects_ass[numberCall-1] = _num_objects_ass;
#ifndef NDEBUG
            //std::cout << "Number associations for first feature " << number_ass << std::endl;
#endif
            int counter = 0;
            int object;
#ifndef NDEBUG
           // std::cout << "Number of considered objects " << _vec_num_objects_ass[numberCall-1]<< std::endl;
#endif
            //for (int i=0 ; i < _objects_ass[numberCall-1]; i++)
            for (int i=0 ; i < _vec_num_objects_ass[numberCall-1]; i++)
            {
#ifndef NDEBUG
               //std::cout << "Number of considered objects " << _vec_num_objects_ass[numberCall-1]<< std::endl;
#endif
               // for (int teller=0 ; teller < _vec_num_objects_ass[numberCall-1]; teller++)
               // {
#ifndef NDEBUG
               //     std::cout << _vec_objects_ass[numberCall-1][teller]<< std::endl;
#endif
               // }
               object = _vec_objects_ass[numberCall-1][i];
#ifndef NDEBUG
               //std::cout << "Considering object " << object << std::endl;
#endif
               // delete object from other possibilities if object !=0
               // probs_rest = probabilities for features if object for feature
               // 1 is chosen
#ifndef NDEBUG
               //std::cout << "probs " << probs << std::endl;
#endif
#ifndef NDEBUG
               //std::cout << "_vec_probs[numberCall-1] " << _vec_probs[numberCall-1] << std::endl;
#endif
#ifndef NDEBUG
               //std::cout << "numberCall " << numberCall << std::endl;
#endif
               for (int teller_objects_probs = 1 ; teller_objects_probs <=number_objects ; teller_objects_probs ++ )
               {
                   for (int teller_features_probs = 1 ; teller_features_probs <=number_features-1 ; teller_features_probs ++ )
                   {
                      //_probs_rest(teller_objects_probs,teller_features_probs) = probs(teller_objects_probs,teller_features_probs+1);
                      _probs_rest(teller_objects_probs,teller_features_probs) = _vec_probs[numberCall-1](teller_objects_probs,teller_features_probs+1);

                   }
                   
               }   
#ifndef NDEBUG
               //std::cout << "probs_rest " << _probs_rest << std::endl;
#endif
#ifndef NDEBUG
               //std::cout << "calculated from probs  " << probs << std::endl;
#endif
#ifndef NDEBUG
               //std::cout << "probs " << probs << std::endl;
#endif
#ifndef NDEBUG
               //std::cout << "_vec_probs[numberCall-1] " << _vec_probs[numberCall-1] << std::endl;
#endif
#ifndef NDEBUG
               std::cout << "numberCall " << numberCall << std::endl;
#endif
               // object can only be associated with one feature => put
               // _probs_rest to zero
               if (object != 0)
               {
                   for (int k = 1 ; k<=number_features - 1 ; k++)
                   {
                       _probs_rest(object,k) = -1;
                   }
               }
               // will fill in _vec_ass[numberGetAssociationsCalls]
#ifndef NDEBUG
               std::cout << "_numberGetAssociationsCalls  " << _numberGetAssociationsCalls << std::endl;
#endif
               int call_send = _numberGetAssociationsCalls + 1;
               if(call_send > _maxCalls -1)
                  cerr << "DATA ASSOCIATION ERROR: increase the number of _maxNumberCalls" << endl;
#ifndef NDEBUG
               //std::cout << "call_send  " << call_send << " for probs " << _probs_rest<< std::endl;
#endif
#ifndef NDEBUG
               //std::cout << "calculated from probs  " << probs << std::endl;
#endif
               //std::cout << "_vec_probs[numberCall-1] " << _vec_probs[numberCall-1] << std::endl;
#ifndef NDEBUG
               std::cout << "numberCall " << numberCall << std::endl;
#endif
#ifndef NDEBUG
               std::cout << "call get associations for _probs_rest" << std::endl;
#endif
               GetAssociations(_probs_rest,number_features-1);
#ifndef NDEBUG
               std::cout << "gotassociations" << std::endl;
#endif
               //int num_associations_other = _vec_number_associations[numberCall];
               int num_associations_other = _vec_number_associations[call_send-1];
#ifndef NDEBUG
               std::cout << "num_associations_other " << num_associations_other << std::endl;
               std::cout << "associations for numberCall " << numberCall << std::endl;
               std::cout << "vec_ass.size() " << _vec_ass.size() << std::endl;
#endif
               for (int teller_other = 0 ; teller_other < num_associations_other ; teller_other ++)
               {
#ifndef NDEBUG
                   cout << "teller_other " << teller_other << endl;
                   cout << "counter + teller_other " << counter + teller_other << endl;
                   cout << "_maxAssociations " << _maxAssociations << endl;
#endif
                   if(counter+teller_other > _maxAssociations-1)
                       cerr << "DATA ASSOCIATION ERROR: increase the number of maximumAssociations" << endl;
                   _vec_ass[numberCall-1][counter + teller_other][0] = object;
                   if(numberCall == 1)
                   {
                       _associations[counter + teller_other][0] = object;
                   }
                   _associations[counter + teller_other][0] = object;
#ifndef NDEBUG
                   //cout << "associations[" << counter+teller_other << "][ 0] " << _vec_ass[numberCall-1][counter + teller_other][0] << endl;
#endif
                   for (int teller_features = 1 ; teller_features <= number_features-1 ; teller_features++ )
                   {
#ifndef NDEBUG
                       cout << "teller_features " << teller_features << endl;
#endif
                       _vec_ass[numberCall-1][counter + teller_other][teller_features] = _vec_ass[call_send-1][teller_other][teller_features-1];
#ifndef NDEBUG
                       //cout << "associations[" << counter+teller_other << " ][ " << teller_features << "] " << _vec_ass[numberCall-1][counter + teller_other][teller_features] << endl;
#endif
                       if(numberCall == 1)
                       {
                           _associations[counter + teller_other][teller_features] = _vec_ass[call_send-1][teller_other][teller_features-1];
                       }
                   }
               }
               counter = counter + num_associations_other;
            }
#ifndef NDEBUG
            //std::cout << "Calculated " << counter << " associations " << std::endl;
            //std::cout << "numberCall " << numberCall << std::endl;
#endif
            //for (int i =0 ; i < counter ; i++)
            //{
#ifndef NDEBUG
            //     std::cout << "Association " << i+1 << std::endl;
#endif
            //     for (int j=0 ; j < number_features ; j++)
#ifndef NDEBUG
            //          std::cout << _vec_ass[numberCall-1][i][j] << std::endl;
#endif
            //}
            _vec_number_associations[numberCall-1] = counter;
            if(numberCall == 1)
            {
                _num_ass = counter;
#ifndef NDEBUG
                //cout << "RESETTING" << endl;
#endif
                _numberGetAssociationsCalls = 0;
            }
#ifndef NDEBUG
            //cout << "number associations " << _num_ass << endl;
            //std::cout << "END *****************************"  <<std::endl;
#endif
            //return _associations;
        }
    }

} // End namespace BFL

#endif // __DATA_ASSOCIATION_FILTER__
