// $Id: dataAssociationFilterMCPdf.h 
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

#ifndef __DATA_ASSOCIATION_FILTER_MCPDF__
#define __DATA_ASSOCIATION_FILTER_MCPDF__

#include "dataAssociationFilter.h"
#include "../model/systemmodel.h"
#include "../model/measurementmodel.h"
#include "../pdf/pdf.h"
#include "../pdf/mcpdf.h"
#include "../filter/filter.h"
#include <vector>
#include "../wrappers/matrix/matrix_wrapper.h"

namespace BFL
{
  using namespace std;
  using namespace MatrixWrapper;

  /// Class for data association filters in which the underlying filters are
  //particle filters
  /** This is a class that defines the data association filters
      These filters are all related to an underlying set of filters
  */
  template <typename StateVar, typename MeasVar> class DataAssociationFilterMCPdf : public DataAssociationFilter<StateVar, MeasVar>
    {

    protected:
      // vector of pointers to posteriors of filters
      std::vector< MCPdf<StateVar>* > _posts; 
      // iterator for vector of pointers to posteriors of filters
      typename std::vector<MCPdf<StateVar>* >::iterator _iter_posts; 
      // vector containing prob(z|particle) for all particles
      vector<Matrix> _prob_particles;
      bool  _measProbsCalculated; 
      vector<MeasVar> _measMeasProbsCalculated;

      int _maxParticles; 



      /// Implementation of Update
      // calls update on the underlying filters
      /** @param sysmodel pointer to the used system model
	  @param u input param for proposal density
	  @param measmodel pointer to the used measurementmodel
	  @param z measurement param for proposal density
	  @param s sensor param for proposal density
      */
      bool UpdateInternal(SystemModel<StateVar>* const sysmodel,
				  const StateVar& u,
				  MeasurementModel<MeasVar,StateVar>* const measmodel,
				  const vector<MeasVar>& z,
				  const StateVar& s);

    public:
      /// Constructor
      /** @pre you created the prior
	  @param prior pointer to the prior Pdf
      */
      DataAssociationFilterMCPdf (vector<Filter<StateVar,MeasVar> *> filters, double gamma=0.0, double treshold=0.0);

      /// copy constructor
      DataAssociationFilterMCPdf (const DataAssociationFilterMCPdf<StateVar,MeasVar> & filters);

      /// destructor
      virtual ~DataAssociationFilterMCPdf();

      // Get the probabilities of the measurements for each of the
      // filters TODO: set protected again
      MatrixWrapper::Matrix GetMeasProbs(MeasurementModel<MeasVar,StateVar>* const measmodel , const vector<MeasVar>& z, const StateVar& s);

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

    };


////////////// IMPLEMENTATION/////////////////////////////////////

    #define StateVar SVar
    #define MeasVar MVar

    // Constructor
    template<typename SVar, typename MVar> 
    DataAssociationFilterMCPdf<SVar,MVar>::DataAssociationFilterMCPdf(vector< Filter<SVar,MVar> * > filters, double gamma, double treshold)
      : DataAssociationFilter<SVar,MVar>(filters,gamma,treshold)
    {
      typename std::vector< Filter<StateVar,MeasVar>* >::iterator iter_filters; 
      _posts.resize(this->_maxFilters);
      _iter_posts = _posts.begin();
      for(iter_filters = filters.begin() ; iter_filters != filters.end() ; iter_filters++) 
      {
        (*_iter_posts) = (MCPdf<SVar>*)((*iter_filters)->PostGet());
        _iter_posts++;
      }
      _iter_posts = _posts.begin();

      _maxParticles = 5000; //maximum number of features
     _prob_particles.resize(this->_maxFilters);
     _prob_particles.assign(this->_maxFilters,Matrix(_maxParticles,this->_maxFeatures));
    }
     
    // Copy Constructor
    template<typename SVar, typename MVar> 
    DataAssociationFilterMCPdf<SVar,MVar>::DataAssociationFilterMCPdf(const DataAssociationFilterMCPdf& filters)
      : DataAssociationFilter<SVar,MVar>(filters)   
    {
      _maxParticles = 5000; //maximum number of features
      _prob_particles(this->_maxFilters,Matrix(_maxParticles,this->_maxFeatures));
      for(this->_iter_filters = this->_filters.begin() ; this->_iter_filters != this->_filters.end() ; this->_iter_filters++) 
      {
        *_iter_posts = (*this->_iter_filters)->PostGet();
        _iter_posts++;
      }
      _iter_posts = _posts.begin();
    }

    template<typename SVar, typename MVar> 
    DataAssociationFilterMCPdf<SVar,MVar>::~DataAssociationFilterMCPdf()
    {}
    

    template<typename SVar, typename MVar>  void
    DataAssociationFilterMCPdf<SVar,MVar>::AddFilter(Filter<SVar,MVar>* filter)
    {
        DataAssociationFilter<SVar,MVar>::AddFilter(filter);
        //_prob_particles.resize(this->NumFiltersGet());
        //_posts.resize(this->NumFiltersGet());
        _posts[this->NumFiltersGet()-1] = (MCPdf<SVar>*)filter->PostGet();
        _measProbsCalculated = false;
    }

    template<typename SVar, typename MVar>  bool
    DataAssociationFilterMCPdf<SVar,MVar>::RemoveFilter(int index)
    {
        DataAssociationFilter<SVar,MVar>::RemoveFilter(index);
        int size = this->NumFiltersGet();
        if(index< 0 || index > size-1)
        {
            return false;
        }
        else
        {
            //_prob_particles.resize(this->NumFiltersGet());
            // get iterator to the index'th element
            _iter_posts = _posts.begin();
            for(int i = 0 ; i<index ; i++ )
            {
                _iter_posts++;
            }
            _posts.erase(_iter_posts);
            _posts.resize(this->_maxFilters);
            return true;
        }
        _measProbsCalculated = false;
    }

    template<typename SVar, typename MVar> MatrixWrapper::Matrix 
    DataAssociationFilterMCPdf<SVar,MVar>::GetMeasProbs(MeasurementModel<MeasVar,StateVar>* const measmodel , const vector<MeasVar>& z, const StateVar& s)
    {
      bool measSame = true;
      if (! _measProbsCalculated)
      {
            // test if measurements same as before
            measSame = (measSame && (_measMeasProbsCalculated.size() == z.size()));
            for(int i = 0 ; i < _measMeasProbsCalculated.size(); i++)
                measSame = (measSame && (_measMeasProbsCalculated[i] == z[i]));
      }
      if(!(_measProbsCalculated && measSame )) 
      {
        // probabilities not calculated yet
        //cout << "probabilies not calculated yet" << endl;

        //cout << "enter GetMeasProbs" << endl;
        int number_features = z.size();
        int number_objects = this->NumFiltersGet();
        //cout << "number objects " << number_objects << endl;
        //cout << "number features " << number_features << endl;
        //cout << "size posts " << _posts.size() << endl;

        //if (!( (_probs.columns()== number_features) && (_probs.rows()==number_objects) ))
        //{ 
        //      _probs.resize(number_objects,number_features);
        //      cout << "PROBS RESIZED" << endl;
        //}

        int num_particles = 0;
        vector< WeightedSample< SVar > > los;
    
        double inter = 0.0;
        Probability prob;

        this->_iter_filters = this->_filters.begin();
        _iter_posts = _posts.begin();
        for(int object = 1 ; object != this->NumFiltersGet() + 1 ; object++) 
        {
          //cout << "object " << object << endl;
          // loop over the different filters (i.e. different objects)
          num_particles =  (*_iter_posts)->NumSamplesGet();
          //cout << "number particles  " << num_particles << endl;
          los = (*_iter_posts)->ListOfSamplesGet() ;

          // do only assignement if necessary
          //if (!(num_particles == _prob_particles[object-1].rows() && number_features == _prob_particles[object-1].columns()))
          //{
          //    _prob_particles[object-1].resize(num_particles,number_features); 
          //    cout << "PROB_PARTICLES RESIZED" << endl;
          //}
          //
          for (int feature = 1 ; feature<=number_features ; feature ++ )
          {
              inter = 0.0;
              //cout << "feature " << z[feature-1] << endl;
              for (int teller_part = 0; teller_part < num_particles; teller_part++)
              {
                  //cout << "particle " << los[teller_part].ValueGet()<< endl;
                  //cout << "weight " << los[teller_part].WeightGet()<< endl;
                  if (measmodel->SystemWithoutSensorParams())
                      prob =  measmodel->ProbabilityGet(z[feature-1], los[teller_part].ValueGet());
                  else
                      prob =  measmodel->ProbabilityGet(z[feature-1], los[teller_part].ValueGet(), s);
                  _prob_particles[object-1](teller_part+1,feature) = prob;
                  //cout << "prob " << prob << endl;
                  inter = inter + (double)prob * los[teller_part].WeightGet();
                  //cout << "inter " << inter << endl;
              }
              //cout << "inter " << inter << endl;
              this->_probs(object,feature) = inter ;
              //cout << "prob " << _probs(object,feature) << endl;
          }
          //cout << "prob_particles " << _prob_particles[object-1] << endl;
          _iter_posts++;
          this->_iter_filters++;
        }
        _measProbsCalculated = true; 
        _measMeasProbsCalculated = z;
        //cout << "PROBS " << probs << endl;
      }
      else
      {
        //cout << "probabilies were already calculated" << endl;
      }
      //return this->_probs.sub(1,number_objects,1,number_features);
      return this->_probs;
    }

    template<typename SVar, typename MVar> bool
    DataAssociationFilterMCPdf<SVar,MVar>::UpdateInternal(SystemModel<SVar>* const sysmodel,
    			  const SVar& u,
    			  MeasurementModel<MVar,SVar>* const measmodel,
    			  const vector<MVar>& z,
    			  const SVar& s)
    {
      //cout << "entering update internal" << endl;
      bool return_bool = true;

      if(sysmodel!=NULL)
      {
         // system update
         _iter_posts = _posts.begin();
         this->_iter_filters = this->_filters.begin();
         for( int object = 0 ; object < this->NumFiltersGet() ; object++)
         {
           return_bool = (return_bool && (*this->_iter_filters)->Update(sysmodel,u) ); 
           *_iter_posts = (MCPdf<SVar>*)( (*this->_iter_filters)->PostGet() );
           _iter_posts++;
           this->_iter_filters++;
         }
      }

      if(measmodel!=NULL)
      {
        // cout << "entering meas update" << endl;
        // measurement update
        vector<vector<double> > association_probs = this->GetAssociationProbs(measmodel , z, s);
        
        int object = 0;
        double temp = 0.0;
        _iter_posts = _posts.begin();
        this->_iter_filters = this->_filters.begin();
        for( int object = 0 ; object < this->NumFiltersGet() ; object++)
        {
          //cout << "expected value" <<(*DataAssociationFilter<SVar,MVar>::_iter_filters)->PostGet()->ExpectedValueGet() << endl;
          //cout << "covariance" <<(*DataAssociationFilter<SVar,MVar>::_iter_filters)->PostGet()->CovarianceGet()<< endl;
          vector<WeightedSample<SVar> > particles = (*_iter_posts)->ListOfSamplesGet();
          for (int particle = 0; particle < (*_iter_posts)->NumSamplesGet() ; particle++)
          {
            //cout << "object " << object << endl;
             temp = 0.0;
             for(int feature = 0; feature < z.size(); feature ++)
             {
                  // only select measurements that have a minimum probability
                  if (association_probs[feature][object]>0.00005)
                  {
                      temp = temp + association_probs[feature][object] * _prob_particles[object](particle+1,feature+1) * particles[particle].WeightGet();
                  }
             }
            if (temp>0.00005)
                particles[particle].WeightSet(temp);
          }
          return_bool = ( return_bool && ((MCPdf<SVar>*)( (*(this->_iter_filters))->PostGet()))->ListOfSamplesUpdate(particles) );
          // listOfSamplesSet takes care of normalization
          *_iter_posts = (MCPdf<SVar>*) ((*this->_iter_filters)->PostGet() );
          // sets posterior too
          _iter_posts++;
          this->_iter_filters++;
        }
      }
    //cout << "leaving update internal" << endl;
    _measProbsCalculated = false; 
    return return_bool;
    }

} // End namespace BFL

#endif // __DATA_ASSOCIATION_FILTER_MCPDF__
