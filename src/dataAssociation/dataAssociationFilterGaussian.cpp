#include "dataAssociationFilterGaussian.h"

namespace BFL
{
  using namespace std;
  using namespace MatrixWrapper;

    // Constructor
    DataAssociationFilterGaussian::DataAssociationFilterGaussian(vector< Filter<ColumnVector,ColumnVector> * > filters, double gamma, double treshold)
      : DataAssociationFilter<ColumnVector,ColumnVector>(filters,gamma,treshold)
    {
     std::vector< Filter<ColumnVector,ColumnVector>* >::iterator iter_filters; 
      _posts.resize(this->_maxFilters);
      _iter_posts = _posts.begin();
      for(iter_filters = filters.begin() ; iter_filters != filters.end() ; iter_filters++) 
      {
        (*_iter_posts) = (Gaussian*)((*iter_filters)->PostGet());
        _iter_posts++;
      }
      _iter_posts = _posts.begin();
      if(this->NumFiltersGet()==0)
      {
      }
      else
      {
        int dimension = (*_iter_posts)->DimensionGet();
        _x.resize(dimension);
        _Mu_new.resize(dimension);
        _Sigma_temp.resize(dimension,dimension);
        _Sigma_temp2.resize(dimension,dimension);
        _Sigma_temp_par.resize(dimension,dimension);
        _Sigma_new.resize(dimension);
      }
    }
     
//    // Copy Constructor
//    DataAssociationFilterGaussian::DataAssociationFilterGaussian(const DataAssociationFilterGaussian& filters)
//      : DataAssociationFilter<ColumnVector,ColumnVector>(filters)   
//    {
//      for(this->_iter_filters = this->_filters.begin() ; this->_iter_filters != this->_filters.end() ; this->_iter_filters++) 
//      {
//        *_iter_posts = (Gaussian*)(*this->_iter_filters)->PostGet();
//        _iter_posts++;
//      }
//      _iter_posts = _posts.begin();
//    }

    DataAssociationFilterGaussian::~DataAssociationFilterGaussian()
    {}
    

    void
    DataAssociationFilterGaussian::AddFilter(Filter<ColumnVector,ColumnVector>* filter)
    {
        //cout << "entered AddFilter" << endl;
        DataAssociationFilter<ColumnVector,ColumnVector>::AddFilter(filter);
        //cout << "baseclass add called" << endl;
        _posts[this->NumFiltersGet()-1] = (Gaussian*)filter->PostGet();
        //cout << "posts set" << endl;
        _measProbsCalculated = false;
        int dimension = _posts[0]->DimensionGet();
        //cout << "dimension " << dimension  << endl;
        _x.resize(dimension);
        //cout << "_x resized " << endl;
        _Mu_new.resize(dimension);
        //cout << "_Mu_new resized " << endl;
        _Sigma_temp.resize(dimension,dimension);
        //cout << "_Sigma_temp resized " << endl;
        _Sigma_temp2.resize(dimension,dimension);
        //cout << "_Sigma_temp2 resized " << endl;
        _Sigma_temp_par.resize(dimension,dimension);
        //cout << "_Sigma_temp_par resized " << endl;
        //_Sigma_new.resize(dimension);
        //cout << " _measProbsCalculated " << _measProbsCalculated << endl; 
        //cout << "finished AddFilter" << endl;
    }

    bool
    DataAssociationFilterGaussian::RemoveFilter(int index)
    {
        DataAssociationFilter<ColumnVector,ColumnVector>::RemoveFilter(index);
        int size = this->NumFiltersGet();
        if(index< 0 || index > size-1)
        {
            return false;
        }
        else
        {
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

    MatrixWrapper::Matrix 
    DataAssociationFilterGaussian::GetMeasProbs(MeasurementModel<ColumnVector,ColumnVector>* const measmodel , const vector<ColumnVector>& z, const ColumnVector& s)
    {
      //cout << "measProbsCalculated " << _measProbsCalculated << endl;
      bool measSame = true;
      if (! _measProbsCalculated)
      {
            // test if measurements same as before
            measSame = (measSame && (_measMeasProbsCalculated.size() == z.size()));
            for(int i = 0 ; i < _measMeasProbsCalculated.size(); i++)
                measSame = (measSame && (_measMeasProbsCalculated[i] == z[i]));
      }
      //cout << "measSame " << measSame << endl;
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

        int num_particles = 0;
    
        double inter = 0.0;
        Probability prob;

        AllocateMeasModel(z[0].rows());
        (_mapMeasUpdateVariables_it->second)._Hmeas = ((LinearAnalyticMeasurementModelGaussianUncertainty*)measmodel)->HGet();
        //cout << "_Hmeas.rows" <<  _Hmeas.rows() << endl;
        //cout << "_Hmeas.columns" <<  _Hmeas.columns() << endl;
        (_mapMeasUpdateVariables_it->second)._Covmeas = ((LinearAnalyticMeasurementModelGaussianUncertainty*)measmodel)->MeasurementPdfGet()->CovarianceGet();

        this->_iter_filters = this->_filters.begin();
        _iter_posts = _posts.begin();
        //cout << "this->NumFiltersGet()" << this->NumFiltersGet() << endl;
        for(int object = 1 ; object != this->NumFiltersGet() + 1 ; object++) 
        {
          //cout << "object " << object << endl;
          // loop over the different filters (i.e. different objects)
          //
          for (int feature = 1 ; feature<=number_features ; feature ++ )
          {
              // create new Gaussian
              (_mapMeasUpdateVariables_it->second)._pyMu = (_mapMeasUpdateVariables_it->second)._Hmeas * ((*_iter_posts)->ExpectedValueGet());
              Matrix temp = (_mapMeasUpdateVariables_it->second)._Hmeas * (Matrix)((*_iter_posts)->CovarianceGet());
              temp = temp * (_mapMeasUpdateVariables_it->second)._Hmeas.transpose();
              temp.convertToSymmetricMatrix((_mapMeasUpdateVariables_it->second)._pySigma);
              //_pySigma = _Hmeas * (Matrix)((*_iter_posts)->CovarianceGet()) * _Hmeas.transpose();
              (_mapMeasUpdateVariables_it->second)._pySigma += (_mapMeasUpdateVariables_it->second)._Covmeas;
              _py.DimensionSet(z[0].rows());
              //cout << "_pyMu.rows" <<  _pyMu.rows() << endl;
              _py.ExpectedValueSet((_mapMeasUpdateVariables_it->second)._pyMu);
              _py.CovarianceSet((_mapMeasUpdateVariables_it->second)._pySigma);

              prob =  _py.ProbabilityGet(z[feature-1]);

              this->_probs(object,feature) = prob ;
              //cout << "prob " << _probs(object,feature) << endl;
          }
          //cout << "prob_particles " << _prob_particles[object-1] << endl;
          _iter_posts++;
          this->_iter_filters++;
        }
        _measProbsCalculated = true; 
        _measMeasProbsCalculated = z;
      }
      else
      {
        //cout << "probabilies were already calculated" << endl;
      }
      //return this->_probs.sub(1,number_objects,1,number_features);
      //cout << "PROBS " << _probs << endl;
      return this->_probs;
    }

    bool
    DataAssociationFilterGaussian::UpdateInternal(SystemModel<ColumnVector>* const sysmodel,
    			  const ColumnVector& u,
    			  MeasurementModel<ColumnVector,ColumnVector>* const measmodel,
    			  const vector<ColumnVector>& z,
    			  const ColumnVector& s)
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
           *_iter_posts = (Gaussian*)( (*this->_iter_filters)->PostGet() );
           _iter_posts++;
           this->_iter_filters++;
         }
      }

      if(measmodel!=NULL)
      {
        //cout << "entering meas update" << endl;
        // measurement update
        //vector<vector<double> > association_probs = this->GetAssociationProbs(measmodel , z, s);
        this->GetAssociationProbs(measmodel , z, s);
        // writes _association_probs
        //cout << "got association probs" << endl;
        
        int object = 0;
        double sum_association_probs = 0.0;
        _iter_posts = _posts.begin();
        this->_iter_filters = this->_filters.begin();

        for( int object = 0 ; object < this->NumFiltersGet() ; object++)
        {
           //cout << "object " << object << endl;
           // allocate measurement for z.rows() if needed
           AllocateMeasModel(z[0].rows());

           (_mapMeasUpdateVariables_it->second)._innov = 0.0;
           sum_association_probs = 0.0;
           _x = (*_iter_posts)->ExpectedValueGet();
           //cout << "_x " << _x << endl;

           (_mapMeasUpdateVariables_it->second)._zpred = ((AnalyticMeasurementModelGaussianUncertainty*)measmodel)->PredictionGet(s,_x);
           // cout << "_zpred " << (_mapMeasUpdateVariables_it->second)._zpred << endl;
           // calculate weighted innovation nu_i = sum_j beta_ji * (zj-zjpred)
           for (int feature = 0 ; feature<z.size() ; feature ++ )
           {
                //cout << "feature " << feature << endl;
                //cout << "_z " << z[feature]<< endl;
                (_mapMeasUpdateVariables_it->second)._innov += (z[feature]-(_mapMeasUpdateVariables_it->second)._zpred) * this->_association_probs[feature][object];
                //cout << "_innov " << (_mapMeasUpdateVariables_it->second)._innov << endl;
                sum_association_probs+=this->_association_probs[feature][object];
           }
           //cout << "sum_association_probs " << sum_association_probs << endl;
           // check if there are really measurements associated with this object
           if (sum_association_probs > 0.00005)
           {
                (_mapMeasUpdateVariables_it->second)._innov/=sum_association_probs;
                //cout << "_innov " << (_mapMeasUpdateVariables_it->second)._innov << endl;
    
                //update post with the weighted innovation
                (_mapMeasUpdateVariables_it->second)._postHT = (Matrix)((*_iter_posts)->CovarianceGet()) * (_mapMeasUpdateVariables_it->second)._Hmeas.transpose() ;
                (_mapMeasUpdateVariables_it->second)._S =  (_mapMeasUpdateVariables_it->second)._Hmeas  * (_mapMeasUpdateVariables_it->second)._postHT;
                (_mapMeasUpdateVariables_it->second)._S += (Matrix)(_mapMeasUpdateVariables_it->second)._Covmeas;
        
                 // _K = covariance * H' * S(-1)
                 (_mapMeasUpdateVariables_it->second)._K = (_mapMeasUpdateVariables_it->second)._postHT * ((_mapMeasUpdateVariables_it->second)._S.inverse());
                  
                // calcutate new state gaussian
                // Mu = expectedValue + K*(z-Z)
                 _Mu_new  =  (_mapMeasUpdateVariables_it->second)._K * (_mapMeasUpdateVariables_it->second)._innov  ;
                 _Mu_new  +=  (*_iter_posts)->ExpectedValueGet() ;
                // Sigma = post - K*H*post
                _Sigma_temp = (Matrix) ( (*_iter_posts)->CovarianceGet());
                _Sigma_temp_par = (_mapMeasUpdateVariables_it->second)._K * (_mapMeasUpdateVariables_it->second)._Hmeas ;
                _Sigma_temp2 = _Sigma_temp_par * _Sigma_temp;
                _Sigma_temp -=  _Sigma_temp2 ;
                // convert to symmetric matrix
                _Sigma_temp.convertToSymmetricMatrix(_Sigma_new);
                
                // set new state gaussian
                (*_iter_posts)->ExpectedValueSet( _Mu_new );
                (*_iter_posts)->CovarianceSet( _Sigma_new );
            }
 
            // sets posterior too
            _iter_posts++;
            this->_iter_filters++;
        }
      }
    //cout << "leaving update internal" << endl;
    _measProbsCalculated = false; 
    return return_bool;
    }

    void
    DataAssociationFilterGaussian::AllocateMeasModel(const vector<unsigned int>& meas_dimensions)
    {
      unsigned int meas_dimension;
      for(int i = 0 ; i< meas_dimensions.size(); i++)
      {
          // find if variables with size meas_sizes[i] are already allocated
          meas_dimension = meas_dimensions[i];
          _mapMeasUpdateVariables_it =  _mapMeasUpdateVariables.find(meas_dimension);
          if( _mapMeasUpdateVariables_it == _mapMeasUpdateVariables.end())
          {
              //variables with size z.rows() not allocated yet
              _mapMeasUpdateVariables_it = (_mapMeasUpdateVariables.insert
                  (std::pair<unsigned int, MeasUpdateVariables>( meas_dimension,MeasUpdateVariables(meas_dimension,_Mu_new.rows()) ))).first;
           }
       }
    }
    
    void
    DataAssociationFilterGaussian::AllocateMeasModel(const unsigned int& meas_dimension)
    {
       // find if variables with size meas_sizes[i] are already allocated
       _mapMeasUpdateVariables_it =  _mapMeasUpdateVariables.find(meas_dimension);
       if( _mapMeasUpdateVariables_it == _mapMeasUpdateVariables.end())
       {
           //variables with size z.rows() not allocated yet
           _mapMeasUpdateVariables_it = (_mapMeasUpdateVariables.insert
               (std::pair<unsigned int, MeasUpdateVariables>( meas_dimension,MeasUpdateVariables(meas_dimension,_Mu_new.rows()) ))).first;
        }
    }

} // End namespace BFL
