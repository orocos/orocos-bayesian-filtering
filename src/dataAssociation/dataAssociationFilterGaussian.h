// $Id: dataAssociationFilterGaussian.h 
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

#ifndef __DATA_ASSOCIATION_FILTER_GAUSSIAN__
#define __DATA_ASSOCIATION_FILTER_GAUSSIAN__

#include "dataAssociationFilter.h"
#include "../model/systemmodel.h"
#include "../model/measurementmodel.h"
#include "../model/linearanalyticmeasurementmodel_gaussianuncertainty.h"
#include "../pdf/pdf.h"
#include "../pdf/gaussian.h"
#include "../filter/filter.h"
#include <vector>
#include <map>
#include "../wrappers/matrix/matrix_wrapper.h"

namespace BFL
{
  using namespace std;
  using namespace MatrixWrapper;
  struct MeasUpdateVariables
  {
     Matrix _S;
     Matrix _postHT;
     Matrix _K;
     Matrix _Hmeas;
     ColumnVector _innov;
     ColumnVector _zpred;
     SymmetricMatrix _Covmeas;
     ColumnVector          _pyMu;
     SymmetricMatrix       _pySigma;
     MeasUpdateVariables() {};
     MeasUpdateVariables(unsigned int meas_dimension, unsigned int state_dimension):
       _S(meas_dimension,meas_dimension)
     , _K(state_dimension,meas_dimension)
     , _innov(meas_dimension)
     , _zpred(meas_dimension)
     , _pyMu(meas_dimension)
     , _postHT(state_dimension,meas_dimension)
     , _Hmeas(meas_dimension,state_dimension)
     , _Covmeas(meas_dimension)
     , _pySigma(meas_dimension)
        {};
  }; //struct

  /// Class for data association filters in which the underlying filters are
  //Kalman filters
  /** This is a class that defines the data association filters
      These filters are all related to an underlying set of filters
  */
  class DataAssociationFilterGaussian : public DataAssociationFilter<ColumnVector, ColumnVector>
    {

    protected:
      // vector of pointers to posteriors of filters
      std::vector< Gaussian* > _posts; 
      // iterator for vector of pointers to posteriors of filters
      std::vector<Gaussian* >::iterator _iter_posts; 
      bool                  _measProbsCalculated; 
      vector<ColumnVector>  _measMeasProbsCalculated;
      Gaussian              _py;
      ColumnVector          _x;
      ColumnVector          _Mu_new;
      Matrix                _Sigma_temp;
      Matrix                _Sigma_temp2;
      Matrix                _Sigma_temp_par;
      SymmetricMatrix       _Sigma_new;

      std::map<unsigned int, MeasUpdateVariables> _mapMeasUpdateVariables;
      std::map<unsigned int, MeasUpdateVariables>::iterator _mapMeasUpdateVariables_it;



      /// Implementation of Update
      // calls update on the underlying filters
      /** @param sysmodel pointer to the used system model
	  @param u input param for proposal density
	  @param measmodel pointer to the used measurementmodel
	  @param z measurement param for proposal density
	  @param s sensor param for proposal density
      */
      bool UpdateInternal(SystemModel<ColumnVector>* const sysmodel,
				  const ColumnVector& u,
				  MeasurementModel<ColumnVector,ColumnVector>* const measmodel,
				  const vector<ColumnVector>& z,
				  const ColumnVector& s);

      void AllocateMeasModel(const vector<unsigned int>& meas_dimensions);
      void AllocateMeasModel(const unsigned int& meas_dimension);

    public:
      /// Constructor
      /** @pre you created the prior
	  @param prior pointer to the prior Pdf
      */
      DataAssociationFilterGaussian (vector<Filter<ColumnVector,ColumnVector> *> filters, double gamma=0.0, double treshold=0.0);

      /// copy constructor
      //DataAssociationFilterGaussian (const DataAssociationFilterGaussian & filters);

      /// destructor
      virtual ~DataAssociationFilterGaussian();

      /// Get the probabilities of the measurements for each of the filters 
      /** @param sysmodel pointer to the used system model
	  @param measmodel pointer to the used measurementmodel
	  @param z vector of measurements 
	  @param s sensor param for proposal density
      @return Matrix containing the probability of each of the measurements given each of the filters
      */
      MatrixWrapper::Matrix GetMeasProbs(MeasurementModel<ColumnVector,ColumnVector>* const measmodel , const vector<ColumnVector>& z, const ColumnVector& s);

      /// Add a filter
      /** @param filter pointer to the filter that will be added
      */
      void AddFilter(Filter<ColumnVector,ColumnVector>* filter);

      /// Remove a filter
      /** @param index index of the filter that will be removed. This index
         should be between 0 and NumFitlersGet()-1, else false is returned
          @return bool indicating of the removal was succesfull or not
      */
      bool RemoveFilter(int index);

    };


} // End namespace BFL

#endif // __DATA_ASSOCIATION_FILTER_GAUSSIAN__
