// $Id$
// Copyright (C) 2003 Klaas Gadeyne <first dot last at gmail dot com>
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

#include "filterproposaldensity.h"

namespace BFL
{
  using namespace MatrixWrapper;


#define FilterPropDens FilterProposalDensity

  FilterPropDens::FilterPropDens(AnalyticSystemModelGaussianUncertainty * SysModel,
				 AnalyticMeasurementModelGaussianUncertainty * MeasModel)
    : AnalyticConditionalGaussian(),
      _sysmodel(SysModel),
      _measmodel(MeasModel)
  {
    if (SysModel != NULL)
      {
	_TmpPrior = new Gaussian(SysModel->StateSizeGet());
	_sample_cov.resize(SysModel->StateSizeGet());

	this->DimensionSet(SysModel->StateSizeGet());
	if (MeasModel != NULL)
	  this->NumConditionalArgumentsSet(SysModel->SystemPdfGet()->NumConditionalArgumentsGet()
					   + MeasModel->MeasurementPdfGet()->NumConditionalArgumentsGet());
      }
    else
      {
	_TmpPrior = new Gaussian();
      }

    _sysmodel = SysModel;
    _measmodel = MeasModel;

    /*  _filter remains uninitialised here!!!
	This is Only an interface class!  Or should we make an
	enumerate or something else here, like some nameserved stuff?
    */
  }

  FilterPropDens::~FilterPropDens(){}

  // BUG Copy constructor not implemented yet
  FilterPropDens::FilterPropDens(const FilterPropDens & fpd){}

  void
  FilterPropDens::SystemModelSet(AnalyticSystemModelGaussianUncertainty * SysModel)
  {
    assert ( SysModel != NULL );
    assert ( (this->DimensionGet() == 0) || (this->DimensionGet() == (unsigned int)SysModel->StateSizeGet()) );
    if ((this->DimensionGet() == 0))
      {
	_TmpPrior->DimensionSet(SysModel->StateSizeGet());
	_sample_cov.resize(SysModel->StateSizeGet());
      }
    this->DimensionSet(SysModel->StateSizeGet());
    if (_measmodel != NULL)
      this->NumConditionalArgumentsSet(SysModel->SystemPdfGet()->NumConditionalArgumentsGet()
				       + _measmodel->MeasurementPdfGet()->NumConditionalArgumentsGet());
    _sysmodel = SysModel;
  }


  void
  FilterPropDens::MeasurementModelSet(AnalyticMeasurementModelGaussianUncertainty * MeasModel)
  {
    assert ( MeasModel != NULL );
    if (_sysmodel != NULL)
      this->NumConditionalArgumentsSet(_sysmodel->SystemPdfGet()->NumConditionalArgumentsGet()
				       + MeasModel->MeasurementPdfGet()->NumConditionalArgumentsGet());
    _measmodel = MeasModel;
  }

  void
  FilterPropDens::SampleCovSet(SymmetricMatrix & cov)
  {
    assert (cov.rows() == this->DimensionGet());
    _sample_cov = cov;
  }

  ColumnVector
  FilterPropDens::ExpectedValueGet() const
  {
    this->FilterStep();
    return (this->_filter->PostGet()->ExpectedValueGet());
  }

  SymmetricMatrix
  FilterPropDens::CovarianceGet() const
  {
    this->FilterStep();
    return (_filter->PostGet()->CovarianceGet());
  }

  void
  FilterPropDens::FilterStep() const
  {
    this->_TmpPrior->ExpectedValueSet(this->ConditionalArgumentGet(0));
    // See above, last argument is previous covariance Matrix
    this->_TmpPrior->CovarianceSet(_sample_cov);
    this->_filter->Reset(_TmpPrior);
    // Updatestep
    if ( _sysmodel == NULL )
      {
	if (_measmodel->SystemWithoutSensorParams() == false){ // 2 conditional arguments
	  _filter->Update(_measmodel,this->ConditionalArgumentGet(1),this->ConditionalArgumentGet(2));}
	else{
	  _filter->Update(_measmodel,this->ConditionalArgumentGet(1));}
      }
    else if ( _measmodel == NULL )
      {
	if (_sysmodel->SystemWithoutInputs() == false){// 2 conditional arguments
	  _filter->Update(_sysmodel,this->ConditionalArgumentGet(1));}
	else _filter->Update(_sysmodel);
      }
    else
      {
	if ( ( _sysmodel ->SystemWithoutInputs() == false) && (_measmodel->SystemWithoutSensorParams() == false) ){
	  _filter->Update(_sysmodel,this->ConditionalArgumentGet(1),
			  _measmodel,this->ConditionalArgumentGet(2),this->ConditionalArgumentGet(3));}
	else if ( ( _sysmodel->SystemWithoutInputs() == true) && (_measmodel->SystemWithoutSensorParams() == false) ){
	  _filter->Update(_sysmodel,_measmodel,this->ConditionalArgumentGet(1),this->ConditionalArgumentGet(2));}
	else if ( (_sysmodel->SystemWithoutInputs() == false) && (_measmodel->SystemWithoutSensorParams() == true) ){
	  _filter->Update(_sysmodel,this->ConditionalArgumentGet(1),
			  _measmodel,this->ConditionalArgumentGet(2));}
	else // No inputs, no sensor parameters
	  _filter->Update(_sysmodel,_measmodel,this->ConditionalArgumentGet(1));
      }
  }

  Matrix
  FilterPropDens::dfGet(unsigned int i) const
  {
    cerr << "FilterPropDens::dfGet() never necessary?" << endl;
    exit(-BFL_ERRMISUSE);
  }

} // End namespace BFL
