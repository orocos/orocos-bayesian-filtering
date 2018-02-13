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
// Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
//

#include "EKparticlefilter.h"

#define MeasModel MeasurementModel
#define CV ColumnVector
namespace BFL
{

  EKParticleFilter::EKParticleFilter(MCPdf<CV> * prior,
				     int resampleperiod,
				     double resamplethreshold,
				     int resamplescheme)
    : ParticleFilter<ColumnVector,ColumnVector>(prior,NULL,
						resampleperiod,
						resamplethreshold,
						resamplescheme)
    , _dimension(prior->DimensionGet())
    , _num_samples(prior->NumSamplesGet())
  {
    this->_proposal_depends_on_meas = true;
    _proposal = new EKFProposalDensity(NULL,NULL);

    _sampleCov.assign(_num_samples,prior->CovarianceGet());
    _tmpCov.assign(_num_samples,prior->CovarianceGet());
    _old_samples.assign(_num_samples,WeightedSample<ColumnVector>(_dimension));
    _result_samples.assign(_num_samples,WeightedSample<ColumnVector>(_dimension));
    _unif_samples.assign(_num_samples,0.0);
    _CumPDF.assign(_num_samples,0.0);
    _x_old.resize(_num_samples);
    _sample.DimensionSet(_dimension);
  }

  EKParticleFilter::~EKParticleFilter(){
    delete _proposal;
  }

  // KG: 20040702:  This method will be exactly the same for.  Extra
  // class between?
  bool
  EKParticleFilter::UpdateInternal(SystemModel<CV>* const sysmodel,
				   const CV& u,
				   MeasModel<CV,CV>* const measmodel,
				   const CV& z,
				   const CV& s)
  {
    bool result = true;
    dynamic_cast<FilterProposalDensity *>(_proposal)->SystemModelSet(dynamic_cast<AnalyticSystemModelGaussianUncertainty *>(sysmodel));
    dynamic_cast<FilterProposalDensity *>(_proposal)->MeasurementModelSet(dynamic_cast<AnalyticMeasurementModelGaussianUncertainty *>(measmodel));
    // Proposalstep is redefined here...
    this->StaticResampleStep();
    result = this->ProposalStepInternal(sysmodel,u,measmodel,z,s) && result;
    result = this->UpdateWeightsInternal(sysmodel,u,measmodel,z,s) && result;
    this->DynamicResampleStep();

    return result;
  }

  bool
  EKParticleFilter::Resample()
  {
    // They're should be a cleaner solution then doubling the code
    // from mcpdf.h, doesn't it??
    // ONLY RIPLEY SAMPLING SUPPORTED FOR NOW!
    _old_samples = (dynamic_cast<MCPdf<CV> *>(this->_post))->ListOfSamplesGet();
    int numsamples = _old_samples.size();
    for ( int i = 0; i < numsamples ; i++) _unif_samples[i] = runif();
    _unif_samples[numsamples-1] = pow(_unif_samples[numsamples-1], double (1.0/numsamples));
    for ( int i = numsamples-2; i >= 0 ; i--){
		_unif_samples[i] = pow(_unif_samples[i],double (1.0/(i+1))) * _unif_samples[i+1];}

    unsigned int index = 0;
    _oit = _old_samples.begin();
    _CumPDF = (dynamic_cast<MCPdf<CV> *>(this->_post))->CumulativePDFGet();
    _CumPDFit = _CumPDF.begin();
    _rit = _result_samples.begin();
    _cov_it = _sampleCov.begin(); _tmpCovit = _tmpCov.begin();

    for ( int i = 0; i < numsamples ; i++)
      {
	while ( _unif_samples[i] > *_CumPDFit )
	  {
	    assert(index <= (unsigned int)numsamples);
	    index++; _oit++; _CumPDFit++; _cov_it++;

	  }
	_oit--; _cov_it--;
	*(_rit) = *(_oit);
    *_tmpCovit = *_cov_it;
	_oit++; _cov_it++;

	_rit++; _tmpCovit++;
      }

    // Update lists
    _sampleCov = _tmpCov;
    return (dynamic_cast<MCPdf<CV> *>(this->_post))->ListOfSamplesUpdate(_result_samples);
  }

  bool
  EKParticleFilter::ProposalStepInternal(SystemModel<CV>* const sysmodel,
					 const CV& u,
					 MeasModel<CV,CV>* const measmodel,
					 const CV& z,
					 const CV& s)
  {
    _old_samples = (dynamic_cast<MCPdf<CV> *>(this->_post))->ListOfSamplesGet();

    _ns_it = _new_samples.begin();
    _cov_it = _sampleCov.begin();

    for ( _os_it=_old_samples.begin(); _os_it != _old_samples.end() ; _os_it++)
      {
	_x_old = _os_it->ValueGet();

	// Set sample Covariance
	dynamic_cast<FilterProposalDensity *>(this->_proposal)->SampleCovSet(*_cov_it);

	_proposal->ConditionalArgumentSet(0,_x_old);

	if (!sysmodel->SystemWithoutInputs())
	  {
	    _proposal->ConditionalArgumentSet(1,u);
	    if (this->_proposal_depends_on_meas)
	      {
		_proposal->ConditionalArgumentSet(2,z);
		if (!measmodel->SystemWithoutSensorParams())
		  _proposal->ConditionalArgumentSet(3,s);
	      }
	  }
	else // System without inputs
	  {
	    if (this->_proposal_depends_on_meas)
	      {
		_proposal->ConditionalArgumentSet(1,z);
		if (!measmodel->SystemWithoutSensorParams())
		  _proposal->ConditionalArgumentSet(2,s);
	      }
	  }
	// Bug, make sampling method a parameter!
	_proposal->SampleFrom(_sample, SampleMthd::DEFAULT,NULL);

	_ns_it->ValueSet(_sample.ValueGet());
	_ns_it->WeightSet(_os_it->WeightGet());
	_ns_it++;

	// Update Covariances here
	*_cov_it = _proposal->CovarianceGet();
	_cov_it++;

      }

    (this->_timestep)++;
    // Update the list of samples
    return (dynamic_cast<MCPdf<CV> *>(this->_post))->ListOfSamplesUpdate(_new_samples);

  }
}
