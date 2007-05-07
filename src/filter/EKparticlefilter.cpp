// Copyright (C) 2003 Klaas Gadeyne <klaas dot gadeyne at mech dot
// kuleuven dot ac dot be>
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
// $Id$

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
  {
    this->_proposal_depends_on_meas = true;
    _proposal = new EKFProposalDensity(NULL,NULL);

    // Is this correct?
    SymmetricMatrix tmp(prior->DimensionGet()); tmp = prior->CovarianceGet();
    _sampleCov.resize(prior->NumSamplesGet(), tmp);
    _tmpCov.resize(prior->NumSamplesGet(), tmp);
  }

  EKParticleFilter::~EKParticleFilter(){}

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
    vector<WeightedSample<CV> > old_samples = (dynamic_cast<MCPdf<CV> *>(this->_post))->ListOfSamplesGet();
    int numsamples = old_samples.size();
    vector<Sample<CV> > result(numsamples);
    std::vector<double> unif_samples(numsamples);

    for ( int i = 0; i < numsamples ; i++) unif_samples[i] = runif();
    unif_samples[numsamples-1] = pow(unif_samples[numsamples-1], double (1.0/numsamples));
    for ( int i = numsamples-2; i >= 0 ; i--){
		unif_samples[i] = pow(unif_samples[i],double (1.0/(i+1))) * unif_samples[i+1];}

    unsigned int index = 0;
    vector<WeightedSample<CV> >::const_iterator it = old_samples.begin();
    vector<double> CumPDF = (dynamic_cast<MCPdf<CV> *>(this->_post))->CumulativePDFGet();
    vector<double>::const_iterator CumPDFit = CumPDF.begin();
    vector<Sample<CV> >::iterator sit = result.begin();
    _cov_it = _sampleCov.begin(); _tmpCovit = _tmpCov.begin();
    
    for ( int i = 0; i < numsamples ; i++)
      {
	while ( unif_samples[i] > *CumPDFit )
	  {
	    assert(index <= (unsigned int)numsamples);
	    index++; it++; CumPDFit++; _cov_it++;
	  }
	it--; _cov_it--;
	*sit = *it; *_tmpCovit = *_cov_it;
	it++; _cov_it++;
	
	sit++; _tmpCovit++;
      }
    
    // Update lists
    _sampleCov = _tmpCov;
    return (dynamic_cast<MCPdf<CV> *>(this->_post))->ListOfSamplesUpdate(result);
  }

  bool
  EKParticleFilter::ProposalStepInternal(SystemModel<CV>* const sysmodel,
					 const CV& u,
					 MeasModel<CV,CV>* const measmodel,
					 const CV& z,
					 const CV& s)
  {
    _old_samples = (dynamic_cast<MCPdf<CV> *>(this->_post))->ListOfSamplesGet();
    CV x_old;
    Sample<CV> sample;

    _ns_it = _new_samples.begin();
    _cov_it = _sampleCov.begin();
  
    for ( _os_it=_old_samples.begin(); _os_it != _old_samples.end() ; _os_it++)
      {
	x_old = _os_it->ValueGet();

	// Set sample Covariance
	dynamic_cast<FilterProposalDensity *>(this->_proposal)->SampleCovSet(*_cov_it);
      
	_proposal->ConditionalArgumentSet(0,x_old);

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
	_proposal->SampleFrom(sample, DEFAULT,NULL);

	_ns_it->ValueSet(sample.ValueGet());
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
