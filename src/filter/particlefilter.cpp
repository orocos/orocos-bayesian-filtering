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
// $Id$

#include "particlefilter.h"
#include "../pdf/mcpdf.h"

#define SV StateVar
#define MV MeasVar
#define MeasModel MeasurementModel

#define STATE_AND_MEAS_VAR_DIFFERENT

template <typename SV, typename MV> 
ParticleFilter<SV,MV>::ParticleFilter(MCPdf<SV> * prior, 
				      ConditionalPdf<SV,SV> * proposal, 
				      int resampleperiod,
				      double resamplethreshold,
				      int resamplescheme)
  : Filter<SV,MV>(prior),
    _proposal(proposal),
    _resampleScheme(resamplescheme),
    _created_post(true)
{
  /* Initialize Post, at time = 0, post = prior
     To be more clean, this should be done in the filter base class,
     but this is impossible because of the pure virtuals...
  */
  this->_post = new MCPdf<SV>(prior->NumSamplesGet(),prior->DimensionGet());
  // Post is equal to prior at timetep 1 
  /* Note: Dirty cast should be avoided by not demanding an MCPdf as
     prior and just sample from the prior instead  :-( 
  */
  bool ret = (dynamic_cast<MCPdf<SV> *>(this->_post))->ListOfSamplesSet(prior->ListOfSamplesGet());
  assert(ret);

  // Initialise lists of samples
  _old_samples = (prior->ListOfSamplesGet());
  _new_samples = _old_samples;
  

  // You have to choose for dynamic resampling by specifying a threshold != 0 OR give me a fixed resample period != 0
  assert(!(resampleperiod == 0 && resamplethreshold == 0));
  assert(!(resampleperiod != 0 && resamplethreshold != 0));

  // dynamic resampling
  if (resampleperiod == 0)  
   _dynamicResampling = true;
  // fixed period resampling
  else
    _dynamicResampling = false;
  _resamplePeriod = resampleperiod;
  _resampleThreshold = resamplethreshold;
}



template <typename SV, typename MV> 
ParticleFilter<SV,MV>::ParticleFilter(MCPdf<SV> * prior, 
				      MCPdf<SV> * post,
				      ConditionalPdf<SV,SV> * proposal, 
				      int resampleperiod,
				      double resamplethreshold,
				      int resamplescheme)
  : Filter<SV,MV>(prior),
    _proposal(proposal),
    _resampleScheme(resamplescheme),
    _created_post(false)
{
  this->_post = post;
  // Post is equal to prior at timestep 1 
  /* Note: Dirty cast should be avoided by not demanding an MCPdf as
     prior and just sample from the prior instead  :-( 
  */
  bool ret = (dynamic_cast<MCPdf<SV> *>(this->_post))->ListOfSamplesSet(prior->ListOfSamplesGet());
  assert(ret);

  // Initialise lists of samples
  _old_samples = (prior->ListOfSamplesGet());
  _new_samples = _old_samples;
  
  // You have to choose for dynamic resampling by specifying a threshold != 0 OR give me a fixed resample period != 0
  assert(!(resampleperiod == 0 && resamplethreshold == 0));
  assert(!(resampleperiod != 0 && resamplethreshold != 0));

  // dynamic resampling
  if (resampleperiod == 0)  
   _dynamicResampling = true;
  // fixed period resampling
  else
    _dynamicResampling = false;

  _resamplePeriod = resampleperiod;
  _resampleThreshold = resamplethreshold;
}




template <typename SV, typename MV> 
ParticleFilter<SV,MV>::~ParticleFilter()
{
  if (_created_post)
    delete this->_post;
}

template <typename SV, typename MV> 
ParticleFilter<SV,MV>::ParticleFilter(const ParticleFilter<SV,MV> & filter)
  : Filter<SV,MV>(filter),
    _created_post(true)
{
  // Copy constructor of MCPdf
  // Probably a bug...
  this->_post = new MCPdf<SV>(dynamic_cast<MCPdf<SV> *>(filter._post));
}
  
template <typename SV, typename MV> void
ParticleFilter<SV,MV>::ProposalSet(ConditionalPdf<SV,SV> * const cpdf)
{
  _proposal = cpdf;
}

template <typename SV, typename MV> ConditionalPdf<SV,SV> *
ParticleFilter<SV,MV>::ProposalGet()
{
  return _proposal;
}

template <typename SV, typename MV> bool
ParticleFilter<SV,MV>::ProposalStepInternal(SystemModel<SV> * const sysmodel,
					    const SV & u,
					    MeasurementModel<MV,SV> * const measmodel,
					    const MV & z,
					    const SV & s)
{
  // Get all samples from the current post through proposal density
  _old_samples = (dynamic_cast<MCPdf<SV> *>(this->_post))->ListOfSamplesGet();
  Sample<SV> sample;

  _ns_it = _new_samples.begin();
  for ( _os_it=_old_samples.begin(); _os_it != _old_samples.end() ; _os_it++)
    {
      const SV& x_old = _os_it->ValueGet();
      _proposal->ConditionalArgumentSet(0,x_old);

      if (!sysmodel->SystemWithoutInputs())
	{
	  _proposal->ConditionalArgumentSet(1,u);
	  if (this->_proposal_depends_on_meas)
	    {
          #ifndef STATE_AND_MEAS_VAR_DIFFERENT
	      _proposal->ConditionalArgumentSet(2,z);
	      if (!measmodel->SystemWithoutSensorParams())
		_proposal->ConditionalArgumentSet(3,s);
              #endif
	    }
 
	}
      else // System without inputs
	{
	  if (this->_proposal_depends_on_meas)
	    {
          #ifndef STATE_AND_MEAS_VAR_DIFFERENT
	      _proposal->ConditionalArgumentSet(1,z);
	      if (!measmodel->SystemWithoutSensorParams())
       		_proposal->ConditionalArgumentSet(2,s);
              #endif

	    }
	}
      // Bug, make sampling method a parameter!
      _proposal->SampleFrom(sample, DEFAULT,NULL);
      _ns_it->ValueSet(sample.ValueGet());
      _ns_it->WeightSet(_os_it->WeightGet());
      _ns_it++;
    }
  
  (this->_timestep)++;

  // Update the list of samples
  return (dynamic_cast<MCPdf<SV> *>(this->_post))->ListOfSamplesUpdate(_new_samples);

}


template <typename SV, typename MV> bool
ParticleFilter<SV,MV>::UpdateWeightsInternal(SystemModel<SV> * const sysmodel,
					     const SV & u,
					     MeasurementModel<MV,SV> * const measmodel,
					     const MV & z, 
					     const SV & s)
{
  Probability weightfactor = 1;
  
  // Update the weights
  // Same remarks as for the system update!
  _new_samples = (dynamic_cast<MCPdf<SV> *>(this->_post))->ListOfSamplesGet();
  _os_it = _old_samples.begin();

  for ( _ns_it=_new_samples.begin(); _ns_it != _new_samples.end() ; _ns_it++)
    {
      const SV& x_new = _ns_it->ValueGet();
      const SV& x_old = _os_it->ValueGet();
      
      if (sysmodel == NULL)
	{
	  if (measmodel->SystemWithoutSensorParams() == true)
	    weightfactor = measmodel->ProbabilityGet(z,x_new);
	  else
	    weightfactor = measmodel->ProbabilityGet(z,x_new,s);
	}
      else // We do have a system model
	{
      _proposal->ConditionalArgumentSet(0,x_old);
	  if (measmodel->SystemWithoutSensorParams() == true)
	    {
	      weightfactor = measmodel->ProbabilityGet(z,x_new);
	      if (sysmodel->SystemWithoutInputs() == false)
		{
		  _proposal->ConditionalArgumentSet(1,u);
		  if (this->_proposal_depends_on_meas){
                    #ifndef STATE_AND_MEAS_VAR_DIFFERENT
                    _proposal->ConditionalArgumentSet(2,z);
                    #endif
		  }

		  if (_proposal->ProbabilityGet(x_new) != 0)
		    weightfactor = weightfactor * ( sysmodel->ProbabilityGet(x_new,x_old,u) / _proposal->ProbabilityGet(x_new) );
		  else weightfactor = 0;
		}
	      else // we do have a system without inputs
		{
		  if (this->_proposal_depends_on_meas){
                    #ifndef STATE_AND_MEAS_VAR_DIFFERENT
		    _proposal->ConditionalArgumentSet(1,z);
                    #endif
		  }
		  if ( _proposal->ProbabilityGet(x_new) != 0)
		    weightfactor = weightfactor * ( sysmodel->ProbabilityGet(x_new,_os_it->ValueGet()) / _proposal->ProbabilityGet(x_new) );
		  else weightfactor = 0;
		}
	    }
	  else // System with sensor Parameters
	    {
	      weightfactor = measmodel->ProbabilityGet(z,x_new,s);
	    }
	}
      _ns_it->WeightSet(_ns_it->WeightGet() * weightfactor);

      _os_it++;
    }
  // Update the sample list of post the SumofWeights of the pdf
  return (dynamic_cast<MCPdf<SV> *>(this->_post))->ListOfSamplesUpdate(_new_samples);

}

template <typename SV, typename MV> bool
ParticleFilter<SV,MV>::DynamicResampleStep()
{
  // Resampling?
  bool resampling = false;
  double sum_sq_weigths = 0.0;

  // Resampling if necessary
  if ( this->_dynamicResampling)
    { 
      // Check if sum of 1 / \sum{(wi_normalised)^2} < threshold
      // This is the criterion proposed by Liu
      // BUG  foresee other methods of approximation/calculating
      // effective sample size.  Let the user implement this in !
      _new_samples = (dynamic_cast<MCPdf<SV> *>(this->_post))->ListOfSamplesGet();
      for ( _ns_it=_new_samples.begin(); _ns_it != _new_samples.end() ; _ns_it++)
	{
	  sum_sq_weigths += pow(_ns_it->WeightGet(),2);
	}
      if ((1.0 / sum_sq_weigths) < _resampleThreshold)
	{
	  // #define __RESAMPLE_DEBUG__
#ifdef __RESAMPLE_DEBUG__
	  cout << "resampling now: " << this->_timestep 
	       << "\tN_eff: " << (1.0 / sum_sq_weigths) << endl;
#endif // __RESAMPLE_DEBUG__
	  resampling = true;
	}
    }
  if (resampling == true)
    return this->Resample();
  else
    return true;
}


template <typename SV, typename MV> bool
ParticleFilter<SV,MV>::StaticResampleStep()
{
  // Resampling if necessary
  if ( (!this->_dynamicResampling) &&  (((this->_timestep) % _resamplePeriod) == 0) && (this->_timestep != 0))
    return this->Resample();
  return true;
}


template <typename SV, typename MV> bool
ParticleFilter<SV,MV>::UpdateInternal(SystemModel<StateVar>* const sysmodel,
				      const StateVar& u,
				      MeasurementModel<MeasVar,StateVar>* const measmodel,
				      const MeasVar& z,
				      const StateVar& s)
{
  bool result = true;

  // Only makes sense if there is a system model?
  // Bug, not completely true, but should do for now...
  if (sysmodel != NULL)
    {
      result = result && this->StaticResampleStep();
      result = result && this->ProposalStepInternal(sysmodel,u,measmodel,z,s);

    }
  // Updating the weights only makes sense using a measurement model
  if (measmodel != NULL)
    {
      result = result && this->UpdateWeightsInternal(sysmodel,u,measmodel,z,s);
      result = result && this->DynamicResampleStep();
    }

  return result;
}

template <typename SV, typename MV> bool
ParticleFilter<SV,MV>::Resample()
{
  int NumSamples = (dynamic_cast<MCPdf<SV> *>(this->_post))->NumSamplesGet();
  vector<Sample<SV> > new_samples(NumSamples);
  // #define __PARTICLEFILTER_DEBUG__
#ifdef __PARTICLEFILTER_DEBUG__
  cout << "PARTICLEFILTER: resampling now" << endl;
  _new_samples= (dynamic_cast<MCPdf<SV> *>(this->_post))->ListOfSamplesGet();
  for ( _ns_it=_new_samples.begin(); _ns_it != _new_samples.end() ; _ns_it++)
    {
      cout << "PF: Old samples:\n";
      cout << _ns_it->ValueGet() << _ns_it->WeightGet() << endl;
    }
#endif // __PARTICLEFILTER_DEBUG
  switch(_resampleScheme)
    {
    case MULTINOMIAL_RS:
      {
	(dynamic_cast<MCPdf<SV> *>(this->_post))->SampleFrom(new_samples, NumSamples,RIPLEY,NULL);
	break;
      }
    case SYSTEMATIC_RS:{break;}
    case STRATIFIED_RS:{break;}
    case RESIDUAL_RS:{break;}
    default:
      {
	cerr << "Sampling method not supported" << endl;
	break;
      }
    }
  bool result = (dynamic_cast<MCPdf<SV> *>(this->_post))->ListOfSamplesUpdate(new_samples);
#ifdef __PARTICLEFILTER_DEBUG__
  cout << "PARTICLEFILTER: after resampling" << endl;
  _new_samples= (dynamic_cast<MCPdf<SV> *>(this->_post))->ListOfSamplesGet();
  for ( _ns_it=_new_samples.begin(); _ns_it != _new_samples.end() ; _ns_it++)
    {
      cout << "PF: New samples:\n";
      cout << _ns_it->ValueGet() << _ns_it->WeightGet() << endl;
    }
#endif // __PARTICLEFILTER_DEBUG

  return result;
}


template<typename SV, typename MV> MCPdf<SV> *
ParticleFilter<SV,MV>::PostGet()
{
  return (MCPdf<SV>*)Filter<SV,MV>::PostGet();
}
