// Copyright (C) 2009 Tinne De Laet <first dot last at gmail dot com>
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
// $Id: mixtureParticleFilter.cpp 2009-02-02 tdelaet $

#include "mixtureParticleFilter.h"
#include "../pdf/mixture.h"

#define SV StateVar
#define MV MeasVar
#define MeasModel MeasurementModel

#define STATE_AND_MEAS_VAR_DIFFERENT

template <typename SV, typename MV>
MixtureParticleFilter<SV,MV>::MixtureParticleFilter(Mixture<SV> * prior,
				      ConditionalPdf<SV,SV> * proposal,
				      int resampleperiod,
				      double resamplethreshold,
				      int resamplescheme,
                      int maintainMixturePeriod)
  : Filter<SV,MV>(prior)
  , _proposal(proposal)
  , _sample(WeightedSample<SV>(prior->DimensionGet()))
  , _resampleScheme(resamplescheme)
  , _created_post(true)
  , _newMixtureWeights(prior->NumComponentsGet())
  , _sumWeights(prior->NumComponentsGet())
  , _old_samplesVec(prior->NumComponentsGet())
  , _new_samplesVec(prior->NumComponentsGet())
  , _new_samples_unweightedVec(prior->NumComponentsGet())
  , _maintainMixturePeriod(maintainMixturePeriod)
{
  /* Initialize Post, at time = 0, post = prior
     To be more clean, this should be done in the filter base class,
     but this is impossible because of the pure virtuals...
  */
  // Post is equal to prior at timetep 1
  this->_post = prior->Clone();

  // Initialise vector of lists of samples
  for(int i = 0 ; i< dynamic_cast<Mixture<SV> *>(this->_post)->NumComponentsGet(); i++)
  {
    _old_samplesVec[i] = (dynamic_cast<MCPdf<SV> *>(prior->ComponentGet(i))->ListOfSamplesGet());
  }
  _new_samplesVec = _old_samplesVec;



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
MixtureParticleFilter<SV,MV>::MixtureParticleFilter(Mixture<SV> * prior,
				      Mixture<SV> * post,
				      ConditionalPdf<SV,SV> * proposal,
				      int resampleperiod,
				      double resamplethreshold,
				      int resamplescheme,
                      int maintainMixturePeriod)
  : Filter<SV,MV>(prior)
  ,  _proposal(proposal)
  ,  _resampleScheme(resamplescheme)
  ,  _created_post(false)
  , _newMixtureWeights(prior->NumComponentsGet())
  , _sumWeights(prior->NumComponentsGet())
  , _old_samplesVec(prior->NumComponentsGet())
  , _new_samplesVec(prior->NumComponentsGet())
  , _new_samples_unweightedVec(prior->NumComponentsGet())
  , _maintainMixturePeriod(maintainMixturePeriod)
{
  this->_post = post;
  // Post is equal to prior at timestep 1
  /* Note: Dirty cast should be avoided by not demanding an MCPdf as
     prior and just sample from the prior instead  :-(
  */
  bool ret = (dynamic_cast<MCPdf<SV> *>(this->_post))->ListOfSamplesSet(prior->ListOfSamplesGet());
  assert(ret);
  for(int i =0 ; i < prior.NumComponentsGet() ; i++)
  {
    bool ret = (dynamic_cast<MCPdf<SV> *>(this->_post->ComponentGet(i)))->ListOfSamplesSet(prior->ComponentGet(i)->ListOfSamplesGet());
  }

  // Initialise vector of lists of samples
  for(int i =0 ; i < prior.NumComponentsGet() ; i++)
  {
    _old_samplesVec[i] = (prior.ComponentGet(i)->ListOfSamplesGet());
  }
  _new_samplesVec = _old_samplesVec;

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
MixtureParticleFilter<SV,MV>::~MixtureParticleFilter()
{
  if (_created_post)
    delete this->_post;
}

template <typename SV, typename MV>
MixtureParticleFilter<SV,MV>::MixtureParticleFilter(const MixtureParticleFilter<SV,MV> & filter)
  : Filter<SV,MV>(filter)
  , _created_post(true)
{
  // Clone the Mixture posterior of filter 
  this->_post = filter.PostGet().Clone();
  this->_newMixtureWeights.resize(this->_post->NumComponentsGet());
  this->_sumWeights.resize(this->_post->NumComponentsGet());
}

template <typename SV, typename MV> void
MixtureParticleFilter<SV,MV>::Reset(Mixture<SV> * prior)
{
    this->_prior = prior;
    delete this->_post;
    this->_post = prior->Clone();
    this->_newMixtureWeights.resize(dynamic_cast<Mixture<SV> *>(this->_post)->NumComponentsGet());
    this->_sumWeights.resize(dynamic_cast<Mixture<SV> *>(this->_post)->NumComponentsGet());
}

template <typename SV, typename MV> void
MixtureParticleFilter<SV,MV>::ProposalSet(ConditionalPdf<SV,SV> * const cpdf)
{
  _proposal = cpdf;
}

template <typename SV, typename MV> ConditionalPdf<SV,SV> *
MixtureParticleFilter<SV,MV>::ProposalGet()
{
  return _proposal;
}

// Proposal step can be executed for each component in the mixture seperately
template <typename SV, typename MV> bool
MixtureParticleFilter<SV,MV>::ProposalStepInternal(SystemModel<SV> * const sysmodel,
					    const SV & u,
					    MeasurementModel<MV,SV> * const measmodel,
					    const MV & z,
					    const SV & s)
{
    bool result = true;
    for(int i = 0 ; i< dynamic_cast<Mixture<SV> *>(this->_post)->NumComponentsGet(); i++)
    {
        result = result && this->ProposalStepInternalOne(i,sysmodel,u,measmodel,z,s);
    }
    return result;
}

// Proposal step for one component
template <typename SV, typename MV> bool
MixtureParticleFilter<SV,MV>::ProposalStepInternalOne(int component, SystemModel<SV> * const sysmodel,
					    const SV & u,
					    MeasurementModel<MV,SV> * const measmodel,
					    const MV & z,
					    const SV & s)
{
  // Get all samples from the current post through proposal density
  _old_samplesVec[component]= (dynamic_cast<MCPdf<SV> *>(dynamic_cast<Mixture<SV> *>(this->_post)->ComponentGet(component)))->ListOfSamplesGet();

  _ns_it = _new_samplesVec[component].begin();
  for ( _os_it=_old_samplesVec[component].begin(); _os_it != _old_samplesVec[component].end() ; _os_it++)
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
      _proposal->SampleFrom(_sample, DEFAULT,NULL);
      _ns_it->ValueSet(_sample.ValueGet());
      _ns_it->WeightSet(_os_it->WeightGet());
      _ns_it++;
    }

  (this->_timestep)++;

  // Update the list of samples
  return (dynamic_cast<MCPdf<SV> *>(dynamic_cast<Mixture<SV> *>(this->_post)->ComponentGet(component)))->ListOfSamplesUpdate(_new_samplesVec[component]);

}

template <typename SV, typename MV> bool
MixtureParticleFilter<SV,MV>::UpdateWeightsInternal(SystemModel<SV> * const sysmodel,
					     const SV & u,
					     MeasurementModel<MV,SV> * const measmodel,
					     const MV & z,
					     const SV & s)
{
  // TODO: first calculate new weights of mixture components
    bool result = true;
    for(int i = 0 ; i< dynamic_cast<Mixture<SV> *>(this->_post)->NumComponentsGet(); i++)
    {
        this->UpdateWeightsInternalOne(i,sysmodel,u,measmodel,z,s);
    }
    for(int i = 0 ; i< dynamic_cast<Mixture<SV> *>(this->_post)->NumComponentsGet(); i++)
    {
        // calculate the new unnormalized mixture weights:
        // new_weight = old_weight * _sumWeights[i]
        // _sumWeights[i] is the sum of the particle weights of the i'th component
        // MPCdf and is an approximation of the i-th component likelihood
        _newMixtureWeights[i] = dynamic_cast<Mixture<SV> *>(this->_post)->WeightGet(i) * _sumWeights[i];
    }
    // Update the mixture weights
    result = result && dynamic_cast<Mixture<SV> * >(this->_post)->WeightsSet(_newMixtureWeights); // this function automatically takes care of normalization

    return result;
}

template <typename SV, typename MV> bool
MixtureParticleFilter<SV,MV>::UpdateWeightsInternalOne(int component,SystemModel<SV> * const sysmodel,
					     const SV & u,
					     MeasurementModel<MV,SV> * const measmodel,
					     const MV & z,
					     const SV & s)
{
  if(component < 0 || component > dynamic_cast<Mixture<SV> *>(this->_post)->NumComponentsGet())
  {
        cerr<< "MixtureParticleFilter::UpdateWeightsInternalOne called with invalid component number " << endl;
        return false;
  } 
  _sumWeights[component] = 0.0; 
  Probability weightfactor = 1;
  // Update the weights
  // Same remarks as for the system update!
  _new_samplesVec[component] = (dynamic_cast<MCPdf<SV> *>(dynamic_cast<Mixture<SV> *>(this->_post)->ComponentGet(component)))->ListOfSamplesGet();
  _os_it = _old_samplesVec[component].begin();

  for ( _ns_it=_new_samplesVec[component].begin(); _ns_it != _new_samplesVec[component].end() ; _ns_it++)
    {
      const SV& x_new = _ns_it->ValueGet();
      const SV& x_old = _os_it->ValueGet();

      if (sysmodel == NULL)
	{
	  if (measmodel->SystemWithoutSensorParams() == true)
      {
	    weightfactor = measmodel->ProbabilityGet(z,x_new);
     }
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
      double new_weight = _ns_it->WeightGet() * weightfactor;
      _ns_it->WeightSet(new_weight);
      // add the new weight to the _sumWeights of this component
      _sumWeights[component] = _sumWeights[component] + new_weight; 

      _os_it++;
    }
  // Update the sample list of post the SumofWeights of the pdf
  // Update the mixture PDfs
  return (dynamic_cast<MCPdf<SV> *>(dynamic_cast<Mixture<SV> *>(this->_post)->ComponentGet(component)))->ListOfSamplesUpdate(_new_samplesVec[component]);
}

template <typename SV, typename MV> bool
MixtureParticleFilter<SV,MV>::DynamicResampleStep()
{
    bool result = true;
    // Independent resampling for different components
    for(int i = 0 ; i< dynamic_cast<Mixture<SV> *>(this->_post)->NumComponentsGet(); i++)
    {
        result == result && this->DynamicResampleStepOne(i);
    }
    return result;
}

template <typename SV, typename MV> bool
MixtureParticleFilter<SV,MV>::DynamicResampleStepOne(int component)
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
      _new_samplesVec[component] = (dynamic_cast<MCPdf<SV> *>(dynamic_cast<Mixture<SV> *>(this->_post)->ComponentGet(component)))->ListOfSamplesGet();
      for ( _ns_it=_new_samplesVec[component].begin(); _ns_it != _new_samplesVec[component].end() ; _ns_it++)
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
    return this->ResampleOne(component);
  else
    return true;
}


template <typename SV, typename MV> bool
MixtureParticleFilter<SV,MV>::StaticResampleStep()
{
  // Resampling if necessary
  if ( (!this->_dynamicResampling) &&  (((this->_timestep) % _resamplePeriod) == 0) && (this->_timestep != 0))
    return this->Resample();
  return true;
}


template <typename SV, typename MV> bool
MixtureParticleFilter<SV,MV>::UpdateInternal(SystemModel<StateVar>* const sysmodel,
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

  // Mixture Computation: recompute mixture representation to take into account
  // possibly varying number of modes.
  result = result && this->MaintainMixture();

  return result;
}

template <typename SV, typename MV> bool
MixtureParticleFilter<SV,MV>::Resample()
{
    bool result = true;
    // Independent resampling for different components
    for(int i = 0 ; i< dynamic_cast<Mixture<SV> *>(this->_post)->NumComponentsGet(); i++)
    {
        result == result && this->ResampleOne(i);
    }
    return result;
}

template <typename SV, typename MV> bool
MixtureParticleFilter<SV,MV>::ResampleOne(int component)
{
  int NumSamples = (dynamic_cast<MCPdf<SV> *>(dynamic_cast<Mixture<SV> *>(this->_post)->ComponentGet(component)))->NumSamplesGet();
  // #define __MixtureParticleFilter_DEBUG__
#ifdef __MixtureParticleFilter_DEBUG__
  cout << "MixtureParticleFilter: resampling now" << endl;
  _new_samplesVec[component]= (dynamic_cast<MCPdf<SV> *>(dynamic_cast<Mixture<SV> *>(this->_post)->ComponentGet(component)))->ListOfSamplesGet();
  for ( _ns_it=_new_samplesVec[component].begin(); _ns_it != _new_samplesVec[component].end() ; _ns_it++)
    {
      cout << "PF: Old samples:\n";
      cout << _ns_it->ValueGet() << _ns_it->WeightGet() << endl;
    }
#endif // __MixtureParticleFilter_DEBUG
  switch(_resampleScheme)
    {
    case MULTINOMIAL_RS:
      {
	(dynamic_cast<MCPdf<SV> *>(dynamic_cast<Mixture<SV> *>(this->_post)->ComponentGet(component)))->SampleFrom(_new_samples_unweightedVec[component], NumSamples,RIPLEY,NULL);
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
  bool result = (dynamic_cast<MCPdf<SV> *>(dynamic_cast<Mixture<SV> *>(this->_post)->ComponentGet(component)))->ListOfSamplesUpdate(_new_samples_unweightedVec[component]);
#ifdef __MixtureParticleFilter_DEBUG__
  cout << "MixtureParticleFilter: after resampling" << endl;
  _new_samplesVec[component]= (dynamic_cast<MCPdf<SV> *>(dynamic_cast<Mixture<SV> *>(this->_post)->ComponentGet(component)))->ListOfSamplesGet();
  for ( _ns_it=_new_samplesVec[component].begin(); _ns_it != _new_samplesVec[component].end() ; _ns_it++)
    {
      cout << "PF: New samples:\n";
      cout << _ns_it->ValueGet() << _ns_it->WeightGet() << endl;
    }
#endif // __MixtureParticleFilter_DEBUG

  return result;
}


template<typename SV, typename MV> Mixture<SV> *
MixtureParticleFilter<SV,MV>::PostGet()
{
  return (Mixture<SV>*)Filter<SV,MV>::PostGet();
}

template <typename SV, typename MV> bool
MixtureParticleFilter<SV,MV>::MaintainMixtureStep()
{
  // Resampling if necessary
  if ( (((this->_timestep) % _maintainMixturePeriod) == 0) && (this->_timestep != 0))
    return this->MaintainMixture();
  return true;
}

template <typename SV, typename MV> bool
MixtureParticleFilter<SV,MV>::MaintainMixture()
{
  // Default method doesn't take care of Mixture Maintainance
  return true;
}
