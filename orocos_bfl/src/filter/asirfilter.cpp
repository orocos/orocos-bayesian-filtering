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

#include "asirfilter.h"
#include "../sample/weightedsample.h"

#define StateVar SVar
#define MeasVar MVar

template <typename SVar, typename MVar>
ASIRFilter<SVar,MVar>::ASIRFilter(MCPdf<SVar> * prior,
				  int resampleperiod,
				  double resamplethreshold,
				  int resamplescheme)
  : ParticleFilter<SVar,MVar>(prior,NULL,resampleperiod,
			      resamplethreshold, resamplescheme)
{
  /* Well euhm, technically speaking it does, but we solve this
     software-matically in another way (by reusing the
     ProposalStepInternal from particlefilter.cpp).  This is because
     otherwise we would have to create a hybrid proposal density,
     which requires far more programming.
  */
  this->_proposal_depends_on_meas = false;
}


template <typename SVar, typename MVar>
ASIRFilter<SVar,MVar>::~ASIRFilter(){}


template <typename SVar, typename MVar> void
ASIRFilter<SVar,MVar>::UpdateInternal(SystemModel<SVar>* const sysmodel,
				      const SVar& u,
				      MeasurementModel<MVar,SVar>* const measmodel,
				      const MVar& z,
				      const SVar& s)
{
  if (sysmodel != NULL)
    {
      this->ProposalSet(sysmodel->SystemPdfGet());
    }

  /* The following code differs from standard SIR filter */
  this->_old_samples = (dynamic_cast<MCPdf<SV> *>(this->_post))->ListOfSamplesGet();
  int NumSamples = (dynamic_cast<MCPdf<SV> *>(this->_post))->NumSamplesGet();

  static SV x_old;
  static vector<SV> mu(NumSamples);
  typename vector<SV>::iterator mu_it = mu.begin();
  static  Probability weightfactor;
  static Sample<SV> sample;

  // Only if there's a system model (for now, probably not general
  // enough)
  if (sysmodel != NULL)
    {
      // Step 1:  Calculate beta values and adapt weights
      for ( this->_os_it=this->_old_samples.begin();
	    this->_os_it != this->_old_samples.end() ;
	    this->_os_it++)
	{
	  // Since the proposal is equal to the SystemPdf of the model
	  // here, simulating the systemmodel will do fine
	  x_old = this->_os_it->ValueGet();
	  this->_proposal->ConditionalArgumentSet(0,x_old);

	  if (!sysmodel->SystemWithoutInputs())
	    {
	      this->_proposal->ConditionalArgumentSet(1,u);
	      if (this->_proposal_depends_on_meas)
		{
		  this->_proposal->ConditionalArgumentSet(2,z);
		  if (!measmodel->SystemWithoutSensorParams())
		    this->_proposal->ConditionalArgumentSet(3,s);
		}
	    }
	  else // System without inputs
	    {
	      if (this->_proposal_depends_on_meas)
		{
		  this->_proposal->ConditionalArgumentSet(1,z);
		  if (!measmodel->SystemWithoutSensorParams())
		    this->_proposal->ConditionalArgumentSet(2,s);
		}
	    }
	  // Get Expected Value of of this particles
	  *mu_it = this->_proposal->ExpectedValueGet();
	  // Calculate likelihood of this particle
	  if (!measmodel->SystemWithoutSensorParams())
	    weightfactor = measmodel->ProbabilityGet(z,*mu_it,s);
	  else
	    weightfactor = measmodel->ProbabilityGet(z,*mu_it);

	  // Set new weight
	  this->_os_it->WeightSet(this->_os_it->WeightGet() * weightfactor);
	  mu_it++;
	}
      // Update list of samples
      (dynamic_cast<MCPdf<SV> *>(this->_post))->ListOfSamplesUpdate(this->_old_samples);

      // Step 2: Sample discrete index k (in O (N) ops)
      this->Resample();

      // Step 3: Proposal step
      this->ProposalStepInternal(sysmodel,u,measmodel,z,s);
    }

  // Step 4: Update the weights
  if (measmodel != NULL)
    {
      this->_new_samples = (dynamic_cast<MCPdf<SV> *>(this->_post))->ListOfSamplesGet();

      static SV x_new;
      mu_it=mu.begin();
      for ( this->_ns_it=this->_new_samples.begin();
	    this->_ns_it != this->_new_samples.end() ;
	    this->_ns_it++)
	{
	  x_new = this->_ns_it->ValueGet();
	  if (measmodel->SystemWithoutSensorParams() == true)
	    {
	      if (measmodel->ProbabilityGet(z,*mu_it) != 0)
		weightfactor = measmodel->ProbabilityGet(z,x_new) / measmodel->ProbabilityGet(z,*mu_it);
	      else weightfactor = 0;
	    }
	  else // System with sensor Parameters
	    if (measmodel->ProbabilityGet(z,*mu_it,s) != 0)
	      weightfactor = measmodel->ProbabilityGet(z,x_new,s) / measmodel->ProbabilityGet(z,*mu_it,s);
	    else weightfactor = 0;
	  this->_ns_it->WeightSet(this->_ns_it->WeightGet() * weightfactor);
	  mu_it++;
	}
      // Update the sample list of post the SumofWeights of the pdf
      (dynamic_cast<MCPdf<SV> *>(this->_post))->ListOfSamplesUpdate(this->_new_samples);

      // Resample
      this->ResampleStep();
    }
}
