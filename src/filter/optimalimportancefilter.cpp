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

#include "optimalimportancefilter.h"

#define SV StateVar
#define MV MeasVar
#define MeasModel MeasurementModel
#define OptImpFilter OptimalImportanceFilter

template <typename SV, typename MV> void
OptimalImportanceFilter::ConstructProposal(SystemModel<SV>* const sysmodel,
				MeasModel<MV,SV>* const measmodel)
{
  
  // DIRTY IMPLEMENTATION
  if (_proposal == NULL) // First time, no memory allocated yet
    {
      // Allocate necessary memory
      _proposal = new ConditionalGaussian(sysmodel->SystemPdfGet());
    }
  else // Check if we really need to allocate memory
    if (_proposal->DimensionGet())
      ;
}

template <typename SV, typename MV> 
OptimalImportanceFilter<SV,MV>::OptimalImportanceFilter(MCPdf<SV> * prior,
				  int resampleperiod,
				  double resamplethreshold,
				  int resamplescheme)
  : ParticleFilter<SV,MV>(prior, NULL, resampleperiod, 
			  resamplethreshold,
			  resamplescheme)
{};


template <typename SV, typename MV> 
OptimalImportanceFilter<SV,MV>::~OptimalImportanceFilter()
{
  delete _proposal;
};

template <typename SV, typename MV> 
OptimalImportanceFilter<SV,MV>::OptimalImportanceFilter(const OptimalImportanceFilter & filter) 
  : ParticleFilter<SV,MV>(filter){}

template <typename SV, typename MV> void 
OptimalImportanceFilter<SV,MV>::Update(SystemModel<SV>* const sysmodel,
			    const SV& u,
			    MeasModel<MV,SV>* const measmodel,
			    const MV& z,
			    const SV& s)
{
  // Change here
  this->ProposalSet();
  this->ProposalStepInput(sysmodel,u);
  this->MeasUpdate(measmodel,z,s);
};

template <typename SV, typename MV> void
OptimalImportanceFilter<SV,MV>::Update(SystemModel<SV>* const sysmodel,
			    MeasModel<MV,SV>* const measmodel,
			    const MV& z,
			    const SV& s)
{
  this->ProposalSet();
  this->ProposalStep(sysmodel);
  this->MeasUpdate(measmodel,z,s);
};

template <typename SV, typename MV> void
OptimalImportanceFilter<SV,MV>::Update(SystemModel<SV>* const sysmodel,
			    MeasModel<MV,SV>* const measmodel,
			    const MV& z)
{
  this->ProposalSet();
  this->ProposalStep(sysmodel);
  this->MeasUpdate(measmodel,z);
};

template <typename SV, typename MV> void
OptimalImportanceFilter<SV,MV>::Update(SystemModel<SV>* const sysmodel,
			    const SV& u,
			    MeasModel<MV,SV>* const measmodel,
			    const MV& z)
{
  this->ProposalSet();
  this->ProposalStepInput(sysmodel,u);
  this->MeasUpdate(measmodel,z);
};

/// Only sysupdate
template <typename SV, typename MV> void
OptimalImportanceFilter<SV,MV>::Update(SystemModel<SV>* const sysmodel,
			    const SV& u)
{
  this->ProposalSet();
  this->ProposalStepInput(sysmodel,u);
}

template <typename SV, typename MV> void
OptimalImportanceFilter<SV,MV>::Update(SystemModel<SV>* const sysmodel)
{
  this->ProposalSet();
  this->ProposalStep(sysmodel);
}

/// Only measupdate
template <typename SV, typename MV> void
OptimalImportanceFilter<SV,MV>::Update(MeasModel<MV,SV>* const measmodel,
			    const MV& z,
			    const SV& s)
{
  this->MeasUpdate(measmodel,z,s);
}

template <typename SV, typename MV> void
OptimalImportanceFilter<SV,MV>::Update(MeasModel<MV,SV>* const measmodel,
			    const MV& z)
{
  this->MeasUpdate(measmodel,z);
}

