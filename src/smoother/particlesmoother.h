// $Id: particlesmoother.h 6736 2006-12-21 11:24:42Z tdelaet $
// Copyright (C) 2006 Tinne De Laet <first dot last at mech dot kuleuven dot be>
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


// implementation based on "A smoothing filter for condensation", by Isard and
// Blake
// (http://www.springerlink.com/content/?k=a+smoothing+filter+for+condensation)
// The bacward stage version
// @bug still needs extra testing!

#ifndef __PARTICLE_SMOOTHER__
#define __PARTICLE_SMOOTHER__

#include "backwardfilter.h"
#include "../pdf/conditionalpdf.h"
#include "../pdf/mcpdf.h"

namespace BFL
{

  /// Class representing a particle backward filter 
  template <typename StateVar> class ParticleSmoother 
    : public BackwardFilter<StateVar>
    {
    protected:
      virtual bool UpdateInternal(SystemModel<StateVar>* const sysmodel,
	    			  const StateVar& u, Pdf<StateVar>* const filtered_post);

      virtual void SysUpdate(SystemModel<StateVar>* const sysmodel, const StateVar& u , Pdf<StateVar>* const filtered_post);

      /// While updating store list of old samples
      vector<WeightedSample<StateVar> > _old_samples;
      /// While updating store list of new samples
      vector<WeightedSample<StateVar> > _new_samples;
      /// While updating store list of filtered samples
      vector<WeightedSample<StateVar> > _filtered_samples;
      /// Iterator for old list of samples
      typename vector<WeightedSample<StateVar> >::iterator _os_it;
      /// Iterator for new list of samples
      typename vector<WeightedSample<StateVar> >::iterator _ns_it;
      /// Iterator for list of filtered samples
      typename vector<WeightedSample<StateVar> >::iterator _fs_it;

    public:
      /// Constructor
      ParticleSmoother(MCPdf<StateVar> * prior); 

      /// Destructor
      virtual ~ParticleSmoother();
  
    };
#include "particlesmoother.cpp"

} // End namespace BFL

#endif // __PARTICLE_FILTER__
