// $Id: particlesmoother.h 6736 2006-12-21 11:24:42Z tdelaet $
// Copyright (C) 2006 Tinne De Laet <first dot last at mech dot kuleuven dot be>
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
