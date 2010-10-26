// $ Id: $
// Copyright (C) 2002 Klaas Gadeyne <first dot last at gmail dot com>
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

#ifndef __EK_PARTICLE_FILTER__
#define __EK_PARTICLE_FILTER__

#include "particlefilter.h"
#include "../pdf/EKF_proposaldensity.h"
#include "../pdf/mcpdf.h"
#include <list>

namespace BFL
{
  /// Particle filter using EKF for proposal step
  /** NOTE: Only applicable to continuous problems with additive
      gaussian noise, so the models you specify ...
      specify should by
  */
  class EKParticleFilter
    : public ParticleFilter<ColumnVector,ColumnVector>
    {
    protected:
      /// Sample Covariances for use with EKF Proposal density
      std::vector<SymmetricMatrix> _sampleCov;
      std::vector<SymmetricMatrix>::iterator _cov_it;

      std::vector<SymmetricMatrix> _tmpCov;
      std::vector<SymmetricMatrix>::iterator _tmpCovit;

      // helper variables for resampleing to prevent memory allocation on the heap
      const int _dimension;
      const int _num_samples;
      std::vector<WeightedSample<ColumnVector> > _old_samples;
      std::vector<WeightedSample<ColumnVector> >::iterator _oit;
      std::vector<WeightedSample<ColumnVector> > _result_samples;
      std::vector<WeightedSample<ColumnVector> >::iterator _rit;
      std::vector<double> _unif_samples;
      std::vector<double> _CumPDF;
      std::vector<double>::const_iterator _CumPDFit;
      ColumnVector _x_old;
      Sample<ColumnVector> _sample;

      virtual bool UpdateInternal(SystemModel<ColumnVector>* const sysmodel,
				  const ColumnVector& u,
				  MeasurementModel<ColumnVector,ColumnVector>* const measmodel,
				  const ColumnVector& z,
				  const ColumnVector& s);

      /// Proposalstep is redefined here since we have to take into
      /// account the sample covariances here!
      virtual bool ProposalStepInternal(SystemModel<ColumnVector> * const sysmodel,
	                                const ColumnVector & u,
	                                MeasurementModel<ColumnVector,ColumnVector> * const measmodel,
					const ColumnVector & z,
					const ColumnVector & s);

      /// Resample also redefined for the same reasons...
      virtual bool Resample();

    public:
      /// Constructor
      /** @pre you created the necessary models and the prior
	  @param prior pointer to the Monte Carlo Pdf prior density
	  @param resampleperiod fixed resampling period (if desired)
	  @param resamplethreshold threshold used when dynamic resampling
	  @param resamplescheme resampling scheme, see header file for
	  different defines and their meaning
	  @bug prior should be of type pdf and not mcpdf.  See also
	  notes with implementation
      */
      EKParticleFilter(MCPdf<ColumnVector> * prior,
		       int resampleperiod = 0,
		       double resamplethreshold = 0,
		       int resamplescheme = DEFAULT_RS);

      /// Destructor
      virtual ~EKParticleFilter();
    };

} // End namespace BFL

#endif // __EK_PARTICLE_FILTER__
