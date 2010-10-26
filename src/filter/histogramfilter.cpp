// $Id: histogramfilter.cpp 14935 2007-12-17  $
// Copyright (C) 2007 Tinne De Laet  <tinne dot delaet at mech dot kuleuven dot be>
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
#include "histogramfilter.h"
#include <cmath>
#include <iostream>


  template <typename MeasVar>
  HistogramFilter<MeasVar>::HistogramFilter(DiscretePdf * prior)
    : Filter<int,MeasVar>(prior)
  {
    // create posterior density
    this->_post = new DiscretePdf(*prior);
    // Initialise lists of probabilities
    _old_prob = (prior->ProbabilitiesGet());
    _new_prob = _old_prob;
  }

  template <typename MeasVar>
  HistogramFilter<MeasVar>::~HistogramFilter()
  {
    delete this->_post;
  }

  template <typename MeasVar>
  void
  HistogramFilter<MeasVar>::SysUpdate(SystemModel<int>* const sysmodel, const int& u)
  {
    bool without_inputs = sysmodel->SystemWithoutInputs();
    int num_states = ( (DiscretePdf*)this->_post )->NumStatesGet();
    _old_prob = ( (DiscretePdf*)this->_post )->ProbabilitiesGet();

    for (int to_state = 0; to_state< num_states ; to_state++)
    {
        double temp = 0;
        if (without_inputs)
        {
            for (int from_state = 0; from_state< num_states ; from_state++)
            {
                temp = temp + (double)sysmodel->ProbabilityGet(to_state,from_state) * _old_prob[from_state];
            }
        }
        else
        {
            for (int from_state = 0; from_state< num_states ; from_state++)
            {
                temp = temp + (double)sysmodel->ProbabilityGet(to_state,from_state,u) * _old_prob[from_state];
            }
        }
       _new_prob[to_state] = temp;
    }
    ( (DiscretePdf*)this->_post )->ProbabilitiesSet(_new_prob);
  }

  template <typename MeasVar>
  void
  HistogramFilter<MeasVar>::MeasUpdate(MeasurementModel<MeasVar,int>* const measmodel, const MeasVar& z, const int& s)
  {
      int num_states = ( (DiscretePdf*)this->_post )->NumStatesGet();
      _old_prob = ( (DiscretePdf*)this->_post )->ProbabilitiesGet();
      for (int state = 0; state< num_states  ; state++)
      {
          if (measmodel->SystemWithoutSensorParams() == true)  _old_prob[state] = _old_prob[state] * measmodel->ProbabilityGet(z,state);
          else _old_prob[state] = _old_prob[state] * measmodel->ProbabilityGet(z,state,s);
          _old_prob[state] = (Probability)((double)( _old_prob[state]) );
      }
      ( (DiscretePdf*)this->_post )->ProbabilitiesSet(_old_prob);
  }

  template <typename MeasVar>
  bool
  HistogramFilter<MeasVar>::UpdateInternal(SystemModel<int>* const sysmodel,
			       const int& u,
			       MeasurementModel<MeasVar,int>* const measmodel,
			       const MeasVar& z, const int& s)
  {
    if (sysmodel != NULL) 	SysUpdate(sysmodel,u);
    if (measmodel != NULL)	MeasUpdate(measmodel,z,s);
    return true;
  }

  template <typename MeasVar>
  DiscretePdf*
  HistogramFilter<MeasVar>::PostGet()
  {
    return (DiscretePdf*)Filter<int,MeasVar>::PostGet();
  }
