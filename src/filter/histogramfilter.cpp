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

namespace BFL
{
  using namespace MatrixWrapper;


  HistogramFilter::HistogramFilter(DiscretePdf * prior)
    : Filter<int,ColumnVector>(prior)
  {
    // create posterior density
    _post = new DiscretePdf(*prior);
  }

  HistogramFilter::~HistogramFilter()
  {
    delete _post;
  }

  void
  HistogramFilter::SysUpdate(SystemModel<int>* const sysmodel, const int& u)
  {
    bool without_inputs = sysmodel->SystemWithoutInputs();
    int num_states = ( (DiscretePdf*)_post )->NumStatesGet();
    vector<Probability> old_prob = ( (DiscretePdf*)_post )->ProbabilitiesGet();
    vector<Probability> new_prob(num_states);
    
    for (int to_state = 0; to_state< num_states ; to_state++)
    { 
        double temp = 0;
        if (without_inputs)
        {
            for (int from_state = 0; from_state< num_states ; from_state++)
            { 
                temp = temp + (double)sysmodel->ProbabilityGet(to_state,from_state) * old_prob[from_state];
            }
        } 
        else
        {
            for (int from_state = 0; from_state< num_states ; from_state++)
            { 
                temp = temp + (double)sysmodel->ProbabilityGet(to_state,from_state,u) * old_prob[from_state];
            }
        }
       new_prob[to_state] = temp;
    }
    ( (DiscretePdf*)_post )->ProbabilitiesSet(new_prob);
  }

  void
  HistogramFilter::MeasUpdate(MeasurementModel<ColumnVector,int>* const measmodel, const ColumnVector& z, const int& s)
  {
      int num_states = ( (DiscretePdf*)_post )->NumStatesGet();
      vector<Probability> prob = ( (DiscretePdf*)_post )->ProbabilitiesGet();
      for (int state = 0; state< num_states  ; state++)
      { 
          if (measmodel->SystemWithoutSensorParams() == true)  prob[state] = prob[state] * measmodel->ProbabilityGet(z,state);
          else prob[state] = prob[state] * measmodel->ProbabilityGet(z,state,s);
          prob[state] = (Probability)((double)( prob[state]) );
      }
      ( (DiscretePdf*)_post )->ProbabilitiesSet(prob);
  }

  bool 
  HistogramFilter::UpdateInternal(SystemModel<int>* const sysmodel,
			       const int& u,
			       MeasurementModel<ColumnVector,int>* const measmodel,
			       const ColumnVector& z, const int& s)
  {
    if (sysmodel != NULL) 	SysUpdate(sysmodel,u); 
    if (measmodel != NULL)	MeasUpdate(measmodel,z,s);
    return true;
  }

  DiscretePdf*
  HistogramFilter::PostGet()
  {
    return (DiscretePdf*)Filter<int,ColumnVector>::PostGet();
  }
  
}
