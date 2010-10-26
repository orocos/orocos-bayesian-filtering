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

#include "particlesmoother.h"
#include "../pdf/mcpdf.h"

#define SV StateVar
#define MV MeasVar

template <typename SV>
ParticleSmoother<SV>::ParticleSmoother(MCPdf<SV> * prior)
  : BackwardFilter<SV>(prior)
{
  this->_post = new MCPdf<SV>(*prior);
  (dynamic_cast<MCPdf<SV> *>(this->_post))->ListOfSamplesSet(prior->ListOfSamplesGet());

  // Initialise lists of samples
  _old_samples = (prior->ListOfSamplesGet());
  _new_samples = _old_samples;
  _filtered_samples = _old_samples;
}


template <typename SV>
ParticleSmoother<SV>::~ParticleSmoother()
{
    delete this->_post;
}

template <typename SV> bool
ParticleSmoother<SV>::UpdateInternal(SystemModel<StateVar>* const sysmodel, const StateVar& u,Pdf<StateVar>* const  filtered_post)
{
  cout << "update started" << endl;
  SysUpdate(sysmodel,u,filtered_post);
  cout << "update fininshed" << endl;
  return true;
}

template <typename SV> void
ParticleSmoother<SV>::SysUpdate(SystemModel<StateVar>* const sysmodel, const StateVar& u,  Pdf<StateVar>* const filtered_post)
{
    // Get all samples from the filtered posterior
    _filtered_samples = (dynamic_cast<MCPdf<SV> *>(filtered_post))->ListOfSamplesGet();
    // Get all samples from the current post through proposal density
    _old_samples = (dynamic_cast<MCPdf<SV> *>(this->_post))->ListOfSamplesGet();
    int N  = (dynamic_cast<MCPdf<SV> *>(filtered_post))->NumSamplesGet();
    int M = (dynamic_cast<MCPdf<SV> *>(this->_post))->NumSamplesGet();

    vector<double> gamma(M);
    Matrix prob(M,N);
    int i = 1;
    for (_os_it = _old_samples.begin(); _os_it !=_old_samples.end() ; _os_it++)
      {
          const SV& x_smo = _os_it->ValueGet();
          double gamma_loc = 0.0; //correction factor for weight of new particle
          int j=1;
          for ( _fs_it=_filtered_samples.begin(); _fs_it != _filtered_samples.end() ; _fs_it++)
          {
              const SV& x_fil = _fs_it->ValueGet();
              // calculate prediction probabilities
              if (!sysmodel->SystemWithoutInputs()){
                prob(i,j) = sysmodel->ProbabilityGet(x_smo,x_fil,u);
              }
              else{
                prob(i,j) = sysmodel->ProbabilityGet(x_smo,x_fil);
              }
              gamma_loc = gamma_loc + _fs_it->WeightGet() * prob(i,j);
              // calculate correction factors
              gamma[j-1] = gamma_loc;
              j++;
             }
      i++;
      }
      //cout << "probabilities " << prob << endl;
      // make copy of filtered sample list => new smoothed samples
      _new_samples = _filtered_samples;
      // new weights for filtered samples
      i = 1 ;
      double tot_weight; // sum of weights needed to normalize
      vector<double> weight (_new_samples.size()); // vector of weights of particles
      for ( _ns_it=_new_samples.begin(); _ns_it != _new_samples.end() ; _ns_it++)
        {
          double delta = 0.0;
          int j=1;
          for (_os_it = _old_samples.begin(); _os_it !=_old_samples.end() ; _os_it++)
          {
              if (gamma[j-1]!=0){
                  delta = delta + _os_it->WeightGet() * prob(i,j) / gamma[j-1];
              }
              else{}
              j++;
          }
           // cout << "probability for i " << i << " is " << prob.columnCopy(i) << endl;
           // cout << "Weight: " << _ns_it->WeightGet() << " * " << delta << endl;
          weight[i-1] =_ns_it->WeightGet()*delta;
          tot_weight = tot_weight + weight[i-1];
          i++;
        }
    // Normalize the weights
      i = 1;
      for ( _ns_it=_new_samples.begin(); _ns_it != _new_samples.end() ; _ns_it++)
        {
          _ns_it->WeightSet(weight[i-1]/tot_weight);
          i++;
        }

    // Update the list of samples
    (dynamic_cast<MCPdf<SV> *>(this->_post))->ListOfSamplesSet(_new_samples);

}
