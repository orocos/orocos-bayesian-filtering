// $Id$
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

#include "filter.h"

#define StateVar SVar
#define MeasVar MVar

// Constructor
template<typename SVar, typename MVar>
Filter<SVar,MVar>::Filter(Pdf<SVar>* prior)
  : _prior(prior),
    _timestep(0)
{}

template<typename SVar, typename MVar>
Filter<SVar,MVar>::~Filter(){}

/// @todo Check if we should make a copy of the pdf's too?
/// @bug we should make a copy of the pdf's too
template<typename SVar, typename MVar>
Filter<SVar,MVar>::Filter(const Filter<SVar,MVar>& filt)
{}

template<typename SVar, typename MVar>  void
Filter<SVar,MVar>::Reset(Pdf<SVar> * prior)
{
  _prior = prior;
  _post = prior;
  // cout << "Filter::Reset() Post = " << _post->ExpectedValueGet() << endl;
}

template<typename SVar, typename MVar> int
Filter<SVar,MVar>::TimeStepGet() const
{
  return _timestep;
}

template<typename SVar, typename MVar> bool
Filter<SVar,MVar>::Update(SystemModel<SVar>* const sysmodel,
			  const SVar& u,
			  MeasurementModel<MVar,SVar>* const measmodel,
			  const MVar& z,
			  const SVar& s)
{
  return this->UpdateInternal(sysmodel,u,measmodel,z,s);
}

template<typename SVar, typename MVar> bool
Filter<SVar,MVar>::Update(SystemModel<SVar>* const sysmodel,
			  MeasurementModel<MVar,SVar>* const measmodel,
			  const MVar& z,
			  const SVar& s)
{
  SVar u;
  return this->UpdateInternal(sysmodel,u,measmodel,z,s);
}

template<typename SVar, typename MVar> bool
Filter<SVar,MVar>::Update(SystemModel<SVar>* const sysmodel,
			  MeasurementModel<MVar,SVar>* const measmodel,
			  const MVar& z)
{
  SVar s; SVar u;
  return this->UpdateInternal(sysmodel,u,measmodel,z,s);
}

template<typename SVar, typename MVar> bool
Filter<SVar,MVar>::Update(SystemModel<SVar>* const sysmodel,
			  const SVar& u,
			  MeasurementModel<MVar,SVar>* const measmodel,
			  const MVar& z)
{
  SVar s;
  return this->UpdateInternal(sysmodel,u,measmodel,z,s);
}

template<typename SVar, typename MVar> bool
Filter<SVar,MVar>::Update(SystemModel<SVar>* const sysmodel,
			  const SVar& u)
{
  SVar s; MVar z;
  return this->UpdateInternal(sysmodel,u,NULL,z,s);
}

template<typename SVar, typename MVar> bool
Filter<SVar,MVar>::Update(SystemModel<SVar>* const sysmodel)
{
  SVar s; MVar z; SVar u;
  return this->UpdateInternal(sysmodel,u,NULL,z,s);
}

template<typename SVar, typename MVar> bool
Filter<SVar,MVar>::Update(MeasurementModel<MVar,SVar>* const measmodel,
			  const MVar& z,
			  const SVar& s)
{
  SVar u;
  return this->UpdateInternal(NULL,u,measmodel,z,s);
}

template<typename SVar, typename MVar> bool
Filter<SVar,MVar>:: Update(MeasurementModel<MVar,SVar>* const measmodel,
			   const MVar& z)
{
  SVar u; SVar s;
  return this->UpdateInternal(NULL,u,measmodel,z,s);
}

template<typename SVar, typename MVar> Pdf<SVar> *
Filter<SVar,MVar>::PostGet()
{
  return _post;
}
