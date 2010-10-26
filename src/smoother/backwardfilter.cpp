// $Id: backwardfilter.cpp 6736 2006-12-22 11:24:42Z tdelaet $
// Copyright (C) 2006  Tinne De Laet <first dot last at mech dot kuleuven dot be>
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

#include "backwardfilter.h"

#define StateVar SVar

// Constructor
template<typename SVar>
BackwardFilter<SVar>::BackwardFilter(Pdf<SVar>* prior)
  : _prior(prior),
    _timestep(0)
{}

template<typename SVar>
BackwardFilter<SVar>::~BackwardFilter(){}

template<typename SVar>
BackwardFilter<SVar>::BackwardFilter(const BackwardFilter<SVar>& backwardfilter)
{}

template<typename SVar>  void
BackwardFilter<SVar>::Reset(Pdf<SVar> * prior)
{
  _prior = prior;
  _post = prior;
}

template<typename SVar> int
BackwardFilter<SVar>::TimeStepGet() const
{
  return _timestep;
}

template<typename SVar> bool
BackwardFilter<SVar>::Update(SystemModel<SVar>* const sysmodel, const SVar& u, Pdf<SVar>* const filtered_post)
{
  return this->UpdateInternal(sysmodel,u,filtered_post);
}

template<typename SVar> bool
BackwardFilter<SVar>::Update(SystemModel<SVar>* const sysmodel, Pdf<SVar>* const filtered_post)
{
  SVar u;
  return this->UpdateInternal(sysmodel,u,filtered_post);
}

template<typename SVar> Pdf<SVar> *
BackwardFilter<SVar>::PostGet()
{
  return _post;
}
