// $Id$
// Copyright (C) 2002 Klaas Gadeyne <first dot last at gmail dot com>
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
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
//
#include "discretesystemmodel.h"

namespace BFL
{

  DiscreteSystemModel::DiscreteSystemModel(DiscreteConditionalPdf * systempdf) :
    SystemModel<int>((ConditionalPdf<int, int> *) systempdf){}

  DiscreteSystemModel::~DiscreteSystemModel(){}

  DiscreteSystemModel::DiscreteSystemModel(const DiscreteSystemModel & model) :
    SystemModel<int>(model){}

  unsigned int DiscreteSystemModel::NumStatesGet()const
  {
    return ( (DiscreteConditionalPdf*)_SystemPdf )->NumStatesGet();
  }

}
