// $Id: innoationCheck.cpp tdelaet $
// Copyright (C) 2007
//                    Tinne De Laet <tinne dot delaet at mech dot kuleuven dot ac dot be>
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
#include "innovationCheck.h"
#include "../wrappers/rng/rng.h" // Wrapper around several rng libraries

namespace BFL
{
  using namespace MatrixWrapper;


  InnovationCheck::InnovationCheck(double min_inno)
    : min_innovation(min_inno)
  {}

  InnovationCheck::~InnovationCheck()
  {}

  bool
  InnovationCheck::check(ColumnVector innovation)
  {
        // basic implementation which checks if the norm of the innovation is
        // higher than the minimum innovation specified
        return (innovation.transpose() * innovation  >= min_innovation);
  }
}
