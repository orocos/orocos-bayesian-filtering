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

// This file only contains template specialisation code
#include "sample.h"
#include <cassert>

namespace BFL
{
  using namespace MatrixWrapper;


  // Template Specialisation for T =ColumnVector
  template <> inline
  Sample<ColumnVector>::Sample (unsigned int dimension)
    : Value(dimension)
    {};


  template <> inline unsigned int
  Sample<ColumnVector>::DimensionGet() const
  {
    return Value.rows();
  };

  template <> inline void
  Sample<ColumnVector>::DimensionSet(unsigned int dim)
  {
    return Value.resize(dim);
  };

  // Template Specialisation for T = double
  template <> inline unsigned int
  Sample<double>::DimensionGet() const
  {
    return 1;
  };

  template <> inline void
  Sample<double>::DimensionSet(unsigned int dim)
  {
    assert(dim == 1);
  };

  // Template Specialisation for T = int
  template <> inline unsigned int
  Sample<int>::DimensionGet() const
  {
    return 1;
  };

  template <> inline void
  Sample<int>::DimensionSet(unsigned int dim)
  {
    assert(dim == 1);
  };
}





