// Copyright (C) 2001-2006 Klaas Gadeyne <first dot last at gmail dot com>
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
// Wrapper around several Pseudo RNG-libraries
#ifndef __ORO_PSEUDORNG__
#define __ORO_PSEUDORNG__

#include "../../sample/sample.h"

namespace BFL
{
  // Sample from univariate normal distribution with mu and sigma
  // Maybe this should become of type sample in the future!
  double rnorm (const double & mu, const double & sigma);
  double runif ();
  double runif (const double & min, const double & max);
}

#endif // __ORO_PSEUDORNG
