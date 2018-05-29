// Copyright (C) 2002-2006 Klaas Gadeyne <first dot last at gmail dot com>
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
#include "rng.h"

#include <../config.h>
#ifdef __RNGWRAPPER_BOOST__  // BOOST RANDOM LIBRARY
// THE BOOST RANDOM NUMBER GENERATION LIBRARY

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_01.hpp>
#include <boost/random/uniform_real.hpp>

static boost::mt19937 Boost_Rng; // Source of randomness

static boost::uniform_real<double> Uniform_Distribution; // Uniform distribution
static boost::variate_generator<boost::mt19937&,boost::uniform_real<double> > roll(Boost_Rng,Uniform_Distribution);

double BFL::rnorm(const double& mu,const double& sigma)
{
  boost::normal_distribution<double> TestDist(mu,sigma);
  boost::variate_generator <boost::mt19937 &,boost::normal_distribution<double> > TestGen(Boost_Rng,TestDist);
  return TestGen();
}

double BFL::runif()
{
  return roll();
}

double BFL::runif(const double &min, const double& max)
{
  boost::uniform_real<double> Uniform_DistributionMinMax(min,max); // Uniform distribution
  boost::variate_generator<boost::mt19937&,boost::uniform_real<double> > roll(Boost_Rng,Uniform_DistributionMinMax);
  return roll();
}

#endif // __RNGWRAPPER_BOOST__






