// $Id: uniform.cpp tdelaet$
// Copyright (C) 2007 Tinne De Laet <first dot last at gmail dot com>
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
#include "uniform.h"

#include "../wrappers/rng/rng.h" // Wrapper around several rng libraries

#include <cmath> 
#include <cassert>

namespace BFL
{
  using namespace MatrixWrapper;

  Uniform::Uniform (const ColumnVector& center, const ColumnVector& width) 
    : Pdf<ColumnVector> ( center.rows() )
    , _samples(DimensionGet())
  {
    // check if inputs are consistent
    assert (center.rows() == width.rows());
    _Lower = center - width/2.0;
    _Higher = center + width/2.0;
    _Height = 1;
    for (int i=1 ; i < width.rows()+1 ; i++ )
    {
        _Height = _Height / width(i);
    }
  }

  Uniform::Uniform (int dimension)
    : Pdf<ColumnVector>(dimension)
    , _samples(dimension)
  {
    _Lower.resize(dimension);
    _Higher.resize(dimension);
  }

  Uniform::~Uniform(){}

  std::ostream& operator<< (std::ostream& os, const Uniform& u)
  {
    os << "\nCenter: \n"    << u.CenterGet()
       << "\nWidth: \n" << u.WidthGet() << endl;
    return os;
  }

  Probability Uniform::ProbabilityGet(const ColumnVector& input) const
  {
    // test if input is located in area of Uniform distribution
    for (int i=1; i<input.rows()+1; i++) 
    {
        if ( ( input(i)>_Higher(i) ) || ( input(i) < _Lower(i) ) ) return 0;
    }    
    return _Height;
  }

  bool
  Uniform::SampleFrom (vector<Sample<ColumnVector> >& list_samples, const int num_samples, int method, void * args) const
  {
    // Perform memory allocation
    list_samples.resize(num_samples); // will break real-timeness if list_samples.size()!=num_samples
    vector<Sample<ColumnVector> >::iterator rit = list_samples.begin();
    switch(method)
      {
         case DEFAULT: 
         {
        	  while (rit != list_samples.end())
        	   {
                 for (unsigned int j=1; j < DimensionGet()+1; j++) _samples(j) = runif(_Lower(j) , _Higher(j) ); 
        	     rit->ValueSet(_samples);
        	     rit++;
        	   }
        	  return true;
          }
          default:
             return false;
       }
  }

  bool
  Uniform::SampleFrom (Sample<ColumnVector>& one_sample, int method, void * args) const
  {
    switch(method)
    {
      case DEFAULT: 
      {
         for (unsigned int j=1; j < DimensionGet()+1; j++) _samples(j) = runif(_Lower(j) , _Higher(j) ); 
	     one_sample.ValueSet(_samples);
	     return true;
      }
      default:
    	return false;
      }
  }

  ColumnVector
  Uniform::CenterGet (  ) const 
  { 
    return (_Higher+_Lower)/2.0;
  }

  ColumnVector
  Uniform::WidthGet () const
  {
    return (_Higher-_Lower);
  }

  void 
  Uniform::UniformSet (const ColumnVector& center,const ColumnVector& width)
  { 
    assert(center.rows() == width.rows());
    _Lower = center - width/2.0;
    _Higher = center + width/2.0;
    _Height = 1;
    for (int i=1 ; i < width.rows()+1 ; i++ )
    {
        _Height = _Height / width(i);
    }
    if (this->DimensionGet() == 0) this->DimensionSet(center.rows());
    assert(this->DimensionGet() == center.rows());
  }


} // End namespace BFL
