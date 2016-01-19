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
#ifndef __BFL_CONSTANTS_H__
#define __BFL_CONSTANTS_H__

#define NUMERIC_PRECISION 0.000000001

#ifndef M_PI
#define M_PI 3.141592653589793284626433832795
#endif

// Check if double is probability value (maybe this should be solved
// by a type probability_t or something.  Needs thinking.
#include <iostream>
#include <cmath>
#include <cassert>

namespace BFL
{
  using namespace std;
  /// Class representing a probability (a double between 0 and 1)
  class Probability
  {
    public:
      Probability(){};
      Probability(double p)
    {
          #ifndef _MSC_VER
       assert(std::isfinite(p) != 0);
         #endif
      assert( p >= 0 );
      _prob = p;
    };
    virtual ~Probability(){};
/*
    ostream & operator<< (ostream & stream)
    {
      stream << this->getValue() << endl;
      return stream;
    };

  istream & operator>> (istream & stream)
    {
        double value;
        stream >> value;
        _prob = Probability(value);
        return stream; 
    };
*/
    friend ostream & operator<< (ostream & stream,Probability& prob);

    friend istream & operator>> (istream & stream,Probability& prob);

    double getValue() const { return _prob;} ;
    double& getValue() { return _prob;} ;

      operator double(){return _prob;};
      Probability operator *(Probability p)
    { return ((Probability) (this->_prob * (double) p));};
      Probability operator /(Probability p)
    { return ((Probability) (this->_prob / (double) p));};


  private:
    double _prob;
  };

} // End namespace

/*!\mainpage BFL
 *
 * <img
 * src="http://www.orocos.org/files/images/particles.img_assist_custom.png"
 * alt="Pallet localization with a Sick Scanner and BFL particle filter" width="550">
 *
 * \section Introduction
 *
 * Please see <a href="http://www.orocos.org/bfl">this
 * page</a> for more information on this software, as well as
 * download/installation/getting started instructions...
 *
 */

#endif


