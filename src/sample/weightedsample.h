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
// Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
//  

#ifndef WEIGHTEDSAMPLE_H
#define WEIGHTEDSAMPLE_H

#include "sample.h"
#include <assert.h>


namespace BFL
{

  /** Template Class representing a sample of a continuous or
      discrete pdf, together with a weight denoting the relative
      importance of that sample.  Inheritance is virtual (only important
      for a particular class hybridweightedsample (see the
      rob/actsens/cubeincorner CVS tree)
  */
  template <typename T> class WeightedSample: virtual public Sample<T> 
    {
    protected:
      /// The weight
      double Weight;
  
    public:
      /// Constructor
      /**
	 @param dimension of the ColumnVector for the continuous samples,
	 number of discrete states for the discrete case
      */ 
      WeightedSample (int dimension = 0 );
      /// Destructor
      virtual ~WeightedSample();
      /// Copy constructor
      WeightedSample ( const WeightedSample<T> & my_weighted_sample );

      /// Get the weight
      /** @return the weight
       */
      double WeightGet (  ) const;

      /// Set the weight
      /** @param weight the weight :-)
	  @return true if weight is a reasonable value
       */
      void WeightSet ( double weight );
  
      /// Print a weighted sample
      /** @param stream the stream to be returned
	  @param mws the weighted sample to be printed
	  @return the stream :-)
      */
      template <typename S> friend ostream & operator<< (ostream & stream, 
							 WeightedSample<S> & mws);

      /// Operator =
      WeightedSample<T> & operator= (const WeightedSample<T> & my_sample);

      /// Turn sample into weighted one (weight = 1)
      WeightedSample<T> & operator= (const Sample<T> & my_sample);
    };


  template <typename T> WeightedSample<T>::WeightedSample(int dimension) 
    : Sample<T>(dimension){}

  template <typename T> WeightedSample<T>::~WeightedSample(){}

  template <typename T> WeightedSample<T>::WeightedSample (const WeightedSample<T> & mws) 
    : Sample<T>(mws)
    {  
      Weight = mws.Weight;
    }

  template <typename T> double WeightedSample<T>::WeightGet (  ) const 
    { 
      return Weight;
    }

  template <typename T> void WeightedSample<T>::WeightSet ( double weight )
    { 
      assert(weight >= 0);
      
      Weight = weight;
    }

  template <typename S> ostream & operator<< (ostream & stream, 
					      WeightedSample<S> & mws)
    {
      stream << "WeightedSample Value = " << (Sample<S> &) mws
	     << "Weight = " << mws.Weight << endl;
      return stream;
    }

  template <typename T> WeightedSample<T> & WeightedSample<T>::operator= (const WeightedSample<T> & my_sample)
    {
      // TODO: Does this automatically calls the = operator in the
      // baseclass?  NO!!!
      Sample<T> * op1; const Sample<T> * op2;
      op1 = this; op2 = & my_sample;
      *op1 = *op2;
      this->Weight = my_sample.WeightGet();
      return *this;
    }

  // Turn sample into weighted one (weight = 1)
  template <typename T> WeightedSample<T> & WeightedSample<T>::operator= (const Sample<T> & my_sample)
    {
      //: Does this automatically calls the = operator in the baseclass?
      Sample<T> * op1; const Sample<T> * op2;
      op1 = this; op2 = & my_sample;
      *op1 = *op2;
      this->Weight = 1;
      return *this;
    }

} // End namespace BFL

#endif
