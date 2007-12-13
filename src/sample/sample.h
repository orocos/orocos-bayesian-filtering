// $Id$
// Copyright (C) 2002 Klaas Gadeyne <first dot last at gmail dot com>
 /***************************************************************************
 *   This library is free software; you can redistribute it and/or         *
 *   modify it under the terms of the GNU General Public                   *
 *   License as published by the Free Software Foundation;                 *
 *   version 2 of the License.                                             *
 *                                                                         *
 *   As a special exception, you may use this file as part of a free       *
 *   software library without restriction.  Specifically, if other files   *
 *   instantiate templates or use macros or inline functions from this     *
 *   file, or you compile this file and link it with other files to        *
 *   produce an executable, this file does not by itself cause the         *
 *   resulting executable to be covered by the GNU General Public          *
 *   License.  This exception does not however invalidate any other        *
 *   reasons why the executable file might be covered by the GNU General   *
 *   Public License.                                                       *
 *                                                                         *
 *   This library is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU     *
 *   Lesser General Public License for more details.                       *
 *                                                                         *
 *   You should have received a copy of the GNU General Public             *
 *   License along with this library; if not, write to the Free Software   *
 *   Foundation, Inc., 59 Temple Place,                                    *
 *   Suite 330, Boston, MA  02111-1307  USA                                *
 *                                                                         *
 ***************************************************************************/ 
#ifndef SAMPLE_H
#define SAMPLE_H

// Adjust this to absolute path!
// Using <> makes that the g++ standard implementation of vector is chosen!
#include "../bfl_err.h"
#include "../wrappers/matrix/vector_wrapper.h"
#include "../wrappers/matrix/matrix_wrapper.h"
#include <iostream>

namespace BFL
{
  using namespace std;

  /** Template Class representing a basic sample of a continuous or
      discrete pdf
  */
  template <typename T> class Sample 
    {
    protected:
      /// The Sample Value
      T Value;
  
    public:
      /// Constructor
      /**
	 @param dimension of the ColumnVector for the continuous
	 samples.  This parameter is ignored in the discrete case.
      */   
      Sample (unsigned int dimension = 0);

      /// Destructor
      virtual ~Sample();

      /// Copy Constructor
      Sample ( const Sample<T> & my_sample );

      /// Get the value of the Sample
      T& ValueGet (  ) ;

      /// Get the value of the Sample
      const T& ValueGet (  ) const;

      // dimension get
      unsigned int DimensionGet () const;

      // dimension set
      void DimensionSet (unsigned int dim);

      /// Set the value of the Sample
      /** 
	  @param value the value indeed :-)
      */
      void ValueSet ( const T& value );
  
      /// Print a sample
      /** @param stream the stream to be returned
	  @param my_sample the sample to be printed
	  @return the stream :-)
      */
      template <typename S> friend ostream & operator<< (ostream & stream, 
							 Sample<S> & my_sample);
      /// Operator = 
      Sample & operator= (const Sample & my_sample);
  
    };








  // constructor
  template <typename T> Sample<T>::Sample (unsigned int dimension)
    {};


  // destructor
  template <typename T> Sample<T>::~Sample( )
    {};


  // copy constructor
  template <typename T> Sample<T>::Sample ( const Sample & my_sample )
    {
      Value     = my_sample.ValueGet();
    }


  // set value
  template <typename T> void Sample<T>::ValueSet (const T& value)
    {
       Value = value;
    }


  // get value
  template <typename T> T& Sample<T>::ValueGet (  ) 
    {
      return Value;
    }


  // get value
  template <typename T> const T& Sample<T>::ValueGet (  )  const
    {
      return Value;
    }

  // get dimension
  template <typename T> unsigned int Sample<T>::DimensionGet (  )  const
    {
      return 0;
    }

  // set dimension
  template <typename T> void Sample<T>::DimensionSet (unsigned int dim)  
    {}

  // stream
  template <typename S> ostream & operator<< (ostream & stream, Sample<S> & my_sample)
    {
      stream << my_sample.ValueGet() << endl;
      return stream;
    }

  // operator =
  template <typename T> Sample<T> & Sample<T>::operator= ( const Sample<T> & my_sample)
    {
      Value     = my_sample.ValueGet();
      return *this;
    }



} // End namespace BFL

#include "sample.cpp"

#endif
