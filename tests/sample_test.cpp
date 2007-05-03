// Copyright (C) 2007 Klaas Gadeyne <first dot last at gmail dot com>
//  
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//  
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//  
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
//  
 
#include "sample_test.hpp"

// Registers the fixture into the 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION( SampleTest );
using namespace BFL;

void 
SampleTest::setUp()
{
  a = ColumnVector(3);
  a(1) = 1; a(2) = 2; a(3) = 1;
  a_sample_cv.ValueSet(a);
  a_sample_int.ValueSet(2);
  a_sample_double.ValueSet(3.1);

  b_sample_cv = a_sample_cv;

}


void 
SampleTest::tearDown()
{
}

void 
SampleTest::testSampleValue()
{
    // Test one way
    CPPUNIT_ASSERT_EQUAL( a, a_sample_cv.ValueGet());
    CPPUNIT_ASSERT_EQUAL( 2.0  , (a_sample_cv.ValueGet())(2));
    //    CPPUNIT_ASSERT_EQUAL( a_sample_cv, b_sample_cv);
    CPPUNIT_ASSERT_EQUAL( a_sample_cv.ValueGet(), b_sample_cv.ValueGet());

    CPPUNIT_ASSERT_EQUAL( 2, a_sample_int.ValueGet());
    CPPUNIT_ASSERT_EQUAL( 3.1, a_sample_double.ValueGet());
}

void
SampleTest::testSampleDimension()
{
  CPPUNIT_ASSERT_EQUAL((unsigned int) 3,a_sample_cv.DimensionGet());
  CPPUNIT_ASSERT_EQUAL((unsigned int) 3,b_sample_cv.DimensionGet());
  CPPUNIT_ASSERT_EQUAL((unsigned int) 1,a_sample_int.DimensionGet());
  CPPUNIT_ASSERT_EQUAL((unsigned int) 1,a_sample_double.DimensionGet());
}

void
SampleTest::testWeightedSample()
{
  double weight = 0.75;
  a_weighted_sample_cv.ValueSet(a);
  CPPUNIT_ASSERT_EQUAL( true,a_weighted_sample_cv.WeightSet(weight));

  b_weighted_sample_cv = a_weighted_sample_cv;
  CPPUNIT_ASSERT_EQUAL( a_weighted_sample_cv.ValueGet(), b_weighted_sample_cv.ValueGet());
  CPPUNIT_ASSERT_EQUAL( a_weighted_sample_cv.WeightGet(), b_weighted_sample_cv.WeightGet());
}

// void 
// SampleTest::testTimeProgress()
// {
//     CPPUNIT_ASSERT( t !=  hbg->getTicks() );
// }

    

