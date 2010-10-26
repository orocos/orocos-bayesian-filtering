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

#ifndef SAMPLE_TEST_HPP
#define SAMPLE_TEST_HPP

#include <cppunit/extensions/HelperMacros.h>
#include <sample/sample.h>
#include <sample/weightedsample.h>
#include <string>

using namespace std;
using namespace BFL;
using namespace MatrixWrapper;

class SampleTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE( SampleTest );
  CPPUNIT_TEST( testSampleValue );
  CPPUNIT_TEST( testSampleDimension );
  CPPUNIT_TEST( testWeightedSample );
  CPPUNIT_TEST_SUITE_END();

  ColumnVector a;
  Sample<ColumnVector> a_sample_cv;
  Sample<ColumnVector> b_sample_cv;
  Sample<int> a_sample_int;
  Sample<double> a_sample_double;

  WeightedSample<ColumnVector> a_weighted_sample_cv;
  WeightedSample<ColumnVector> b_weighted_sample_cv;

public:
  void setUp();
  void tearDown();

  void testSampleValue();
  void testSampleDimension();
  void testWeightedSample();

};

#endif  // SAMPLE_TEST_HPP
