// Copyright (C) 2007 Tinne De Laet <first dot last at mech dot kuleuven dot be>
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
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
//

#ifndef MODEL_TEST_HPP
#define MODEL_TEST_HPP

#include <cppunit/extensions/HelperMacros.h>
#include <model/discretesystemmodel.h>
#include <model/linearanalyticsystemmodel_gaussianuncertainty.h>
#include <model/linearanalyticmeasurementmodel_gaussianuncertainty.h>
#include <wrappers/matrix/matrix_wrapper.h>

using namespace std;
using namespace BFL;
using namespace MatrixWrapper;

class ModelTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE( ModelTest );
  CPPUNIT_TEST( testDiscreteSystemModel );
  CPPUNIT_TEST( testLinearAnalyticSystemModelGaussianUncertainty );
  CPPUNIT_TEST( testLinearAnalyticMeasurementModelGaussianUncertainty );
  CPPUNIT_TEST_SUITE_END();


public:
  void setUp();
  void tearDown();

  void testDiscreteSystemModel();
  void testLinearAnalyticSystemModelGaussianUncertainty();
  void testLinearAnalyticMeasurementModelGaussianUncertainty();

};

#endif  // MODEL_TEST_HPP
