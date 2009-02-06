// Copyright (C) 2007 Wim Meeussen <wim DOT meeussen AT mech DOT kuleuven DOT be>
// Copyright (C) 2008 Tinne De Laet <tinne DOT delaet AT mech DOT kuleuven DOT be>
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

#ifndef COMPLETE_FILTER_TEST_HPP
#define COMPLETE_FILTER_TEST_HPP

#include <cppunit/extensions/HelperMacros.h>
#include <wrappers/matrix/matrix_wrapper.h>
#include <wrappers/matrix/vector_wrapper.h>

#include <filter/bootstrapfilter.h>
#include <filter/extendedkalmanfilter.h>
#include <filter/iteratedextendedkalmanfilter.h>
#include <filter/asirfilter.h>
#include <filter/EKparticlefilter.h>
#include <filter/histogramfilter.h>

#include <model/analyticsystemmodel_gaussianuncertainty.h>
#include <model/linearanalyticmeasurementmodel_gaussianuncertainty.h>
#include <pdf/analyticconditionalgaussian.h>
#include <pdf/discreteconditionalpdf.h>
#include <pdf/discretepdf.h>
#include <pdf/gaussian.h>

#include "../examples/mobile_robot_wall_cts.h"
#include "../examples/compare_filters/nonlinearanalyticconditionalgaussianmobile.h"
#include "../examples/mobile_robot.h"
#include "../examples/discrete_filter/conditionalUniformMeasPdf1d.h"

#include <iostream>
#include <fstream>
#include <string>


using namespace std;
using namespace BFL;
using namespace MatrixWrapper;

class Complete_FilterTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE( Complete_FilterTest );
  CPPUNIT_TEST( testComplete_FilterValue_Cont );
  CPPUNIT_TEST( testComplete_FilterValue_Discr );
  CPPUNIT_TEST_SUITE_END();

public:
  void setUp();
  void tearDown();

  void testComplete_FilterValue_Cont();
  void testComplete_FilterValue_Discr();

};

#endif  // COMPLETE_FILTER_TEST_HPP
