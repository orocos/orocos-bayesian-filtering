// Copyright (C) 2008 Willow Garage inc.
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

#ifndef EKF_TEST_HPP
#define EKF_TEST_HPP

#include <cppunit/extensions/HelperMacros.h>
#include <filter/extendedkalmanfilter.h>
#include <model/linearanalyticsystemmodel_gaussianuncertainty.h>
#include <model/linearanalyticmeasurementmodel_gaussianuncertainty.h>
#include <pdf/linearanalyticconditionalgaussian.h>


namespace BFL
{

class EKFTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE( EKFTest );
  CPPUNIT_TEST( testMeasUpdate );
  CPPUNIT_TEST_SUITE_END();

public:
  void setUp();
  void tearDown();

  void testMeasUpdate();

private:
  double epsilon_;
  LinearAnalyticMeasurementModelGaussianUncertainty* meas_model_;
  LinearAnalyticConditionalGaussian*                 meas_pdf_;
  Gaussian*                                          prior_;
  ExtendedKalmanFilter*                              filter_;

};
}

#endif  // EKF_TEST_HPP
