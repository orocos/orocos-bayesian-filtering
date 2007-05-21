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
// Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
//  

#ifndef SMOOTHER_TEST_HPP
#define SMOOTHER_TEST_HPP

#include <cppunit/extensions/HelperMacros.h>
#include <wrappers/matrix/matrix_wrapper.h>
#include <smoother/rauchtungstriebel.h>
#include <smoother/particlesmoother.h>

using namespace std;
using namespace BFL;
using namespace MatrixWrapper;

class SmootherTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE( SmootherTest );
  CPPUNIT_TEST( testKalmanSmoother );
  CPPUNIT_TEST( testParticleSmoother );
  CPPUNIT_TEST_SUITE_END();

public:
  void setUp();
  void tearDown();
  
  void testKalmanSmoother();
  void testParticleSmoother();

private:
  double epsilon;
  
};

#endif  // SMOOTHER_TEST_HPP
