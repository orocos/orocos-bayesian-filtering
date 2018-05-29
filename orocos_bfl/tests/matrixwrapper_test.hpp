// Copyright (C) 2007 Wim Meeussen <wim.meeussen@mech.kuleuven.be>
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

#ifndef MATRIXWRAPPER_TEST_HPP
#define MATRIXWRAPPER_TEST_HPP

#include <cppunit/extensions/HelperMacros.h>
#include <wrappers/matrix/matrix_wrapper.h>
#include <wrappers/matrix/vector_wrapper.h>
#include <string>

using namespace std;
using namespace MatrixWrapper;

class MatrixwrapperTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE( MatrixwrapperTest );
  CPPUNIT_TEST( testMatrixwrapperValue );
  CPPUNIT_TEST_SUITE_END();

public:
  void setUp();
  void tearDown();

  void testMatrixwrapperValue();

};

#endif  // MATRIXWRAPPER_TEST_HPP
