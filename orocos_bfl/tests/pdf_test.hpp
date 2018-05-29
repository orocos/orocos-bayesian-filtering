// Copyright (C) 2007 Klaas Gadeyne <first dot last at gmail dot com>
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

#ifndef PDF_TEST_HPP
#define PDF_TEST_HPP

#include <cppunit/extensions/HelperMacros.h>
#include <pdf/pdf.h>
#include <pdf/gaussian.h>
#include <pdf/uniform.h>
#include <pdf/discretepdf.h>
#include <pdf/linearanalyticconditionalgaussian.h>
#include <pdf/discreteconditionalpdf.h>
#include <pdf/mcpdf.h>
#include <pdf/mixture.h>
#include <wrappers/matrix/matrix_wrapper.h>

#include <iostream>
using namespace std;
using namespace BFL;
using namespace MatrixWrapper;

class PdfTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE( PdfTest );
  CPPUNIT_TEST( testMcpdfType );
  CPPUNIT_TEST( testUniform );
  CPPUNIT_TEST( testDiscretePdf );
  CPPUNIT_TEST( testLinearAnalyticConditionalGaussian );
  CPPUNIT_TEST( testDiscreteConditionalPdf );
  CPPUNIT_TEST( testMcpdf );
  CPPUNIT_TEST( testMixture );
  CPPUNIT_TEST( testGaussian );
  CPPUNIT_TEST_SUITE_END();

  ColumnVector _mu;
  SymmetricMatrix _sigma;
  ColumnVector _width;
  void testMcpdfType();

public:
  void setUp();
  void tearDown();

  void testGaussian();
  void testUniform();
  void testDiscretePdf();
  void testLinearAnalyticConditionalGaussian();
  void testDiscreteConditionalPdf();
  void testMcpdf();
  void testMixture();

private:
  double epsilon;

};

#endif  // PDF_TEST_HPP
