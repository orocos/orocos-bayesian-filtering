// $Id: conditionalUniformMeasPdf1d.h tdelaet $
// Copyright (C) 2007  Tinne De Laet <first dot last at mech dot kuleuven dot be>
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

#ifndef __CONDITIONAL_UNIFORM_PDF__
#define __CONDITIONAL_UNIFORM_PDF__

#include <pdf/conditionalpdf.h>
#include <pdf/gaussian.h>
#include <pdf/uniform.h>

namespace BFL
{
  /** Conditional Uniform Measurement Pdf for 1d mobile robot example
   * Specific class especially created for the measurement pdf for the 1d mobile robot localisation
   * example with a histogram filter.
   * The measurement is direct measurement (ultrasonic sensor) of the robot's 1d
   * position (=2*robot position).
   * The conditional distribution takes into account some extra Gaussian
   * measurement noise
   */
  class ConditionalUniformMeasPdf1d : public ConditionalPdf<MatrixWrapper::ColumnVector, int>
    {
    public:
      /// Constructor
      /** 
	  @param measNoise additiveNoise Pdf representing the extra additive noise
      */
      ConditionalUniformMeasPdf1d( const Gaussian& measNoise);

      /// Destructor
      virtual ~ConditionalUniformMeasPdf1d();

      virtual Probability ProbabilityGet(const MatrixWrapper::ColumnVector& measurement) const;

    private:
      Gaussian _measNoise;

    };

} // End namespace BFL
 
#endif //  
