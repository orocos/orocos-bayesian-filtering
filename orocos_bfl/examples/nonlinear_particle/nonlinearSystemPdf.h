// $Id: nonlinearanalyticconditionalgaussianmobile.h 5374 2005-05-06 14:57:05Z TDeLaet $
// Copyright (C) 2006  Tinne De Laet <first dot last at mech dot kuleuven dot be>
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


#ifndef __NON_LINEAR_SYSTEM_MOBILE__
#define __NON_LINEAR_SYSTEM_MOBILE__

#include <pdf/conditionalpdf.h>
#include <pdf/gaussian.h>

namespace BFL
{
  /// Non Linear Conditional Gaussian
  class NonlinearSystemPdf : public ConditionalPdf<MatrixWrapper::ColumnVector, MatrixWrapper::ColumnVector>
    {
    public:
      /// Constructor
      /** @param additiveNoise Pdf representing the additive Gaussian uncertainty
      */
      NonlinearSystemPdf( const Gaussian& additiveNoise);

      /// Destructor
      virtual ~NonlinearSystemPdf();

      // implement this virtual function for system model of a particle filter
      virtual bool SampleFrom (Sample<MatrixWrapper::ColumnVector>& one_sample, const SampleMthd method=SampleMthd::DEFAULT, void * args=NULL) const;

    private:
      Gaussian _additiveNoise;

    };

} // End namespace BFL

#endif //
