// $Id: innovationCheck.h  tdelaet $
// Copyright (C) 2007 Tinne De Laet <first dot last at mech dot kuleuven dot be>
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

#ifndef __INNOVATION_CHECK__
#define __INNOVATION_CHECK__

#include <vector>
#include "../wrappers/matrix/vector_wrapper.h"
#include "../wrappers/matrix/matrix_wrapper.h"

namespace BFL
{

  /// Class implementing an innovationCheck used in IEKF
  /** This is a class implementing the innovationCheck class used for
      Iterated Extended Kalman Filters innovation check.
  */
  class InnovationCheck
    {

    public:
      /** Constructor
	  @pre
	  @param min_innovation is the minimal innovation desired per iteration of the iteratedExtendedKalmanFilter
      */
      InnovationCheck(double min_innovation = 0.0 );

      /// Destructor
      virtual ~InnovationCheck();

      /// check Innovation
      /** Returns true if the innovation is still considered 'larger' than the
          minimal innovation
          @param innovation innovation to be checked
     */
      bool check(MatrixWrapper::ColumnVector innovation);

    private:
      ///the minimal innovation desired per iteration of the iteratedExtendedKalmanFilter
      double min_innovation;

    };  // class

} // End namespace BFL

#endif // __INNOVATION_CHECK__

