// $Id: mobile_robot.h tdelaet $
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

#ifndef MOBILE_ROBOT_HPP
#define MOBILE_ROBOT_HPP


#include <model/analyticsystemmodel_gaussianuncertainty.h>
#include <model/linearanalyticmeasurementmodel_gaussianuncertainty.h>
#include <pdf/gaussian.h>
#include <wrappers/matrix/matrix_wrapper.h>
#include <wrappers/matrix/vector_wrapper.h>


namespace BFL{

  class MobileRobot
    {
    public:
      // Constructor
      MobileRobot();
      ~MobileRobot();

      void Move(MatrixWrapper::ColumnVector inputs);

      MatrixWrapper::ColumnVector Measure();
      MatrixWrapper::ColumnVector GetState(); //method only for simulation purposes

    private:
      Gaussian* _sys_noise;
      Gaussian* _meas_noise;
      MatrixWrapper::Matrix _meas_model;
      MatrixWrapper::ColumnVector _state;
    };
}

#endif
