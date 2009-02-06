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

#include "mobile_robot_wall_cts.h"
#include "nonlinearanalyticconditionalgaussianmobile.h"


namespace BFL{

/// This is a class simulating a mobile robot 
/** The state of the mobile robot is represented with a ColumnVector of three
* elements: the x and y position and the orientation. 
* The inputs of the robot are the linear velocity and the angular velocity.
* The mobile robot is equipped with a ultrasonic sensor returning the distance
* to a wall.
* The initial position of the mobile robot is read from mobile_robot_wall_cts.h
* During construction time the measurement model and system model are
* constructed and their properties are read from mobile_robot_wall_cts.h
*/
     

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
      Gaussian* _system_Uncertainty;
      NonLinearAnalyticConditionalGaussianMobile* _sys_pdf;
      AnalyticSystemModelGaussianUncertainty* _sys_model;
      Gaussian* _measurement_Uncertainty;
      LinearAnalyticConditionalGaussian* _meas_pdf;
      LinearAnalyticMeasurementModelGaussianUncertainty* _meas_model;
      MatrixWrapper::ColumnVector _state;
    };
}

#endif
