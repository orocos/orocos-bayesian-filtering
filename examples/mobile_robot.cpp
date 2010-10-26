// $Id: mobile_robot.cpp tdelaet $
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

#include "mobile_robot.h"

using namespace MatrixWrapper;

namespace BFL
{

    MobileRobot::MobileRobot():
        _state(STATE_SIZE)
      {
	// initial state
	_state(1) = X_0;
	_state(2) = Y_0;
	_state(3) = THETA_0;

	// sys noise
	ColumnVector sys_noise_Mu(STATE_SIZE);
	sys_noise_Mu(1) = MU_SYSTEM_NOISE_X_ROB;
	sys_noise_Mu(2) = MU_SYSTEM_NOISE_Y_ROB;
	sys_noise_Mu(3) = MU_SYSTEM_NOISE_THETA_ROB;
	SymmetricMatrix sys_noise_Cov(STATE_SIZE);
    sys_noise_Cov = 0.0;
	sys_noise_Cov(1,1) = SIGMA_SYSTEM_NOISE_X_ROB;
	sys_noise_Cov(1,2) = 0.0;
	sys_noise_Cov(1,3) = 0.0;
	sys_noise_Cov(2,1) = 0.0;
	sys_noise_Cov(2,2) = SIGMA_SYSTEM_NOISE_Y_ROB;
	sys_noise_Cov(2,3) = 0.0;
	sys_noise_Cov(3,1) = 0.0;
	sys_noise_Cov(3,2) = 0.0;
	sys_noise_Cov(3,3) = SIGMA_SYSTEM_NOISE_THETA_ROB;
    _system_Uncertainty = new Gaussian(sys_noise_Mu, sys_noise_Cov);

    // create the model
    _sys_pdf = new NonLinearAnalyticConditionalGaussianMobile(*_system_Uncertainty);
    _sys_model = new AnalyticSystemModelGaussianUncertainty(_sys_pdf);

	// meas noise
	SymmetricMatrix meas_noise_Cov(MEAS_SIZE);
	meas_noise_Cov(1,1) = SIGMA_MEAS_NOISE_ROB;
	ColumnVector meas_noise_Mu(MEAS_SIZE);
	meas_noise_Mu(1) = MU_MEAS_NOISE_ROB;
	_measurement_Uncertainty = new Gaussian(meas_noise_Mu, meas_noise_Cov);

    // create matrix _meas_model for linear measurement model
    double wall_ct = 2/(sqrt(pow(RICO_WALL,2.0) + 1));
	Matrix H(MEAS_SIZE,STATE_SIZE);
    H = 0.0;
    H(1,1) = wall_ct * RICO_WALL;
    H(1,2) = 0 - wall_ct;
    H(1,3) = 0.0;

    // create the measurement model
    _meas_pdf = new LinearAnalyticConditionalGaussian(H, *_measurement_Uncertainty);
    _meas_model = new LinearAnalyticMeasurementModelGaussianUncertainty(_meas_pdf);

      }

    MobileRobot::~MobileRobot()
      {
    delete _system_Uncertainty;
    delete _sys_pdf;
	delete _sys_model;
	delete _measurement_Uncertainty;
	delete _meas_pdf;
	delete _meas_model;
      }

    void
    MobileRobot::Move(ColumnVector inputs)
    {
      _state = _sys_model->Simulate(_state,inputs);
    }


    ColumnVector
    MobileRobot::Measure()
    {
      return _meas_model->Simulate(_state);
    }


    ColumnVector
    MobileRobot::GetState()
    {
      return _state;
    }
}

