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

    MobileRobot::MobileRobot( )
      {
	// sys noise
	SymmetricMatrix sys_noise_Cov(3);
	sys_noise_Cov(1,1) = pow(0.001,2); 
	sys_noise_Cov(1,2) = 0.0; 
	sys_noise_Cov(1,3) = 0.0; 
	sys_noise_Cov(2,1) = 0.0; 
	sys_noise_Cov(2,2) = pow(0.001,2); 
	sys_noise_Cov(2,3) = 0.0; 
	sys_noise_Cov(3,1) = 0.0; 
	sys_noise_Cov(3,2) = 0.0; 
	sys_noise_Cov(3,3) = pow(0.1*3.14/180,2); 
	ColumnVector sys_noise_Mu(3);
	sys_noise_Mu = 0.0;
	_sys_noise = new Gaussian(sys_noise_Mu, sys_noise_Cov);

	// meas noise
	SymmetricMatrix meas_noise_Cov(1);
	meas_noise_Cov(1,1) = pow(0.05,2);
	ColumnVector meas_noise_Mu(1);
	meas_noise_Mu(1) = 0.0;
	_meas_noise = new Gaussian(meas_noise_Mu, meas_noise_Cov);

	_meas_model = Matrix(1,3);
	_meas_model(1,1) = 0.0;
	_meas_model(1,2) = 2.0;
	_meas_model(1,3) = 0.0;

	// initial state
	_state = ColumnVector(3);  
	_state(1) = 0.0;
	_state(2) = 0.0;
	_state(3) = 0.8;
      }
    
    MobileRobot::~MobileRobot()
      {
	delete _sys_noise;
	delete _meas_noise;
      }
    
    void
    MobileRobot::Move(ColumnVector inputs)
    {
      // exact position
      _state(1) += cos(_state(3)) * inputs(1);
      _state(2) += sin(_state(3)) * inputs(1);
      _state(3) += inputs(2);

      // add uncertainty
      Sample<ColumnVector> sample;
      _sys_noise->SampleFrom(sample);
      _state += sample.ValueGet();
    }
    

    ColumnVector
    MobileRobot::Measure()
    {
      // exact model
      ColumnVector meas = _meas_model * _state;

      // add uncertainty
      Sample<ColumnVector> sample;
      _meas_noise->SampleFrom(sample);

      return meas + sample.ValueGet();
    }


    ColumnVector
    MobileRobot::GetState()
    {
      return _state;
    }
}

