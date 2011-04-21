// $Id$
// Copyright (C) 2002 Klaas Gadeyne <first dot last at gmail dot com>
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

#include "measurementmodel.h"
#include <iostream>

using namespace std;

// Constructor
template<typename MeasVar, typename StateVar> 
MeasurementModel<MeasVar,StateVar>::MeasurementModel(ConditionalPdf<MeasVar,StateVar>* measurementpdf)
{
#ifdef __CONSTRUCTOR__
  cout << "MeasurementModel::Constructor" << endl;
#endif // __CONSTRUCTOR__
  if (measurementpdf != NULL)
    {
      switch(measurementpdf->NumConditionalArgumentsGet())
	{
	case 1: 
	  {
	    _systemWithoutSensorParams = true;
	    _MeasurementPdf  = measurementpdf;
	    break;
	  }
	case 2:
	  {
	    _systemWithoutSensorParams = false;
	    _MeasurementPdf  = measurementpdf;
	    break;
	  }
	default:{
	  cerr << "MeasurementModel::Constructor : MeasPdf can only have 1 or 2 conditional Arguments (x and u, in that order!))" << endl;
	  exit(-BFL_ERRMISUSE);
	}
	}
    }
}

// Destructor
template<typename MeasVar, typename StateVar>
 MeasurementModel<MeasVar,StateVar>::~MeasurementModel()
{
#ifdef __DESTRUCTOR__
  cout << "MeasurementModel::Destructor" << endl;
#endif // __DESTRUCTOR__
  /* KG: Probably a memory leak
     Who should clean this up? Sometimes the user will have created
     this Pdf, sometimes not (eg. by copy constructor).  If we allways
     delete it here.
     There has to be a cleaner way to implement this!
  */
  // delete _MeasurementPdf;
}

// BUG: Should have copy constructor here?

// Get Measurement Size
template<typename MeasVar, typename StateVar> int 
MeasurementModel<MeasVar,StateVar>::MeasurementSizeGet() const 
{ 
  return _MeasurementPdf->DimensionGet();
}

template<typename MeasVar, typename StateVar> bool 
MeasurementModel<MeasVar,StateVar>::SystemWithoutSensorParams() const 
{ 
  return _systemWithoutSensorParams;
}

// Get MeasurementPdf
template<typename MeasVar, typename StateVar> ConditionalPdf<MeasVar,StateVar>*
MeasurementModel<MeasVar,StateVar>::MeasurementPdfGet()
{ 
  return _MeasurementPdf;
}

// Set MeasurementPdf
template<typename MeasVar, typename StateVar> void 
MeasurementModel<MeasVar,StateVar>::MeasurementPdfSet(ConditionalPdf<MeasVar,StateVar> * pdf)
{ 
  assert(pdf != NULL);
  switch(pdf->NumConditionalArgumentsGet())
    {
    case 1: 
      {
	_systemWithoutSensorParams = true;
	_MeasurementPdf  = pdf;
	break;
      }
    case 2:
      {
	_systemWithoutSensorParams = false;
	_MeasurementPdf  = pdf;
	break;
      }
    default:
      {
	cerr << "MeasurementModel::Constructor : MeasPdf can only have 1 or 2 conditional Arguments (x and u, in that order!))" << endl;
	exit(-BFL_ERRMISUSE);
      }
    }
}

template<typename MeasVar, typename StateVar> MeasVar
MeasurementModel<MeasVar,StateVar>::Simulate (const StateVar& x, 
					      const StateVar& s, 
					      int sampling_method, 
					      void * sampling_args)
{
  assert(_systemWithoutSensorParams == false);
  _MeasurementPdf->ConditionalArgumentSet(0,x);
  _MeasurementPdf->ConditionalArgumentSet(1,s);
  Sample<StateVar> Simulated(MeasurementSizeGet());
  _MeasurementPdf->SampleFrom(Simulated, sampling_method,sampling_args);
  MeasVar result = Simulated.ValueGet();
  return result;
}


template<typename MeasVar, typename StateVar> MeasVar
MeasurementModel<MeasVar,StateVar>::Simulate(const StateVar& x, 
					     int sampling_method, 
					     void * sampling_args)
{
  assert(_systemWithoutSensorParams == true);
  _MeasurementPdf->ConditionalArgumentSet(0,x);
  Sample<StateVar> Simulated(MeasurementSizeGet());
  _MeasurementPdf->SampleFrom(Simulated, sampling_method,sampling_args);
  MatrixWrapper::ColumnVector result = Simulated.ValueGet();
  return result;
}

template <typename MeasVar, typename StateVar> Probability 
MeasurementModel<MeasVar,StateVar>::ProbabilityGet (const MeasVar& z, 
						    const StateVar& x, 
						    const StateVar& s)
{
  assert(_systemWithoutSensorParams == false);
  _MeasurementPdf->ConditionalArgumentSet(0,x);
  _MeasurementPdf->ConditionalArgumentSet(1,s);
  return _MeasurementPdf->ProbabilityGet(z);
}

template <typename MeasVar, typename StateVar> Probability
MeasurementModel<MeasVar,StateVar>::ProbabilityGet (const MeasVar& z, 
						    const StateVar& x)
{
  assert(_systemWithoutSensorParams == true);
  _MeasurementPdf->ConditionalArgumentSet(0,x);
  return _MeasurementPdf->ProbabilityGet(z);
}
