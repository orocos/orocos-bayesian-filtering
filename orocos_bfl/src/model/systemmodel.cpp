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
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
//
#include "systemmodel.h"

// Constructor
template<typename T> SystemModel<T>::SystemModel(ConditionalPdf<T,T>* systempdf)
{
#ifdef __CONSTRUCTOR__
  cout << "SystemModel::Constructor" << endl;
#endif // __CONSTRUCTOR__
  if (systempdf != NULL)
    {
      switch(systempdf->NumConditionalArgumentsGet())
	{
	case 1:
	  {
	    _systemWithoutInputs = true;
	    _SystemPdf  = systempdf;
	    break;
	  }
	case 2:
	  {
	  _systemWithoutInputs = false;
	  _SystemPdf  = systempdf;
	  break;
	  }
	default:
	  {
	    cerr << "SystemModel::Constructor : SystemPdf can only have 1 or 2 conditional Arguments (x and u, in that order!))" << endl;
	    exit(-BFL_ERRMISUSE);
	  }
	}
    }
}

// Destructor
template<typename T>
SystemModel<T>::~SystemModel()
{
#ifdef __DESTRUCTOR__
  cout << "SystemModel::Destructor" << endl;
#endif // __DESTRUCTOR__
  /* KG: Probably a memory leak
     Who should clean this up? Sometimes the user will have created
     this Pdf, sometimes not (eg. by copy constructor).  If we allways
     delete it here.
     There has to be a cleaner way to implement this!
  */
  // delete SystemPdf;
}

// Copy constructor
/*
template<typename T>
SystemModel<T>::SystemModel(const SystemModel<T>& model)
{
  SystemPdf  = &(model.SystemPdfGet());
}
*/

// Get State Size
template<typename T> int
SystemModel<T>::StateSizeGet() const
{
  return _SystemPdf->DimensionGet();
}

template<typename T> bool
SystemModel<T>::SystemWithoutInputs() const
{
  return _systemWithoutInputs;
}

// Get SystemPdf
template<typename T> ConditionalPdf<T,T>*
SystemModel<T>::SystemPdfGet()
{
  return _SystemPdf;
}

// Set SystemPdf
template<typename T> void
SystemModel<T>::SystemPdfSet(ConditionalPdf<T,T>* pdf)
{
  assert(pdf != NULL);
  switch(pdf->NumConditionalArgumentsGet())
    {
    case 1:
      {
	_systemWithoutInputs = true;
	_SystemPdf  = pdf;
	break;
      }
    case 2:
      {
	_systemWithoutInputs = false;
	_SystemPdf  = pdf;
	break;
      }
    default:
      {
	cerr << "SystemModel::SystemPdfSet() : SystemPdf can only have 1 or 2 conditional Arguments (x and u, in that order!))" << endl;
	exit(-BFL_ERRMISUSE);
      }
    }
}

// Simulate from the system model
template<typename T> T
SystemModel<T>::Simulate (const T& x, const T& u, const SampleMthd sampling_method,
		          void * sampling_args)
{
  assert(_systemWithoutInputs == false);
  _SystemPdf->ConditionalArgumentSet(0,x);
  _SystemPdf->ConditionalArgumentSet(1,u);
  Sample<T> Simulated(StateSizeGet());
  _SystemPdf->SampleFrom(Simulated, sampling_method,sampling_args);
  T result = Simulated.ValueGet();
  return result;
}

template<typename T> T
SystemModel<T>::Simulate (const T& x, const SampleMthd sampling_method,
			  void * sampling_args)
{
  assert(_systemWithoutInputs == true);
  _SystemPdf->ConditionalArgumentSet(0,x);
  Sample<T> Simulated(StateSizeGet());
  _SystemPdf->SampleFrom(Simulated, sampling_method,sampling_args);
  T result = Simulated.ValueGet();
  return result;
}

template <typename T> Probability
SystemModel<T>::ProbabilityGet (const T& x_k, const T& x_kminusone,
                                const T& u)
{
  assert(_systemWithoutInputs == false);
  _SystemPdf->ConditionalArgumentSet(0,x_kminusone);
  _SystemPdf->ConditionalArgumentSet(1,u);
  return _SystemPdf->ProbabilityGet(x_k);
}

template <typename T> Probability
SystemModel<T>::ProbabilityGet (const T& x_k, const T& x_kminusone)
{
  assert(_systemWithoutInputs == true);
  _SystemPdf->ConditionalArgumentSet(0,x_kminusone);
  return _SystemPdf->ProbabilityGet(x_k);
}
