// $Id: mixture.cpp 2009-01-22 tdelaet $
// Copyright (C) 2009 Tinne De Laet <first dot last at gmail dot com>
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

// This file only contains template specialisation code
#include "mixture.h"

namespace BFL
{
  using namespace MatrixWrapper;


//TODO: this still contains allocation
  // Template Specialisation for T =ColumnVector
  template <> inline
  ColumnVector Mixture<ColumnVector>::ExpectedValueGet (  ) const
  {
    TestNotInit();
    ColumnVector expectedValue(DimensionGet());
    expectedValue = 0.0;
    for (int i=0; i<NumComponentsGet();i++)
        expectedValue = expectedValue + ( (*_componentPdfs)[i]->ExpectedValueGet() * (double)((*_componentWeights)[i]) );
    return expectedValue;
  }


  // Template Specialisation for T =unsigned int
  // THIS RETURNS THE MOST PROBABLE STATE
  // TODO: only works for DiscretePdf
  template <> inline
  unsigned int Mixture<unsigned int>::ExpectedValueGet (  ) const
  {
    TestNotInit();
    unsigned int expectedValue =0;
    double probState = 0.0;
    double mostProbState = -1.0;
    int numStates = ((DiscretePdf*)((*_componentPdfs)[0]))->NumStatesGet();
    for (int i=0; i<NumComponentsGet();i++) //loop over all components
    {
        if(numStates != ((DiscretePdf*)((*_componentPdfs)[i]))->NumStatesGet())
        {
            cerr << "Mixture::ExpectedValueGet failed since the different components in the mixture don't have the same number of states" << endl; 
            assert(0);
        }
    }
    for (int j=0; j<numStates ;j++ ) //loop over all states
    {
        probState = 0.0;
        for (int i=0; i<NumComponentsGet();i++) //loop over all components
        {
            probState  +=   ( ((double)(*_componentPdfs)[i]->ProbabilityGet(j)) * (double)((*_componentWeights)[i]) );
        }
        if (probState > mostProbState)
        {
            expectedValue = (unsigned int)j;
            mostProbState = probState;
        }
    }
    return expectedValue;
  }


  // Template Specialisation for T =int
  // THIS RETURNS THE MOST PROBABLE STATE
  // TODO: only works for DiscretePdf
  template <> inline
  int Mixture<int>::ExpectedValueGet (  ) const
  {
    TestNotInit();
    int expectedValue =0;
    double probState = 0.0;
    double mostProbState = -1.0;
    int numStates = ((DiscretePdf*)((*_componentPdfs)[0]))->NumStatesGet();
    for (int i=0; i<NumComponentsGet();i++) //loop over all components
    {
        if(numStates != ((DiscretePdf*)((*_componentPdfs)[i]))->NumStatesGet())
        {
            cerr << "Mixture::ExpectedValueGet failed since the different components in the mixture don't have the same number of states" << endl; 
            assert(0);
        }
    }
    for (int j=0; j< numStates;j++ ) //loop over all states
    {
        probState = 0.0;
        for (int i=0; i<NumComponentsGet();i++) //loop over all components
        {
            probState  +=   ( ((double)(*_componentPdfs)[i]->ProbabilityGet(j)) * (double)((*_componentWeights)[i]) );
        }
        if (probState > mostProbState)
        {
            expectedValue = j;
            mostProbState = probState;
        }
    }
    return expectedValue;
  }


  // Template Specialisation for T =double
  template <> inline
  double Mixture<double>::ExpectedValueGet (  ) const
  {
    TestNotInit();
    double expectedValue;
    expectedValue = 0.0;
    for (int i=0; i<NumComponentsGet();i++)
        expectedValue = expectedValue + ( (*_componentPdfs)[i]->ExpectedValueGet() * (double)((*_componentWeights)[i]) );
    return expectedValue;
  }
}

