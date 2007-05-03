// $Id$
// Copyright (C) 2003 Klaas Gadeyne <first dot last at gmail dot com>
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
#include "mcpdf.h"

namespace BFL
{
  using namespace MatrixWrapper;
  

  // Template Specialisation for T =ColumnVector
  template <> inline
  ColumnVector MCPdf<ColumnVector>::ExpectedValueGet (  ) const
  {
    ColumnVector CumSum(DimensionGet());
    CumSum=0.0;
    vector<WeightedSample<ColumnVector> > los = _listOfSamples;
    vector<WeightedSample<ColumnVector> >::iterator it;
    for ( it = los.begin() ; it != los.end() ; it++ )
      CumSum += ( it->ValueGet() * it->WeightGet() );
    return CumSum/_SumWeights;
  }


  template <> inline
  SymmetricMatrix MCPdf<ColumnVector>::CovarianceGet (  ) const
  {
    ColumnVector mean(this->ExpectedValueGet());
    ColumnVector diff(DimensionGet()); // Temporary storage
    SymmetricMatrix covariance(DimensionGet());
    Matrix diffsum(DimensionGet(), DimensionGet());
    vector<WeightedSample<ColumnVector> > los = _listOfSamples;
    diffsum = 0.0;
    // Actual calculation
    static vector<WeightedSample<ColumnVector> >::iterator it;
    for (it = los.begin(); it != los.end(); it++)
      {
	diff = (it->ValueGet() - mean);
	diffsum += diff * (diff.transpose() * it->WeightGet());
      }
    // Biased estimator!! (unbiased possible with weighted samples??)
    (diffsum/_SumWeights).convertToSymmetricMatrix(covariance);
    return covariance;
  }
 





  // Template Specialisation for T =unsigned int
  template <> inline
  unsigned int MCPdf<unsigned int>::ExpectedValueGet (  ) const
  {
    unsigned int result;
    double CumSum = 0;

    double current_weight;
    vector<WeightedSample<unsigned int> > los = _listOfSamples;
    vector<WeightedSample<unsigned int> >::iterator it;
    for ( it = los.begin() ; it != los.end() ; it++ )
      {
	current_weight = it->WeightGet();
	CumSum += ( ((double)it->ValueGet()) * current_weight );
      }
    result = (unsigned int)((CumSum/_SumWeights) + 0.5);
    return result;
  }




  template <> inline
  SymmetricMatrix MCPdf<unsigned int>::CovarianceGet (  ) const
  {
    unsigned int mean = this->ExpectedValueGet();
    unsigned int diff;
    double diffsum, current_weight;
    vector<WeightedSample<unsigned int> > los = _listOfSamples;
    diffsum = 0.0;
    // Actual calculation
    static vector<WeightedSample<unsigned int> >::iterator it;
    for (it = los.begin(); it != los.end(); it++)
      {
	current_weight = it->WeightGet();
	diff = (it->ValueGet() - mean);
	diffsum += (((double)(diff * diff)) * current_weight);
      }

    // Biased estimator!! (unbiased possible with weighted samples??)
    SymmetricMatrix Covariance(1);
    Covariance(1,1) = (diffsum / _SumWeights);
    return Covariance;
  }

}
