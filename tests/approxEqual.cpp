// Copyright (C) 2007 Wim Meeussen
//  
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//  
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//  
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
//  

#include "approxEqual.hpp"

bool approxEqual(double a, double b, double epsilon)
{
  return (a <= b+epsilon) && (a >= b-epsilon);
}


bool approxEqual(const MatrixWrapper::Matrix& a, const MatrixWrapper::Matrix& b, double epsilon)
{
  if (a.rows() != b.rows()) return false;
  if (a.columns() != b.columns()) return false;
  
  for (unsigned int r=0; r<a.rows(); r++)
    for (unsigned int c=0; c<a.columns(); c++)
      if (!approxEqual(a(r+1, c+1), b(r+1,c+1),epsilon)) return false;

  return true;
}


bool approxEqual(const MatrixWrapper::SymmetricMatrix& a, const MatrixWrapper::SymmetricMatrix& b, double epsilon)
{
  if (a.rows() != b.rows()) return false;
  if (a.columns() != b.columns()) return false;
  
  for (unsigned int r=0; r<a.rows(); r++)
    for (unsigned int c=0; c<a.columns(); c++)
      if (!approxEqual(a(r+1, c+1), b(r+1,c+1),epsilon)) return false;

  return true;
}

bool approxEqual(const MatrixWrapper::ColumnVector& a, const MatrixWrapper::ColumnVector& b, double epsilon)
{
  if (a.rows() != b.rows()) return false;
  
  for (unsigned int r=0; r<a.rows(); r++)
    if (!approxEqual(a(r+1), b(r+1),epsilon)) return false;

  return true;
}
