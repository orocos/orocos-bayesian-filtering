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
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
//

#ifndef APPROX_EQUAL_HPP
#define APPROX_EQUAL_HPP

#include <wrappers/matrix/matrix_wrapper.h>
#include <wrappers/matrix/vector_wrapper.h>

bool approxEqual(double a, double b, double epsilon);
bool approxEqual(const MatrixWrapper::Matrix& a, const MatrixWrapper::Matrix& b, double epsilon);
bool approxEqual(const MatrixWrapper::SymmetricMatrix& a, const MatrixWrapper::SymmetricMatrix& b, double epsilon);
bool approxEqual(const MatrixWrapper::ColumnVector& a, const MatrixWrapper::ColumnVector& b, double epsilon);

#endif
