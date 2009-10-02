// Copyright  (C)  2009  Tinne De Laet <tinne dot delaet at mech dot kuleuven dot be>

// Author: Tinne De Laet 
// Maintainer: Tinne De Laet

// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.

// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.

// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

#ifndef BFL_TOOLKIT_HPP
#define BFL_TOOLKIT_HPP

#include <rtt/ToolkitPlugin.hpp>

//using namespace MatrixWrapper;
namespace BFL{

class bflToolkitPlugin: public RTT::ToolkitPlugin {
public:
    virtual std::string getName();

    virtual bool loadTypes();
    virtual bool loadConstructors();
    virtual bool loadOperators();
};

  extern bflToolkitPlugin bflToolkit;

    template<typename T>
    struct VectorAssignChecker
        : public std::binary_function<T, T, bool>
    {
        bool operator()(const T& v1, const T& v2) const
        {
            return v1.size()==v2.size();
        }
    };

    template<typename T>
    struct MatrixAssignChecker
        : public std::binary_function<T, T, bool>
    {
        bool operator()(const T& m1, const T& m2) const
        {
            return (m1.rows()==m2.rows())&&(m1.columns()==m2.columns());
        }
    };

    template<typename T>
    struct MatrixIndexChecker
        : public std::binary_function< T, unsigned int, bool>
    {
        bool operator()(const T& m, unsigned int i ) const
        {
            return (i > 0) && (i < m.rows());
        }
    };


}
#endif

