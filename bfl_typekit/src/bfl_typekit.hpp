// Copyright  (C)  2010  Tinne De Laet <tinne dot delaet at mech dot kuleuven dot be>

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

#include <rtt/types/StructTypeInfo.hpp>
#include <rtt/types/VectorTemplateComposition.hpp>
#include <rtt/types/SequenceTypeInfo.hpp>
#include <rtt/types/TemplateConstructor.hpp>
#include <rtt/types/TypeInfoRepository.hpp>    
#include <rtt/types/TypekitPlugin.hpp>
#include <rtt/Property.hpp>
#include <rtt/PropertyBag.hpp>
#include <rtt/types/SequenceConstructor.hpp>
#include <rtt/types/TemplateConstructor.hpp>

#include <bfl/wrappers/matrix/vector_wrapper.h>
#include <bfl/wrappers/matrix/matrix_wrapper.h>
#include <bfl/bfl_constants.h>
#include <bfl/sample/sample.h>
#include <bfl/sample/weightedsample.h>

namespace BFL{

class bflTypekitPlugin: public RTT::types::TypekitPlugin {
    public:
        virtual std::string getName();
    
        virtual bool loadTypes();
        virtual bool loadConstructors();
        virtual bool loadOperators();
    };
}
#endif

