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

#include "bfl_typekit.hpp"

using namespace MatrixWrapper;
namespace BFL{ 

    using namespace RTT;
    using namespace RTT::detail;

    bflTypekitPlugin bflTypekit;

    void loadProbabilityTypes();
    void loadProbabilitysTypes();
    void loadVectorTypes();
    void loadVectorsTypes();
    void loadMatrixTypes();
    void loadMatrixsTypes();
    void loadSampleTypes();
    void loadWeightedSampleTypes();
    void loadSamplesTypes();
    void loadWeightedSamplesTypes();

    std::string bflTypekitPlugin::getName()
    {
        return "BFL_Typekit";
    }

    bool bflTypekitPlugin::loadTypes()
    {
        // load probability types  
        loadProbabilityTypes();
        // load probabilitys types  
        loadProbabilitysTypes();
        // load vector types
        loadVectorTypes();
        // load vectors types
        loadVectorsTypes();
        // load matrix types
        loadMatrixTypes();
        // load matrixs types
        loadMatrixsTypes();
        // load Sample types
        loadSampleTypes();
        // load Samples types
        loadSamplesTypes();
        // load WeightedSample types
        loadWeightedSampleTypes();
        // load WeightedSamples types
        loadWeightedSamplesTypes();

        return true;
    }

    bool bflTypekitPlugin::loadOperators()
    {
        return true;
    }
}

ORO_TYPEKIT_PLUGIN(BFL::bflTypekitPlugin)
