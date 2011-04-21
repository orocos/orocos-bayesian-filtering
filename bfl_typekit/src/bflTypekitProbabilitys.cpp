#include "bfl_typekit.hpp"

namespace BFL{

    using namespace std;
    using namespace RTT;
    using namespace RTT::detail;

    struct ProbabilityVectorTypeInfo: public SequenceTypeInfo<std::vector<Probability>,false > 
    {
         ProbabilityVectorTypeInfo():SequenceTypeInfo< std::vector<Probability>,false > ("Probabilitys")
        {
        }
    };

    void loadProbabilitysTypes(){
        RTT::types::Types()->addType( new ProbabilityVectorTypeInfo() );
    };

}
