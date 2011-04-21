#include "bfl_typekit.hpp"

namespace BFL{
    using namespace std;
    using namespace RTT;
    using namespace RTT::detail;
    struct IntVectorTypeInfo: public SequenceTypeInfo<std::vector<int>,false > 
    {
         IntVectorTypeInfo():SequenceTypeInfo< std::vector<int>,false > ("ints")
        {
        }
    };

    struct VectorTypeInfo : public SequenceTypeInfo<ColumnVector,true>
    {
        VectorTypeInfo():SequenceTypeInfo<ColumnVector, true >("ColumnVector")
        {
        }
    };

    struct RVectorTypeInfo : public SequenceTypeInfo<RowVector,true>
    {
        RVectorTypeInfo():SequenceTypeInfo<RowVector, true >("RowVector")
        {
        }
    };

    void loadVectorTypes(){

          RTT::types::Types()->addType( new IntVectorTypeInfo() );
          RTT::types::Types()->addType( new VectorTypeInfo() );
          RTT::types::Types()->addType( new RVectorTypeInfo() );
    };
}
