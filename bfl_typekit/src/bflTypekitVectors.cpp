#include "bfl_typekit.hpp"

namespace BFL{
    using namespace std;
    using namespace RTT;
    using namespace RTT::detail;

    struct ColumnVectorsTypeInfo: public SequenceTypeInfo<std::vector<ColumnVector>,false > 
    {
         ColumnVectorsTypeInfo():SequenceTypeInfo< std::vector<ColumnVector>,false > ("ColumnVectors")
        {
        }
    };

    struct RowVectorsTypeInfo: public SequenceTypeInfo<std::vector<RowVector>,false > 
    {
         RowVectorsTypeInfo():SequenceTypeInfo< std::vector<RowVector>,false > ("RowVectors")
        {
        }
    };

    void loadVectorsTypes(){
          RTT::types::Types()->addType( new ColumnVectorsTypeInfo() );
          RTT::types::Types()->addType( new RowVectorsTypeInfo() );
    };
}
