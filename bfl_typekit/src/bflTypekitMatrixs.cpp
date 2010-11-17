#include "bfl_typekit.hpp"

namespace BFL{
    using namespace std;
    using namespace RTT;
    using namespace RTT::detail;
    struct MatrixsTypeInfo: public SequenceTypeInfo<std::vector<Matrix>,false > 
    {
         MatrixsTypeInfo():SequenceTypeInfo< std::vector<Matrix>,false > ("Matrixs")
        {
        }
    };

    struct SymmetricMatrixsTypeInfo: public SequenceTypeInfo<std::vector<SymmetricMatrix>,false > 
    {
         SymmetricMatrixsTypeInfo():SequenceTypeInfo< std::vector<SymmetricMatrix>,false > ("SymmetricMatrixs")
        {
        }
    };
    void loadMatrixsTypes(){
          RTT::types::Types()->addType( new MatrixsTypeInfo() );
          RTT::types::Types()->addType( new SymmetricMatrixsTypeInfo() );
    };
}
