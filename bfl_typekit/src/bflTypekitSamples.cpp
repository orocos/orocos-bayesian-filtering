#include "bfl_typekit.hpp"

namespace BFL{
    using namespace RTT;
    using namespace RTT::detail;
    using namespace MatrixWrapper;
    class PropertyIntrospection;
/*****************************************************************************
 * SAMPLES
 * **************************************************************************/

    template<class T>
    struct SamplesTypeInfo: public SequenceTypeInfo<std::vector<Sample<T> >,false > 
    {
         SamplesTypeInfo<T>(std::string name):SequenceTypeInfo< std::vector<Sample<T> >,false > (name)
        {
        };
    };

    void loadSamplesTypes(){
          RTT::types::Types()->addType( new SamplesTypeInfo<int >("SampleInts") );
          RTT::types::Types()->addType( new SamplesTypeInfo<double >("SampleDoubles") );
          RTT::types::Types()->addType( new SamplesTypeInfo<ColumnVector >("SampleColumnVectors") );
    };
}
