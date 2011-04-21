#include "bfl_typekit.hpp"

namespace BFL{
    using namespace std;
    using namespace RTT;
    using namespace RTT::detail;
/*****************************************************************************
 * WEIGHTED SAMPLES
 * **************************************************************************/

    template<class T>
    struct WeightedSamplesTypeInfo: public SequenceTypeInfo<std::vector<WeightedSample<T> >,false > 
    {
         WeightedSamplesTypeInfo<T>(std::string name):SequenceTypeInfo< std::vector<WeightedSample<T> >,false > (name)
        {
        };
    };
    

  void loadWeightedSamplesTypes(){

        RTT::types::Types()->addType( new WeightedSamplesTypeInfo<int>("WeightedSampleInts") );
        RTT::types::Types()->addType( new WeightedSamplesTypeInfo<double> ("WeightedSampleDoubles") );
        RTT::types::Types()->addType( new WeightedSamplesTypeInfo<ColumnVector> ("WeightedSampleColumnVectors") );
  }
}
