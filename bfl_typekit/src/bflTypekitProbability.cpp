#include "bfl_typekit.hpp"
#include <bfl/bfl_constants.h>

#ifndef __PROB_STREAM__
#define __PROB_STREAM__
namespace BFL
{
  using namespace std;
  ostream & operator<< (ostream & stream, Probability& prob)
    {
      stream << prob.getValue() << endl;
      return stream;
    }

  istream & operator>> (istream & stream, Probability& prob)
    {
        double value;
        stream >> value;
        prob = Probability(value);
        return stream; 
    }
};
#endif //__PROBSTREAM__

namespace BFL{

    using namespace std;
    using namespace RTT;
    using namespace RTT::detail;


    template<class Archive>
    void serialize(Archive & a, Probability & prob, unsigned int){
        using boost::serialization::make_nvp;   
        a & make_nvp("Probability", prob.getValue() );
    }

    struct ProbabilityTypeInfo : public StructTypeInfo<Probability,false>
    {
        ProbabilityTypeInfo():StructTypeInfo<Probability,false>("Probability")
        {
        };
    };

    void loadProbabilityTypes(){
        RTT::types::Types()->addType( new ProbabilityTypeInfo() );
    };

}
