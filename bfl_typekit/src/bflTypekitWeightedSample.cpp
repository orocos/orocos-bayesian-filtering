#include "bfl_typekit.hpp"

namespace BFL{
    using namespace std;
    using namespace RTT;
    using namespace RTT::detail;
/*****************************************************************************
 * WEIGHTED SAMPLE
 * **************************************************************************/

    /**
     * A composeProperty method for composing a property of a vector<T>
     * The dimension of the vector must be less than 100.
     */
    template<class T>
    bool composeProperty(const PropertyBag& bag, WeightedSample<T>& weightedSample)
    {
        std::string tname = internal::DataSourceTypeInfo<T>::getType();
        if ( bag.getType() == std::string("WeightedSample") ) {
            // Get values of sample
            Property<PropertyBag>* el_bag =  bag.getPropertyType<PropertyBag>("WeightedSampleValue");


            if(el_bag==NULL){
                // Works for properties in WeightedSample
                PropertyBase* element = bag.getItem( 0 );
                Property<T> my_property_t (element->getName(),element->getDescription());
                if(my_property_t.getType()!=element->getType())
                {
                    log(Error)<< "Type of "<< element->getName() << " does not match type of WeightedSample"<< "OR "<<"Could not read WeightedSample Value "<<endlog();
                    return false;
                }
                else{
                    my_property_t.getTypeInfo()->composeType(element->getDataSource(),my_property_t.getDataSource());
                    weightedSample.ValueSet( my_property_t.get());
                }
            }
            else{
                // Works for propertybags in WeightedSample
                const std::string el_bagType = el_bag->getType();
                Property<T > el_p(el_bag->getName(),el_bag->getDescription());
                //if(!(el_p.getDataSource()->composeType(el_bag->getDataSource()))){
                if(!(el_p.getTypeInfo()->composeType(el_bag->getDataSource(),el_p.getDataSource()))){
                    log(Error)<<"Could not compose WeightedSampleValue "<<endlog();
                    return false;
                }
                if(el_p.ready()){
                    weightedSample.ValueSet( el_p.get());
                }else{
                    log(Error)<<"Property of WeightedSampleValue was not ready for use"<<endlog();
                    return false;
                }
            }
            // Get weight of sample
            Property<double>* weightProp =  bag.getPropertyType<double>("WeightedSampleWeight");

            if(!weightProp)
            {
                log(Error)<< "Error reading weight of WeightedSample"<<endlog();
                return false;
            }
            else{
                weightedSample.WeightSet( weightProp->get());
            }
        }
        else {
            Logger::log() << Logger::Error << "Composing Property< WeightedSample<T> > :"
                          << " type mismatch, got type '"<< bag.getType()
                          << "', expected type "<<tname<<"."<<Logger::endl;
            return false;
        }
        return true;
    };

    template <typename T>
    struct WeightedSampleTypeInfo
        : public TemplateTypeInfo<WeightedSample<T>, false>
    {
        WeightedSampleTypeInfo<T>(std::string name)
            : TemplateTypeInfo<WeightedSample<T>, false >(name)
        {
        };

        bool composeTypeImpl(const PropertyBag& bag, WeightedSample<T>& result) const
        {
            return composeProperty<T>( bag, result );
        }

        std::vector<std::string> getMemberNames() const
        {
            std::vector<std::string> result;
            result.push_back("WeightedSampleValue");
            result.push_back("WeightedSampleWeight");
            return result;
        }           
        
        base::DataSourceBase::shared_ptr getMember(base::DataSourceBase::shared_ptr source, const std::string& name) const{
            typename internal::DataSource<WeightedSample<T> >::shared_ptr ds = internal::DataSource<WeightedSample<T> >::narrow( source.get() );
            if(name=="WeightedSampleValue"){
                return new internal::ValueDataSource<T>(ds->get().ValueGet());
            }
            if(name=="WeightedSampleWeight"){
                return new internal::ValueDataSource<double>(ds->get().WeightGet());
            }
            return base::DataSourceBase::shared_ptr();
        }

        bool composeTypeImpl(const PropertyBag& bag, Sample<T>& result) const
        {
            return composeProperty<T>( bag, result );
        }
    };

    WeightedSample<ColumnVector> createWeightedSampleColumnVector(int dimension)
    {
        return WeightedSample<ColumnVector>(dimension);
    };
    WeightedSample<ColumnVector> createWeightedSampleInt(int dimension)
    {
        return WeightedSample<ColumnVector>(dimension);
    };

  void loadWeightedSampleTypes(){

        RTT::types::Types()->addType( new WeightedSampleTypeInfo<int>("WeightedSampleInt") );
        RTT::types::Types()->addType( new WeightedSampleTypeInfo<double> ("WeightedSampleDouble") );
        RTT::types::Types()->addType( new WeightedSampleTypeInfo<ColumnVector> ("WeightedSampleColumnVector") );
  }
}
