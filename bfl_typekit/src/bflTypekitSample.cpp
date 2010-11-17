#include "bfl_typekit.hpp"

namespace BFL{
    using namespace RTT;
    using namespace RTT::detail;
    using namespace MatrixWrapper;
    class PropertyIntrospection;
/*****************************************************************************
 * SAMPLE
 * **************************************************************************/

    /**
     * A composeProperty method for composing a property of a Sample<T>
     */
    template<class T>
    bool composeProperty(const PropertyBag& bag, Sample<T>& sample)
    {
        std::string tname = internal::DataSourceTypeInfo<T>::getType();
        if ( bag.getType() == std::string("Sample") ) {
            // Get values
            Property<PropertyBag>* el_bag  =  bag.getPropertyType<PropertyBag>("SampleValue");

            if(el_bag==NULL){
                // Works for properties in Sample
                PropertyBase* element = bag.getItem( 0 );
                Property<T> my_property_t (element->getName(),element->getDescription());
                if(my_property_t.getType()!=element->getType())
                {
                    log(Error)<< "Type of "<< element->getName() << " does not match type of Sample"<< "OR "<<"Could not read Sample Value "<<endlog();
                    return false;
                }
                else{
                    my_property_t.getTypeInfo()->composeType(element->getDataSource(),my_property_t.getDataSource());
                    sample.ValueSet( my_property_t.get());
                }
            }
            else{
                // Works for propertybags in Sample
                const std::string el_bagType = el_bag->getType();
                Property<T > el_p(el_bag->getName(),el_bag->getDescription());
                if(!(el_p.getTypeInfo()->composeType(el_bag->getDataSource(),el_p.getDataSource()))){
                    log(Error)<<"Could not compose SampleValue "<<endlog();
                    return false;
                }
                if(el_p.ready()){
                    sample.ValueSet( el_p.get());
                }else{
                    log(Error)<<"Property of SampleValue was not ready for use"<<endlog();
                    return false;
                }
            }
        }
        else {
            Logger::log() << Logger::Error << "Composing Property< Sample<T> > :"
                          << " type mismatch, got type '"<< bag.getType()
                          << "', expected type "<<tname<<"."<<Logger::endl;
            return false;
        }
        return true;
    };

    template <typename T>
    struct SampleTypeInfo
        : public TemplateTypeInfo<Sample<T>, true>
    {
        SampleTypeInfo<T>(std::string name)
            : TemplateTypeInfo<Sample<T>, true >(name)
        {
        };

        std::vector<std::string> getMemberNames() const
        {
            std::vector<std::string> result;
            result.push_back("SampleValue");
            return result;
        }           
        
        base::DataSourceBase::shared_ptr getMember(base::DataSourceBase::shared_ptr source, const std::string& name) const{
            typename internal::DataSource<Sample<T> >::shared_ptr ds = internal::DataSource<Sample<T> >::narrow( source.get() );
            if(name=="SampleValue"){
                return new internal::ValueDataSource<T>(ds->get().ValueGet());
            }
            return base::DataSourceBase::shared_ptr();
        }

        bool composeTypeImpl(const PropertyBag& bag, Sample<T>& result) const
        {
            return composeProperty<T>( bag, result );
        }


    };

    template<typename T>
    struct Sample_ctor
        : public std::unary_function<T, const Sample<T>&>
    {
        typedef const Sample<T>& (Signature)( T );
        mutable boost::shared_ptr< Sample<T> > ptr;
        Sample_ctor()
            : ptr( new Sample<T>() ) {}
        const Sample<T>& operator()( T value ) const
        {
            ptr->ValueSet( value );
            return *(ptr);
        }
    };

    Sample<ColumnVector> createSampleColumnVector(int dimension)
    {
        return Sample<ColumnVector>(dimension);
    }

    void loadSampleTypes(){

          RTT::types::Types()->addType( new SampleTypeInfo<int >("SampleInt") );
          RTT::types::Types()->addType( new SampleTypeInfo<double >("SampleDouble") );
          RTT::types::Types()->addType( new SampleTypeInfo<ColumnVector >("SampleColumnVector") );
    }
}
