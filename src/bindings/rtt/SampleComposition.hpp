/***************************************************************************
  tag: Tinne De Laet 2009  SampleComposition.hpp

 ***************************************************************************
 *   This library is free software; you can redistribute it and/or         *
 *   modify it under the terms of the GNU General Public                   *
 *   License as published by the Free Software Foundation;                 *
 *   version 2 of the License.                                             *
 *                                                                         *
 *   As a special exception, you may use this file as part of a free       *
 *   software library without restriction.  Specifically, if other files   *
 *   instantiate templates or use macros or inline functions from this     *
 *   file, or you compile this file and link it with other files to        *
 *   produce an executable, this file does not by itself cause the         *
 *   resulting executable to be covered by the GNU General Public          *
 *   License.  This exception does not however invalidate any other        *
 *   reasons why the executable file might be covered by the GNU General   *
 *   Public License.                                                       *
 *                                                                         *
 *   This library is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU     *
 *   Lesser General Public License for more details.                       *
 *                                                                         *
 *   You should have received a copy of the GNU General Public             *
 *   License along with this library; if not, write to the Free Software   *
 *   Foundation, Inc., 59 Temple Place,                                    *
 *   Suite 330, Boston, MA  02111-1307  USA                                *
 *                                                                         *
 ***************************************************************************/

#ifndef SAMPLE_COMPOSITION_HPP
#define SAMPLE_COMPOSITION_HPP

#include <rtt/Property.hpp>
#include <rtt/PropertyBag.hpp>
#include <rtt/TemplateTypeInfo.hpp>
#include <rtt/Types.hpp>
#include <rtt/Logger.hpp>
#include <rtt/DataSources.hpp>
#include <ostream>
#include <sstream>
#include <vector>

#include "../../sample/sample.h"
#include "../../sample/weightedsample.h"
namespace BFL
{
    using namespace RTT;
    using namespace RTT::detail;
    using namespace MatrixWrapper;
    class PropertyIntrospection;
/*****************************************************************************
 * SAMPLE
 * **************************************************************************/

    /**
     * A decomposePropertyBag method for decomposing a sample<T>
     * into a PropertyBag with Property<T>'s.
     */
    template<class T>
    void decomposeProperty(const Sample<T>& sample, PropertyBag& targetbag)
    {
        std::string tname = detail::DataSourceTypeInfo<T>::getType();
        targetbag.setType("Sample");
        //std::string str;

        assert( targetbag.empty() );

        bool result =true;
        //std::stringstream out;
        //out << i+1;
        //str = out.str();

        Property<PropertyBag>* el_bag = new Property<PropertyBag>("SampleValue", "Sample Value"); 
        Property<T> el("SampleValue" , "Sample value ",sample.ValueGet())  ;
        if(    el.getTypeInfo()->decomposeType(el.getDataSource(),el_bag->value()) )
        {
            //log(Debug)<<"Element type "<<el.getType()<<" is a bag"<<endlog();
            targetbag.add( el_bag ); // Put variables in the bag
        }
        else
        {
            //log(Debug)<<"Element type "<<el.getType()<<" is not a bag"<<endlog();
            //For Property
            targetbag.add( new Property<T>("SampleValue" ,"Sample Value",sample.ValueGet() )); // Put variables in the bag
        }
    };

    /**
     * A composeProperty method for composing a property of a vector<T>
     * The dimension of the vector must be less than 100.
     */
    template<class T>
    bool composeProperty(const PropertyBag& bag, Sample<T>& sample)
    {
        //log(Debug) << "composeProperty of sample " << endlog();
        std::string tname = detail::DataSourceTypeInfo<T>::getType();
        if ( bag.getType() == std::string("Sample") ) {
            // Get values
            Property<PropertyBag>* el_bag =  bag.getProperty<PropertyBag>("SampleValue");

            if(el_bag==NULL){
                // Works for properties in Sample
                PropertyBase* element = bag.getItem( 0 );
                //log(Debug)<<element->getName()<<", "<< element->getDescription()<<endlog();
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
                if(!(el_p.getDataSource()->composeType(el_bag->getDataSource()))){
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

        bool decomposeTypeImpl(const Sample<T>& sample, PropertyBag& targetbag) const
        {
            decomposeProperty<T>( sample, targetbag );
            return true;
        };

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

/*****************************************************************************
 * WEIGHTED SAMPLE
 * **************************************************************************/
    /**
     * A decomposePropertyBag method for decomposing a sample<T>
     * into a PropertyBag with Property<T>'s.
     */
    template<class T>
    void decomposeProperty(const WeightedSample<T>& weightedSample, PropertyBag& targetbag)
    {
        std::string tname = detail::DataSourceTypeInfo<T>::getType();
        targetbag.setType("WeightedSample");
        //std::string str;

        assert( targetbag.empty() );

        bool result =true;
        //std::stringstream out;
        //out << i+1;
        //str = out.str();

        // add value
        Property<PropertyBag>* el_bag = new Property<PropertyBag>("WeightedSampleValue", "WeightedSample Value"); 
        Property<T> el("WeightedSampleValue" , "WeightedSample value ",weightedSample.ValueGet())  ;
        if(    el.getTypeInfo()->decomposeType(el.getDataSource(),el_bag->value()) )
        {
            //log(Debug)<<"Element type "<<el.getType()<<" is a bag"<<endlog();
            targetbag.add( el_bag ); // Put variables in the bag
        }
        else
        {
            //log(Debug)<<"Element type "<<el.getType()<<" is not a bag"<<endlog();
            //For Property
            targetbag.add( new Property<T>("WeightedSampleValue" ,"WeightedSample Value",weightedSample.ValueGet() )); // Put variables in the bag
        }
        // add weight
        targetbag.add( new Property<double>("WeightedSampleWeight" ,"WeightedSample Weight",weightedSample.WeightGet() )); // Put variables in the bag
    };

    /**
     * A composeProperty method for composing a property of a vector<T>
     * The dimension of the vector must be less than 100.
     */
    template<class T>
    bool composeProperty(const PropertyBag& bag, WeightedSample<T>& weightedSample)
    {
        //log(Debug) << "composeProperty of WeightedSample " << endlog();
        std::string tname = detail::DataSourceTypeInfo<T>::getType();
        if ( bag.getType() == std::string("WeightedSample") ) {
            // Get values of sample
            Property<PropertyBag>* el_bag =  bag.getProperty<PropertyBag>("WeightedSampleValue");

            if(el_bag==NULL){
                // Works for properties in WeightedSample
                PropertyBase* element = bag.getItem( 0 );
                //log(Debug)<<element->getName()<<", "<< element->getDescription()<<endlog();
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
                if(!(el_p.getDataSource()->composeType(el_bag->getDataSource()))){
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
            Property<double>* weightProp =  bag.getProperty<double>("WeightedSampleWeight");

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
        : public TemplateTypeInfo<WeightedSample<T>, true>
    {
        WeightedSampleTypeInfo<T>(std::string name)
            : TemplateTypeInfo<WeightedSample<T>, true >(name)
        {
        };

        bool decomposeTypeImpl(const WeightedSample<T>& weightedSample, PropertyBag& targetbag) const
        {
            decomposeProperty<T>( weightedSample, targetbag );
            return true;
        };

        bool composeTypeImpl(const PropertyBag& bag, WeightedSample<T>& result) const
        {
            return composeProperty<T>( bag, result );
        }

    };

    template<typename T>
    struct WeightedSample_ctor
        : public std::binary_function<T ,double , const WeightedSample<T>&>
    {
        typedef const WeightedSample<T>& (Signature)( T, double );
        mutable boost::shared_ptr< WeightedSample<T> > ptr;
        WeightedSample_ctor()
            : ptr( new WeightedSample<T>() ) {}
        const WeightedSample<T>& operator()( T value , double weight) const
        {
            ptr->ValueSet( value );
            ptr->WeightSet( weight );
            return *(ptr);
        }
    };
    
    /**
     * See NArityDataSource which requires a function object like
     * this one.
     */
/*
    template<typename T>
    struct Sample_varargs_ctor
    {
        typedef const Sample<T>& result_type;
        typedef T argument_type;
        result_type operator()( const Sample<T>& args ) const
        {
            return args;
        }
    };
*/
    

};
#endif

