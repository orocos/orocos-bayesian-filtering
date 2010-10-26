/***************************************************************************
// Copyright  (C)  2010  Tinne De Laet <tinne dot delaet at mech dot kuleuven dot be>
// Author: Tinne De Laet 
// Maintainer: Tinne De Laet

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
#include <rtt/types/TemplateTypeInfo.hpp>
#include <rtt/Logger.hpp>
#include <ostream>
#include <sstream>
#include <vector>

#include <bfl/sample/sample.h>
#include <bfl/sample/weightedsample.h>
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

    template<class T>
    struct SamplesTypeInfo: public SequenceTypeInfo<std::vector<Sample<T> >,false > 
    {
         SamplesTypeInfo<T>(std::string name):SequenceTypeInfo< std::vector<Sample<T> >,false > (name)
        {
        };
    };

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
        : public TemplateTypeInfo<WeightedSample<T>, true>
    {
        WeightedSampleTypeInfo<T>(std::string name)
            : TemplateTypeInfo<WeightedSample<T>, true >(name)
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

    template<class T>
    struct WeightedSamplesTypeInfo: public SequenceTypeInfo<std::vector<WeightedSample<T> >,false > 
    {
         WeightedSamplesTypeInfo<T>(std::string name):SequenceTypeInfo< std::vector<WeightedSample<T> >,false > (name)
        {
        };
    };
    

};
#endif

