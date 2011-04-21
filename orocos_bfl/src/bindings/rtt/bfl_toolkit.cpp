// Copyright  (C)  2009  Tinne De Laet <tinne dot delaet at mech dot kuleuven dot be>

// Author: Tinne De Laet 
// Maintainer: Tinne De Laet


// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.

// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.

// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

#include "bfl_toolkit.hpp"
#include <rtt/Property.hpp>
#include <rtt/PropertyBag.hpp>
#include <rtt/TemplateTypeInfo.hpp>
#include <rtt/Operators.hpp>
#include <rtt/OperatorTypes.hpp>
#include <rtt/Types.hpp>
#include <rtt/Logger.hpp>
#include <rtt/DataSources.hpp>
#include <rtt/mystd.hpp>
#include <rtt/os/StartStopManager.hpp>
#include <rtt/Toolkit.hpp>

#include "../../wrappers/matrix/vector_wrapper.h"
#include "../../wrappers/matrix/matrix_wrapper.h"
#include "../../bfl_constants.h"
#include "../../sample/sample.h"

#include <rtt/VectorTemplateComposition.hpp>
#include "SampleComposition.hpp"


/*
std::ostream& operator<<(std::ostream& os, const Probability p) {
 return os << (double)p;
}

std::istream& operator>>(std::istream& is, Probability p) {
 char c;
 double p;
 return os >> p ;
}
*/



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


using namespace MatrixWrapper;
namespace BFL{ 

    using namespace RTT;
    using namespace RTT::detail;

    bflToolkitPlugin bflToolkit;

    struct ProbabilityTypeInfo : public TemplateTypeInfo<Probability,true>
    {
        ProbabilityTypeInfo():TemplateTypeInfo<Probability,true>("Probability")
        {
        };
        virtual bool decomposeTypeImpl(const Probability& prob, PropertyBag& targetbag) const
        {
                targetbag.setType("Probability");
                double probDouble = prob.getValue();
                targetbag.add( new Property<double>("Probability","Probability value", probDouble ) ); // Put variables in the bag
                return true;
        };
        virtual bool composeTypeImpl(const PropertyBag& bag, Probability& prob) const
        {
                Property<double>* probability = bag.getProperty<double>("Probability");
                if (!probability)
                    return false;
                prob = (probability)->get();
                return true;
        };
    };


    
    struct VectorTypeInfo : public TemplateContainerTypeInfo<ColumnVector,int,double,ArrayIndexChecker<ColumnVector >, VectorAssignChecker<ColumnVector >,true>
    {
        VectorTypeInfo():TemplateContainerTypeInfo<ColumnVector,int,double,ArrayIndexChecker<ColumnVector>, VectorAssignChecker<ColumnVector >, true >("ColumnVector")
        {
        };
        
        bool decomposeTypeImpl(const ColumnVector& vec, PropertyBag& targetbag) const
        {
            targetbag.setType("ColumnVector");
            int dimension = vec.size();
            std::string str;
            
            for ( int i=1; i <= dimension ; i++){
                std::stringstream out;
                out << i;
                str = out.str();
                targetbag.add( new Property<double>(str, str +"th element of vector",vec(i)) ); // Put variables in the bag
            }

            return true;
        };

        bool composeTypeImpl(const PropertyBag& bag, ColumnVector& result) const{

            if ( bag.getType() == "ColumnVector" ) {
                int dimension = bag.size();
                result.resize( dimension );
                
                // Get values
                for (int i = 1; i <= dimension ; i++) {
                    std::stringstream out;
                    out << i;
                    Property<double>* elem = bag.getProperty<double>(out.str());
                    if(elem->ready())
                        result(i) = elem->get();
                    else{
                        log(Error)<<"Could not read element "<<i<<endlog();
                        return false;
                    }
                }
            }else{
                log(Error) << "Composing Property< ColumnVector > :"
                           << " type mismatch, got type '"<< bag.getType()
                           << "', expected type "<<"ColumnVector."<<endlog();
                return false;
            }
            return true;
        };
    };

    struct RVectorTypeInfo : public TemplateContainerTypeInfo<RowVector,int,double,ArrayIndexChecker<RowVector >, VectorAssignChecker<RowVector >,true>
    {
        RVectorTypeInfo():TemplateContainerTypeInfo<RowVector,int,double,ArrayIndexChecker<RowVector>, VectorAssignChecker<RowVector >, true >("RowVector")
        {
        };

        bool decomposeTypeImpl(const RowVector& vec, PropertyBag& targetbag) const
        {
            targetbag.setType("RowVector");
            int dimension = vec.size();
            std::string str;


            for ( int i=1; i <= dimension ; i++){
                std::stringstream out;
                out << i;
                str = out.str();
                targetbag.add( new Property<double>(str, str +"th element of vector",vec(i)) ); // Put variables in the bag
            }

            return true;
        };

        bool composeTypeImpl(const PropertyBag& bag, RowVector& result) const{

            if ( bag.getType() == "RowVector" ) {
                int dimension = bag.size();
                result.resize( dimension );

                // Get values
                for (int i = 1; i <= dimension ; i++) {
                    std::stringstream out;
                    out << i;
                    Property<double>* elem = bag.getProperty<double>(out.str());
                    if(elem->ready())
                        result(i) = elem->get();
                    else{
                        log(Error)<<"Could not read element "<<i<<endlog();
                        return false;
                    }
                }
            }else{
                log(Error) << "Composing Property< RowVector > :"
                           << " type mismatch, got type '"<< bag.getType()
                           << "', expected type "<<"RowVector."<<endlog();
                return false;
            }
            return true;
        };
    };

    struct MatrixTypeInfo : public TemplateContainerTypeInfo<Matrix,int,RowVector,MatrixIndexChecker<Matrix>, MatrixAssignChecker<Matrix>, true>
    {
        MatrixTypeInfo():TemplateContainerTypeInfo<Matrix,int,RowVector,MatrixIndexChecker<Matrix>, MatrixAssignChecker<Matrix>,true>("Matrix"){
        };

        bool decomposeTypeImpl(const Matrix& mat, PropertyBag& targetbag) const{
            targetbag.setType("Matrix");
            unsigned int dimension = mat.rows();

            for ( unsigned int i=1; i <= dimension ; i++){
                std::stringstream out;
                out << i;
                Property<PropertyBag>* row_bag = new Property<PropertyBag>(out.str(), out.str() +"th row of matrix"); 
                Property<RowVector> row(out.str(), out.str() +"th row of matrix",mat.rowCopy(i))  ;
                row.getTypeInfo()->decomposeType(row.getDataSource(),row_bag->value());
                targetbag.add( row_bag ); // Put variables in the bag
            }

            return true;
        };

        bool composeTypeImpl(const PropertyBag& bag, Matrix& result) const{
            if ( bag.getType() == "Matrix" ) {
                unsigned int rows = bag.size();
                unsigned int cols = 0;
                // Get values
                for (unsigned int i = 1; i <= rows ; i++) {
                    std::stringstream out;
                    out << i;
                    Property<PropertyBag>* row_bag =  bag.getProperty<PropertyBag>(out.str());
                    if(row_bag==NULL){
                        log(Error)<<"Could not read row "<<i<<endlog();
                        return false;
                    }
                    Property<RowVector> row_p(row_bag->getName(),row_bag->getDescription());
                    if(!(row_p.getDataSource()->composeType(row_bag->getDataSource()))){
                        log(Error)<<"Could not decompose row "<<i<<endlog();
                        return false;
                    }
                    if(row_p.ready()){
                        if(i==1){
                            cols = row_p.get().size();
                            result.resize(rows,cols);
                        } else
                            if(row_p.get().size()!=cols){
                                log(Error)<<"Row "<<i+1<<" size does not match matrix columns"<<endlog();
                                return false;
                            }
                        for ( unsigned int j=1; j <= row_p.get().size() ; j++){
                            result(i,j)=row_p.get()(j);
                        }
                    }else{
                        log(Error)<<"Property of Row "<<i<<"was not ready for use"<<endlog();
                        return false;
                    }
                }
            }else {
                log(Error) << "Composing Property< Matrix > :"
                           << " type mismatch, got type '"<< bag.getType()
                           << "', expected type "<<"Matrix."<<endlog();
                return false;
            }
            return true;
        };
    };

    struct SymmetricMatrixTypeInfo : public TemplateContainerTypeInfo<SymmetricMatrix,int,RowVector,MatrixIndexChecker<SymmetricMatrix>, MatrixAssignChecker<SymmetricMatrix>, true>
    {
        SymmetricMatrixTypeInfo():TemplateContainerTypeInfo<SymmetricMatrix,int,RowVector,MatrixIndexChecker<SymmetricMatrix>, MatrixAssignChecker<SymmetricMatrix>,true>("SymmetricMatrix"){
        };

        bool decomposeTypeImpl(const SymmetricMatrix& mat, PropertyBag& targetbag) const{
            targetbag.setType("SymmetricMatrix");
            unsigned int dimension = mat.rows();

            for ( unsigned int i=1; i <= dimension ; i++){
                std::stringstream out;
                out << i;
                //targetbag.add( new Property<RowVector >(out.str(), out.str() +"th row of matrix",((Matrix)mat).rowCopy(i) )); // Put variables in the bag
                Property<PropertyBag>* row_bag = new Property<PropertyBag>(out.str(), out.str() +"th row of matrix"); 
                Property<RowVector> row(out.str(), out.str() +"th row of matrix",((Matrix)mat).rowCopy(i))  ;
                row.getTypeInfo()->decomposeType(row.getDataSource(),row_bag->value());
                targetbag.add( row_bag ); // Put variables in the bag
            }

            return true;
        };

        bool composeTypeImpl(const PropertyBag& bag, SymmetricMatrix& result) const{
            Matrix matrix;
            if ( bag.getType() == "SymmetricMatrix" ) {
                unsigned int rows = bag.size();
                unsigned int cols = 0;
                // Get values
                for (unsigned int i = 1; i <= rows ; i++) {
                    std::stringstream out;
                    out << i;
                    Property<PropertyBag>* row_bag =  bag.getProperty<PropertyBag>(out.str());
                    if(row_bag==NULL){
                        log(Error)<<"Could not read row "<<i<<endlog();
                        return false;
                    }
                    Property<RowVector > row_p(row_bag->getName(),row_bag->getDescription());
                    if(!(row_p.getDataSource()->composeType(row_bag->getDataSource()))){
                        log(Error)<<"Could not decompose row "<<i<<endlog();
                        return false;
                    }
                    if(row_p.ready()){
                        if(i==1){
                            cols = row_p.get().size();
                            matrix.resize(rows,cols);
                        } else
                            if(row_p.get().size()!=cols){
                                log(Error)<<"Row "<<i+1<<" size does not match matrix columns"<<endlog();
                                return false;
                            }
                        for ( unsigned int j=1; j <= row_p.get().size() ; j++){
                            matrix(i,j)=row_p.get()(j);
                        }
                    }else{
                        log(Error)<<"Property of Row "<<i<<"was not ready for use"<<endlog();
                        return false;
                    }
                }
            }else {
                log(Error) << "Composing Property< SymmetricMatrix > :"
                           << " type mismatch, got type '"<< bag.getType()
                           << "', expected type "<<"SymmetricMatrix."<<endlog();
                return false;
            }
            matrix.convertToSymmetricMatrix(result);
            return true;
        };
    };

    struct vector_index
        : public std::binary_function<const ColumnVector&, int, double>
    {
        double operator()(const ColumnVector& v, int index) const
        {
            if ( index > (int)(v.size()) || index < 1)
                return NAN;
            return v(index);
        }
    };

    struct rvector_index
        : public std::binary_function<const RowVector&, int, double>
    {
        double operator()(const RowVector& v, int index) const
        {
            if ( index > (int)(v.size()) || index < 1)
                return NAN;
            return v(index);
        }
    };

    struct get_size
        : public std::unary_function<const ColumnVector&, int>
    {
        int operator()(const ColumnVector& cont ) const
        {
            return cont.size();
        }
    };

    struct rget_size
        : public std::unary_function<const RowVector&, int>
    {
        int operator()(const RowVector& cont ) const
        {
            return cont.size();
        }
    };

    struct vector_index_constructor
        : public std::unary_function<int,const ColumnVector&>
    {
        typedef const ColumnVector& (Signature)( int );
        mutable boost::shared_ptr< ColumnVector > ptr;
        vector_index_constructor() :
            ptr( new ColumnVector() ){}
        const ColumnVector& operator()(int size ) const
        {
            ptr->resize(size);
            return *(ptr);
        }
    };

    struct rvector_index_constructor
        : public std::unary_function<int,const RowVector&>
    {
        typedef const RowVector& (Signature)( int );
        mutable boost::shared_ptr< RowVector > ptr;
        rvector_index_constructor() :
            ptr( new RowVector() ){}
        const RowVector& operator()(int size ) const
        {
            ptr->resize(size);
            return *(ptr);
        }
    };

    struct matrix_index
        : public std::binary_function<const Matrix&, int, const RowVector&>
    {
        const RowVector& operator()(const Matrix& m, int index) const
        {
            if ( index > (int)(m.rows()) || index < 1)
                {
                    log(Error) << "index error" << endlog();
                    return RowVector(0);
                }
            return m.rowCopy(index);
        }
    };

    //struct matrix_index
    //    : public std::ternary_function<const Matrix&, int, int, double>
    //{
    //    double operator()(const Matrix& m, int i, int j) const{
    //        if ( i > (int)(m.rows()) || i < 1 || j<1 || j> (int)(m.columns()))
    //            return NAN;
    //        return m(i,j);
    //    }
    //};

    struct matrix_i_j_constructor
        : public std::binary_function<int,int,const Matrix&>
    {
        typedef const Matrix& (Signature)( int, int );
        mutable boost::shared_ptr< Matrix > ptr;
        matrix_i_j_constructor() :
            ptr( new Matrix() ){}
        const Matrix& operator()(int size1,int size2) const
        {
            ptr->resize(size1,size2);
            return *(ptr);
        }
    };

    struct symmetricMatrix_index_constructor
        : public std::unary_function<int,const SymmetricMatrix&>
    {
        typedef const SymmetricMatrix& (Signature)( int );
        mutable boost::shared_ptr< SymmetricMatrix > ptr;
        symmetricMatrix_index_constructor() :
            ptr( new SymmetricMatrix() ){}
        const SymmetricMatrix& operator()(int size) const
        {
            ptr->resize(size);
            return *(ptr);
        }
    };

    struct Probability_ctor
        : public std::unary_function<double, const Probability&>
    {
        typedef const Probability& (Signature)( double );
        mutable boost::shared_ptr< Probability > ptr;
        Probability_ctor()
            : ptr( new Probability() ) {}
        const Probability& operator()( double value ) const
        {
            //ptr = new Probability(value);
            //return *(ptr);
            return *(new Probability(value));
        }
    };

    std::string bflToolkitPlugin::getName()
    {
        return "bfl_toolkit";
    }

    bool bflToolkitPlugin::loadTypes()
    {
        RTT::TypeInfoRepository::Instance()->addType( new VectorTypeInfo() );
        RTT::TypeInfoRepository::Instance()->addType( new RVectorTypeInfo() );
        RTT::TypeInfoRepository::Instance()->addType( new MatrixTypeInfo() );
        RTT::TypeInfoRepository::Instance()->addType( new SymmetricMatrixTypeInfo() );
        RTT::TypeInfoRepository::Instance()->addType( new SampleTypeInfo<int>("SampleInt") );
        RTT::TypeInfoRepository::Instance()->addType( new SampleTypeInfo<double>("SampleDouble") );
        RTT::TypeInfoRepository::Instance()->addType( new SampleTypeInfo<ColumnVector>("SampleColumnVector") );
        RTT::TypeInfoRepository::Instance()->addType( new WeightedSampleTypeInfo<int>("WeightedSampleInt") );
        RTT::TypeInfoRepository::Instance()->addType( new WeightedSampleTypeInfo<double>("WeightedSampleDouble") );
        RTT::TypeInfoRepository::Instance()->addType( new WeightedSampleTypeInfo<ColumnVector>("WeightedSampleColumnVector") );
        RTT::TypeInfoRepository::Instance()->addType( new ProbabilityTypeInfo() );

        RTT::TypeInfoRepository::Instance()->addType( new StdVectorTemplateTypeInfo<ColumnVector,true>("ColumnVectors") );
        RTT::TypeInfoRepository::Instance()->addType( new StdVectorTemplateTypeInfo<RowVector,true>("RowVectors") );
        RTT::TypeInfoRepository::Instance()->addType( new StdVectorTemplateTypeInfo<Matrix,true>("Matrixs") );
        RTT::TypeInfoRepository::Instance()->addType( new StdVectorTemplateTypeInfo<SymmetricMatrix,true>("SymmetricMatrixs") );
        RTT::TypeInfoRepository::Instance()->addType( new StdVectorTemplateTypeInfo<Sample<int>,false>("SampleInts") );
        RTT::TypeInfoRepository::Instance()->addType( new StdVectorTemplateTypeInfo<Sample<double>,false>("SampleDoubles") );
        RTT::TypeInfoRepository::Instance()->addType( new StdVectorTemplateTypeInfo<Sample<ColumnVector>,false>("SampleColumnVectors") );
        RTT::TypeInfoRepository::Instance()->addType( new StdVectorTemplateTypeInfo<WeightedSample<int>,false>("WeightedSampleInts") );
        RTT::TypeInfoRepository::Instance()->addType( new StdVectorTemplateTypeInfo<WeightedSample<double>,false>("WeightedSampleDoubles") );
        RTT::TypeInfoRepository::Instance()->addType( new StdVectorTemplateTypeInfo<WeightedSample<ColumnVector>,false>("WeightedSampleColumnVectors") );
        RTT::TypeInfoRepository::Instance()->addType( new StdVectorTemplateTypeInfo<Probability,false>("Probabilitys") );

        RTT::TypeInfoRepository::Instance()->addType( new StdVectorTemplateTypeInfo<vector<Sample<int> >,false>("VecSampleInts") );
        RTT::TypeInfoRepository::Instance()->addType( new StdVectorTemplateTypeInfo<vector<Sample<double> >,false>("VecSampleDoubles") );
        RTT::TypeInfoRepository::Instance()->addType( new StdVectorTemplateTypeInfo<vector<Sample<ColumnVector> >,false>("VecSampleColumnVectors") );
        RTT::TypeInfoRepository::Instance()->addType( new StdVectorTemplateTypeInfo<vector<WeightedSample<int> >,false>("VecWeightedSampleInts") );
        RTT::TypeInfoRepository::Instance()->addType( new StdVectorTemplateTypeInfo<vector<WeightedSample<double> >,false>("VecWeightedSampleDoubles") );
        RTT::TypeInfoRepository::Instance()->addType( new StdVectorTemplateTypeInfo<vector<WeightedSample<ColumnVector> >,false>("VecWeightedSampleColumnVectors") );
        return true;
    }

    bool bflToolkitPlugin::loadConstructors()
    {
        RTT::TypeInfoRepository::Instance()->type("ColumnVector")->addConstructor(newConstructor(vector_index_constructor()));
        RTT::TypeInfoRepository::Instance()->type("RowVector")->addConstructor(newConstructor(rvector_index_constructor()));
        RTT::TypeInfoRepository::Instance()->type("Matrix")->addConstructor(newConstructor(matrix_i_j_constructor()));
        RTT::TypeInfoRepository::Instance()->type("SymmetricMatrix")->addConstructor(newConstructor(symmetricMatrix_index_constructor()));
        //RTT::TypeInfoRepository::Instance()->type("ColumnVector")->addConstructor(newConstructor(vector_index_value_constructor()));
        //
        RTT::TypeInfoRepository::Instance()->type("ColumnVectors")->addConstructor(newConstructor(stdvector_ctor<ColumnVector>() ) );
        RTT::TypeInfoRepository::Instance()->type("ColumnVectors")->addConstructor(newConstructor(stdvector_ctor2<ColumnVector>() ) );
        RTT::TypeInfoRepository::Instance()->type("ColumnVectors")->addConstructor(new StdVectorBuilder<ColumnVector>() );

        RTT::TypeInfoRepository::Instance()->type("RowVectors")->addConstructor(newConstructor(stdvector_ctor<RowVector>() ) );
        RTT::TypeInfoRepository::Instance()->type("RowVectors")->addConstructor(newConstructor(stdvector_ctor2<RowVector>() ) );
        RTT::TypeInfoRepository::Instance()->type("RowVectors")->addConstructor(new StdVectorBuilder<RowVector>() );

        RTT::TypeInfoRepository::Instance()->type("Matrixs")->addConstructor(newConstructor(stdvector_ctor<Matrix>() ) );
        RTT::TypeInfoRepository::Instance()->type("Matrixs")->addConstructor(newConstructor(stdvector_ctor2<Matrix>() ) );
        RTT::TypeInfoRepository::Instance()->type("Matrixs")->addConstructor(new StdVectorBuilder<Matrix>() );

        RTT::TypeInfoRepository::Instance()->type("SymmetricMatrixs")->addConstructor(newConstructor(stdvector_ctor<SymmetricMatrix>() ) );
        RTT::TypeInfoRepository::Instance()->type("SymmetricMatrixs")->addConstructor(newConstructor(stdvector_ctor2<SymmetricMatrix>() ) );
        RTT::TypeInfoRepository::Instance()->type("SymmetricMatrixs")->addConstructor(new StdVectorBuilder<SymmetricMatrix>() );

        RTT::TypeInfoRepository::Instance()->type("SampleInts")->addConstructor(newConstructor(stdvector_ctor<Sample<int> >() ) );
        RTT::TypeInfoRepository::Instance()->type("SampleInts")->addConstructor(newConstructor(stdvector_ctor2<Sample<int> >() ) );
        RTT::TypeInfoRepository::Instance()->type("SampleInts")->addConstructor(new StdVectorBuilder<Sample<int> >() );

        RTT::TypeInfoRepository::Instance()->type("SampleDoubles")->addConstructor(newConstructor(stdvector_ctor<Sample<double> >() ) );
        RTT::TypeInfoRepository::Instance()->type("SampleDoubles")->addConstructor(newConstructor(stdvector_ctor2<Sample<double> >() ) );
        RTT::TypeInfoRepository::Instance()->type("SampleDoubles")->addConstructor(new StdVectorBuilder<Sample<double> >() );

        RTT::TypeInfoRepository::Instance()->type("SampleColumnVectors")->addConstructor(newConstructor(stdvector_ctor<Sample<ColumnVector> >() ) );
        RTT::TypeInfoRepository::Instance()->type("SampleColumnVectors")->addConstructor(newConstructor(stdvector_ctor2<Sample<ColumnVector> >() ) );
        RTT::TypeInfoRepository::Instance()->type("SampleColumnVectors")->addConstructor(new StdVectorBuilder<Sample<ColumnVector> >() );


        RTT::TypeInfoRepository::Instance()->type("WeightedSampleInts")->addConstructor(newConstructor(stdvector_ctor<WeightedSample<int> >() ) );
        RTT::TypeInfoRepository::Instance()->type("WeightedSampleInts")->addConstructor(newConstructor(stdvector_ctor2<WeightedSample<int> >() ) );
        RTT::TypeInfoRepository::Instance()->type("WeightedSampleInts")->addConstructor(new StdVectorBuilder<WeightedSample<int> >() );

        RTT::TypeInfoRepository::Instance()->type("WeightedSampleDoubles")->addConstructor(newConstructor(stdvector_ctor<WeightedSample<double> >() ) );
        RTT::TypeInfoRepository::Instance()->type("WeightedSampleDoubles")->addConstructor(newConstructor(stdvector_ctor2<WeightedSample<double> >() ) );
        RTT::TypeInfoRepository::Instance()->type("WeightedSampleDoubles")->addConstructor(new StdVectorBuilder<WeightedSample<double> >() );

        RTT::TypeInfoRepository::Instance()->type("WeightedSampleColumnVectors")->addConstructor(newConstructor(stdvector_ctor<WeightedSample<ColumnVector> >() ) );
        RTT::TypeInfoRepository::Instance()->type("WeightedSampleColumnVectors")->addConstructor(newConstructor(stdvector_ctor2<WeightedSample<ColumnVector> >() ) );
        RTT::TypeInfoRepository::Instance()->type("WeightedSampleColumnVectors")->addConstructor(new StdVectorBuilder<WeightedSample<ColumnVector> >() );

        RTT::TypeInfoRepository::Instance()->type("VecSampleInts")->addConstructor(newConstructor(stdvector_ctor<vector<Sample<int> > >() ) );
        RTT::TypeInfoRepository::Instance()->type("VecSampleInts")->addConstructor(newConstructor(stdvector_ctor2<vector<Sample<int> > >() ) );
        RTT::TypeInfoRepository::Instance()->type("VecSampleInts")->addConstructor(new StdVectorBuilder<vector<Sample<int> > >() );

        RTT::TypeInfoRepository::Instance()->type("VecSampleDoubles")->addConstructor(newConstructor(stdvector_ctor<vector<Sample<double> > >() ) );
        RTT::TypeInfoRepository::Instance()->type("VecSampleDoubles")->addConstructor(newConstructor(stdvector_ctor2<vector<Sample<double> > >() ) );
        RTT::TypeInfoRepository::Instance()->type("VecSampleDoubles")->addConstructor(new StdVectorBuilder<vector<Sample<double>  > >() );

        RTT::TypeInfoRepository::Instance()->type("VecSampleColumnVectors")->addConstructor(newConstructor(stdvector_ctor<vector<Sample<ColumnVector> > >() ) );
        RTT::TypeInfoRepository::Instance()->type("VecSampleColumnVectors")->addConstructor(newConstructor(stdvector_ctor2<vector<Sample<ColumnVector> > >() ) );
        RTT::TypeInfoRepository::Instance()->type("VecSampleColumnVectors")->addConstructor(new StdVectorBuilder<vector<Sample<ColumnVector> > >() );

        RTT::TypeInfoRepository::Instance()->type("VecWeightedSampleInts")->addConstructor(newConstructor(stdvector_ctor<vector<WeightedSample<int> > >() ) );
        RTT::TypeInfoRepository::Instance()->type("VecWeightedSampleInts")->addConstructor(newConstructor(stdvector_ctor2<vector<WeightedSample<int> > >() ) );
        RTT::TypeInfoRepository::Instance()->type("VecWeightedSampleInts")->addConstructor(new StdVectorBuilder<vector<WeightedSample<int> > >() );

        RTT::TypeInfoRepository::Instance()->type("VecWeightedSampleDoubles")->addConstructor(newConstructor(stdvector_ctor<vector<WeightedSample<double> > >() ) );
        RTT::TypeInfoRepository::Instance()->type("VecWeightedSampleDoubles")->addConstructor(newConstructor(stdvector_ctor2<vector<WeightedSample<double> > >() ) );
        RTT::TypeInfoRepository::Instance()->type("VecWeightedSampleDoubles")->addConstructor(new StdVectorBuilder<vector<WeightedSample<double>  > >() );

        RTT::TypeInfoRepository::Instance()->type("VecWeightedSampleColumnVectors")->addConstructor(newConstructor(stdvector_ctor<vector<WeightedSample<ColumnVector> > >() ) );
        RTT::TypeInfoRepository::Instance()->type("VecWeightedSampleColumnVectors")->addConstructor(newConstructor(stdvector_ctor2<vector<WeightedSample<ColumnVector> > >() ) );
        RTT::TypeInfoRepository::Instance()->type("VecWeightedSampleColumnVectors")->addConstructor(new StdVectorBuilder<vector<WeightedSample<ColumnVector> > >() );

        RTT::TypeInfoRepository::Instance()->type("Probabilitys")->addConstructor(newConstructor(stdvector_ctor<Probability>() ) );
        RTT::TypeInfoRepository::Instance()->type("Probabilitys")->addConstructor(newConstructor(stdvector_ctor2<Probability>() ) );
        RTT::TypeInfoRepository::Instance()->type("Probabilitys")->addConstructor(new StdVectorBuilder<Probability>() );

        RTT::OperatorRepository::Instance()->add( newBinaryOperator( "[]", stdvector_index<ColumnVector>() ) );
        RTT::OperatorRepository::Instance()->add( newDotOperator( "size", RTT::get_size<const std::vector<ColumnVector>&>() ) );

        RTT::TypeInfoRepository::Instance()->type("SampleInt")->addConstructor(newConstructor(Sample_ctor<int>() ) );
        RTT::TypeInfoRepository::Instance()->type("SampleDouble")->addConstructor(newConstructor(Sample_ctor<double>() ) );
        RTT::TypeInfoRepository::Instance()->type("SampleColumnVector")->addConstructor(newConstructor(Sample_ctor<ColumnVector>() ) );
        RTT::TypeInfoRepository::Instance()->type("WeightedSampleInt")->addConstructor(newConstructor(WeightedSample_ctor<int>() ) );
        RTT::TypeInfoRepository::Instance()->type("WeightedSampleDouble")->addConstructor(newConstructor(WeightedSample_ctor<double>() ) );
        RTT::TypeInfoRepository::Instance()->type("WeightedSampleColumnVector")->addConstructor(newConstructor(WeightedSample_ctor<ColumnVector>() ) );

        RTT::TypeInfoRepository::Instance()->type("Probability")->addConstructor(newConstructor(Probability_ctor() ) );

        return true;
    }

    bool bflToolkitPlugin::loadOperators()
    {
        RTT::OperatorRepository::Instance()->add( newBinaryOperator( "[]", vector_index() ) );
        RTT::OperatorRepository::Instance()->add( newBinaryOperator( "[]", rvector_index() ) );
        RTT::OperatorRepository::Instance()->add( newDotOperator( "size", get_size() ) );
        RTT::OperatorRepository::Instance()->add( newDotOperator( "size", rget_size() ) );
        RTT::OperatorRepository::Instance()->add( newBinaryOperator( "+", std::plus<ColumnVector>() ) );
        RTT::OperatorRepository::Instance()->add( newBinaryOperator( "+", std::plus<RowVector>() ) );
        //RTT::OperatorRepository::Instance()->add( newTernaryOperator( "[,]", matrix_index() ) );
        //RTT::OperatorRepository::Instance()->add( newBinaryOperator( "[]", matrix_index() ) );
        return true;
    }
}

ORO_TOOLKIT_PLUGIN(BFL::bflToolkit)
