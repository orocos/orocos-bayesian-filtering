// Copyright  (C)  2010  Tinne De Laet <tinne dot delaet at mech dot kuleuven dot be>

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

#include "bfl_typekit.hpp"
#include "SampleComposition.hpp"


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

    bflTypekitPlugin bflTypekit;

    struct IntVectorTypeInfo: public SequenceTypeInfo<std::vector<int>,false > 
    {
         IntVectorTypeInfo():SequenceTypeInfo< std::vector<int>,false > ("ints")
        {
        }
    };

    template<class Archive>
    void serialize(Archive & a, Probability & prob, unsigned int){
        using boost::serialization::make_nvp;   
        a & make_nvp("Probability", prob.getValue() );
    }

    struct ProbabilityTypeInfo : public StructTypeInfo<Probability,true>
    {
        ProbabilityTypeInfo():StructTypeInfo<Probability,true>("Probability")
        {
        };
    };

    struct ProbabilityVectorTypeInfo: public SequenceTypeInfo<std::vector<Probability>,false > 
    {
         ProbabilityVectorTypeInfo():SequenceTypeInfo< std::vector<Probability>,false > ("Probabilitys")
        {
        }
    };

    struct VectorTypeInfo : public SequenceTypeInfo<ColumnVector,true>
    {
        VectorTypeInfo():SequenceTypeInfo<ColumnVector, true >("ColumnVector")
        {
        }
    };

    struct ColumnVectorsTypeInfo: public SequenceTypeInfo<std::vector<ColumnVector>,false > 
    {
         ColumnVectorsTypeInfo():SequenceTypeInfo< std::vector<ColumnVector>,false > ("ColumnVectors")
        {
        }
    };

    struct RVectorTypeInfo : public SequenceTypeInfo<RowVector,true>
    {
        RVectorTypeInfo():SequenceTypeInfo<RowVector, true >("RowVector")
        {
        }
    };

    struct RowVectorsTypeInfo: public SequenceTypeInfo<std::vector<RowVector>,false > 
    {
         RowVectorsTypeInfo():SequenceTypeInfo< std::vector<RowVector>,false > ("RowVectors")
        {
        }
    };


    struct MatrixTypeInfo : public TemplateTypeInfo<Matrix, true>
    {
        MatrixTypeInfo():TemplateTypeInfo<Matrix,true>("Matrix"){
        };
        

        base::DataSourceBase::shared_ptr convertType(base::DataSourceBase::shared_ptr source) const
         {
           log(Debug)<<"Converting Matrix to PropertyBag"<<endlog();
           internal::DataSource<Matrix>::shared_ptr ds = internal::DataSource<Matrix>::narrow( source.get() );
           PropertyBag targetbag;
           targetbag.setType("Matrix");
           unsigned int dimension = ds->get().rows(); //ds->get() returns matrix

           for ( unsigned int i=1; i <= dimension ; i++)
            {
                std::stringstream out;
                out << i;
                Property<PropertyBag>* row_bag = new Property<PropertyBag>(out.str(), out.str() +"th row of matrix"); 
                Property<RowVector> row(out.str(), out.str() +"th row of matrix",ds->get().rowCopy(i))  ;
                typeDecomposition( row.getDataSource(), row_bag->value());
                targetbag.add(row_bag);
            }
           internal::ValueDataSource<PropertyBag>::shared_ptr targetbag_ptr = new internal::ValueDataSource<PropertyBag>(targetbag);
           return targetbag_ptr;
         }


        bool composeTypeImpl(const PropertyBag& bag, Matrix& result) const{
            if ( bag.getType() == "Matrix" ) {
                unsigned int rows = bag.size();
                unsigned int cols = 0;
                // Get values
                for (unsigned int i = 1; i <= rows ; i++) {
                    std::stringstream out;
                    out << i;
                    Property<PropertyBag>* row_bag ;
                    row_bag=  bag.getPropertyType<PropertyBag>(out.str());
                    if(row_bag==NULL){
                        log(Error)<<"Could not read row "<<i<<endlog();
                        return false;
                    }
                    Property<RowVector> row_p(row_bag->getName(),row_bag->getDescription());
                    if(!(row_p.getTypeInfo()->composeType(row_bag->getDataSource(),row_p.getDataSource()))){
                        log(Error)<<"Could not compose row "<<i<<endlog();
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

    struct MatrixsTypeInfo: public SequenceTypeInfo<std::vector<Matrix>,false > 
    {
         MatrixsTypeInfo():SequenceTypeInfo< std::vector<Matrix>,false > ("Matrixs")
        {
        }
    };

    struct SymmetricMatrixTypeInfo : public TemplateTypeInfo<SymmetricMatrix,true>
    {
        SymmetricMatrixTypeInfo():TemplateTypeInfo<SymmetricMatrix,true>("SymmetricMatrix"){
        };

        base::DataSourceBase::shared_ptr convertType(base::DataSourceBase::shared_ptr source) const
         {
           log(Debug)<<"Converting SymmetricMatrix to PropertyBag"<<endlog();
           internal::DataSource<SymmetricMatrix>::shared_ptr ds = internal::DataSource<SymmetricMatrix>::narrow( source.get() );
           PropertyBag targetbag;
           targetbag.setType("SymmetricMatrix");
           unsigned int dimension = ds->get().rows(); //ds->get() returns matrix

           for ( unsigned int i=1; i <= dimension ; i++)
            {
                std::stringstream out;
                out << i;
                Property<PropertyBag>* row_bag = new Property<PropertyBag>(out.str(), out.str() +"th row of matrix"); 
                Property<RowVector> row(out.str(), out.str() +"th row of matrix",ds->get().rowCopy(i))  ;
                typeDecomposition( row.getDataSource(), row_bag->value());
                targetbag.add(row_bag);
            }
           internal::ValueDataSource<PropertyBag>::shared_ptr targetbag_ptr = new internal::ValueDataSource<PropertyBag>(targetbag);
           return targetbag_ptr;
         }

        bool composeTypeImpl(const PropertyBag& bag, SymmetricMatrix& result) const{
            Matrix matrix;
            if ( bag.getType() == "SymmetricMatrix" ) {
                unsigned int rows = bag.size();
                unsigned int cols = 0;
                // Get values
                for (unsigned int i = 1; i <= rows ; i++) {
                    std::stringstream out;
                    out << i;
                    Property<PropertyBag>* row_bag =  bag.getPropertyType<PropertyBag>(out.str());
                    if(row_bag==NULL){
                        log(Error)<<"Could not read row "<<i<<endlog();
                        return false;
                    }
                    Property<RowVector > row_p(row_bag->getName(),row_bag->getDescription());
                    if(!(row_p.getTypeInfo()->composeType(row_bag->getDataSource(),row_p.getDataSource()))){
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

    struct SymmetricMatrixsTypeInfo: public SequenceTypeInfo<std::vector<SymmetricMatrix>,false > 
    {
         SymmetricMatrixsTypeInfo():SequenceTypeInfo< std::vector<SymmetricMatrix>,false > ("SymmetricMatrixs")
        {
        }
    };

    std::string bflTypekitPlugin::getName()
    {
        return "BFL_Typekit";
    }

    bool bflTypekitPlugin::loadTypes()
    {
        RTT::types::Types()->addType( new IntVectorTypeInfo() );

        RTT::types::Types()->addType( new ProbabilityTypeInfo() );

        RTT::types::Types()->addType( new VectorTypeInfo() );
        RTT::types::Types()->addType( new RVectorTypeInfo() );

        RTT::types::Types()->addType( new MatrixTypeInfo() );
        RTT::types::Types()->addType( new SymmetricMatrixTypeInfo() );

        RTT::types::Types()->addType( new SampleTypeInfo<int >("SampleInt") );
        RTT::types::Types()->addType( new SampleTypeInfo<double >("SampleDouble") );
        RTT::types::Types()->addType( new SampleTypeInfo<ColumnVector >("SampleColumnVector") );
        RTT::types::Types()->addType( new WeightedSampleTypeInfo<int>("WeightedSampleInt") );
        RTT::types::Types()->addType( new WeightedSampleTypeInfo<double> ("WeightedSampleDouble") );
        RTT::types::Types()->addType( new WeightedSampleTypeInfo<ColumnVector> ("WeightedSampleColumnVector") );

        RTT::types::Types()->addType( new ProbabilityVectorTypeInfo() );
        RTT::types::Types()->addType( new ColumnVectorsTypeInfo() );
        RTT::types::Types()->addType( new RowVectorsTypeInfo() );
        RTT::types::Types()->addType( new MatrixsTypeInfo() );
        RTT::types::Types()->addType( new SymmetricMatrixsTypeInfo() );
        RTT::types::Types()->addType( new SamplesTypeInfo<int >("SampleInts") );
        RTT::types::Types()->addType( new SamplesTypeInfo<double >("SampleDoubles") );
        RTT::types::Types()->addType( new SamplesTypeInfo<ColumnVector >("SampleColumnVectors") );
        RTT::types::Types()->addType( new WeightedSamplesTypeInfo<int>("WeightedSampleInts") );
        RTT::types::Types()->addType( new WeightedSamplesTypeInfo<double> ("WeightedSampleDoubles") );
        RTT::types::Types()->addType( new WeightedSamplesTypeInfo<ColumnVector> ("WeightedSampleColumnVectors") );

        return true;
    }

    /*****************
    CONSTRUCTORS
    *******************/
    Probability createProbability(double p)
    {
        return Probability(p);
    }

    ColumnVector createColumnVector1(int rows)
    {
        return ColumnVector(rows);
    }

    ColumnVector createColumnVector2(int rows, double value)
    {
        return ColumnVector(rows,value);
    }

    RowVector createRowVector1(int cols)
    {
        return RowVector(cols);
    }

    RowVector createRowVector2(int cols, double value)
    {
        return RowVector(cols,value);
    }

    Matrix createMatrix1(int rows, int cols)
    {
        return Matrix(rows,cols);
    }

    SymmetricMatrix createSymmetricMatrix(int dimension )
    {
        return SymmetricMatrix(dimension);
    }


    bool bflTypekitPlugin::loadConstructors()
    {
        RTT::types::Types()->type("Probability")->addConstructor(newConstructor(&createProbability));

        RTT::types::Types()->type("ColumnVector")->addConstructor(newConstructor(&createColumnVector1));
        RTT::types::Types()->type("ColumnVector")->addConstructor(newConstructor(&createColumnVector2));

        RTT::types::Types()->type("RowVector")->addConstructor(newConstructor(&createRowVector1));
        RTT::types::Types()->type("RowVector")->addConstructor(newConstructor(&createRowVector2));

        RTT::types::Types()->type("Matrix")->addConstructor(newConstructor(&createMatrix1));
        RTT::types::Types()->type("SymmetricMatrix")->addConstructor(newConstructor(&createSymmetricMatrix));

        RTT::types::Types()->type("SampleColumnVector")->addConstructor(newConstructor(&createSampleColumnVector));
        RTT::types::Types()->type("WeightedSampleColumnVector")->addConstructor(newConstructor(&createWeightedSampleColumnVector));


        return true;

    }

    bool bflTypekitPlugin::loadOperators()
    {
        return true;
    }
}

ORO_TYPEKIT_PLUGIN(BFL::bflTypekitPlugin)
