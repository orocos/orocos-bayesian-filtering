#include "bfl_typekit.hpp"

namespace BFL{
    using namespace std;
    using namespace RTT;
    using namespace RTT::detail;

    /*****************************************************************************
     * MATRIX
     * **************************************************************************/
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

    /*****************************************************************************
     * SYMMETRIXMATRIX
     * **************************************************************************/
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

    void loadMatrixTypes(){
          RTT::types::Types()->addType( new MatrixTypeInfo() );
          RTT::types::Types()->addType( new SymmetricMatrixTypeInfo() );
    }
}
