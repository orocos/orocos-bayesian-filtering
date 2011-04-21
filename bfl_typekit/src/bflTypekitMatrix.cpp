#include "bfl_typekit.hpp"

namespace BFL{
    using namespace std;
    using namespace RTT;
    using namespace RTT::detail;

    template<class T>
    int get_rows(T const& cont)
    {
        return cont.rows();
    }

    template<class T>
    int get_columns(T const& cont)
    {
        return cont.columns();
    }

    template<class T>
    RowVector  get_matrix_item(T & cont, int index)
    {
        if (index >= (int) (cont.rows()) || index < 0)
            return internal::NA<RowVector>::na();
        return cont[index];
    }
     
    template<class T>
    RowVector get_matrix_item_copy(const T & cont, int index)
    {
        if (index >= (int) (cont.rows()) || index < 0)
            return internal::NA<RowVector>::na();
        return cont[index];
    }

    /*****************************************************************************
     * MATRIX
     * **************************************************************************/
    struct MatrixTypeInfo : public TemplateTypeInfo<Matrix, true>
    {
        MatrixTypeInfo():TemplateTypeInfo<Matrix,true>("Matrix"){
        };

        virtual base::DataSourceBase::shared_ptr getMember(base::DataSourceBase::shared_ptr item, const std::string& name) const {
            // the only thing we do is to check for an integer in name, otherwise, assume a part (size/capacity) is accessed:
            try {
                unsigned int indx = boost::lexical_cast<unsigned int>(name);
                // @todo could also return a direct reference to item indx using another DS type that respects updated().
                return getMember( item, new internal::ConstantDataSource<int>(indx));
            } catch(...) {}

            return getMember( item, new internal::ConstantDataSource<std::string>(name) );
        }

        virtual base::DataSourceBase::shared_ptr getMember(base::DataSourceBase::shared_ptr item,
                                                         base::DataSourceBase::shared_ptr id) const {
            // discover if user gave us a part name or index:
            internal::DataSource<int>::shared_ptr id_indx = internal::DataSource<int>::narrow( internal::DataSourceTypeInfo<int>::getTypeInfo()->convert(id).get() );
            internal::DataSource<string>::shared_ptr id_name = internal::DataSource<string>::narrow( id.get() );
            if ( id_name ) {
                if ( id_name->get() == "rows" ) {
                    try {
                        return internal::newFunctorDataSource(&get_rows<Matrix>, internal::GenerateDataSource()(item.get()) );
                    } catch(...) {}
                }
                if ( id_name->get() == "columns" ) {
                    try {
                        return internal::newFunctorDataSource(&get_columns<Matrix>, internal::GenerateDataSource()(item.get()) );
                    } catch(...) {}
                }
            }

            if ( id_indx ) {
                try {
                    if ( item->isAssignable() )
                            return internal::newFunctorDataSource(&get_matrix_item<Matrix>, internal::GenerateDataSource()(item.get(), id_indx.get() ) );
                        else
                            return internal::newFunctorDataSource(&get_matrix_item_copy<Matrix>, internal::GenerateDataSource()(item.get(), id_indx.get() ) );
                } catch(...) {}
            }
            if (id_name) {
                log(Error) << "SequenceTypeInfo: No such member : " << id_name->get() << endlog();
            }
            if (id_indx) {
                log(Error) << "SequenceTypeInfo: Invalid index : " << id_indx->get() <<":"<< id_indx->getTypeName() << endlog();
            }
            if ( !id_name && ! id_indx)
                log(Error) << "SequenceTypeInfo: Not a member or index : " << id <<":"<< id->getTypeName() << endlog();
            return base::DataSourceBase::shared_ptr();
        };
        

        bool decomposeTypeImpl(const Matrix& mat, PropertyBag& targetbag) const{
            targetbag.setType("Matrix");
            unsigned int dimension = mat.rows();

            for ( unsigned int i=1; i <= dimension ; i++){
                std::stringstream out;
                out << i;
                Property<PropertyBag>* row_bag = new Property<PropertyBag>(out.str(), out.str() +"th row of matrix"); 
                Property<RowVector> row(out.str(), out.str() +"th row of matrix",mat.rowCopy(i))  ;
                typeDecomposition( row.getDataSource(), row_bag->value());
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

        virtual base::DataSourceBase::shared_ptr getMember(base::DataSourceBase::shared_ptr item, const std::string& name) const {
            // the only thing we do is to check for an integer in name, otherwise, assume a part (size/capacity) is accessed:
            try {
                unsigned int indx = boost::lexical_cast<unsigned int>(name);
                // @todo could also return a direct reference to item indx using another DS type that respects updated().
                return getMember( item, new internal::ConstantDataSource<int>(indx));
            } catch(...) {}

            return getMember( item, new internal::ConstantDataSource<std::string>(name) );
        }

        virtual base::DataSourceBase::shared_ptr getMember(base::DataSourceBase::shared_ptr item,
                                                         base::DataSourceBase::shared_ptr id) const {
            // discover if user gave us a part name or index:
            internal::DataSource<int>::shared_ptr id_indx = internal::DataSource<int>::narrow( internal::DataSourceTypeInfo<int>::getTypeInfo()->convert(id).get() );
            internal::DataSource<string>::shared_ptr id_name = internal::DataSource<string>::narrow( id.get() );
            if ( id_name ) {
                if ( id_name->get() == "rows" ) {
                    try {
                        return internal::newFunctorDataSource(&get_rows<SymmetricMatrix>, internal::GenerateDataSource()(item.get()) );
                    } catch(...) {}
                }
                if ( id_name->get() == "columns" ) {
                    try {
                        return internal::newFunctorDataSource(&get_columns<SymmetricMatrix>, internal::GenerateDataSource()(item.get()) );
                    } catch(...) {}
                }
            }

            if ( id_indx ) {
                try {
                    if ( item->isAssignable() )
                            return internal::newFunctorDataSource(&get_matrix_item<SymmetricMatrix>, internal::GenerateDataSource()(item.get(), id_indx.get() ) );
                        else
                            return internal::newFunctorDataSource(&get_matrix_item_copy<SymmetricMatrix>, internal::GenerateDataSource()(item.get(), id_indx.get() ) );
                } catch(...) {}
            }
            if (id_name) {
                log(Error) << "SequenceTypeInfo: No such member : " << id_name->get() << endlog();
            }
            if (id_indx) {
                log(Error) << "SequenceTypeInfo: Invalid index : " << id_indx->get() <<":"<< id_indx->getTypeName() << endlog();
            }
            if ( !id_name && ! id_indx)
                log(Error) << "SequenceTypeInfo: Not a member or index : " << id <<":"<< id->getTypeName() << endlog();
            return base::DataSourceBase::shared_ptr();
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
                typeDecomposition( row.getDataSource(), row_bag->value());
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
                    Property<PropertyBag>* row_bag =  bag.getPropertyType<PropertyBag>(out.str());
                    if(row_bag==NULL){
                        log(Error)<<"Could not read row "<<i<<endlog();
                        return false;
                    }
                    Property<RowVector > row_p(row_bag->getName(),row_bag->getDescription());
                    if(!(row_p.getTypeInfo()->composeType(row_bag->getDataSource(),row_p.getDataSource()))){
                        log(Error)<<"Could not compose row "<<i<<endlog();
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
