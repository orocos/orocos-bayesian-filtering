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

using namespace MatrixWrapper;
namespace BFL{ 

    using namespace RTT;
    using namespace RTT::detail;

    bflTypekitPlugin bflTypekit;

    void loadProbabilityTypes();
    void loadProbabilitysTypes();
    void loadVectorTypes();
    void loadVectorsTypes();
    void loadMatrixTypes();
    void loadMatrixsTypes();
    void loadSampleTypes();
    void loadWeightedSampleTypes();
    void loadSamplesTypes();
    void loadWeightedSamplesTypes();

    std::string bflTypekitPlugin::getName()
    {
        return "BFL_Typekit";
    }

    bool bflTypekitPlugin::loadTypes()
    {
        // load probability types  
        loadProbabilityTypes();
        // load probabilitys types  
        loadProbabilitysTypes();
        // load vector types
        loadVectorTypes();
        // load vectors types
        loadVectorsTypes();
        // load matrix types
        loadMatrixTypes();
        // load matrixs types
        loadMatrixsTypes();
        // load Sample types
        loadSampleTypes();
        // load Samples types
        loadSamplesTypes();
        // load WeightedSample types
        loadWeightedSampleTypes();
        // load WeightedSamples types
        loadWeightedSamplesTypes();

        return true;
    }

    // resize matrix
    void resize_matrix(Matrix& m, int rows, int cols)
    {
         m.resize(rows,cols);
    }

    // resize SymmetricMatrix
    void resize_symmetricMatrix(SymmetricMatrix& m, int size)
    {
         m.resize(size);
    }

    // get element of matrix  
    template<class T>                                                                                                                                        
    double & get_matrix_el(T& m, int rows, int cols)
    {
         return m(rows,cols);
    }

    // set element of matrix  
    template<class T>                                                                                                                                        
    void set_matrix_el(T& m, int rows, int cols, double value)
    {
         m(rows,cols)= value;
    }

    // resize ColumnVector
    template<class T>                                                                                                                                        
    void resize_vector(T& v, int size)
    {
         v.resize(size);
    }


    //does not work: ask Peter/Ruben
    template<class T>
    struct set_T_d
        : public std::binary_function<const T, const double , T>
    {
        T operator()(T t , double d ) const
        {
            return t=d;
        }
    };


    //struct multiply_m_sm
    //    : public std::binary_function<const Matrix, const SymmetricMatrix, Matrix>
    //{
    //    Matrix operator()(Matrix m, SymmetricMatrix sm ) const
    //    {
    //        return m*sm;
    //    }
    //};

    struct multiply_sm_m
        : public std::binary_function<const SymmetricMatrix, const Matrix, Matrix>
    {
        Matrix operator()(SymmetricMatrix sm, Matrix m ) const
        {
            return sm*m;
        }
    };

    struct add_sm_m
        : public std::binary_function<const SymmetricMatrix, const Matrix, Matrix>
    {
        Matrix operator()(SymmetricMatrix sm, Matrix m ) const
        {
            return sm+m;
        }
    };

    struct subtract_sm_m
        : public std::binary_function<const SymmetricMatrix, const Matrix, Matrix>
    {
        Matrix operator()(SymmetricMatrix sm, Matrix m ) const
        {
            return sm-m;
        }
    };

    struct multiply_sm_cv
        : public std::binary_function<const SymmetricMatrix, const ColumnVector, ColumnVector>
    {
        ColumnVector operator()(SymmetricMatrix sm, ColumnVector cv ) const
        {
            return sm*cv;
        }
    };

    struct multiply_m_cv
        : public std::binary_function<const Matrix, const ColumnVector, ColumnVector>
    {
        ColumnVector operator()(Matrix m, ColumnVector cv ) const
        {
            return m*cv;
        }
    };


    struct multiply_cv_rv
        : public std::binary_function<const ColumnVector, const RowVector, Matrix>                                                                                                                 
    {                                                                                                                                                        
        Matrix operator()(ColumnVector cv, RowVector rv ) const                                                                                                                        
        {                                                                                                                                                    
            return cv*rv;                                                                                                                              
        }                                                                                                                                                    
    };                                                                                                                                                       

    struct multiply_rv_cv
        : public std::binary_function<const RowVector, const ColumnVector, double>                                                                                                                 
    {                                                                                                                                                        
        double operator()(RowVector rv, ColumnVector cv ) const                                                                                                                        
        {                                                                                                                                                    
            return rv*cv;                                                                                                                              
        }                                                                                                                                                    
    };                                                                                                                                                       

    template<class T>
    struct multiply_T_double
        : public std::binary_function<const T, const double, T>
    {
        T operator()(T t, double d ) const
        {
            return t*d;
        }
    };

    template<class T>
    struct divide_T_double
        : public std::binary_function<const T, const double, T>
    {
        T operator()(T t, double d ) const
        {
            return t/d;
        }
    };

    template<class T>
    struct add_T_double
        : public std::binary_function<const T, const double, T>
    {
        T operator()(T t, double d ) const
        {
            return t+d;
        }
    };

    template<class T>
    struct subtract_T_double
        : public std::binary_function<const T, const double, T>
    {
        T operator()(T t, double d ) const
        {
            return t-d;
        }
    };

    bool bflTypekitPlugin::loadOperators()
    {
        RTT::types::OperatorRepository::shared_ptr oreg = RTT::types::OperatorRepository::Instance();
        //oreg->add( RTT::types::newBinaryOperator( "[]", matrix_row() ) );
        oreg->add( RTT::types::newBinaryOperator( "+", std::plus<ColumnVector>() ) );
        oreg->add( RTT::types::newBinaryOperator( "+", std::plus<RowVector>() ) );
        oreg->add( RTT::types::newBinaryOperator( "+", std::plus<Matrix>() ) );
        oreg->add( RTT::types::newBinaryOperator( "+", std::plus<SymmetricMatrix>() ) );
        oreg->add( RTT::types::newBinaryOperator( "-", std::minus<ColumnVector>() ) );
        oreg->add( RTT::types::newBinaryOperator( "-", std::minus<RowVector>() ) );
        oreg->add( RTT::types::newBinaryOperator( "-", std::minus<Matrix>() ) );
        oreg->add( RTT::types::newBinaryOperator( "-", std::minus<SymmetricMatrix>() ) );
        oreg->add( RTT::types::newBinaryOperator( "*", std::multiplies<Matrix>() ) );
        oreg->add( RTT::types::newBinaryOperator( "==", std::equal_to<ColumnVector>() ) );
        oreg->add( RTT::types::newBinaryOperator( "==", std::equal_to<RowVector>() ) );
        oreg->add( RTT::types::newBinaryOperator( "==", std::equal_to<Matrix>() ) );
        oreg->add( RTT::types::newBinaryOperator( "==", std::equal_to<SymmetricMatrix>() ) );
        oreg->add( RTT::types::newBinaryOperator( "*" ,multiply_cv_rv() ) );
        oreg->add( RTT::types::newBinaryOperator( "*" ,multiply_rv_cv() ) );
        oreg->add( RTT::types::newBinaryOperator( "=" ,set_T_d<Matrix>() ) );
        oreg->add( RTT::types::newBinaryOperator( "=" ,set_T_d<SymmetricMatrix>() ) );
        oreg->add( RTT::types::newBinaryOperator( "=" ,set_T_d<ColumnVector>() ) );
        oreg->add( RTT::types::newBinaryOperator( "=" ,set_T_d<RowVector>() ) );
        oreg->add( RTT::types::newBinaryOperator( "*" ,multiply_sm_m() ) );
        oreg->add( RTT::types::newBinaryOperator( "*" ,add_sm_m() ) );
        oreg->add( RTT::types::newBinaryOperator( "*" ,subtract_sm_m() ) );
        oreg->add( RTT::types::newBinaryOperator( "*" ,multiply_m_cv() ) );
        oreg->add( RTT::types::newBinaryOperator( "*" ,multiply_sm_cv() ) );
        //oreg->add( RTT::types::newBinaryOperator( "*" ,multiply_m_sm() ) );
        oreg->add( RTT::types::newBinaryOperator( "*" ,multiply_T_double<ColumnVector>() ) );
        oreg->add( RTT::types::newBinaryOperator( "*" ,multiply_T_double<RowVector>() ) );
        oreg->add( RTT::types::newBinaryOperator( "*" ,multiply_T_double<Matrix>() ) );
        oreg->add( RTT::types::newBinaryOperator( "*" ,multiply_T_double<SymmetricMatrix>() ) );
        oreg->add( RTT::types::newBinaryOperator( "/" ,divide_T_double<ColumnVector>() ) );
        oreg->add( RTT::types::newBinaryOperator( "/" ,divide_T_double<RowVector>() ) );
        oreg->add( RTT::types::newBinaryOperator( "/" ,divide_T_double<Matrix>() ) );
        oreg->add( RTT::types::newBinaryOperator( "/" ,divide_T_double<SymmetricMatrix>() ) );
        oreg->add( RTT::types::newBinaryOperator( "+" ,add_T_double<ColumnVector>() ) );
        oreg->add( RTT::types::newBinaryOperator( "+" ,add_T_double<RowVector>() ) );
        oreg->add( RTT::types::newBinaryOperator( "+" ,add_T_double<Matrix>() ) );
        oreg->add( RTT::types::newBinaryOperator( "+" ,add_T_double<SymmetricMatrix>() ) );
        oreg->add( RTT::types::newBinaryOperator( "-" ,subtract_T_double<ColumnVector>() ) );
        oreg->add( RTT::types::newBinaryOperator( "-" ,subtract_T_double<RowVector>() ) );
        oreg->add( RTT::types::newBinaryOperator( "-" ,subtract_T_double<Matrix>() ) );
        oreg->add( RTT::types::newBinaryOperator( "-" ,subtract_T_double<SymmetricMatrix>() ) );
        //oreg->add( RTT::types::newBinaryOperator( "!=", std::not_equal_to<ColumnVector>() ) );
        //oreg->add( RTT::types::newBinaryOperator( "!=", std::not_equal_to<RowVector>() ) );
        //oreg->add( RTT::types::newBinaryOperator( "!=", std::not_equal_to<Matrix>() ) );
        //oreg->add( RTT::types::newBinaryOperator( "!=", std::not_equal_to<SymmetricMatrix>() ) );
        //oreg->add( RTT::types::newUnaryOperator( "-", std::negate<ColumnVector>() ) );
        //oreg->add( RTT::types::newUnaryOperator( "-", std::negate<RowVector>() ) );

        RTT::internal::GlobalService::Instance()->provides("bfl")->addOperation("resize_Matrix",&resize_matrix).doc("resize matrix");
        RTT::internal::GlobalService::Instance()->provides("bfl")->addOperation("resize_SymmetricMatrix",&resize_matrix).doc("resize SymmetricMatrix");
        RTT::internal::GlobalService::Instance()->provides("bfl")->addOperation("resize_ColumnVector",&resize_vector<ColumnVector>).doc("resize ColumnVector");
        RTT::internal::GlobalService::Instance()->provides("bfl")->addOperation("resize_RowVector",&resize_vector<RowVector>).doc("resize RowVector");
        RTT::internal::GlobalService::Instance()->provides("bfl")->addOperation("get_element_Matrix",&get_matrix_el<Matrix>).doc("get element of Matrix");
        RTT::internal::GlobalService::Instance()->provides("bfl")->addOperation("get_element_SymmetricMatrix",&get_matrix_el<SymmetricMatrix>).doc("get element of SymmetricMatrix");
        RTT::internal::GlobalService::Instance()->provides("bfl")->addOperation("set_element_Matrix",&set_matrix_el<Matrix>).doc("set element of Matrix");
        RTT::internal::GlobalService::Instance()->provides("bfl")->addOperation("set_element_SymmetricMatrix",&set_matrix_el<SymmetricMatrix>).doc("set element of SymmetricMatrix");

        return true;
    }
}

ORO_TYPEKIT_PLUGIN(BFL::bflTypekitPlugin)
