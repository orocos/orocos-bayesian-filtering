
#include "bfl_typekit.hpp"
using namespace MatrixWrapper;
namespace BFL{

    /*****************
    CONSTRUCTORS
    *******************/
    Probability createProbability(double p)
    {
        return Probability(p);
    }

    ColumnVector createColumnVector2(int rows, double value)
    {
        return ColumnVector(rows,value);
    }

    RowVector createRowVector2(int cols, double value)
    {
        return RowVector(cols,value);
    }

    Matrix createMatrix1(int rows, int cols)
    {
        return Matrix(rows,cols);
    }

    Matrix createMatrix2(int rows, RowVector row)
    {
        return Matrix(rows,row);
    }

    SymmetricMatrix createSymmetricMatrix(int dimension )
    {
        return SymmetricMatrix(dimension);
    }

    bool bflTypekitPlugin::loadConstructors()
    {
        RTT::types::TypeInfoRepository::shared_ptr ti = RTT::types::TypeInfoRepository::Instance(); 
        ti->type("Probability")->addConstructor( RTT::types::newConstructor(&createProbability) );
        ti->type("ColumnVector")->addConstructor( RTT::types::newConstructor(&createColumnVector2) );
        ti->type("RowVector")->addConstructor( RTT::types::newConstructor(&createRowVector2) );

        ti->type("Matrix")->addConstructor( RTT::types::newConstructor(&createMatrix1) );
        ti->type("Matrix")->addConstructor( RTT::types::newConstructor(&createMatrix2) );
        ti->type("SymmetricMatrix")->addConstructor( RTT::types::newConstructor(&createSymmetricMatrix) );

        //RTT::types::Types()->type("SampleColumnVector")->addConstructor(newConstructor(&createSampleColumnVector));
        //RTT::types::Types()->type("WeightedSampleColumnVector")->addConstructor(newConstructor(&createWeightedSampleColumnVector));
        return true;
    }

}
