
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
        //RTT::types::Types()->type("Probability")->addConstructor(newConstructor(&createProbability));

        //RTT::types::Types()->type("ColumnVector")->addConstructor(newConstructor(&createColumnVector1));
        //RTT::types::Types()->type("ColumnVector")->addConstructor(newConstructor(&createColumnVector2));

        //RTT::types::Types()->type("RowVector")->addConstructor(newConstructor(&createRowVector1));
        //RTT::types::Types()->type("RowVector")->addConstructor(newConstructor(&createRowVector2));

        //RTT::types::Types()->type("Matrix")->addConstructor(newConstructor(&createMatrix1));
        //RTT::types::Types()->type("SymmetricMatrix")->addConstructor(newConstructor(&createSymmetricMatrix));

        //RTT::types::Types()->type("SampleColumnVector")->addConstructor(newConstructor(&createSampleColumnVector));
        //RTT::types::Types()->type("WeightedSampleColumnVector")->addConstructor(newConstructor(&createWeightedSampleColumnVector));


        return true;

    }

}
