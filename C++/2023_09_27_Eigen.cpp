#include<iostream>
#include<Eigen/Dense>

using namespace std;
using namespace Eigen;

int main()
{
    //explicit declaration
    Matrix <float, 3, 3> matrixA;
    matrixA.setZero();
    cout << matrixA << endl;

    // typedef declaration
    Matrix3f matrixA1;
    matrixA1.setZero();
    cout << "\n" << matrixA1 << endl;
    // define a dynamic matrix, explixit declaration
    Matrix <float, Dynamic, Dynamic> matrixB;

    //define a dynamic matrix, typedef declaration
    MatrixXf matrixB1;

    //constructor, allocate memory, but do not initialize
    MatrixXd matrixC(10,10);

    //assigning matrix entries
    MatrixXd matrixC1(2,2);
    matrixC1(0,0) = 1;
    matrixC1(0,1) = 2;
    matrixC1(1,0) = 3;
    matrixC1(1,1) = 4;

    cout << endl << matrixC1 << endl;

    // fill-in the matrix entries using comma seperated values and print the matrix
    Matrix4f matrixD;
    matrixD << 1,2,3,4,
        5,6,7,8,
        9,10,11,12,
        13,14,15,16;

    cout << endl << matrixD << endl;

    // setting matrix entries - two approaches
    int rowNumber = 5;
    int columnNumber = 5;

    // matrix of zeros
    MatrixXf matrix1zeros;
    matrix1zeros = MatrixXf::Zero(rowNumber, columnNumber);
    cout << "\n \n" << matrix1zeros << endl;
    // another option
    MatrixXf matrix1zeros1;
    matrix1zeros1.setZero(rowNumber, columnNumber);
    cout << "\n \n" << matrix1zeros1 << endl;

    // matrix of ones
    MatrixXf matrix1ones;
    matrix1ones = MatrixXf::Ones(rowNumber, columnNumber);
    cout << "\n \n" << matrix1ones << endl;
    // another option
    MatrixXf matrix1ones1;
    matrix1ones1.setOnes(rowNumber,columnNumber);
    cout << "\n \n" << matrix1ones1 << endl;

    // matrix of constants
    float value = 1.1;
    MatrixXf matrix1const;
    matrix1const = MatrixXf::Constant(rowNumber, columnNumber, value);
    cout << endl << matrix1const << endl;
    // another option
    MatrixXf matrix1const1;
    matrix1const1.setConstant(rowNumber, columnNumber, value);
    cout << endl << matrix1const1 << endl;
    // identity matrix, two approaches
    
    rowNumber = 10;
    columnNumber = 10;
    
    MatrixXd matrixIdentity;
    matrixIdentity = MatrixXd::Identity(rowNumber, columnNumber);
    cout << endl << matrixIdentity << endl;

    //accessing matrix blocks
    MatrixXd matrixV(4,4);
    matrixV << 101, 102, 103, 104,
            105, 106, 107,108,
            109, 110, 111, 112,
            113, 114, 115, 116;
    //access the matrix composed of 1:2 rows and 1:2 columns of matrixV
    MatrixXd matrixVpartition = matrixV.block(0,0,2,2);
    cout << endl << matrixVpartition << endl;

    MatrixXd matrixVpartition2 = matrixV.block(1,1,2,2);
    cout << endl << matrixVpartition2 << endl;

    //accessing columns and rows of a matrix
    cout << endl << "Row1 of matrixV is \n" << matrixV.row(0);
    cout << endl << "Column1 of matrixV is \n" << matrixV.col(0);

    Matrix <double, 3, 1> vector1;
    vector1 << 1,2,3;
    MatrixXd diagonalMatrix;
    diagonalMatrix = vector1.asDiagonal();
    cout << "Diagonal matrix is \n\n" << diagonalMatrix;

    MatrixXd A1(2,2);
    MatrixXd B1(2,2);

    A1 << 1,2,
        3,4;
    B1 << 3,4,
        5,6;
    
    //additionn and subtraction
    MatrixXd C1 = A1 + B1;
    cout << C1 << endl;
    
    //matrix multiplication
    MatrixXd D1 = A1 * B1;
    cout << D1 << endl;

    //scalar multiplication
    int s1 = 2;
    MatrixXd F1;
    F1 = s1 * A1;
    cout << F1 << endl; 

    //matrix transpose
    MatrixXd At;
    MatrixXd R1;

    At = A1.transpose();
    cout << "Original matrix A1 \n" << A1;
    cout << "Its transpose \n" << At;

    R1 = A1.transpose() + B1;

    //safeway to do the matrix transpose
    A1.transposeInPlace();
    
    //matrix inverse
    MatrixXd G1;
    G1 = A1.inverse();

    cout << "The inverse of the matrix A1 is" << G1;
    cout << "Double Check" << G1 * A1;
}