#include <iostream>
#include <cmath> // You included <math.h>, but cmath is preferred in C++
#include <vector>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

int main() {
    MatrixXcd Innerprod(2,2);  // Initialize matrix size
    MatrixXcd Transpose(2,2); // Initialize matrix size

    Innerprod << complex<double>(1, 1), complex<double>(0, 0),
                 complex<double>(0, 0), complex<double>(1, -1);

    Transpose = Innerprod.adjoint();

    cout << Innerprod * Transpose << endl;

    return 0;
}