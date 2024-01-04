#include<iostream>
#include<Eigen/Dense>
using namespace std;
using namespace Eigen;

int main()
{
    MatrixXd A(3,3);

    A << 1, 2, 3,
        4, 5, 6,
        7, 8, 9;

    SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver(A);
    
    if (eigensolver.info() == Success)
    {
        VectorXd eigenvalues = eigensolver.eigenvalues();
        MatrixXd eigenvectors = eigensolver.eigenvectors();

        cout << "Eigenvaules: " << eigenvalues.transpose() << endl;
        cout << "Eigenvectors: \n" << eigenvectors << endl;
    }
    else
    {
        cerr << "Eigenvalue decomposition failed!" << endl;    
    }
    

}