#include<iostream>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <vector>
#include <cmath>
#include <OCA_Function_def.hpp>

using namespace std;
using namespace Eigen;

MAIN_DEF MD;

int MAIN_DEF::k = MD.tau_grid.size();
vector<double> k_mode(100, 1);
double g_ma = 1;


vector<double> MAIN_DEF::green(vector<double> tau)
{
    double T = 273;
    vector<int> one_vec(k, 1);
    vector<double> bose_dist(k);

    for (int i = 0; i < k; i++)
    {
        bose_dist[i] = one_vec[i] / (exp(tau_grid[k - 1] * k_mode[i]) - 1);
    }

    vector<double> Test_green(k);

    for (int j = 0; j < k; j++)
    {
        Test_green[j] = ((bose_dist[j] + 1) * exp(-1 * k_mode[j] * tau[j]) + (bose_dist[j]) * exp(k_mode[j] * tau[j]));
    }

    return Test_green;
}

vector<double> MAIN_DEF::coupling(double v, double g, double W)
{
    vector<double> v_array(k_mode.size(), v);
    vector<double> g_array(k_mode.size(), g);
    vector<double> W_array(k_mode.size(), W);
    vector<double> coupling_array(k_mode.size());

    for (int i = 0; i < k_mode.size(); i++)
    {
        coupling_array[i] = g_array[i] * sqrt(abs(k_mode[i]) * v_array[i] / (1 + pow((abs(k_mode[i]) * v_array[i] / W_array[i]), 2)));
    }

    return coupling_array;
}
////////////////////////////////////////////////////////////////////////////////////

vector<double> MAIN_DEF::Interact_V(vector<double>coupling, vector<double> tau, double omega)
{
    double coupling_const = coupling[0];

    vector<double> hpcos(tau.size(), 0);
    vector<double> hpsin(tau.size(), 0);
    vector<double> coupling_arr(tau.size(), coupling_const * coupling_const);
    vector<double> V_arr(tau.size(), 0);

    for (int i = 0; i < tau.size(); i++)
    {
        hpcos[i] = cosh(tau[i] - tau[tau.size() - 1] / 2) * omega;
        hpsin[i] = sinh(tau[tau.size() - 1] * omega / 2);
        V_arr[i] = (coupling_arr[i] * hpcos[i] / hpsin[i]);

        //cout << "this is V_arr " << V_arr[i] << endl;
    }

    return V_arr;
}

////////////////////////////////////////////////////////////////////////////////////

MatrixXd MAIN_DEF::Eigenvector_Even()
{
    MatrixXd a;

    SelfAdjointEigenSolver<MatrixXd> es(Matrix_Even(3, g_ma));
    a = es.eigenvectors();

    return a;
}

MatrixXd MAIN_DEF::Eigenvalue_Even()
{
    MatrixXd b;

    SelfAdjointEigenSolver<MatrixXd> es(Matrix_Even(3, g_ma));
    b = es.eigenvalues();

    return b;
}

MatrixXd MAIN_DEF::Eigenvector_Odd()
{
    MatrixXd a;

    SelfAdjointEigenSolver<MatrixXd> es(Matrix_Odd(3, g_ma));
    a = es.eigenvectors();

    return a;
}

MatrixXd MAIN_DEF::Eigenvalue_Odd()
{
    MatrixXd b;

    SelfAdjointEigenSolver<MatrixXd> es(Matrix_Odd(3, g_ma));
    b = es.eigenvalues();

    return b;
}

///////////////////////////////////////////////////////////////////////


MatrixXd MAIN_DEF::Hamiltonian_N(MatrixXd even, MatrixXd odd, double g)
{
    MatrixXd odd_eigenvec;
    MatrixXd even_eigenvec;

    odd_eigenvec = odd.transpose();
    even_eigenvec = even;

    MatrixXd c;
    c = odd_eigenvec * even_eigenvec;

    MatrixXd d = MatrixXd::Zero(3, 3);

    d(0, 1) = g * c(0, 0);
    d(1, 0) = g * c(0, 0);
    d(1, 2) = g * c(0, 1);
    d(2, 1) = g * c(0, 1);

    return d;
}

vector<MatrixXd> MAIN_DEF::Hamiltonian_exp(MatrixXd a, MatrixXd b)
{
    //g_0 
    MatrixXd Even = a;
    MatrixXd Odd = b;

    double zeroth = exp(Even(0));
    double first = exp(Odd(0));
    double second = exp(Even(1));

    vector<MatrixXd> array_with_Matrix(k);

    MatrixXd Hamiltonian_exp;

    for (int i = 0; i < k; i++)
    {
        Hamiltonian_exp = MatrixXd::Zero(3, 3);

        Hamiltonian_exp(0, 0) = tau_grid[i] * zeroth;
        Hamiltonian_exp(1, 1) = tau_grid[i] * first;
        Hamiltonian_exp(2, 2) = tau_grid[i] * second;

        array_with_Matrix[i] = Hamiltonian_exp;
    }

    return array_with_Matrix;
}



MatrixXd MAIN_DEF::Hamiltonian_loc(MatrixXd a, MatrixXd b)
{
    MatrixXd Hamiltonian = MatrixXd::Zero(3, 3);

    Hamiltonian(0, 0) = a(0);
    Hamiltonian(1, 1) = b(0);
    Hamiltonian(2, 2) = a(1);

    return Hamiltonian;
}

////////////////////////////////////////////////////////////////////////////////

