#include<iostream>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <vector>
#include <cmath>
#include <OCA_bath.hpp>
#include <chrono>

using namespace std;
using namespace Eigen;

double g_ma = 1;
int siz = 0;
int sys = 0;

//////////////////////////////////////////////////////////////////////////////////////

MD_OC::MD_OC(double beta, int grid)
     : tau_grid(linspace(0,beta,grid)) , t(grid-1)
{
    mode_grid = linspace(1,30000,30000);

    Delta_t = tau_grid[1] - tau_grid[0];

    M = mode_grid.size();
    t = tau_grid.size();
    H_N = MatrixXd::Zero(3,3);
    H_loc = MatrixXd::Zero(3,3);
    
    coup_Arr.resize(M);
    omega_Arr.resize(M);
    INT_Arr.resize(t);

}

MD_OC::~MD_OC()
{
    //blank;
}

////////////////////////////////////////////////////////////////////////////////////

MatrixXd MD_OC::Eigenvector_Even()
{
    SelfAdjointEigenSolver<MatrixXd> es(Matrix_Even(siz, g_ma));
    return es.eigenvectors();
}

MatrixXd MD_OC::Eigenvalue_Even()
{
    SelfAdjointEigenSolver<MatrixXd> es(Matrix_Even(siz, g_ma));
    return es.eigenvalues();
}

MatrixXd MD_OC::Eigenvector_Odd()
{
    SelfAdjointEigenSolver<MatrixXd> es(Matrix_Odd(siz, g_ma));
    return es.eigenvectors();
}

MatrixXd MD_OC::Eigenvalue_Odd()
{
    SelfAdjointEigenSolver<MatrixXd> es(Matrix_Odd(siz, g_ma));
    return es.eigenvalues();
}

int main()
{
    double beta = 1;
    int grid = 1;

    MD_OC MD(beta,grid);

    std::ofstream outputFile ("./");

    int& size = siz;
    int& syst = sys;
    double& gamm = g_ma;

    cout << " * Set System size : ";
    cin >> syst;

    cout << " * Set calculation size : ";
    cin >> size;

    if (size > syst){
        cout << "**************** Program will shutdown *******************" << endl;
        exit(1);
    }

    vector<double> g_ma_arr(5,0);
    for (int i = 0; i < 5 ; i++)
    {
        if (i==0)
        {
            g_ma_arr[i] = 0;
        }
        if (i!=0)
        {
            g_ma_arr[i] = g_ma_arr[i-1] + 0.2;
        }
    }

    for (int ga=0; ga<g_ma_arr.size();ga++)
    {
        gamm = g_ma_arr[ga];

        string Eig_name = "EIGENVEC_EVE_GAM_";

        std::stringstream gam;
        std::stringstream sy;
        std::stringstream si;

        gam << g_ma;
        sy << sys;
        si << siz;

        Eig_name += gam.str();
        Eig_name += "_SYS_";
        Eig_name += sy.str();
        Eig_name += "_SIZ_";
        Eig_name += si.str();
        Eig_name += ".txt";

        outputFile.open(Eig_name);
        for (int j=0; j < siz; j++){
            for (int i = 0; i < siz; i++)
            {
                outputFile << MD.Eigenvector_Even()(i,j) << "\t";
            }
            outputFile << "\n";
        }
        outputFile.close();

        string Oig_name = "EIGENVEC_ODD_GAM_";

        Oig_name += gam.str();
        Oig_name += "_SYS_";
        Oig_name += sy.str();
        Oig_name += "_SIZ_";
        Oig_name += si.str();
        Oig_name += ".txt";

        outputFile.open(Oig_name);
        for (int j=0; j < siz; j++){
            for (int i = 0; i < siz; i++)
            {
                outputFile << MD.Eigenvector_Odd()(i,j) << "\t";
            }
            outputFile << "\n";
        }
        outputFile.close();
    }



}