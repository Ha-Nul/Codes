#include <iostream>
#include <fstream>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <vector>
#include <cmath>
#include <OCA_bath.hpp>
#include <string>
#include <sstream>

using namespace std;
using namespace Eigen;

///////////////////////////////////////////////////////////////////////

double g_ma;
double alpha;
int sys;

///////////////////////////////////////////////////////////////////////

MD_OC::MD_OC(double beta, int grid)
     : tau_grid(linspace(0,beta,grid)) , t(grid-1)
{
    mode_grid = linspace(1,30000,30000);

    Delta_t = tau_grid[1] - tau_grid[0];

    M = mode_grid.size();
    t = tau_grid.size();
    
    coup_Arr.resize(M);
    omega_Arr.resize(M);
    INT_Arr.resize(t);
    //vector<double> k_mode(100,1);

}

MD_OC::~MD_OC()
{
    //blank
}

/*
vector<double> MD_OC::green()
{
    vector<double> bose_dist(t);

    for (int i = 0; i < t; i++)
    {
        bose_dist[i] = 1 / (exp(tau_grid[t - 1] * k_mode[i]) - 1);
    }

    vector<double> green(t);

    for (int j = 0; j < t; j++)
    {
        green[j] = ((bose_dist[j] + 1) * exp(-1 * k_mode[j] * tau_grid[j]) + (bose_dist[j]) * exp(k_mode[j] * tau_grid[j]));
    }

    return green;
}
*/

void MD_OC::Tilde_g_calculation_function(double alpha, double k_cutoff)
{
    double nu = pi * k_cutoff / alpha;
    //Initializing Array
    for (int j = 0; j < M; j++)
    {
        omega_Arr[j] = 0;
        coup_Arr[j] = 0;
    }

    for (int i = 0; i < M; i++)
    {
        //omega_Arr[i] = k_cutoff * (mode_grid[i]/mode_grid[M-1]);
        //coup_Arr[i] = sqrt((2 * k_cutoff / (alpha * M)) * (omega_Arr[i] / (1 + pow(nu * omega_Arr[i] / k_cutoff,2))));

        //simpson formulae
        omega_Arr[0] = -0.05;
        omega_Arr[i] = (mode_grid[i] / mode_grid[M - 1]); // fix to x to adjust simpson's rule
        coup_Arr[i] = sqrt((2 * k_cutoff / (alpha)) * (k_cutoff * omega_Arr[i] / (1 + pow(nu * omega_Arr[i], 2)))); // fix to adjust simpson's rule
    }

    if (alpha == 0)
    {
        for (int i = 0; i < M; i++)
        {
            coup_Arr[i] = 0;
        }
    }
}

////////////////////////////////////////////////////////////////////////////////////

void MD_OC::Interact_V(double k_cutoff)
{
    //Initializing Interaction array
    for (int i = 0; i < t; i++)
    {
        INT_Arr[i] = 0;
    }

    for (int i = 0; i < t; i++)
    {
        for (int j = 0; j < M; j++)
        {
            //INT_Arr[i] += -pow(coup_Arr[j],2) * cosh((tau_grid[i] - tau_grid[t - 1] / 2) * omega_Arr[j])/sinh(tau_grid[t - 1] * omega_Arr[j] / 2); //caution for sign

            //simpson formulae

            if (j == 0 || j == M - 1)
            {
                INT_Arr[i] += -(1.0 / (3 * M)) * pow(coup_Arr[j], 2) * cosh((tau_grid[i] - tau_grid[t - 1] / 2) * k_cutoff * omega_Arr[j]) / sinh(tau_grid[t - 1] * k_cutoff * omega_Arr[j] / 2);
            }

            else if (j % 2 != 0)
            {
                INT_Arr[i] += -(1.0 / (3 * M)) * 4 * pow(coup_Arr[j], 2) * cosh((tau_grid[i] - tau_grid[t - 1] / 2) * k_cutoff * omega_Arr[j]) / sinh(tau_grid[t - 1] * k_cutoff * omega_Arr[j] / 2);
            }

            else if (j % 2 == 0)
            {
                INT_Arr[i] += -(1.0 / (3 * M)) * 2 * pow(coup_Arr[j], 2) * cosh((tau_grid[i] - tau_grid[t - 1] / 2) * k_cutoff * omega_Arr[j]) / sinh(tau_grid[t - 1] * k_cutoff * omega_Arr[j] / 2);
            }

        }
        //INT_Arr[i] += -0.05;
    }
}


void MD_OC::CAL_COUP_INT_with_g_arr(double alpha, double k_cutoff)
{
    Tilde_g_calculation_function(alpha,k_cutoff);
    Interact_V(k_cutoff);
}

///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////

void MD_OC::Dataoutput(double g_ma, double alpha)
{
    std::ofstream outputFile("./");
    std::stringstream gam;
    std::stringstream alp;

    gam << g_ma;
    alp << alpha;

    string INT= "INT_Arr_g";

    INT += gam.str();
    INT += "_a";
    INT += alp.str();
    INT += ".dat";

    outputFile.open(INT);
    for (int i = 0; i < t; i++)
    {
        outputFile << tau_grid[i] << "\t" << INT_Arr[i] << endl;
    }
    outputFile.close();

}
/////////////////////////////////////////////////////////////////////////////////

int main()
{
    double beta = 10;
    int grid = 101;

    MD_OC OC(beta,grid);

    double& ref_g_ma = g_ma;
    double& alp = alpha;
    int& syst = sys;
    double k_cutoff = 20;

    syst = 21;

    string input;
    getline(cin, input); // 한 줄을 입력받음

    istringstream iss(input);
    iss >> g_ma >> alpha;
    cout << " Value of g_ma : " << g_ma << ", alpha : " << alpha << endl;

    OC.CAL_COUP_INT_with_g_arr(alpha,k_cutoff);
    OC.Dataoutput(g_ma,alpha);

    return 0;
}

