#include<iostream>
#include <fstream>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <vector>
#include <cmath>
#include <iomanip>
#include <string>
#include <OCA_Function_def.hpp>

using namespace std;
using namespace Eigen;

int main()
{
    MAIN_DEF MD;

    vector<double> g_array(25, 0);
    for (int j = 1; j < 25; ++j)
    {
        if (j < 21)
        {
            g_array[j] = (g_array[j - 1] + 0.05);
        }

        else
        {
            g_array[j] = g_array[j - 1] + 1;
        }
    }

    for (int m = 0; m < 21; m++)
    {
        g_array[m] = g_array[m] * g_array[m];
    }
    
    for (int i = 0; i < 1; i++)
    {
        MD.CAL_COUP_INT_with_g_arr(1);
        vector<MatrixXd> ITER = MD.Iteration(3);
        vector<double> a = MD.Chi_sp_Function(ITER);

        std::ofstream outputFile;

        //string name = "20240111_Trap_beta_0_4_g_";
        string name = "WINDOWTEST_100_ite3";
        //std::stringstream back;
        //back << g_array[k];

        //name += back.str();
        name += ".txt";

        outputFile.open(name);
        //vector<double> a = test.Interact_V(test.coupling(velocity,g_array[k],cutoff),test.grid,omega);

        for (int j = 0; j < a.size(); j++)
        {
            cout << a[j] << endl;
            outputFile << MD.tau_grid[j] << "\t" << a[j] << endl;
        }
        outputFile.close();
    }
    

    return 0;

}
