#include<iostream>
#include <fstream>
#include <Eigen/Dense>
#include <vector>
#include <cmath>
#include <iomanip>
#include <string>
#include <const.h>

using namespace std;
using namespace Eigen;
using namespace dlib;

vector<double> linspace(const double& min, const double& max, int n)
    {
        vector<double> result;
        // vector iterator
        int iterator = 0;

        for (int i = 0; i <= n - 2; i++)
        {
            double temp = min + i * (max - min) / (floor((double)n) - 1);
            result.insert(result.begin() + iterator, temp);
            iterator += 1;
        }

        //iterator += 1;

        result.insert(result.begin() + iterator, max);
        return result;
    }

vector<double> omega_arr = linspace(1,10,10);
vector<double> tilde_g_arr = linspace(1,10,10);
vector<double> tau_arr = linspace(0,1,100);
double re_planck_cst = planck_cst/(2*pi);

void tilde_g_calculation_function(double alpha, double k_cutoff)
{
    double nu = pi * k_cutoff / alpha;

    for (int i=0; i<omega_arr.size(); i++)
    {
        omega_arr[i] = k_cutoff * (omega_arr[i]/omega_arr[omega_arr.size()-1]);
        tilde_g_arr[i] = sqrt((2 * k_cutoff / (alpha * omega_arr.size())) * (omega_arr[i] / (1 + pow(nu * omega_arr[i] / k_cutoff,2))));
        //tilde_g_arr[i] = sqrt( (omega_arr[i] / (1 + pow(nu * omega_arr[i] / k_cutoff,2))));
        //tilde_g_arr[i] = sqrt((2 * k_cutoff / (alpha * omega_arr.size())) * (re_planck_cst * omega_arr[i] / (1 + pow(nu * re_planck_cst * omega_arr[i] / k_cutoff,2))));
    }
}

vector<double> Interact_V(vector<double> tau)
{
    vector<double> V_arr(tau.size(), 0);

    for (int i = 0; i < tau.size(); i++)
    {
        for (int j = 0; j < tilde_g_arr.size();j++)
        {
            V_arr[i] += pow(tilde_g_arr[j],2) * cosh((tau[i] - tau[tau.size() - 1] / 2) * omega_arr[j])/sinh(tau[tau.size() - 1] * omega_arr[j] / 2);
            //cout << "\t" << j <<" V_arr : " << V_arr[i] << " with tau-beta/2 : " << tau[i] - tau[tau.size()-1]/2 << endl;
        }
        cout << "this is V_arr " << V_arr[i] << endl;
        cout << setprecision(16);
    }

    return V_arr;
}

int main()
{
    double k_cutoff = 20;
    for (int i=0; i<omega_arr.size(); i++)
    {
        cout << omega_arr[i] << endl;
    }

    tilde_g_calculation_function(0.5,k_cutoff);
    vector<double> a = Interact_V(tau_arr);

    std::ofstream outputFile;
    string name = "coupTEST";
    name += ".txt";

    outputFile.open(name);

    for (int j = 0; j < tau_arr.size(); j++)
    {
        cout << a[j] << endl;
        outputFile << tau_arr[j] << "\t" << a[j] << endl;
    }
    

    outputFile.close();
/*    for (int i=0; i<omega_arr.size(); i++)
    {
        cout << g_arr[i] << endl;
    }
*/

    return 0;
}
