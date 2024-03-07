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

vector<double> g_array = linspace(0,10,10);

void g_calculation_function(double alpha, double k_cutoff, double mode)
{
    double nu = pi * k_cutoff / alpha;

    for (int i=0; i<g_array.size(); i++)
    {
        g_array[i] = (1/planck_cst) * sqrt((2 * planck_cst * k_cutoff / alpha * mode) * (planck_cst * g_array[i] / 1 + pow(nu * planck_cst * g_array[i] / k_cutoff,2)));
    }
}

int main()
{
    g_calculation_function(1,10,10);

    for(int i=0; i<g_array.size(); i++)
    {
        cout << g_array[i] << endl;
    }
}