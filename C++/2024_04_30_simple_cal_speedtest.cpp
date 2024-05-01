#include<iostream>
#include <fstream>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <vector>
#include <cmath>
#include <iomanip>
#include <string>
#include <chrono>

using namespace std;

vector<double> linspace(const double &min,const double &max, int n)
        {
            vector<double> result;
            // vector iterator
            int iterator = 0;

            for (int i = 0; i <= n-2; i++)	
            {
                double temp = min + i*(max-min)/(floor((double)n) - 1);
                result.insert(result.begin() + iterator, temp);
                iterator += 1;
            }

            //iterator += 1;

            result.insert(result.begin() + iterator, max);
            return result;
        }

int main()
{
    vector<double> sample_arr = linspace(0,1,10);
    vector<double> asmple_arr = linspace(0,1,10);

    for (int i = 0; i < sample_arr.size(); i++)
    {
        std::chrono::system_clock::time_point start= std::chrono::system_clock::now();
        cout << "Iteration " << i << " Starts" << endl;

        cout << sample_arr[i] * asmple_arr[i] << endl;

        std::chrono::system_clock::time_point sec = std::chrono::system_clock::now();
        std::chrono::duration<double> nanoseconds = std::chrono::duration_cast<std::chrono::nanoseconds>(sec-start);
        cout << "Process ends in : " << nanoseconds.count() << "[sec]" << endl;
        cout << "-----------------------------" << endl;
    }

    return 0;

}