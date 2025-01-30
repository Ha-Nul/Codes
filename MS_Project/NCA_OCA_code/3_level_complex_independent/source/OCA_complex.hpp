#include <iostream>
#include <fstream>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <vector>
#include <cmath>
#include <iomanip>
#include <string>
#include <complex>  // Include for complex numbers
#include <const.h>

using namespace std;
using namespace Eigen;

class MD_OC //MAIN_DEF_OCA
{
public:

    MD_OC(double beta, int grid);
    ~MD_OC();

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

    double Limit;
    void SetLimit(double value);
    void Setgrid();

    ////////////////////////////////////////////////////////////////////////////////

    vector<double> alp_arr;
    vector<double> gam_arr;

    ////////////////////////////////////////////////////////////////////////////////

    vector<double> tau_grid;
    vector<double> mode_grid;

    int M;
    int t;
    //int siz;
    //int sys;
    //double g_ma;

    double Delta_t;

    double pi = dlib::pi;
    //vector<double> green();   
    void Interact_V(double k_cutoff);

    MatrixXcd Eigenvector_Even();
    MatrixXcd Eigenvalue_Even();
    MatrixXcd Eigenvector_Odd();
    MatrixXcd Eigenvalue_Odd();


    void Hamiltonian_N(MatrixXcd even, MatrixXcd odd);
    void Hamiltonian_loc(MatrixXcd a, MatrixXcd b);

    void CAL_COUP_INT_with_g_arr(double alp, double cutoff);
    void Tilde_g_calculation_function(double alpha, double k_cutoff);

    void Dataoutput(double gamma, double alpha);


    ///////////////////////////////////////////////////////////////////////////////////////////////

    vector<double> coup_Arr;
    vector<double> omega_Arr;
    vector<double> INT_Arr;
    //vector<double> k_mode;

    MatrixXcd H_N;

    ///////////////////////////////////////////////////////////////////////////////////////////////

    void Ordercal(MatrixXcd even, MatrixXcd odd);
    MatrixXcd Order_param;
    MatrixXcd H_loc;

    vector<complex<double> > Chi_Arr;
    vector<vector<MatrixXcd> > T;
    vector<vector<MatrixXcd> > Chi_st;
    vector<MatrixXcd> SELF_E;
    vector<MatrixXcd> Prop;
    vector<MatrixXcd> GELL;

    ////////////////////////////////////////////////////////////////////////////////////////////////

    void readVfunc();
    MatrixXcd Make_N_Matrix(int lineNumber);
    MatrixXcd Make_Loc_Matrix(int lineNumber);
    void data_store(int lineNumber);

    ////////////////////////////////////////////////////////////////////////////////////////////////

    void NCA_self();
    void OCA_self();
    void OCA_T();
    void SELF_Energy();

   MatrixXcd round_propagator_ite(const MatrixXcd& loc, const vector<MatrixXcd>& sigma, const vector<MatrixXcd>& ite, int n, int boolean);
    vector<MatrixXcd> Propagator(const vector<MatrixXcd>& array, const MatrixXcd& loc);

    complex<double> chemical_poten(MatrixXcd prop);

    vector<MatrixXcd> Iteration(const int& iteration);

    //////////////////////////////////////////////////////////////////////////////////////////////

    double temp_minpoint(vector<MatrixXcd>& arr);
    vector<double> temp_itemin(vector<MatrixXcd>& arr, double minpo, int size);

    //////////////////////////////////////////////////////////////////////////////////////////////

    void NCA_Chi_sp(vector<MatrixXcd>& ITER);
    void OCA_store(vector<MatrixXcd>& ITER);
    void OCA_Chi_sp(vector<MatrixXcd>& ITER);
    vector<complex<double> > Chi_sp_Function(vector<MatrixXcd> Iter);
};

/////////////////////////////////////////////////////////////////////////////////////