#include<iostream>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <vector>
#include <cmath>
#include <OCA_Function_def.hpp>

using namespace std;
using namespace Eigen;

///////////////////////////////////////////////////////////////////////////////

double omega = 1;
double velocity = 1;
double cutoff = 1;

int k = 100;

vector<double> INT_Arr(k, 0);
vector<double> Chi_sp(k, 0);
vector<MatrixXd> SELF_E(k, MatrixXd::Zero(3, 3));
MatrixXd MAIN_DEF::H_N;


void MAIN_DEF::CAL_COUP_INT_with_g_arr(double g)
{
    INT_Arr = Interact_V(coupling(velocity, g, cutoff), tau_grid, omega);
    H_N = Hamiltonian_N(Eigenvector_Even(), Eigenvector_Odd(), g);
}


////////////////////////////////////////////////////////////////////////////////


void MAIN_DEF::NCA_self(const MatrixXd& N, const vector<MatrixXd>& Prop, const vector<double>& V)
{
    for (int i = 0; i < k; i++)
    {
        SELF_E[i] = V[i] * (N * Prop[i] * N);
    }
}

void MAIN_DEF::OCA_self(MatrixXd& N, vector<MatrixXd>& Prop, vector<double>& V)
{
    MatrixXd Stmp;

    for (int i = 0; i < k; i++)
    {
        Stmp = MatrixXd::Zero(3, 3);
        for (int n = 0; n <= i; n++) for (int m = 0; m <= n; m++) {
            Stmp += N * Prop[i - n] * N * Prop[n - m] * N * Prop[m] * N * V[i - n] * V[m];
        }
        SELF_E[i] += pow(Delta_t, 2) * Stmp;
    }
}

void MAIN_DEF::SELF_Energy(vector<MatrixXd> Prop)
{
    NCA_self(H_N, Prop, INT_Arr);
    OCA_self(H_N, Prop, INT_Arr);
}


//////////////////////////////////////////////////////////////////////////////


MatrixXd MAIN_DEF::round_propagator_ite(const MatrixXd& loc, const vector<MatrixXd>& sigma, const vector<MatrixXd>& ite, int n, int boolean)
{

    MatrixXd sigsum = MatrixXd::Zero(3, 3);

    if (n == 1)
    {
        sigsum = 0.5 * Delta_t * (sigma[1] * ite[0] + sigma[0] * ite[1]);
    }
    else if (n > 1) {
        for (int i = 0; i < n; i++)
        {
            sigsum += 0.5 * Delta_t * (sigma[n - (i)] * ite[i] + sigma[n - (i + 1)] * ite[i + 1]);

            if (i + 1 == n)
            {
                break;
            }

        }
    }

    //cout << sigsum << endl;

    MatrixXd Bucket = MatrixXd::Zero(3, 3);
    if (boolean == 0)
    {
        Bucket = -loc * ite[n] + sigsum;
    }
    else if (boolean == 1)
    {
        Bucket = sigsum;
    }
    //cout << -loc * ite << endl;
    return Bucket;
}



vector<MatrixXd> MAIN_DEF::Propagator(const vector<MatrixXd>& sig, const MatrixXd& loc)
{
    vector<MatrixXd> P_arr(k, MatrixXd::Zero(3, 3));
    vector<MatrixXd> S_arr(k, MatrixXd::Zero(3, 3));

    P_arr[0] = MatrixXd::Identity(3, 3);
    S_arr[0] = MatrixXd::Identity(3, 3);

    MatrixXd sig_form = MatrixXd::Zero(3, 3);
    MatrixXd sig_late = MatrixXd::Zero(3, 3);

    for (int i = 1; i < k; i++)
    {
        P_arr[1] = P_arr[0];
        sig_late = 0.5 * Delta_t * (0.5 * Delta_t * (sig[1] * P_arr[0] + sig[0] * (P_arr[0] + Delta_t * P_arr[0])));
        P_arr[1] = P_arr[0] - 0.5 * Delta_t * loc * (2 * P_arr[0] + Delta_t * P_arr[0]) + sig_late;
        S_arr[1] = P_arr[1];

        if (i > 1)
        {
            sig_form = round_propagator_ite(loc, sig, P_arr, i - 1, 0);
            S_arr[i] = P_arr[i - 1] + Delta_t * sig_form;

            sig_late = 0.5 * Delta_t * (round_propagator_ite(loc, sig, P_arr, i - 1, 1) + round_propagator_ite(loc, sig, S_arr, i, 1));
            P_arr[i] = P_arr[i - 1] - 0.5 * Delta_t * loc * (2 * P_arr[i - 1] + Delta_t * sig_form) + sig_late;

        }
    }

    return P_arr;
}

/////////////////////////////////////////////////////////////////////////////

double MAIN_DEF::chemical_poten(MatrixXd prop)
{
    double Trace = prop.trace();
    double lambda = -(1 / tau_grid[k - 1]) * log(Trace);

    return lambda;
}

///////////////////////////////////////////////////////////////////////////////

vector<MatrixXd> MAIN_DEF::Iteration(const int& n)
{
    vector<MatrixXd> Prop(k, MatrixXd::Zero(3, 3));

    vector<MatrixXd> H_loc(n + 1, MatrixXd::Zero(3, 3));
    H_loc[0] = Hamiltonian_loc(Eigenvalue_Even(), Eigenvalue_Odd());

    MatrixXd Iden = MatrixXd::Identity(3, 3);

    vector<double> lambda(n + 1, 0);
    double expDtauLambda;
    double factor;

    for (int i = 0; i <= n; i++)
    {
        if (i == 0)
        {
            for (int j = 0; j < k; j++)
            {
                Prop[j](0, 0) = exp(-tau_grid[j] * Hamiltonian_loc(Eigenvalue_Even(), Eigenvalue_Odd())(0, 0));
                Prop[j](1, 1) = exp(-tau_grid[j] * Hamiltonian_loc(Eigenvalue_Even(), Eigenvalue_Odd())(1, 1));
                Prop[j](2, 2) = exp(-tau_grid[j] * Hamiltonian_loc(Eigenvalue_Even(), Eigenvalue_Odd())(2, 2));

                cout << Hamiltonian_loc(Eigenvalue_Even(), Eigenvalue_Odd())(0, 0) << endl;
            }


            lambda[0] = chemical_poten(Prop[k - 1]);
            expDtauLambda = exp((tau_grid[1] - tau_grid[0]) * lambda[0]);
            factor = 1.0;


            for (int j = 0; j < k; j++)
            {
                Prop[j] *= factor;
                factor *= expDtauLambda;
                //cout << Prop[j] << endl;
            }
        }

        else
        {

            H_loc[i] = H_loc[i - 1] - lambda[i - 1] * Iden;
            SELF_Energy(Prop);
            Prop = Propagator(SELF_E, H_loc[i]);

            lambda[i] = chemical_poten(Prop[k - 1]);

            expDtauLambda = exp((tau_grid[1] - tau_grid[0]) * lambda[i]);
            factor = 1.0;

            for (int j = 0; j < k; j++)
            {
                Prop[j] *= factor;
                factor *= expDtauLambda;

                //cout << Prop[j] << endl;
            }

        }

    }

    return Prop;
}

//////////////////////////////////////////////////////////////////////////////

void MAIN_DEF::NCA_Chi_sp(vector<MatrixXd> iter)
{
    MatrixXd GELL_1 = MatrixXd::Zero(3, 3);
    GELL_1(0, 1) = 1;
    GELL_1(1, 0) = 1;

    for (int i = 0; i < k; i++)
    {
        Chi_sp[i] = (iter[k - i - 1] * GELL_1 * iter[i] * GELL_1).trace();
    }
}

void MAIN_DEF::OCA_Chi_sp(vector<MatrixXd> iter)
{
    MatrixXd GELL_1 = MatrixXd::Zero(3, 3);
    GELL_1(0, 1) = 1;
    GELL_1(1, 0) = 1;

    for (int i = 0; i < k; i++)
    {
        MatrixXd Stmp = MatrixXd::Zero(3, 3);

        for (int n = 0; n <= i; n++) for (int m = i; m < k; m++)
        {
            Stmp += INT_Arr[m - n] * iter[k - m - 1] * H_N * iter[m - i] * GELL_1 * iter[i - n] * H_N * iter[n] * GELL_1;
        }

        Chi_sp[i] += pow(Delta_t, 2) * Stmp.trace();
    }
}

vector<double> MAIN_DEF::Chi_sp_Function(vector<MatrixXd> ITE)
{
    NCA_Chi_sp(ITE);
    OCA_Chi_sp(ITE);
    
    return Chi_sp;
    
}
////////////////////////////////////////////////////////////////////////////////////

