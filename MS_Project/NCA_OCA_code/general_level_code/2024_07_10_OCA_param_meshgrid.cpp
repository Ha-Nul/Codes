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

/////////////////////////////////////////////////////////////

MD_OC::MD_OC(double beta, int grid)
    : tau_grid(linspace(0, beta, grid)), t(grid - 1)
{
    mode_grid = linspace(1, 30000, 30000);

    Delta_t = tau_grid[1] - tau_grid[0];

    M = mode_grid.size();
    t = tau_grid.size();
    H_N = MatrixXd::Zero(3, 3);
    H_loc = MatrixXd::Zero(3, 3);

    coup_Arr.resize(M);
    omega_Arr.resize(M);
    INT_Arr.resize(t);

}

MD_OC::~MD_OC()
{
    //blank;
}


//////////////////////////////////////////////////////////////

void MD_OC::Tilde_g_calculation_function(double alpha, double k_cutoff)
{
    //Initializing data
    for (int j = 0; j < M; j++)
    {
        omega_Arr[j] = 0;
        coup_Arr[j] = 0;
    }

    double nu = pi * k_cutoff / alpha;

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



////////////////////////////////////////////////////////////////////////////////////

MatrixXd MD_OC::Eigenvector_Even()
{
    SelfAdjointEigenSolver<MatrixXd> es(Matrix_Even(sys, g_ma));
    return es.eigenvectors();
}

MatrixXd MD_OC::Eigenvalue_Even()
{
    SelfAdjointEigenSolver<MatrixXd> es(Matrix_Even(sys, g_ma));
    return es.eigenvalues();
}

MatrixXd MD_OC::Eigenvector_Odd()
{
    SelfAdjointEigenSolver<MatrixXd> es(Matrix_Odd(sys, g_ma));
    return es.eigenvectors();
}

MatrixXd MD_OC::Eigenvalue_Odd()
{
    SelfAdjointEigenSolver<MatrixXd> es(Matrix_Odd(sys, g_ma));
    return es.eigenvalues();
}

///////////////////////////////////////////////////////////////////////

void MD_OC::Hamiltonian_N(MatrixXd even, MatrixXd odd)
{
    //cout << "input g value :" << g << endl;
    MatrixXd INT_odd = MatrixXd::Zero(sys, sys);
    MatrixXd INT_even = MatrixXd::Zero(sys, sys);

    //H_N initialize
    H_N = MatrixXd::Zero(siz, siz);

    //cout << "initialized check" << endl;
    //cout << H_N << "\n" << endl;

    for (int i = 0; i < sys; i++) for (int j = 0; j < sys; j++)
    {
        INT_even(i, j) = -1 * even(i, j) * i; // -\sum_1^\infty \alpha_i \sin{i\phi}

        if (i < sys - 1)
        {
            INT_odd(i + 1, j) = odd(i, j);
        }

    }

    MatrixXd c = INT_even.transpose() * INT_odd;
    //cout << c << endl;

    //stocks matrix elements
    for (int i = 0; i < siz; i++) for (int j = 0; j < siz; j++)
    {
        if (j % 2 != 0 & i % 2 == 0) {
            H_N(i, j) = c(i / 2, j / 2); // even * diff odd
        }
        else if (j % 2 == 0 & i % 2 != 0) {
            H_N(i, j) = c(j / 2, i / 2); // odd * diff even
        }
    }

    //matching sign
    for (int i = 0; i < siz; i++) for (int j = 0; j < siz; j++)
    {
        if (i > j)
        {
            H_N(i, j) = -H_N(i, j);
        }
    }

    //cout << H_N << endl;
}

void MD_OC::Hamiltonian_loc(MatrixXd a, MatrixXd b)
{
    H_loc = MatrixXd::Zero(siz, siz);

    for (int i = 0; i < siz; i++) for (int j = 0; j < siz; j++)
    {
        if (i == j & i % 2 == 0)
        {
            H_loc(i, j) = a(i / 2);
            /*
            if (a(i / 2) > 30)
            {
                H_loc(i, j) = 30;
            }
            */
        }
        if (i == j & i % 2 != 0)
        {
            H_loc(i, j) = b(i / 2);
            /*
            if (b(i / 2) > 30)
            {
                H_loc(i, j) = 30;
            }
            */
        }
    }
}

void MD_OC::Ordercal(MatrixXd even, MatrixXd odd)
{
    Order_param = MatrixXd::Zero(siz, siz);

    cout << "**EVENMAT**" << endl;
    cout << even << endl;

    //** constructing even matrix
    MatrixXd eve_0 = even;//MatrixXd::Zero(sys,sys);
    for (int i = 0; i < sys; i++) for (int j = 0; j < sys; j++)
    {
        if (i == 0) {
            eve_0(i, j) = (1 / sqrt(2)) * even(i, j);
        }
        else {
            eve_0(i, j) = even(i, j);
        }
    }

    //Even OffDiagonal construct
    MatrixXd eve_1 = MatrixXd::Zero(sys + 1, sys + 1);
    MatrixXd eve_2 = MatrixXd::Zero(sys + 1, sys + 1);

    for (int i = 0; i < sys + 1; i++) for (int j = 0; j < sys + 1; j++)
    {
        if (i > 0 && j < sys) {
            if (i > 1) {
                eve_1(i, j) = 0.5 * eve_0(i - 1, j);
            }
            else {
                eve_1(i, j) = eve_0(i - 1, j);
            }
        }
        if (j < sys && i < sys) {
            eve_2(i, j) = eve_0(i, j);
        }
    }
    //Activate to change the direction of groundstate eigenvector
    for (int i = 0; i < sys + 1; i++) {
        eve_1(i, 0) = -eve_1(i, 0);
        eve_2(i, 0) = -eve_2(i, 0);
    }


    MatrixXd eve_off = eve_1.transpose() * eve_2;
    //cout << "EVEOFF" << endl;
    //cout << eve_off << endl;

    //cout << "\t" << "<Mateven 1>" << endl;
    //cout << eve_1.transpose() << endl;
    /*
    for (int i = 0; i < sys+1; i++) for (int j = 0; j< sys+1; j++)
    {
        if (j<sys && i<sys)
        {
            eve_2(i,j) = eve_0(i,j);
        }
    }
    */
    //cout << "\t" << "<Mateven 2>" << endl;
    //cout << eve_2 << endl;
    MatrixXd eve_ele = MatrixXd::Zero(sys, sys);

    for (int i = 0; i < sys; i++) for (int j = 0; j < sys; j++)
    {
        if ((i != j) && (i % 2 == 0) && (j % 2 == 0)) {
            eve_ele(i, j) = eve_off(i / 2, j / 2) + eve_off(j / 2, i / 2);
            //eve_ele(j,i) = eve_off(i/2,j/2) + eve_off(j/2,i/2);
        }

        if ((i == j) && (i % 2 == 0) && (j % 2 == 0)) {
            eve_ele(i, j) = 2 * eve_off(i / 2, j / 2);
            //eve_ele(j,i) = eve_off(i/2,j/2) + eve_off(j/2,i/2);
        }
    }

    //cout << "*****EVEELE*****" << endl;
    //cout << eve_ele << endl;
    ///////////// even matrix construction complete ///////////////

    //constructing odd matrix

    cout << "**ODD Eigen matrix**" << endl;
    cout << odd << endl;
    //odd matrix calculation structure design
    MatrixXd odd_ele1 = MatrixXd::Zero(sys, sys);
    MatrixXd odd_ele2 = MatrixXd::Zero(sys, sys);
    for (int i = 0; i < sys; i++) for (int j = 0; j < sys; j++)
    {
        if (i != 0) {
            odd_ele1(i, j) = odd(i - 1, j);
        }
    }

    for (int i = 0; i < sys; i++) for (int j = 0; j < sys; j++)
    {
        if (i != 0) {
            odd_ele1(i, j) = odd(i - 1, j);
        }
    }

    for (int i = 0; i < sys; i++) for (int j = 0; j < sys; j++)
    {
        if (i != (sys - 1)) {
            odd_ele2(i, j) = odd(i, j);
        }
    }
    //calculation
    MatrixXd odd_ele = MatrixXd::Zero(sys, sys);
    for (int i = 0; i < sys; i++) for (int j = 0; j < sys; j++)
    {
        if ((i == j) && (i % 2 == 1) && (j % 2 == 1)) {
            odd_ele(i, j) = (odd_ele1.transpose() * odd)(i / 2, j / 2);
        }
        else if ((i < j) && (i % 2 == 1) && (j % 2 == 1)) {
            odd_ele(i, j) = 0.5 * (odd_ele1.transpose() * odd + odd_ele2.transpose() * odd_ele1)(i, j);
            odd_ele(j, i) = odd_ele(i, j);
        }
    }

    /*
    cout << "ODD_element 1 is : " << endl;
    cout << odd_ele1 << endl;

    cout << "\n";

    cout << "ODD_element 2 is : " << endl;
    cout << odd_ele2 << endl;
    */

    cout << "Structure check" << endl;
    cout << odd_ele << endl;
    ///////////// odd matrix construction complete ///////////////

    for (int i = 0; i < siz; i++) for (int j = 0; j < siz; j++)
    {
        if (i % 2 == 0 && j % 2 == 0)
        {
            Order_param(i, j) = eve_ele(i, j);
        }

        else if (i % 2 != 0 && j != 0)
        {
            Order_param(i, j) = odd_ele(i, j);
        }

    }

    cout << "Order : \n" << Order_param << endl;

}


///////////////////////////////////////////////////////////////////////

void MD_OC::CAL_COUP_INT_with_g_arr(double alpha, double k_cutoff)
{
    Tilde_g_calculation_function(alpha, k_cutoff);
    Interact_V(k_cutoff);
    Hamiltonian_N(Eigenvector_Even(), Eigenvector_Odd());
    Hamiltonian_loc(Eigenvalue_Even(), Eigenvalue_Odd());
    Ordercal(Eigenvector_Even(), Eigenvector_Odd());


    cout << "$ H_N value : \n " << H_N << endl;
    cout << "$ H_loc value : \n " << H_loc << endl;
}

////////////////////////////////////////////////////////////////////////////////

void MD_OC::NCA_self()
{
    cout << "\t ** NCA RUN" << endl;
    SELF_E.resize(t);
    for (int i = 0; i < t; i++)
    {
        SELF_E[i] = INT_Arr[i] * (H_N * Prop[i] * H_N);
    }
}

void MD_OC::OCA_T()
{

    cout << "\t ** OCA_T RUN" << endl;
    T.resize(t);
    for (int n = 0; n < t; n++)
    {
        T[n].resize(n + 1);
        for (int m = 0; m <= n; m++)
        {
            T[n][m] = H_N * INT_Arr[n] * Prop[n - m] * H_N * Prop[m] * H_N;
        }
    }
}

void MD_OC::OCA_self()
{
    cout << "\t ** OCA_self RUN" << endl;
    MatrixXd Stmp;
    for (int i = 0; i < t; i++)
    {
        int loop = 0;
        Stmp = MatrixXd::Zero(siz, siz);
        for (int n = 0; n <= i; n++) for (int m = 0; m <= n; m++)
        {
            //std::chrono::system_clock::time_point start= std::chrono::system_clock::now();
            //cout << "\t" << "\t" <<  "For loop count : " << loop  << endl;
            /********************main code**************************/
            Stmp += H_N * Prop[i - n] * T[n][m] * INT_Arr[i - m];
            //cout << "innerloop : " << loop << endl;
            //cout << "(" << n << "," << m << ")" << "th for loop \n" << Stmp << endl;
            /*******************************************************/
            //std::chrono::system_clock::time_point sec = std::chrono::system_clock::now();
            //std::chrono::duration<double> nanoseconds = std::chrono::duration_cast<std::chrono::nanoseconds>(sec-start);
            //cout << "\t" << "\t" << "Calculation ends : " << nanoseconds.count() << "[sec]" << endl;
            //cout << "-----------------------------------------------------" << endl;

        }
        SELF_E[i] += pow(Delta_t, 2) * Stmp;
    }

    //cout << "****************total loop count***************** : " << oloop << endl;
}


void MD_OC::SELF_Energy()
{
    cout << "** SELF_Energy RUN" << endl;
    NCA_self();
    OCA_T();
    std::chrono::system_clock::time_point start = std::chrono::system_clock::now();
    cout << "\t" << "OCA calculation Starts" << endl;
    OCA_self();
    std::chrono::system_clock::time_point sec = std::chrono::system_clock::now();
    std::chrono::duration<double> microseconds = std::chrono::duration_cast<std::chrono::milliseconds>(sec - start);
    cout << "\t" << "Calculation ends : " << microseconds.count() << "[sec]" << endl;
    cout << "-----------------------------" << endl;

    //cout << SELF_E[99] << endl;
}

//////////////////////////////////////////////////////////////////////////////


MatrixXd MD_OC::round_propagator_ite(const MatrixXd& loc, const vector<MatrixXd>& sigma, const vector<MatrixXd>& ite, int n, int boolean)
{

    MatrixXd sigsum = MatrixXd::Zero(siz, siz);

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

    MatrixXd Bucket = MatrixXd::Zero(siz, siz);
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



vector<MatrixXd> MD_OC::Propagator(const vector<MatrixXd>& sig, const MatrixXd& loc)
{
    vector<MatrixXd> P_arr(t, MatrixXd::Zero(siz, siz));
    vector<MatrixXd> S_arr(t, MatrixXd::Zero(siz, siz));

    P_arr[0] = MatrixXd::Identity(siz, siz);
    S_arr[0] = MatrixXd::Identity(siz, siz);

    MatrixXd sig_form = MatrixXd::Zero(siz, siz);
    MatrixXd sig_late = MatrixXd::Zero(siz, siz);

    for (int i = 1; i < t; i++)
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

double MD_OC::chemical_poten(MatrixXd prop)
{
    double Trace = prop.trace();
    double lambda = -(1 / tau_grid[t - 1]) * log(Trace);

    return lambda;
}

///////////////////////////////////////////////////////////////////////////////
vector<double> MD_OC::temp_itemin(vector<MatrixXd>& arrr, double minpo, int size)
{
    vector<double> dist_return(size, 0);

    for (int i = 0; i < size; i++)
    {
        dist_return[i] = arrr[minpo](i, i);
    }

    return dist_return;
}
///////////////////////////////////////////////////////////////////////////////

vector<MatrixXd> MD_OC::Iteration(const int& n)
{
    //cout << "** Iteration RUN " << endl;

    Prop.resize(t, MatrixXd::Zero(siz, siz));
    Prop[0] = MatrixXd::Identity(siz, siz);
    MatrixXd Iden = MatrixXd::Identity(siz, siz);

    vector<double> lambda(n + 1, 0);
    double expDtauLambda;
    double factor;

    ///////////////////////////////////////////////////////////////

    double temp_minpoin = t - 1;
    vector<vector<double> > temp_itemi(2, vector<double>(siz, 0));
    double RELA_ENTROPY;

    ///////////////////////////////////////////////////////////////

    for (int i = 0; i <= n; i++)
    {
        if (i == 0)
        {
            for (int j = 0; j < t; j++)
            {
                for (int k = 0; k < siz; k++)
                {
                    Prop[j](k, k) = exp(-tau_grid[j] * H_loc(k, k));
                }
            }
            //cout << Prop[99] << endl;


            lambda[0] = chemical_poten(Prop[t - 1]);
            expDtauLambda = exp((tau_grid[1] - tau_grid[0]) * lambda[0]);
            factor = 1.0;


            for (int j = 0; j < t; j++)
            {
                Prop[j] *= factor;
                factor *= expDtauLambda;
                //cout << Prop[j] << endl;
            }
            //////////////////////////////////////////////////////////////////////////////

            temp_minpoin = t / 2;

            //////////////////////////////////////////////////////////////////////////////

        }

        else
        {
            //std::chrono::system_clock::time_point start = std::chrono::system_clock::now();
            //cout << "Iteration " << i << " Starts" << endl;
            /////////////////////////////////////////////////////////////////////////////

            temp_itemi[(i - 1) % 2] = temp_itemin(Prop, temp_minpoin, siz); // temporary store for previous iteration data
            RELA_ENTROPY = 0;

            /////////////////////////////////////////////////////////////////////////////
            H_loc = H_loc - lambda[i - 1] * Iden;
            SELF_Energy();
            Prop = Propagator(SELF_E, H_loc);

            lambda[i] = chemical_poten(Prop[t - 1]);

            expDtauLambda = exp((tau_grid[1] - tau_grid[0]) * lambda[i]);
            factor = 1.0;

            for (int j = 0; j < t; j++)
            {
                Prop[j] *= factor;
                factor *= expDtauLambda;

                //cout << Prop[j] << endl;
            }

            /////////////////////////////////////////////////////////////////////////////

            temp_itemi[i % 2] = temp_itemin(Prop, temp_minpoin, siz);
            cout << "\n";

            // Relative entropy calculation

            for (int j = 0; j < siz; j++)
            {
                RELA_ENTROPY += temp_itemi[i % 2][j] * log(temp_itemi[i % 2][j] / temp_itemi[(i - 1) % 2][j]);
            }


            if (i > 1) {
                cout << "\t""\t" << i << " th Iteration stop value : " << fabs(RELA_ENTROPY) << endl;
                if (fabs(RELA_ENTROPY) < 0.00001) {
                    break;
                }
            }

            /////////////////////////////////////////////////////////////////////////////

            //std::chrono::system_clock::time_point sec = std::chrono::system_clock::now();
            //std::chrono::duration<double> microseconds = std::chrono::duration_cast<std::chrono::milliseconds>(sec - start);
            //cout << "Process ends in : " << microseconds.count() << "[sec]" << endl;
            //cout << "-----------------------------" << endl;
        }
    }

    return Prop;
}


//////////////////////////////////////////////////////////////////////////////
/*
void MD_OC::NCA_Chi_sp(vector<MatrixXd>& iter)
{
    Chi_Arr.resize(t);
    MatrixXd GELL = MatrixXd::Zero(siz, siz);
    GELL(0, 1) = 1;
    GELL(1, 0) = 1;

    for (int i = 0; i < t; i++)
    {
        Chi_Arr[i] = (iter[t - i - 1] * GELL * iter[i] * GELL).trace();
    }
}

void MD_OC::OCA_store(vector<MatrixXd>& iter)
{
    MatrixXd GELL = MatrixXd::Zero(siz, siz);
    GELL(0, 1) = 1;
    GELL(1, 0) = 1;

    Chi_st.resize(t);
    for (int n = 0; n < t; n++)
    {
        Chi_st[n].resize(n + 1);
        for (int m = 0; m <= n; m++)
        {
            Chi_st[n][m] = iter[n - m] * H_N * iter[m] * GELL;
            //cout << "pair (n,m) is : " <<  "(" << n << "," << m << ")" << "corresponds with" << "(" << n-m << "," << m << ")" << endl;
        }
    }
}


void MD_OC::OCA_Chi_sp(vector<MatrixXd>& iter)
{
    for (int i = 0; i < t; i++)
    {
        MatrixXd Stmp = MatrixXd::Zero(siz, siz);

        for (int n = 0; n <= i; n++) for (int m = i; m < t; m++)
        {
            Stmp += INT_Arr[m - n] * (Chi_st[t - i - 1][m - i] * Chi_st[i][n]);
            //cout << "pair ("<<n<<","<<m<<") is : " << "(" << k-i-1 << "," << m-i << ")"<< " with " << "(" << i << "," << n << ")" << endl;
        }
        Chi_Arr[i] += pow(Delta_t, 2) * Stmp.trace();
    }
}

vector<double> MD_OC::Chi_sp_Function(vector<MatrixXd> ITE)
{
    NCA_Chi_sp(ITE);
    OCA_store(ITE);
    OCA_Chi_sp(ITE);

    return Chi_Arr;
}
*/
////////////////////////////////////////////////////////////////////////////////////

void MD_OC::NCA_Chi_sp(vector<MatrixXd>& iter)
{
    Chi_Arr.resize(t);
    MatrixXd GELL = MatrixXd::Zero(siz, siz);
    GELL(0, 1) = 1;
    GELL(1, 0) = 1;

    for (int i = 0; i < t; i++)
    {
        Chi_Arr[i] = (iter[t - i - 1] * Order_param * iter[i] * Order_param).trace();
    }
}

void MD_OC::OCA_store(vector<MatrixXd>& iter)
{
    MatrixXd GELL = MatrixXd::Zero(siz, siz);
    GELL(0, 1) = 1;
    GELL(1, 0) = 1;

    Chi_st.resize(t);
    for (int n = 0; n < t; n++)
    {
        Chi_st[n].resize(n + 1);
        for (int m = 0; m <= n; m++)
        {
            Chi_st[n][m] = iter[n - m] * H_N * iter[m] * Order_param;
            //cout << "pair (n,m) is : " <<  "(" << n << "," << m << ")" << "corresponds with" << "(" << n-m << "," << m << ")" << endl;
        }
    }
}

void MD_OC::OCA_Chi_sp(vector<MatrixXd>& iter)
{
    for (int i = 0; i < t; i++)
    {
        MatrixXd Stmp = MatrixXd::Zero(siz, siz);

        for (int n = 0; n <= i; n++) for (int m = i; m < t; m++)
        {
            Stmp += INT_Arr[m - n] * (Chi_st[t - i - 1][m - i] * Chi_st[i][n]);
            //cout << "pair ("<<n<<","<<m<<") is : " << "(" << k-i-1 << "," << m-i << ")"<< " with " << "(" << i << "," << n << ")" << endl;
        }
        Chi_Arr[i] += pow(Delta_t, 2) * Stmp.trace();
    }
}

vector<double> MD_OC::Chi_sp_Function(vector<MatrixXd> ITE)
{
    NCA_Chi_sp(ITE);
    OCA_store(ITE);
    OCA_Chi_sp(ITE);

    return Chi_Arr;
}

///////////////////////////////////////////////////////////////////////////////////////////////

int main()
{
    double beta;
    int grid;

    cout << " * Set beta : ";
    cin >> beta;

    cout << " * Set grid (number of index, not interval count) : ";
    cin >> grid;

    std::chrono::system_clock::time_point P_start = std::chrono::system_clock::now();
    cout << " ## OCA Program begins ##" << endl;
    cout << "-------------------------------" << endl;
    int modeselec = 0;

    MD_OC MD(beta, grid);
    /// Parameter adjustment ////

    double alpha = 0;
    double k_cutoff = 20;
    double& ref_g_ma = g_ma;

    int& size = siz;
    int& syst = sys;

    size = 3;
    syst = 21;


    /////////////////////////////////

    vector<double> alp_arr(21, 0);
    for (int i = 0; i < 21; i++)
    {
        if (i == 0)
        {
            alp_arr[i] = 0;
        }
        if (i != 0)
        {
            alp_arr[i] = alp_arr[i - 1] + 0.1;
        }

    }


    vector<double> g_ma_arr(21, 0);
    for (int i = 0; i < 21; i++)
    {
        if (i == 0)
        {
            g_ma_arr[i] = 0;
        }
        if (i != 0)
        {
            g_ma_arr[i] = g_ma_arr[i - 1] + 0.1;
        }
    }

    std::ofstream outputFile("./");

    string name = "OCA_GENSIZE_";

    std::stringstream cuof;
    std::stringstream bet;
    std::stringstream gri;
    std::stringstream sizz;

    cuof << k_cutoff;
    bet << MD.tau_grid[MD.tau_grid.size() - 1];
    gri << MD.t;
    sizz << siz;

    name += sizz.str();
    name += "_MODE_";
    name += cuof.str();
    name += "_BETA_";
    name += bet.str();
    name += "_GRID_";
    name += gri.str();
    name += ".txt";

    outputFile.open(name);

    for (int ga = 0; ga < g_ma_arr.size(); ga++)
    {
        //ref_g_ma = 2;
        ref_g_ma = g_ma_arr[ga];
        alpha = 1;
        for (int al = 0; al < alp_arr.size(); al++)
        {
            alpha = alp_arr[al];

            /********************\beta * Chi(\beta / 2) Calculation****************************/
            MD.CAL_COUP_INT_with_g_arr(alpha, k_cutoff);
            vector<MatrixXd> ITER = MD.Iteration(25);
            vector<double> a = MD.Chi_sp_Function(ITER);

            outputFile << MD.tau_grid[MD.t - 1] * a[int(MD.t / 2)] << "\t";
            /**************************************************************************/
            ////////////////////DATA OUTPUT ///////////////////
            string Prop_name = "OCA_PROP_GAMMA_";

            std::stringstream gam;
            std::stringstream alp;
            std::stringstream cuof;
            std::stringstream bet;
            std::stringstream gri;
            std::stringstream sizz;
            std::stringstream mod;

            gam << g_ma;
            alp << alpha;
            cuof << k_cutoff;
            mod << MD.mode_grid.size();
            bet << MD.tau_grid[MD.tau_grid.size() - 1];
            gri << MD.t;
            sizz << siz;

            Prop_name += gam.str();
            Prop_name += "_ALPHA_";
            Prop_name += alp.str();
            Prop_name += "_CUTOF_";
            Prop_name += cuof.str();
            Prop_name += "_MODE_";
            Prop_name += mod.str();
            Prop_name += "_BETA_";
            Prop_name += bet.str();
            Prop_name += "_GRID_";
            Prop_name += gri.str();
            Prop_name += "_SIZE_";
            Prop_name += sizz.str();
            Prop_name += "_T.txt";

            outputFile.open(Prop_name);

            for (int k = 0; k < MD.tau_grid.size(); k++) {
                for (int i = 0; i < siz; i++) for (int j = 0; j < siz; j++)
                {
                    outputFile << ITER[k](i, j) << "\t";
                }
                outputFile << "\n";
            }

            outputFile.close();

            string Chi_name = "OCA_CHI_GAMMA_";

            Chi_name += gam.str();
            Chi_name += "_ALPHA_";
            Chi_name += alp.str();
            Chi_name += "_MODE_";
            Chi_name += cuof.str();
            Chi_name += "_BETA_";
            Chi_name += bet.str();
            Chi_name += "_GRID_";
            Chi_name += gri.str();
            Chi_name += "_SIZE_";
            Chi_name += sizz.str();
            Chi_name += "_T.txt";

            outputFile.open(Chi_name);

            for (int j = 0; j < MD.tau_grid.size(); j++)
            {
                outputFile << MD.tau_grid[j] << "\t" << a[j] << endl;
            }

            outputFile.close();
            ////////////////////DATA OUTPUT ///////////////////

        /********************Hybridization Check****************************/
        /*
            std::stringstream gam;
            std::stringstream alp;

            gam << g_ma_arr[ga];
            alp << alp_arr[al];

            string Hyb = "HYB_SIZE_";

            Hyb += sizz.str();
            Hyb += "_MODE_";
            Hyb += cuof.str();
            Hyb += "_BETA_";
            Hyb += bet.str();
            Hyb += "_GRID_";
            Hyb += gri.str();
            Hyb += "_GAM_";
            Hyb += gam.str();
            Hyb += "_ALP_";
            Hyb += alp.str();
            Hyb += ".txt";

            outputFile.open(Hyb);

            for (int j = 0; j < MD.tau_grid.size(); j++)
            {
                outputFile << MD.tau_grid[j] << "\t" << MD.INT_Arr[j] << endl;
            }

            outputFile.close();

        */
        /*************************************************************************/

        }
        outputFile << "\n";

    }

    outputFile.close();


    std::chrono::system_clock::time_point P_sec = std::chrono::system_clock::now();
    std::chrono::duration<double> seconds = std::chrono::duration_cast<std::chrono::seconds>(P_sec - P_start);
    cout << "## Total Process ends with : " << seconds.count() << "[sec] ##" << endl;
    cout << "-----------------------------" << endl;

    return 0;


}
