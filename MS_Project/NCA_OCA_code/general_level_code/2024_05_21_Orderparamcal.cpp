#include<iostream>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <vector>
#include <cmath>
#include <OCA_bath.hpp>
#include <chrono>
#include <const.h>

using namespace std;
using namespace Eigen;

vector<double> k_mode(30000, 1);
double g_ma = 1;
int siz = 0;

///////////////////////////////////////////////////////

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
///////////////////////

void MD_OC::Tilde_g_calculation_function(double alpha, double k_cutoff)
{
    omega_Arr.resize(M);
    coup_Arr.resize(M);
    double nu = pi * k_cutoff / alpha;
    
    for (int i=0; i < M; i++)
    {
        omega_Arr[i] = k_cutoff * (mode_grid[i]/mode_grid[M-1]);
        coup_Arr[i] = sqrt((2 * k_cutoff / (alpha * M)) * (omega_Arr[i] / (1 + pow(nu * omega_Arr[i] / k_cutoff,2))));

        //simpson formulae
        //omega_Arr[i] = (mode_grid[i]/mode_grid[M-1]); // fix to x to adjust simpson's rule
        //coup_Arr[i] = sqrt((2 * k_cutoff / (alpha)) * ( k_cutoff * omega_Arr[i] / (1 + pow(nu * omega_Arr[i],2)))); // fix to adjust simpson's rule
    }

    if (alpha == 0)
    {
        for (int i=0; i < M; i++)
        {
            coup_Arr[i] = 0;
        }
    }
}

////////////////////////////////////////////////////////////////////////////////////

void MD_OC::Interact_V(double k_cutoff)
{
    //Initializing block
    
    for (int i=0; i < t; i++)
    {
        INT_Arr[i] = 0;
    }
    

    for (int i = 0; i < t; i++)
    {
        for (int j = 0; j < M ;j++)
        {
            INT_Arr[i] += -pow(coup_Arr[j],2) * cosh((tau_grid[i] - tau_grid[t - 1] / 2) * omega_Arr[j])/sinh(tau_grid[t - 1] * omega_Arr[j] / 2); //caution for sign
        }
    }

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

////////////////////////////////////////////////////////////////////////////////////

void MD_OC::Hamiltonian_N(MatrixXd even, MatrixXd odd)
{
    //cout << "input g value :" << g << endl;
    MatrixXd INT_odd = MatrixXd::Zero(siz, siz);
    MatrixXd INT_even = MatrixXd::Zero(siz, siz);

    //H_N initialize
    H_N = MatrixXd::Zero(siz, siz);

    //cout << "initialized check" << endl;
    //cout << H_N << "\n" << endl;

    for (int i = 0; i < siz; i++) for (int j = 0; j < siz; j++)
    {
        INT_even(i, j) = -1 * even(i, j) * i; // -\sum_1^\infty \alpha_i \sin{i\phi}

        if (i < siz - 1)
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
///////////////////////////////////////////////////////////////////////

MatrixXd Ordercal(MatrixXd even, MatrixXd odd)
{
    MatrixXd Order_param = MatrixXd::Zero(siz, siz);

    cout << "**EVENMAT**" << endl;
    cout << even << endl;

    //** constructing even matrix
    MatrixXd eve_0 = even;//MatrixXd::Zero(siz,siz);
    for (int i = 0; i < siz; i++) for (int j = 0; j < siz; j++)
    {
        if (i == 0) {
            eve_0(i, j) = (1 / sqrt(2)) * even(i, j);
        }
        else {
            eve_0(i, j) = even(i, j);
        }
    }

    //Even OffDiagonal construct
    MatrixXd eve_1 = MatrixXd::Zero(siz + 1, siz + 1);
    MatrixXd eve_2 = MatrixXd::Zero(siz + 1, siz + 1);

    for (int i = 0; i < siz + 1; i++) for (int j = 0; j < siz + 1; j++)
    {
        if (i > 0 && j < siz) {
            if (i > 1) {
                eve_1(i, j) = 0.5 * eve_0(i - 1, j);
            }
            else {
                eve_1(i, j) = eve_0(i - 1, j);
            }
        }
        if (j < siz && i < siz) {
            eve_2(i, j) = eve_0(i, j);
        }
    }
    //Activate only if size 3
    for (int i = 0; i < siz + 1; i++) {
        eve_1(i, 0) = -eve_1(i, 0);
        eve_2(i, 0) = -eve_2(i, 0);
    }


    MatrixXd eve_off = eve_1.transpose() * eve_2;
    cout << "EVEOFF" << endl;
    cout << eve_off << endl;

    cout << "\t" << "<Mateven 1>" << endl;
    cout << eve_1.transpose() << endl;
    /*
    for (int i = 0; i < siz+1; i++) for (int j = 0; j< siz+1; j++)
    {
        if (j<siz && i<siz)
        {
            eve_2(i,j) = eve_0(i,j);
        }
    }
    */
    cout << "\t" << "<Mateven 2>" << endl;
    cout << eve_2 << endl;
    MatrixXd eve_ele = MatrixXd::Zero(siz, siz);

    for (int i = 0; i < siz; i++) for (int j = 0; j < siz; j++)
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

    cout << "*****EVEELE*****" << endl;
    cout << eve_ele << endl;
    ///////////// even matrix construction complete ///////////////

    //constructing odd matrix

    cout << "**ODD Eigen matrix**" << endl;
    cout << odd << endl;

    MatrixXd odd_ele = MatrixXd::Zero(siz, siz);
    for (int i = 0; i < siz; i++) for (int j = 0; j < siz; j++)
    {
        if ((i == j) && (i % 2 != 0)) {
            odd_ele(i, j) = odd(0, j / 2) * odd(1, j / 2); // Diagonal element analogous with odd basis
        }
        else if ((i < j) && (i % 2 != 0) && (j % 2 != 0)) {
            odd_ele(i, j) = (odd(0, i / 2) * odd(1, j / 2) + odd(1, i / 2) * odd(0, j / 2)) / 2;
            odd_ele(j, i) = (odd(0, i / 2) * odd(1, j / 2) + odd(1, i / 2) * odd(0, j / 2)) / 2;
        }
    }
    cout << "Odd" << endl;
    cout << odd_ele << endl;
    cout << "\n\n";

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

    return Order_param;

}


void MD_OC::Hamiltonian_loc(MatrixXd a, MatrixXd b)
{
    int siz = a.rows();
    H_loc = MatrixXd::Zero(siz,siz);

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

///////////////////////////////////////////////////////////////////////////////


void MD_OC::CAL_COUP_INT_with_g_arr(double alpha, double k_cutoff)
{
    Tilde_g_calculation_function(alpha,k_cutoff);
    Interact_V(k_cutoff);
    Hamiltonian_N(Eigenvector_Even(), Eigenvector_Odd());
    Hamiltonian_loc(Eigenvalue_Even(),Eigenvalue_Odd());


    cout << "$ H_N value : \n " << H_N << endl;
    cout << "$ H_loc value : \n " << H_loc << endl;

}


////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////BLOCK FOR CALCULATE ORDERPARAMETER//////////////////////////////////////////////////////////////
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
    OCA_self();
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
    vector<MatrixXd> P_arr(t, MatrixXd::Zero(siz,siz));
    vector<MatrixXd> S_arr(t, MatrixXd::Zero(siz,siz));

    P_arr[0] = MatrixXd::Identity(siz,siz);
    S_arr[0] = MatrixXd::Identity(siz,siz);

    MatrixXd sig_form = MatrixXd::Zero(siz,siz);
    MatrixXd sig_late = MatrixXd::Zero(siz,siz);

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


vector<MatrixXd> MD_OC::Iteration(const int& n)
{
    cout << "** Iteration RUN " << endl;

    Prop.resize(t, MatrixXd::Zero(siz,siz));
    Prop[0] = MatrixXd::Identity(siz,siz);
    MatrixXd Iden = MatrixXd::Identity(siz,siz);

    vector<double> lambda(n + 1, 0);
    double expDtauLambda;
    double factor;

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
        }

        else
        {
            std::chrono::system_clock::time_point start = std::chrono::system_clock::now();
            cout << "Iteration " << i << " Starts" << endl;
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
            std::chrono::system_clock::time_point sec = std::chrono::system_clock::now();
            std::chrono::duration<double> microseconds = std::chrono::duration_cast<std::chrono::milliseconds>(sec - start);
            cout << "Process ends in : " << microseconds.count() << "[sec]" << endl;
            cout << "-----------------------------" << endl;
        }
    }

    return Prop;
}
/////////////////////////////////////////////////////////BLOCK FOR CALCULATE ORDERPARAMETER//////////////////////////////////////////////////////////////


int main()
{
    std::chrono::system_clock::time_point P_start= std::chrono::system_clock::now();
    cout << " ## Approx Program begins ##" << endl;
    cout << " Program now set in One-crossing mode " << endl;
    cout << "-------------------------------" << endl;
    /*
    double& taulimit = MD.beta;

    cout << " Set BETA values to calculate : ";
    cin >> taulimit;

    cout << "\n" << "Calculation would be done under " << taulimit << " value";
    */

    //while (modeselec != -1)

    //cout << "< Select mode to run >" << "\n"  << " 1. Prop(G), 2. Chi, 3. beta*Chi " << "\n" << "MODE INDEX : ";
    //cin >> modeselec;
    double beta;
    int grid;

    cout << " * Set beta : ";
    cin >> beta;

    cout << " * Set grid (number of index, not interval count) : ";
    cin >> grid;

    MD_OC MD(beta,grid);

    double alpha = 0.5;
    double k_cutoff = 20;
    double& ref_g_ma = g_ma;

    int &size = siz;

    size = 5;
    /*
    vector<double> alp_arr(21,0);
    for (int i = 0; i < 21 ; i++)
    {
        if (i==0)
        {
            alp_arr[i] = 0;
        }
        if (i!=0)
        {
            alp_arr[i] = alp_arr[i-1] + 0.1;
        }
        
    }
    */
    
    vector<double> alp_arr = {0,1,5,20};
    
    
    vector<double> g_ma_arr(11,0);
    for (int i = 0; i < 11 ; i++)
    {
        if (i==0)
        {
            g_ma_arr[i] = 0;
        }
        if (i!=0)
        {
            g_ma_arr[i] = g_ma_arr[i-1] + 0.1;
        }
    }
    
    vector<double> output(g_ma_arr.size(),0);


    for (int al = 0; al < alp_arr.size(); al ++)
    {
        for (int ga = 0; ga < 1 ; ga++)//g_ma_arr.size(); ga++)
        {
            //ref_g_ma = g_ma_arr[ga];
            alpha = 1;
            alpha = alp_arr[al];
            //ref_g_ma = g_ma_arr[ga];
        
                MD.CAL_COUP_INT_with_g_arr(alpha,k_cutoff);
                MatrixXd Order_param = Ordercal(MD.Eigenvector_Even(),MD.Eigenvector_Odd());

                cout << "Value" << endl;
                cout << Order_param << endl;

                vector<MatrixXd> a = MD.Iteration(5);
                double b = (a[(MD.t-1)] * Order_param).trace(); // Trace of Prop(beta)*(Orderparam(cos\phi))
                //cout << "TEST PROP BETA : " << "\n" << a[MD.t-1] << endl;

                output[ga] = b;
                
        }
        
        std::ofstream outputFile("./");

        string name = "OCA_COS_TEST_CAL_";

        std::stringstream gam;
        std::stringstream alp;
        std::stringstream cuof;
        std::stringstream bet;
        std::stringstream gri;
        std::stringstream sizz;

        gam << g_ma;
        alp << alpha;
        cuof << k_cutoff;
        bet << MD.tau_grid[MD.t-1];
        gri << MD.t;
        sizz << siz;

        name += gam.str();
        name += "_ALPHA_";
        name += alp.str();
        name += "_MODE_";
        name += cuof.str();
        name += "_BETA_";
        name += bet.str();
        name += "_GRID_";
        name += gri.str();
        name += "_SIZE_";
        name += sizz.str();
        name += ".txt";

        outputFile.open(name);

        for (int i = 0; i < output.size(); i++)
        {
            outputFile << g_ma_arr[i] << "\t" << output[i] << endl;
        }

        outputFile.close();
    

        }
        
        
        
}
    //}
//}
