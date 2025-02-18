#include <iostream>
#include <fstream>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <vector>
#include <cmath>
#include <iomanip>
#include <string>
#include <const.h>
#include <chrono>

using namespace std;
using namespace Eigen;

class MD_NC
{
private:

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

    vector<MatrixXd> convolve(const vector<MatrixXd>& Signal,
        const vector<MatrixXd>& Kernel, int n, int i)
    {
        size_t SignalLen = i;
        size_t KernelLen = Kernel.size();
        size_t ResultLen = SignalLen + KernelLen - 1;

        vector<MatrixXd> Result(ResultLen, MatrixXd::Zero(n, n));

        for (size_t n = 0; n < ResultLen; ++n)
        {
            size_t kmin = (n >= KernelLen - 1) ? n - (KernelLen - 1) : 0;
            size_t kmax = (n < SignalLen - 1) ? n : SignalLen - 1;

            for (size_t k = kmin; k <= kmax; k++)
            {
                Result[n] += Signal[k] * Kernel[n - k];
            }
        }

        return Result;
    }

    MatrixXd Matrix_Odd(int n, double r)
    {
        MatrixXd Matrix1 = MatrixXd::Zero(n, n);

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++)
            {
                try
                {
                    if (i == j) {
                        Matrix1(i, j) = pow((i + 1), 2);
                    }
                    if (abs(i - j) == 1) {
                        Matrix1(i, j) = -r / 2.0;
                    }
                }
                catch (...) {}
            }
        }
        return Matrix1;
    }

    MatrixXd Matrix_Even(int n, double r)
    {
        MatrixXd Matrix1 = MatrixXd::Zero(n, n);

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++)
            {
                try
                {
                    if (i == j) {
                        Matrix1(i, j) = pow(i, 2);
                    }
                    if (abs(i - j) == 1) {
                        Matrix1(i, j) = -r / 2.0;
                    }
                }
                catch (...) {}
            }
        }
        Matrix1(0, 1) = -r / sqrt(2);
        Matrix1(1, 0) = -r / sqrt(2);

        return Matrix1;
    }



public:
    double pi = dlib::pi;
    double hbar = dlib::planck_cst / 2 * dlib::pi;

    vector<double> tau_grid = linspace(0, 20, 401);
    vector<double> mode_grid = linspace(1, 30000, 30000);
    int beta = tau_grid.size();
    int M = mode_grid.size();
    double Delta_t = tau_grid[1] - tau_grid[0];

    static MatrixXd H_N;

    vector<double> coupling(double v, double g, double W);
    vector<double> Interact_V(vector<double> couplint, vector<double> tau, double omega);


    void Tilde_g_calculation_function(double alpha, double k_cutoff);
    vector<double> Interact_V();


    MatrixXd Eigenvector_Even();
    MatrixXd Eigenvalue_Even();
    MatrixXd Eigenvector_Odd();
    MatrixXd Eigenvalue_Odd();

    void Hamiltonian_N(MatrixXd even, MatrixXd odd);
    vector<MatrixXd> Hamiltonian_exp(MatrixXd a, MatrixXd b);
    MatrixXd Hamiltonian_loc(MatrixXd a, MatrixXd b);
    MatrixXd Hamiltonian_loc_ite(MatrixXd a, MatrixXd b, const double& lambda);

    void CAL_COUP_INT_with_g_arr(double alpha, double k_cutoff);
    void NCA_self(const vector<MatrixXd>& H_exp, const vector<double>& V);

    MatrixXd round_propagater_ite(const MatrixXd& loc, const vector<MatrixXd>& sigma, const vector<MatrixXd>& ite, int n, int boolean);
    vector<MatrixXd> Propagator(const vector<MatrixXd>& array, const MatrixXd& loc);

    double chemical_poten(MatrixXd prop);

    vector<MatrixXd> Iteration(const int& iteration);

    void Chi_sp(int ITE);

};
/////////////////////////////////////////////////////////////////////////////////////

MD_NC MD;

//double gamma = 0;
int siz = 5;
//double nu = MD.pi/0.025;

///////////////////////////////////////////////////////////////

vector<double> G_Arr(MD.M, 0);
vector<double> omega_Arr(MD.M, 0);
//vector<MatrixXd> H_N(MD.M,MatrixXd::Zero(sizsiz);

//////////////////////////////////////////////////////////////

vector<double> INT_Arr(MD.beta, 0);
vector<double> Chi_Arr(MD.beta, 0);

vector<MatrixXd> SELF_E(MD.beta, MatrixXd::Zero(siz, siz));
MatrixXd MD_NC::H_N = MatrixXd::Zero(siz, siz);

//////////////////////////////////////////////////////////////

vector<double> k_mode(100, 1);
double gamma = 1;

double omega = 1;
double velocity = 1;
double cutoff = 1;


/////////////////////////////////////////////////////////////////////////////////////


void MD_NC::Tilde_g_calculation_function(double alpha, double k_cutoff)
{
    double nu = pi * k_cutoff / alpha;

    for (int i = 0; i < M; i++)
    {
        omega_Arr[i] = 0;
        G_Arr[i] = 0;

        //tilde_g_arr[i] = sqrt( (omega_arr[i] / (1 + pow(nu * omega_arr[i] / k_cutoff,2))));
        //tilde_g_arr[i] = sqrt((2 * k_cutoff / (alpha * omega_arr.size())) * (re_planck_cst * omega_arr[i] / (1 + pow(nu * re_planck_cst * omega_arr[i] / k_cutoff,2))));
    }

    for (int i = 0; i < M; i++)
    {
        omega_Arr[i] = k_cutoff * (mode_grid[i] / mode_grid[M - 1]);
        G_Arr[i] = sqrt((2 * k_cutoff / (alpha * M)) * (omega_Arr[i] / (1 + pow(nu * omega_Arr[i] / k_cutoff, 2))));
        //tilde_g_arr[i] = sqrt( (omega_arr[i] / (1 + pow(nu * omega_arr[i] / k_cutoff,2))));
        //tilde_g_arr[i] = sqrt((2 * k_cutoff / (alpha * omega_arr.size())) * (re_planck_cst * omega_arr[i] / (1 + pow(nu * re_planck_cst * omega_arr[i] / k_cutoff,2))));
    }

    if (alpha == 0)
    {
        for (int i = 0; i < M; i++)
        {
            G_Arr[i] = 0;
        }

    }
}


////////////////////////////////////////////////////////////////////////////////////


vector<double> MD_NC::Interact_V()
{

    for (int i = 0; i < beta; i++)
    {
        INT_Arr[i] = 0;
    }


    for (int i = 0; i < beta; i++)
    {
        for (int j = 0; j < M; j++)
        {
            INT_Arr[i] += -pow(G_Arr[j], 2) * cosh((tau_grid[i] - tau_grid[beta - 1] / 2) * omega_Arr[j]) / sinh(tau_grid[beta - 1] * omega_Arr[j] / 2); //caution for sign
            //cout << "\t" << j <<" V_arr : " << V_arr[i] << " with tau-beta/2 : " << tau[i] - tau[tau.size()-1]/2 << endl;
        }
    }

    return INT_Arr;
}


////////////////////////////////////////////////////////////////////////////////////

MatrixXd MD_NC::Eigenvector_Even()
{
    SelfAdjointEigenSolver<MatrixXd> es(Matrix_Even(siz, gamma));
    return es.eigenvectors();
}

MatrixXd MD_NC::Eigenvalue_Even()
{
    SelfAdjointEigenSolver<MatrixXd> es(Matrix_Even(siz, gamma));
    return es.eigenvalues();
}

MatrixXd MD_NC::Eigenvector_Odd()
{
    SelfAdjointEigenSolver<MatrixXd> es(Matrix_Odd(siz, gamma));
    return es.eigenvectors();
}

MatrixXd MD_NC::Eigenvalue_Odd()
{
    SelfAdjointEigenSolver<MatrixXd> es(Matrix_Odd(siz, gamma));
    return es.eigenvalues();
}

///////////////////////////////////////////////////////////////////////


void MD_NC::Hamiltonian_N(MatrixXd even, MatrixXd odd)
{
    int siz = even.rows();

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

    //cout << "ODD" << endl;
    //cout << INT_odd << endl;
    //cout << "EVEN" << endl;
    //cout << INT_even << "\n" << endl;

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

vector<MatrixXd> MD_NC::Hamiltonian_exp(MatrixXd a, MatrixXd b)
{
    //g_0
    MatrixXd Even = a;
    MatrixXd Odd = b;

    double zeroth = exp(Even(0));
    double first = exp(Odd(0));
    double second = exp(Even(1));

    vector<MatrixXd> array_with_Matrix(beta);

    MatrixXd Hamiltonian_exp;

    for (int i = 0; i < beta; i++)
    {
        Hamiltonian_exp = MatrixXd::Zero(siz, siz);

        Hamiltonian_exp(0, 0) = tau_grid[i] * zeroth;
        Hamiltonian_exp(1, 1) = tau_grid[i] * first;
        Hamiltonian_exp(2, 2) = tau_grid[i] * second;

        array_with_Matrix[i] = Hamiltonian_exp;
    }

    return array_with_Matrix;
}



MatrixXd MD_NC::Hamiltonian_loc(MatrixXd a, MatrixXd b)
{
    int siz = a.rows();
    MatrixXd Hamiltonian = MatrixXd::Zero(siz, siz);

    for (int i = 0; i < siz; i++) for (int j = 0; j < siz; j++)
    {
        if (i == j & i % 2 == 0)
        {
            Hamiltonian(i, j) = a(i / 2);
            /*
            if (a(i / 2) > 30)
            {
                Hamiltonian(i, j) = 30;
            }
            */
        }
        if (i == j & i % 2 != 0)
        {
            Hamiltonian(i, j) = b(i / 2);
            /*
            if (b(i / 2) > 30)
            {
                Hamiltonian(i, j) = 30;
            }
            */
        }
    }

    return Hamiltonian;
}

MatrixXd MD_NC::Hamiltonian_loc_ite(MatrixXd a, MatrixXd b, const double& lambda)
{
    int siz = a.rows();
    MatrixXd Hamiltonian = MatrixXd::Zero(siz, siz);

    Hamiltonian(0, 0) = a(0) - lambda;
    Hamiltonian(1, 1) = b(0) - lambda;
    Hamiltonian(2, 2) = a(1) - lambda;

    return Hamiltonian;
}

////////////////////////////////////////////////////////////////////////////////

void MD_NC::CAL_COUP_INT_with_g_arr(double alpha, double k_cutoff)
{
    Tilde_g_calculation_function(alpha, k_cutoff);
    INT_Arr = Interact_V();
    Hamiltonian_N(Eigenvector_Even(), Eigenvector_Odd());
    cout << "** bare H_loc \n" << endl;
    cout << Hamiltonian_loc(Eigenvalue_Even(), Eigenvalue_Odd()) << endl;
    cout << "--------------" << endl;
}


////////////////////////////////////////////////////////////////////////////////


void MD_NC::NCA_self(const vector<MatrixXd>& Prop, const vector<double>& V)
{
    //initializing SELF_Ee
    for (int i = 0; i < beta; i++)
    {
        SELF_E[i] = MatrixXd::Zero(siz, siz);
    }

    for (int i = 0; i < beta; i++)
    {
        SELF_E[i] = V[i] * (H_N * Prop[i] * H_N);
    }
    cout << "     ******* SELF_E (beta) : \n " << SELF_E[beta-1] << endl;
}


//////////////////////////////////////////////////////////////////////////////


MatrixXd MD_NC::round_propagater_ite(const MatrixXd& loc, const vector<MatrixXd>& sigma, const vector<MatrixXd>& ite, int n, int boolean)
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



vector<MatrixXd> MD_NC::Propagator(const vector<MatrixXd>& self_E, const MatrixXd& loc)
{
    vector<MatrixXd> P_arr(beta, MatrixXd::Zero(siz, siz));
    vector<MatrixXd> S_arr(beta, MatrixXd::Zero(siz, siz));

    P_arr[0] = MatrixXd::Identity(siz, siz);
    S_arr[0] = MatrixXd::Identity(siz, siz);

    MatrixXd sig_form = MatrixXd::Zero(siz, siz);
    MatrixXd sig_late = MatrixXd::Zero(siz, siz);

    for (int i = 1; i < beta; i++)
    {
        P_arr[1] = P_arr[0];
        sig_late = 0.5 * Delta_t * (0.5 * Delta_t * (self_E[1] * P_arr[0] + self_E[0] * (P_arr[0] + Delta_t * P_arr[0])));
        P_arr[1] = P_arr[0] - 0.5 * Delta_t * loc * (2 * P_arr[0] + Delta_t * P_arr[0]) + sig_late;
        S_arr[1] = P_arr[1];

        if (i > 1)
        {
            sig_form = round_propagater_ite(loc, self_E, P_arr, i - 1, 0);
            S_arr[i] = P_arr[i - 1] + Delta_t * sig_form;

            sig_late = 0.5 * Delta_t * (round_propagater_ite(loc, self_E, P_arr, i - 1, 1) + round_propagater_ite(loc, self_E, S_arr, i, 1));
            P_arr[i] = P_arr[i - 1] - 0.5 * Delta_t * loc * (2 * P_arr[i - 1] + Delta_t * sig_form) + sig_late;

        }
    }

    return P_arr;
}

/////////////////////////////////////////////////////////////////////////////

double MD_NC::chemical_poten(MatrixXd prop)
{
    double Trace = prop.trace();
    double lambda = -(1 / tau_grid[beta - 1]) * log(Trace);

    cout << " \n lambda checking : " << lambda << endl;

    return lambda;
}

///////////////////////////////////////////////////////////////////////////////

vector<MatrixXd> MD_NC::Iteration(const int& n)
{
    vector<MatrixXd> Prop(beta, MatrixXd::Zero(siz, siz));
    Prop[0] = MatrixXd::Identity(siz, siz);

    vector<MatrixXd> H_loc(n + 1, MatrixXd::Zero(siz, siz));
    H_loc[0] = Hamiltonian_loc(Eigenvalue_Even(), Eigenvalue_Odd());

    MatrixXd Iden = MatrixXd::Identity(siz, siz);

    vector<double> lambda(n + 1, 0);
    double expDtauLambda;
    double factor;

    for (int i = 0; i <= n; i++)
    {
        if (i == 0)
        {
            for (int j = 0; j < beta; j++)
            {
                for (int k = 0; k < siz; k++)
                {
                    Prop[j](k, k) = exp(-tau_grid[j] * H_loc[0](k, k));
                }
            }

            cout << "************* bare Prop at " << beta << " : \n " << Prop[beta-1] << endl;

            lambda[0] = chemical_poten(Prop[beta - 1]);
            expDtauLambda = exp(Delta_t * lambda[0]);
            factor = 1.0;


            for (int j = 0; j < beta; j++)
            {
                Prop[j] *= factor;
                factor *= expDtauLambda;
                //cout << Prop[j].trace() << endl;
            }
        }

        else
        {

            std::chrono::system_clock::time_point start = std::chrono::system_clock::now();
            cout << "Iteration " << i << " Starts" << endl;
            H_loc[i] = H_loc[i - 1] - lambda[i - 1] * Iden;
            cout << " ***  H_loc (BEFORE NORMALIZE)  " << beta << " : \n " << H_loc[i] << endl;
            NCA_self(Prop, INT_Arr);
            Prop = Propagator(SELF_E, H_loc[i]);
            cout << " ***  Prop (BEFORE NORMALIZE)  " << beta << " with Trace " << Prop[beta-1].trace() << ": \n " << Prop[beta-1] << endl;

            //cout << "  Prop at " << beta << " : \n " << Prop[beta-1] << endl;


            lambda[i] = chemical_poten(Prop[beta - 1]);

            expDtauLambda = exp(Delta_t * lambda[i]);
            factor = 1.0;

            for (int j = 0; j < beta; j++)
            {
                Prop[j] *= factor;
                factor *= expDtauLambda;
            }
            cout << " -- Norm -- " << endl;
            cout << i << "th iteration Prop : \n " << Prop[beta - 1] << endl;
            cout << " ----" << endl;
            cout << " ** " << i << "th Trace Prop : " << Prop[beta - 1].trace() << endl;

            std::chrono::system_clock::time_point sec = std::chrono::system_clock::now();
            std::chrono::duration<double> microseconds = std::chrono::duration_cast<std::chrono::milliseconds>(sec - start);
            cout << "Process ends in : " << microseconds.count() << "[sec]" << endl;
            cout << "-----------------------------" << endl;

        }

    }

    return Prop;
}

//////////////////////////////////////////////////////////////////////////////

void MD_NC::Chi_sp(int ITE)
{
    MatrixXd Gellmann_1 = MatrixXd::Zero(siz, siz);
    Gellmann_1(0, 1) = 1;
    Gellmann_1(1, 0) = 1;

    vector<MatrixXd> Ite_ra = Iteration(ITE);

    for (int i = 0; i < beta; i++)
    {
        Chi_Arr[i] = (Ite_ra[beta - i - 1] * Gellmann_1 * Ite_ra[i] * Gellmann_1).trace(); // main code
        //Chi_Arr[i] = (Ite_ra[i] * Gellmann_1).trace();
        cout << setprecision(16);
        //cout << chi_array[i] << endl;
    }
}

int main()
{

    /*
    double& taulimit = MD.bett;

    cout << " Set BETA values to calculate : ";
    cin >> taulimit;


    cout << "\n" << "Calculation would be done under " << MD.bett << " value" << endl;
    */

    std::chrono::system_clock::time_point P_start = std::chrono::system_clock::now();
    double alpha = 1;
    double k_cutoff = 20;
    double& ref_gamma = gamma;

    cout << " ## Program begins ##" << endl;
    cout << "-------------------------------" << endl;
    /*
    double& taulimit = MD.bett;

    cout << " Set BETA values to calculate : ";
    cin >> taulimit;


    cout << "\n" << "Calculation would be done under " << MD.bett << " value" << endl;
    */

    /*while (modeselec != -1)
    {

    cout << "< Select mode to run >" << "\n"  << " 1. Prop(G), 2. Chi, siz beta*Chi " << "\n" << "MODE INDEX : ";
    cin >> modeselec;
    */

    vector<double> gamma_arr(21, 0);
    for (int i = 0; i < 21; i++)
    {
        if (i == 0)
        {
            gamma_arr[i] = 0;
        }
        if (i != 0)
        {
            gamma_arr[i] = gamma_arr[i - 1] + 0.1;
        }

    }
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

    for (int ga = 0; ga < gamma_arr.size(); ga++) for (int al = 0; al < alp_arr.size(); al++)
    {

        //ref_gamma = 1;
        ref_gamma = gamma_arr[ga];
        //alpha = 0.2;
        alpha = alp_arr[al];

        /****************************G(tau) Calcultaion******************************/
        
        for (int i = 0; i < 1; i++)
        {
            std::ofstream outputFile("");

            string name = "NCA_PROP_GAMMA_";

            std::stringstream gam;
            std::stringstream alp;
            std::stringstream cuof;
            std::stringstream bet;
            std::stringstream gri;
            std::stringstream size;
            std::stringstream mod;

            gam << gamma;
            alp << alpha;
            cuof << k_cutoff;
            mod << MD.mode_grid.size();
            bet << MD.tau_grid[MD.tau_grid.size() - 1];
            gri << MD.beta;
            size << siz;

            name += gam.str();
            name += "_ALPHA_";
            name += alp.str();
            name += "_CUTOF_";
            name += cuof.str();
            name += "_MODE_";
            name += mod.str();
            name += "_BETA_";
            name += bet.str();
            name += "_GRID_";
            name += gri.str();
            name += "_SIZE_";
            name += size.str();
            name += ".txt";


            //cout << gamma_arr[ga] << endl;

            outputFile.open(name);
            MD.CAL_COUP_INT_with_g_arr(alpha, k_cutoff);
            vector<MatrixXd> a = MD.Iteration(20);
            /*
            for (int i = 0; i < a.size(); i++)
            {
                //cout << (a[i])[0][0] << (a[i])[0][1] << endl;
                outputFile << MD.tau_grid[i] << "\t" << (a[i])(0, 0) << "\t" << (a[i])(0, 1) << "\t" << (a[i])(0, 2) << "\t"
                    << (a[i])(1, 0) << "\t" << (a[i])(1, 1) << "\t" << (a[i])(1, 2) << "\t"
                    << (a[i])(2, 0) << "\t" << (a[i])(2, 1) << "\t" << (a[i])(2, 2) << "\t" << endl; //변수 a에 값을 할당 후 벡터 각 요소를 반복문으로 불러옴. 이전에는 a 대신 함수를 반복해서 호출하는 방법을 썼는데 그래서 계산 시간이 오래 걸림.
                cout << setprecision(16);
            }
            */

            outputFile << MD.tau_grid[MD.beta - 1] << "\t";

            for (int i = 0; i < siz; i++)
            {
                outputFile << a[MD.beta - 1](i, i) << "\t";
            }

            outputFile.close();


        }
        /****************************************************************************/


        /********************Chi(\tau) Calculation****************************/
        /*
        for (int i = 0; i < 1; i++)
        {
            std::ofstream outputFile("/Users/e2_602_qma/Documents/GitHub/Anaconda/C++_Mac/EXECUTION");

            string name = "NCA_CHI_GAMMA_";

            std::stringstream gam;
            std::stringstream alp;
            std::stringstream cuof;
            std::stringstream bet;
            std::stringstream gri;
            std::stringstream size;

            gam << gamma;
            alp << alpha;
            cuof << k_cutoff;
            bet << MD.tau_grid[MD.tau_grid.size() - 1];
            gri << MD.beta;
            size << siz;

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
            name += size.str();
            name += ".txt";

            outputFile.open(name);
            MD.CAL_COUP_INT_with_g_arr(alpha, k_cutoff);
            MD.Chi_sp(10);

            for (int i = 0; i < MD.beta; i++)
            {
                outputFile << MD.tau_grid[i] << "\t" << Chi_Arr[i] << endl; //변수 a에 값을 할당 후 벡터 각 요소를 반복문으로 불러옴. 이전에는 a 대신 함수를 반복해서 호출하는 방법을 썼는데 그래서 계산 시간이 오래 걸림.
            }
            outputFile.close();

        }
        
        /**************************************************************************/


        /********************\beta * Chi(\beta / 2) Calculation****************************/
        /*
        for (int i = 0; i < 1; i++)
        {
            std::ofstream outputFile("/Users/e2_602_qma/Documents/GitHub/Anaconda/C++_Mac/EXECUTION");

            string name = "NCA_BETATIMES_CHI_GAMMA_";

            std::stringstream gam;
            std::stringstream alp;
            std::stringstream cuof;
            std::stringstream bet;
            std::stringstream gri;


            gam << gamma;
            alp << alpha;
            cuof << k_cutoff;
            bet << MD.tau_grid[MD.tau_grid.size() - 1];
            gri << MD.beta;

            name += gam.str();
            name += "_ALPHA_";
            name += alp.str();
            name += "_MODE_";
            name += cuof.str();
            name += "_BETA_";
            name += bet.str();
            name += "_GRID_";
            name += gri.str();
            name += ".txt";

            outputFile.open(name);
            MD.CAL_COUP_INT_with_g_arr(alpha, k_cutoff);
            MD.Chi_sp(20);

            for (int i = 0; i < MD.beta; i++)
            {
                outputFile << MD.tau_grid[i] << "\t" << MD.tau_grid[(MD.tau_grid.size() - 1) / 2] * Chi_Arr[i] << endl; //변수 a에 값을 할당 후 벡터 각 요소를 반복문으로 불러옴. 이전에는 a 대신 함수를 반복해서 호출하는 방법을 썼는데 그래서 계산 시간이 오래 걸림.
            }
            outputFile.close();

        }
        */
        /**************************************************************************/
    }


    std::chrono::system_clock::time_point P_sec = std::chrono::system_clock::now();
    std::chrono::duration<double> seconds = std::chrono::duration_cast<std::chrono::seconds>(P_sec - P_start);
    cout << "## Total Process ends with : " << seconds.count() << "[sec] ##" << endl;
    cout << "-----------------------------" << endl;

    return 0;

}
