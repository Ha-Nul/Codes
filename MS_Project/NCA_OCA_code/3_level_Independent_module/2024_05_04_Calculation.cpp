#include <iostream>
#include <fstream>
#include <Eigen/Dense>
#include <vector>
#include <cmath>
#include <OCA_bath.hpp>
#include <chrono>

using namespace std;
using namespace Eigen;

MD_OC::MD_OC(double beta, int grid)
    : tau_grid(linspace(0, beta, grid)), t(grid - 1)
{
    Delta_t = tau_grid[1] - tau_grid[0];
    t = tau_grid.size();
}

MD_OC::~MD_OC()
{
    //blank;
}
////////////////////////////////////////////////////////////////////////////////

void MD_OC::readVfunc()
{
    INT_Arr.resize(t);

    ifstream readFile;
    readFile.open("/home/way_ern/Programs/Github/run/LOW_exe/20250120/INT_Arr_g1_a1.dat");

    if (readFile.is_open())
    {
        int i = 0;
        while (!readFile.eof())
        {
            string str;
            getline(readFile, str);
            if (str == "")
            {
                break;
            }
            INT_Arr[i] = stod(str);
            i += 1;
        }
        readFile.close();
    }

}

MatrixXd MD_OC::Make_N_Matrix(int lineNumber){
    ifstream file("/home/way_ern/Programs/Github/Codes/Codes/MS_Project/Matheiu/12_09_Mathieu/Matrix_data/M_N_gam_0to1.txt");
    string line;
    MatrixXd Nmatrix = MatrixXd::Zero(3,3);
    int index = 1;

    if (file.is_open()){
        for (int i = 1; i <= lineNumber; ++i){
            getline(file,line);
            if (i == lineNumber){
                break; // string 자료형으로 받아왔으므로 line[0], line[1] 등에는 공백이나 문자의 형태의 데이터가 저장되어 있음
            }
        }
    }
    file.close();

    vector<double> elements;
    istringstream iss(line);
    double value;

    while (iss >> value) {
        elements.push_back(value);
    }

    for (int i = 0 ; i < 3; i++) for (int j = 0 ; j < 3 ; j++){
        Nmatrix(i,j) = elements[index];
        index += 1;
    }

    return Nmatrix;
}

MatrixXd MD_OC::Make_Loc_Matrix(int lineNumber){
    ifstream file("/home/way_ern/Programs/Github/Codes/Codes/MS_Project/Matheiu/12_09_Mathieu/Matrix_data/M_H_loc_gam_0to1.txt");
    string line;
    MatrixXd Locmatrix = MatrixXd::Zero(3,3);
    int i = 0;

    if (file.is_open()){
        for (int i = 1; i <= lineNumber; ++i){
            getline(file,line);
            if (i == lineNumber){
                break; // string 자료형으로 받아왔으므로 line[0], line[1] 등에는 공백이나 문자의 형태의 데이터가 저장되어 있음
            }
        }
    }
    file.close();

    vector<double> elements;
    istringstream iss(line);
    double value;

    while (iss >> value) {
        elements.push_back(value);
    }

    Locmatrix(0,0) = elements[1];
    Locmatrix(1,1) = elements[2];
    Locmatrix(2,2) = elements[3];

    return Locmatrix;
}

void MD_OC::data_store(int lineNumber){
    readVfunc();
    H_loc = Make_Loc_Matrix(lineNumber);
    H_N = Make_N_Matrix(lineNumber);

    cout << "This is H_loc matrix: " << H_loc << endl;
    cout << "This is a N matrix: " << H_N << endl;
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
        Stmp = MatrixXd::Zero(3, 3);
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



vector<MatrixXd> MD_OC::Propagator(const vector<MatrixXd>& sig, const MatrixXd& loc)
{
    vector<MatrixXd> P_arr(t, MatrixXd::Zero(3, 3));
    vector<MatrixXd> S_arr(t, MatrixXd::Zero(3, 3));

    P_arr[0] = MatrixXd::Identity(3, 3);
    S_arr[0] = MatrixXd::Identity(3, 3);

    MatrixXd sig_form = MatrixXd::Zero(3, 3);
    MatrixXd sig_late = MatrixXd::Zero(3, 3);

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

    Prop.resize(t, MatrixXd::Zero(3, 3));
    Prop[0] = MatrixXd::Identity(3, 3);
    MatrixXd Iden = MatrixXd::Identity(3, 3);

    vector<double> lambda(n + 1, 0);
    double expDtauLambda;
    double factor;

    for (int i = 0; i <= n; i++)
    {
        if (i == 0)
        {
            for (int j = 0; j < t; j++)
            {
                Prop[j](0, 0) = exp(-tau_grid[j] * H_loc(0, 0));
                Prop[j](1, 1) = exp(-tau_grid[j] * H_loc(1, 1));
                Prop[j](2, 2) = exp(-tau_grid[j] * H_loc(2, 2));
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

//////////////////////////////////////////////////////////////////////////////

void MD_OC::NCA_Chi_sp(vector<MatrixXd>& iter)
{
    Chi_Arr.resize(t);
    MatrixXd GELL = MatrixXd::Zero(3, 3);
    GELL(0, 1) = 1;
    GELL(1, 0) = 1;

    for (int i = 0; i < t; i++)
    {
        Chi_Arr[i] = (iter[t - i - 1] * GELL * iter[i] * GELL).trace();
    }
}

void MD_OC::OCA_store(vector<MatrixXd>& iter)
{
    MatrixXd GELL = MatrixXd::Zero(3, 3);
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
        MatrixXd Stmp = MatrixXd::Zero(3, 3);

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

////////////////////////////////////////////////////////////////////////////////////

int main()
{
    double beta = 10;
    int grid = 101;

    MD_OC MD(beta,grid);

    MD.data_store(2);

    vector<MatrixXd> ITER = MD.Iteration(25);
    vector<double> a = MD.Chi_sp_Function(ITER);

    std::ofstream outputFile;
    string name = "Testing.txt";

    outputFile.open(name);

    for (int j = 0; j < MD.tau_grid.size(); j++)
    {
        outputFile << MD.tau_grid[j] << "\t" << MD.tau_grid[MD.t - 1] * a[j] << endl;
    }

    outputFile.close();

}
