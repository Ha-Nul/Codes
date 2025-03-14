#include <iostream>
#include <fstream>
#include <Eigen/Dense>
#include <vector>
#include <cmath>
#include <OCA_complex.hpp>
#include <chrono>
#include <complex>

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

void MD_OC::readVfunc(double gamma, double alpha)
{
    INT_Arr.resize(t);
    
    string filename = "./data/INT_Arr_g";

    std::stringstream gam;
    std::stringstream alp;

    gam << gamma;
    alp << alpha;

    filename += gam.str();
    filename += "_a";
    filename += alp.str();
    filename += ".dat";
    
    ifstream readFile;
    readFile.open(filename);

    if (readFile.is_open())
    {
        cout << "File opened!:" << filename << endl;
        int i = 0;
        string str;
        while (getline(readFile, str))
        {   
            if (str.empty())
            {
                break;
            }
            INT_Arr[i] = stod(str);
            i += 1;
        }
        readFile.close();
    }
    /*
    for (int i =0 ; i<t; i++){
        cout << INT_Arr[i] << endl;
    }
    */

}

MatrixXcd MD_OC::Make_N_Matrix(int lineNumber, double ng) {
    
    string Refilename = "./data/M_c_gam_0to0.4_";
    string Imfilename = "./data/M_c_gam_0to0.4_";

    std::stringstream characval;

    characval << ng;

    Refilename += characval.str();
    Imfilename += characval.str();
    Refilename += "_Re.txt";
    Imfilename += "_Im.txt";
    

    ifstream Refile(Refilename);
    ifstream Imfile(Imfilename);

    if (!Refile.is_open() || !Imfile.is_open()) {
        cerr << "Error opening files!" << endl;
        return MatrixXcd::Zero(3, 3); // Or handle the error in a different way
    }

    string realLine, imagLine;
    MatrixXcd Nmatrix = MatrixXcd::Zero(3, 3);

    // lineNumber 행까지 읽기 (두 파일 모두 lineNumber 행까지 읽어야 함)
    for (int i = 0; i < lineNumber+1; ++i) {
        if (!getline(Refile, realLine) || !getline(Imfile, imagLine)) {
            cerr << "Error reading lines from files!" << endl;
            return MatrixXcd::Zero(3, 3); // Or handle the error in a different way
        }
    }

    // 실수부와 허수부를 공백(탭)으로 구분하여 읽기
    vector<double> realParts, imagParts;
    istringstream realIss(realLine), imagIss(imagLine);
    double realPart, imagPart;

    while (realIss >> realPart) {
        realParts.push_back(realPart);
    }

    while (imagIss >> imagPart) {
        imagParts.push_back(imagPart);
    }

    // 두 벡터의 크기가 같은지 확인
    if (realParts.size() != imagParts.size() || realParts.size() != 9) {
        cerr << "Error: Inconsistent number of elements or not enough elements for 3x3 matrix." << endl;
        return MatrixXcd::Zero(3, 3); // Or handle the error in a different way
    }

    // 복소수 행렬 생성
    int index = 0;
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            Nmatrix(i, j) = complex<double>(realParts[index], imagParts[index]);
            index++;
        }
    }

    return Nmatrix;

}


MatrixXcd MD_OC::Make_Loc_Matrix(int lineNumber, double ng){
    
    string filename = "./data/M_c_H_loc_gam_0to0.4_";

    std::stringstream characval;

    characval << ng;

    filename += characval.str();
    filename += ".txt";
    
    ifstream file(filename);

    string line;
    MatrixXd Locmatrix = MatrixXd::Zero(3,3);
    int i = 0;

    if (file.is_open()){
        cout << "   ** H_loc file opened!: " << filename << endl;
        for (int i = 0; i <= lineNumber; ++i){
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


void MD_OC::data_store(int lineNumber,double ng){
    cout << "** Data store Executed!" <<endl;
    H_loc = Make_Loc_Matrix(lineNumber,ng);
    H_N = Make_N_Matrix(lineNumber,ng);

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
    MatrixXcd Stmp;
    for (int i = 0; i < t; i++)
    {
        int loop = 0;
        Stmp = MatrixXcd::Zero(3, 3);
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


MatrixXcd MD_OC::round_propagator_ite(const MatrixXcd& loc, const vector<MatrixXcd>& sigma, const vector<MatrixXcd>& ite, int n, int boolean)
{

    MatrixXcd sigsum = MatrixXcd::Zero(3, 3);

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

    MatrixXcd Bucket = MatrixXcd::Zero(3, 3);
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



vector<MatrixXcd> MD_OC::Propagator(const vector<MatrixXcd>& sig, const MatrixXcd& loc)
{
    vector<MatrixXcd> P_arr(t, MatrixXcd::Zero(3, 3));
    vector<MatrixXcd> S_arr(t, MatrixXcd::Zero(3, 3));

    P_arr[0] = MatrixXcd::Identity(3, 3);
    S_arr[0] = MatrixXcd::Identity(3, 3);

    MatrixXcd sig_form = MatrixXcd::Zero(3, 3);
    MatrixXcd sig_late = MatrixXcd::Zero(3, 3);

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

complex<double> MD_OC::chemical_poten(MatrixXcd prop)
{
    complex<double> Trace = prop.trace();
    complex<double> lambda = -(1 / tau_grid[t - 1]) * log(Trace);

    return lambda;
}

///////////////////////////////////////////////////////////////////////////////

vector<MatrixXcd> MD_OC::Iteration(const int& n)
{
    cout << "** Iteration RUN " << endl;

    Prop.resize(t, MatrixXcd::Zero(3, 3));
    Prop[0] = MatrixXcd::Identity(3, 3);
    MatrixXcd Iden = MatrixXcd::Identity(3, 3);

    vector<complex<double> > lambda(n + 1, 0);
    complex<double> expDtauLambda;
    complex<double> factor;

    for (int i = 0; i <= n; i++)
    {
        if (i == 0)
        {
            for (int j = 0; j < t; j++)
            {
                Prop[j](0, 0) = std::complex<double>(exp(-tau_grid[j] * H_loc(0, 0)));
                Prop[j](1, 1) = std::complex<double>(exp(-tau_grid[j] * H_loc(1, 1)));
                Prop[j](2, 2) = std::complex<double>(exp(-tau_grid[j] * H_loc(2, 2)));
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

void MD_OC::NCA_Chi_sp(vector<MatrixXcd>& iter)
{
    Chi_Arr.resize(t);
    MatrixXcd GELL = MatrixXcd::Zero(3, 3);
    GELL(0, 1) = 1;
    GELL(1, 0) = 1;

    for (int i = 0; i < t; i++)
    {
        Chi_Arr[i] = (iter[t - i - 1] * GELL * iter[i] * GELL).trace();
    }
}

void MD_OC::OCA_store(vector<MatrixXcd>& iter)
{
    MatrixXcd GELL = MatrixXcd::Zero(3, 3);
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


void MD_OC::OCA_Chi_sp(vector<MatrixXcd>& iter)
{
    for (int i = 0; i < t; i++)
    {
        MatrixXcd Stmp = MatrixXcd::Zero(3, 3);

        for (int n = 0; n <= i; n++) for (int m = i; m < t; m++)
        {
            Stmp += INT_Arr[m - n] * (Chi_st[t - i - 1][m - i] * Chi_st[i][n]);
            //cout << "pair ("<<n<<","<<m<<") is : " << "(" << k-i-1 << "," << m-i << ")"<< " with " << "(" << i << "," << n << ")" << endl;
        }
        Chi_Arr[i] += pow(Delta_t, 2) * Stmp.trace();
    }
}

vector<complex<double> > MD_OC::Chi_sp_Function(vector<MatrixXcd> ITE)
{
    NCA_Chi_sp(ITE);
    OCA_store(ITE);
    OCA_Chi_sp(ITE);

    return Chi_Arr;
}

////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv)
{
    double beta = 10;
    int grid = 401;

    int al = 21;

    MD_OC MD(beta,grid);

    vector<double> gam_arr(al,0);
    vector<double> alp_arr(al,0);
    

    for (int i = 0 ; i < al; i++){
        if (i==0){
            gam_arr[i] == 0;
        }
        else{
            gam_arr[i] += gam_arr[i-1] + 0.02;
        }
    }

    for (int i = 0 ; i < al; i++){
        if (i==0){
            alp_arr[i] == 0;
        }
        else{
            alp_arr[i] += alp_arr[i-1] + 0.05;
        }
    }

    double charval;
    charval = 0;

    //////////////// MPI code activate //////////////////////////////

    MPI_Init(&argc, &argv);

    int world_rank, world_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);


    for (int i = 0; i<al; i++){
        for(int j = 0; j<al; j++){

            MD.readVfunc(gam_arr[i],alp_arr[j]);
            MD.data_store(i,charval);

            vector<MatrixXcd> ITER = MD.Iteration(20);
            vector<complex<double> > a = MD.Chi_sp_Function(ITER);

            std::ofstream outputFile;
            string Chi_name = "M_OCA_CHI_GAMMA_";

            std::stringstream gam;
            std::stringstream alp;
            std::stringstream bet;
            std::stringstream gri;
            std::stringstream mod;
            std::stringstream cv;

            gam << gam_arr[i];
            alp << alp_arr[j];
            mod << 30000;
            bet << 10;
            gri << MD.t;
            cv << charval;


            Chi_name += gam.str();
            Chi_name += "_ALPHA_";
            Chi_name += alp.str();
            Chi_name += "_BETA_";
            Chi_name += bet.str();
            Chi_name += "_NG_";
            Chi_name += cv.str();
            Chi_name += "_GRID_";
            Chi_name += gri.str();
            Chi_name += ".txt";

            outputFile.open(Chi_name);


            if (outputFile.is_open()) { // 파일 열기 확인
                for (int j = 0; j < MD.tau_grid.size(); j++) {
                outputFile << MD.tau_grid[j] << "\t" << a[j].real() << endl; // a[j].real()로 실수부 추출
            }
            outputFile.close();
            } else {
            std::cerr << "Unable to open file: " << Chi_name << endl;
            // 오류 처리 (예: return 1, 예외 발생 등)
            }

            outputFile.close();
        }
    }
}
