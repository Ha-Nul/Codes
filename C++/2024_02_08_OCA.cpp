#include<iostream>
#include <fstream>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <vector>
#include <cmath>
#include <iomanip>
#include <string>

using namespace std;
using namespace Eigen;

vector<double> k_mode(100,1);
double g_ma = 1;
double omega = 1;
double velocity = 1;
double cutoff = 1;

class Testing
{
    private:

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

        vector<MatrixXd> convolve(const vector<MatrixXd>& Signal,
              const vector<MatrixXd>& Kernel, int n, int i)
        {
        size_t SignalLen = i;
        size_t KernelLen = Kernel.size();
        size_t ResultLen = SignalLen + KernelLen - 1;

        vector<MatrixXd> Result(ResultLen,MatrixXd::Zero(n,n));

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
            MatrixXd Matrix1(n,n);

            for (int i=0;i < n;i++){
                for (int j=0;j < n;j++)
                {
                    try
                    {
                        if(i==j){
                            Matrix1(i,j) = pow((i+1),2);
                        }
                        if(abs(i - j) == 1){
                            Matrix1(i,j) = -r/2.0;
                        }
                    }
                    catch (...) {}
                }
            }
            return Matrix1;
        }

        MatrixXd Matrix_Even(int n, double r)
        {
            MatrixXd Matrix1(n,n);

            for (int i=0;i < n;i++){
                for (int j=0;j < n;j++)
                {
                    try
                    {
                        if(i==j){
                            Matrix1(i,j)= pow(i,2);
                        }
                        if(abs(i - j) == 1){
                            Matrix1(i,j) = -r/2.0;
                        }
                    }
                    catch (...) {}
                }
            }
            Matrix1(0,1) = -r/sqrt(2);
            Matrix1(1,0) = -r/sqrt(2);

            return Matrix1;
        }

        vector<double> tau_grid = linspace(0,1,100);
        int k = tau_grid.size();
        double Delta_t = tau_grid[1] - tau_grid[0];

    public:

        const vector<double> grid = linspace(0,1,100);
        vector<double> green(vector<double> tau);
        vector<double> coupling(double v, double g, double W);
        vector<double> Interact(vector<double> coupling, vector<double> tau);
        vector<double> Interact_V(vector<double> couplint, vector<double> tau, double omega);

        MatrixXd Eigenvector_Even();
        MatrixXd Eigenvalue_Even();
        MatrixXd Eigenvector_Odd();
        MatrixXd Eigenvalue_Odd();

        MatrixXd Hamiltonian_N(MatrixXd even, MatrixXd odd, double g);
        vector<MatrixXd> Hamiltonian_exp(MatrixXd a, MatrixXd b);
        MatrixXd Hamiltonian_loc(MatrixXd a, MatrixXd b);
        MatrixXd Hamiltonian_loc_ite(MatrixXd a, MatrixXd b,const double &lambda);

        MatrixXd round_propagator_ite(const MatrixXd &loc, const vector<MatrixXd> &sigma, const vector<MatrixXd> &ite,int weight);
        vector<MatrixXd> NCA_self(const MatrixXd &N,const vector<MatrixXd> &H_exp, const vector<double> &V);
        vector<MatrixXd> OCA_self(MatrixXd &N, vector<MatrixXd> &H_exp, vector<double> &V);
        vector<MatrixXd> Self_E(MatrixXd &N, vector<MatrixXd> &H_exp, vector<double> &V);
        vector<MatrixXd> Propagator(const vector<MatrixXd> &array , const MatrixXd &loc , const double &gvalue);

        double chemical_poten(MatrixXd prop);

        vector<MatrixXd> Iteration(const int &iteration, const double &gvalue);
        vector<double> TestingIteration(const int &n, int testingint);

        vector<double> NCA_Chi_sp(int iteration, const double &gvalue);
        vector<double> OCA_Chi_sp(int iteration, const double &gvalue);
        vector<double> Chi_sp(int ite, const double &g);

};

/////////////////////////////////////////////////////////////////////////////////////

vector<double> Testing::green(vector<double> tau)
{
    double T = 273;
    vector<int> one_vec(k,1); // 원소는 1, 길이는 n 짜리 배열..
    vector<double> bose_dist(k);

    for (int i = 0; i < k; i++)
    {
        bose_dist[i]=one_vec[i]/(exp(tau_grid[tau_grid.size()-1] * k_mode[i])-1);
    }

    vector<double> Test_green(k);

    for (int j = 0; j < tau_grid.size(); j++)
    {
        Test_green[j] = ((bose_dist[j] + 1)*exp(-1 * k_mode[j] * tau[j]) + (bose_dist[j])*exp(k_mode[j] * tau[j]));
    }

    return Test_green;
}

vector<double> Testing::coupling(double v, double g, double W)
{
    vector<double> v_array(k_mode.size(),v);
    vector<double> g_array(k_mode.size(),g);
    vector<double> W_array(k_mode.size(),W);
    vector<double> coupling_array(k_mode.size());

    for (int i = 0; i < k_mode.size() ; i++)
    {
        coupling_array[i] = g_array[i] * sqrt(abs(k_mode[i]) * v_array[i]/(1 + pow((abs(k_mode[i]) * v_array[i]/W_array[i]),2)));
    }
    
    return coupling_array;
}
////////////////////////////////////////////////////////////////////////////////////

vector<double> Testing::Interact_V(vector<double>coupling, vector<double> tau, double omega)
{
    double coupling_const = coupling[0];

    vector<double> hpcos(tau.size(),0);
    vector<double> hpsin(tau.size(),0);
    vector<double> coupling_arr(tau.size(),coupling_const * coupling_const);
    vector<double> V_arr(tau.size(),0);

    for (int i = 0; i < tau.size(); i++)
    {
        hpcos[i] = cosh(tau[i]-tau[tau.size()-1]/2)*omega;
        hpsin[i] = sinh(tau[tau.size()-1] * omega/2);
        V_arr[i] = (coupling_arr[i] * hpcos[i] / hpsin[i]);

        //cout << "this is V_arr " << V_arr[i] << endl;
    }

    return V_arr;
}

////////////////////////////////////////////////////////////////////////////////////

MatrixXd Testing::Eigenvector_Even()
{
	MatrixXd a;

	SelfAdjointEigenSolver<MatrixXd> es(Matrix_Even(3,g_ma));
	a = es.eigenvectors();

	return a;
}

MatrixXd Testing::Eigenvalue_Even()
{
	MatrixXd b;

	SelfAdjointEigenSolver<MatrixXd> es(Matrix_Even(3,g_ma));
	b = es.eigenvalues();

	return b;
}

MatrixXd Testing::Eigenvector_Odd()
{
	MatrixXd a;

	SelfAdjointEigenSolver<MatrixXd> es(Matrix_Odd(3,g_ma));
	a = es.eigenvectors();

	return a;
}

MatrixXd Testing::Eigenvalue_Odd()
{
	MatrixXd b;

	SelfAdjointEigenSolver<MatrixXd> es(Matrix_Odd(3,g_ma));
	b = es.eigenvalues();

	return b;
}

///////////////////////////////////////////////////////////////////////


MatrixXd Testing::Hamiltonian_N(MatrixXd even, MatrixXd odd, double g)
{
    MatrixXd odd_eigenvec;
    MatrixXd even_eigenvec;

    odd_eigenvec = odd.transpose();
    even_eigenvec = even;

    MatrixXd c;
    c = odd_eigenvec * even_eigenvec;

    MatrixXd d = MatrixXd::Zero(3,3);

    d(0,1) = g * c(0,0);
    d(1,0) = g * c(0,0);
    d(1,2) = g * c(0,1);
    d(2,1) = g * c(0,1);

    return d;
}

vector<MatrixXd> Testing::Hamiltonian_exp(MatrixXd a, MatrixXd b)
{
    //g_0 
    MatrixXd Even = a;
    MatrixXd Odd = b;

    double zeroth = exp(Even(0));
    double first = exp(Odd(0));
    double second = exp(Even(1));

    vector<MatrixXd> array_with_Matrix(k);
 
    MatrixXd Hamiltonian_exp;

    for (int i = 0; i < k; i++)
    {
        Hamiltonian_exp = MatrixXd::Zero(3,3);

        Hamiltonian_exp(0,0) = tau_grid[i] * zeroth;
        Hamiltonian_exp(1,1) = tau_grid[i] * first;
        Hamiltonian_exp(2,2) = tau_grid[i] * second;

        array_with_Matrix[i] = Hamiltonian_exp;
    }

    return array_with_Matrix;
}



MatrixXd Testing::Hamiltonian_loc(MatrixXd a, MatrixXd b)
{
    MatrixXd Hamiltonian = MatrixXd::Zero(3,3);

    Hamiltonian(0,0) = a(0);
    Hamiltonian(1,1) = b(0);
    Hamiltonian(2,2) = a(1);

    return Hamiltonian;
}

MatrixXd Testing::Hamiltonian_loc_ite(MatrixXd a, MatrixXd b, const double &lambda)
{
    MatrixXd Hamiltonian = MatrixXd::Zero(3,3);

    Hamiltonian(0,0) = a(0)-lambda;
    Hamiltonian(1,1) = b(0)-lambda;
    Hamiltonian(2,2) = a(1)-lambda;

    return Hamiltonian;
}

////////////////////////////////////////////////////////////////////////////////


vector<MatrixXd> Testing::NCA_self(const MatrixXd &N,const vector<MatrixXd> &Prop, const vector<double> &V)
{

    vector<MatrixXd> Sarray(k);
    
    for (int i=0; i < k ; i++)
    {   
        Sarray[i] = V[i] * (N * Prop[i] * N);
    }
    
    return Sarray;
}

vector<MatrixXd> Testing::OCA_self(MatrixXd &N, vector<MatrixXd> &Prop, vector<double> &V)
{
    vector<MatrixXd> Sarray_tau_0(k,MatrixXd::Zero(3,3));
    vector<MatrixXd> Sarray_tau_1(k,MatrixXd::Zero(3,3));
    vector<MatrixXd> Sarray_tau_2(k,MatrixXd::Zero(3,3));

    for (int i = 0 ; i<k ; i++)
    {
        for (int n=0; n<i; n++)
        {
            for (int m=n ; m<k; m++)
            {
                Sarray_tau_2[m] = N * Prop[k-1-n] * N * Prop[m-n] * N * Prop[m] * N * V[k-1-n] * V[m];
                Sarray_tau_1[n] += Sarray_tau_2[m] * Delta_t;
            }

            Sarray_tau_0[i] += Sarray_tau_1[n] * Delta_t;
        }
    }

    return Sarray_tau_0;
}

vector<MatrixXd> Testing::Self_E(MatrixXd &N, vector<MatrixXd> &Prop, vector<double> &V)
{
    vector<MatrixXd> EArray(k,MatrixXd::Zero(3,3));
    vector<MatrixXd> NCA = NCA_self(N,Prop,V);
    vector<MatrixXd> OCA = OCA_self(N,Prop,V);

    for(int i=0; i<k ; i++)
    {
        EArray[i] = NCA[i] + OCA[i];
    }

    return EArray;
}

//////////////////////////////////////////////////////////////////////////////


MatrixXd Testing::round_propagator_ite(const MatrixXd &loc, const vector<MatrixXd> &Self_Energy, const vector<MatrixXd> &Prop_ite, int n)
{   

    MatrixXd SUM = MatrixXd::Zero(3,3);
    
    if (n == 1)
    {
        SUM = Self_Energy[1]*Prop_ite[0] + Self_Energy[0]*Prop_ite[1];
    }
    else if (n > 1){
        for (int i = 0 ; i < n ; i++)
        {
            SUM += Self_Energy[n-(i)] * Prop_ite[i];

            if (i+1 == n)
            {
                break;
            }

        }
    }

    MatrixXd Bucket = MatrixXd::Zero(3,3);
    Bucket = -loc * Prop_ite[n] + Delta_t * SUM;

    //cout << -loc * ite << endl;
    return Bucket;
}



vector<MatrixXd> Testing::Propagator(const vector<MatrixXd> &Self_E, const MatrixXd &loc, const double &gvalue)
{
    vector<MatrixXd> PArray(k,MatrixXd::Zero(3,3));
    
    PArray[0] = MatrixXd::Identity(3,3);
    
    /*
    vector<double> C = coupling(velocity,gvalue,cutoff);
    vector<double> V = Interact_V(C,tau_grid,omega);
    MatrixXd P_N = Hamiltonian_N(Eigenvector_Even(),Eigenvector_Odd(),gvalue);
    */

    for (int i=1; i < k; i++)
    {
        if (i == 1)
        {
            PArray[1] = PArray[0] + Delta_t * round_propagator_ite(loc,Self_E,PArray,0);
        }

        if (i > 1)
        {
            PArray[i] = PArray[i-1] + Delta_t * round_propagator_ite(loc,Self_E,PArray,i-1);
        }
    }

    return PArray;
}

/////////////////////////////////////////////////////////////////////////////

double Testing::chemical_poten(MatrixXd prop)
{
    double Trace = prop.trace();
    double lambda = -(1/tau_grid[k-1]) * log(Trace);
    
    return lambda;
}

///////////////////////////////////////////////////////////////////////////////

vector<MatrixXd> Testing::Iteration(const int &n, const double &gvalue)
{
    vector<MatrixXd> Sig;
    vector<MatrixXd> Prop;
    vector<MatrixXd> Prop_zeroth(k,MatrixXd::Identity(3,3));
    vector<double> coup = coupling(velocity,gvalue,cutoff);
    vector<double> Int = Interact_V(coup,tau_grid,omega);

    MatrixXd H_loc = Hamiltonian_loc(Eigenvalue_Even(),Eigenvalue_Odd());
    MatrixXd Iden = MatrixXd::Identity(3,3);
    MatrixXd H_N = Hamiltonian_N(Eigenvector_Even(),Eigenvector_Odd(),gvalue);
    vector<MatrixXd> H_e = Hamiltonian_exp(Eigenvalue_Even(),Eigenvalue_Odd());

    double lambda;
    
    for(int i = 0; i <= n; i++)
    {
        if (i==0)
        {   
            Prop = Prop_zeroth;
            for(int j=0; j<k; j++)
            {
                Prop[j](0,0) = exp(-tau_grid[j] * Hamiltonian_loc(Eigenvalue_Even(),Eigenvalue_Odd())(0,0));
                Prop[j](1,1) = exp(-tau_grid[j] * Hamiltonian_loc(Eigenvalue_Even(),Eigenvalue_Odd())(1,1));
                Prop[j](2,2) = exp(-tau_grid[j] * Hamiltonian_loc(Eigenvalue_Even(),Eigenvalue_Odd())(2,2));
            }

            for(int j=0; j<k; j++)
            {
                Prop[j] = Prop[j] * exp(tau_grid[j]*(lambda));
                //cout << Prop[j].trace() << endl;

            }
        }

    
        else
        {   
            H_loc = H_loc - lambda * Iden;

            Sig = Self_E(H_N,Prop,Int);
            Prop = Propagator(Sig,H_loc,gvalue);
            lambda = chemical_poten(Prop[k-1]);
            
            for(int j=0; j<k; j++)
            {
                Prop[j] = Prop[j] * exp(tau_grid[j]*(lambda));
                //cout << Prop[j].trace() << endl;
            }

        }
    
    }

    return Prop;
}

//////////////////////////////////////////////////////////////////////////////

vector<double> Testing::NCA_Chi_sp(int iter, const double &gvalue)
{
    MatrixXd Gellmann_1 = MatrixXd::Zero(3,3);

    Gellmann_1(0,1) = 1;
    Gellmann_1(1,0) = 1;

    vector<double> NCA_chi_array(k,0);
    vector<MatrixXd> Ite_ra = Iteration(iter,gvalue);

    for (int i=0; i<k; i++)
    {
        NCA_chi_array[i] =(Ite_ra[k-i-1] * Gellmann_1 * Ite_ra[i] * Gellmann_1).trace();
        cout << setprecision(16);   
        //cout << chi_array[i] << endl;
    }

    return NCA_chi_array;
}

vector<double> Testing::OCA_Chi_sp(int iter, const double &gvalue)
{
    MatrixXd Gellmann_1 = MatrixXd::Zero(3,3);

    Gellmann_1(0,1) = 1;
    Gellmann_1(1,0) = 1;

    MatrixXd C_N = Hamiltonian_N(Eigenvector_Even(),Eigenvector_Odd(),gvalue);
    vector<double> C_V = Interact_V(coupling(velocity,gvalue,cutoff),tau_grid,omega);

    vector<double> OCA_chi_array0(k,0);
    vector<MatrixXd> OCA_chi_array1(k,MatrixXd::Zero(3,3));
    vector<MatrixXd> OCA_chi_array2(k,MatrixXd::Zero(3,3));

    vector<MatrixXd> Ite = Iteration(iter,gvalue);

    for (int i=0; i<k; i++)
    {
        for (int n=0; n<=i; n++)
        {
            for (int m=i; m<k; m++)
            {
                OCA_chi_array2[m] = Delta_t * C_V[m-n] * Ite[tau_grid[k-1] - m] * C_N * Ite[m - i] * Gellmann_1 * Ite[i-n] * C_N * Ite[n] * Gellmann_1;
                OCA_chi_array1[n] += OCA_chi_array2[m];
            }

            OCA_chi_array0[i] += Delta_t * OCA_chi_array1[n].trace();
            cout << setprecision(16);
        }

    }

    return OCA_chi_array0;

}

vector<double> Testing::Chi_sp(int ite, const double &g)
{
    vector<double> NCA = NCA_Chi_sp(ite,g);
    vector<double> OCA = OCA_Chi_sp(ite,g);

    vector<double> Chi(k,0);

    for (int i=0; i<k ; i++)
    {
        Chi[i] = NCA[i] + OCA[i];
    }

    return Chi;
}


////////////////////////////////////////////////////////////////////////////////////
int main()                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
{

    Testing test;
    /*
    vector<double> g_array(2,0);
    g_array[0] = 0.5;
    g_array[1] = 1.0;
    */

    vector<double> g_array(25,0);
    for (int j=1; j<25; ++j)
    {
        if (j<21)
        {
          g_array[j] = (g_array[j-1] + 0.05);
        }

        else
        {
          g_array[j] = g_array[j-1] + 1;
        }
    }

    for (int m=0; m<21; m++)
    {
        g_array[m] = g_array[m] * g_array[m];
    }
    /*
    for (int n=0; n<25; n++)
    {
        std::ofstream outputFile;

        //string name = "20240111_Trap_beta_0_4_g_";
        string name = "N_matrix_beta_2_g_";
        std::stringstream back;
        back << g_array[n];

        name += back.str();
        name += ".txt";

        outputFile.open(name);

        MatrixXd H_N = test.Hamiltonian_N(test.Eigenvector_Even(),test.Eigenvector_Odd(),g_array[n]);
        outputFile << H_N << endl; //변수 a에 값을 할당 후 벡터 각 요소를 반복문으로 불러옴. 이전에는 a 대신 함수를 반복해서 호출하는 방법을 썼는데 그래서 계산 시간이 오래 걸림.

        outputFile.close();

    }
    */

    /*
    for (int k=0; k<g_array.size(); k++)
    {
        std::ofstream outputFile;

        //string name = "20240111_Trap_beta_0_4_g_";
        string name = "SQ_grid1000_beta_1_g_";
        std::stringstream back;
        back << g_array[k];

        name += back.str();
        name += ".txt";

        outputFile.open(name);
        //vector<double> a = test.Interact_V(test.coupling(velocity,g_array[k],cutoff),test.grid,omega);
        vector<double> a = test.Chi_sp(5,g_array[k]);

        for (int i = 0; i < a.size(); i++)
        {     
            cout << a[i] << endl;
            outputFile << test.grid[i] << "\t" << a[i] << endl; //변수 a에 값을 할당 후 벡터 각 요소를 반복문으로 불러옴. 이전에는 a 대신 함수를 반복해서 호출하는 방법을 썼는데 그래서 계산 시간이 오래 걸림.
        }
        outputFile.close();
    }
    */

   /*
    MatrixXd NN = test.Hamiltonian_N(test.Eigenvector_Even(),test.Eigenvector_Odd(),g_array[10]);
    vector<MatrixXd> He = test.Hamiltonian_exp(test.Eigenvalue_Even(),test.Eigenvalue_Odd());
    MatrixXd Hl = test.Hamiltonian_loc(test.Eigenvalue_Even(),test.Eigenvalue_Odd());
    vector<double> VV = test.Interact_V(test.coupling(velocity,g_array[10],cutoff),test.grid,omega);

    vector<MatrixXd> Self = test.Self_E(NN,He,VV);

    vector<MatrixXd> OCATEST = test.Propagator(Self,Hl,g_array[10]);
    */

    vector<double> ITE = test.Chi_sp(1,g_array[10]);
    
    for (int i = 0; i < test.grid.size() ;  i++)
    {
        //cout << Self[i] << endl;
        cout << ITE[i] << endl;
    }
    
    return 0;

}
