#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <vector>
#include <cmath>

using namespace std;
using namespace Eigen;

vector<double> linspace(double min, double max, int n)
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

//const variable (젼역변수)
const vector<double> const_linspace_test = linspace(1,10,10);
const double const_tau = const_linspace_test[0];

vector<double> Testing_green(vector<double> x, int n, double tau) //vector<double> tau) 
{
    double T = 273;
    vector<int> one_vec(n,1);
    vector<double> bose_dist(n);

    //cout << "Testing_bose" << endl;

    for (int i = 0; i < n; i++)
    {
        bose_dist[i]=one_vec[i]/(exp(2 * x[i])-1);
        //result[i] = one_vec[i] / exp(x[i]);
        //cout << bose_dist[i] << endl;
    }

    vector<double> Test_green(n);

    //cout << "Test_Green\n" << endl;

    for (int j = 0; j < n; j++)
    {
        Test_green[j] = ((bose_dist[j] + 1)*exp(-1 * x[j] * tau) + (bose_dist[j])*exp(x[j] * tau));
        //Test_green[j] = ((bose_dist[j] + 1)*exp(-1 * x[j] * tau[j]) + (bose_dist[j])*exp(x[j] * tau[j]));
        //cout << Test_green[j] << endl;
    }


    return Test_green;
}

vector<double> Testing_coupling(double v, double g, double W, int n)
{
    vector<double> v_array(n,v);
    vector<double> g_array(n,g);
    vector<double> W_array(n,W);
    vector<double> coupling_array(n);

    //cout << "Testing_coupling\n" << endl;

    for (int i = 0; i < n ; i++)
    {
        coupling_array[i] = g_array[i] * sqrt(abs(const_linspace_test[i]) * v_array[i]/(1 + pow((abs(const_linspace_test[i]) * v_array[i]/W_array[i]),2)));
        //cout << coupling_array[i] << endl;
    }
    
    return coupling_array;
    //파이썬 결과랑 경향은 비슷한데, 결과값이 다르게 나옴. 이거 왜 이런다냐
    //coupling이 현재 문제임. 파이썬 결과랑 다름. 나머지는 오차가 소수점 표시하는 마지막에서 반올림하면서 생기는데 얘만 값이 완전히 다름. 왜지?
    //분모에서 제곱하는 term을 빼먹어서 그랬음 ^-^... ㅎㅎ....
}

double Testing_Interact(vector<double> coupling, vector<double> linspace, vector<double> green, int n)
{
    MatrixXd blank_matrix = MatrixXd::Zero(n,n);
    double blank_factor = 0;

    for (int i = 0; i < n; i++)
    {
        blank_matrix(0,i)= (coupling[i] *coupling[i]) * green[i];
        //row of matrix indicates each time interval. 0th row means first time interval, t_0 to t_1.
    }

    //cout << blank_matrix << endl;
    //역시 경향성은 비슷하게 나오는데 값이 다름. 왜 이런거지? Coupling 에서부터 뭐가 달라진 것 같긴 한데
    blank_factor = blank_matrix.sum();

    return blank_factor;
}

MatrixXd Testing_Matrix_Odd(int n, double r, int m)
{
    MatrixXd Matrix1(n,n);
    MatrixXd Matrix2;
    MatrixXd Matrix3;

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

    SelfAdjointEigenSolver<MatrixXd> es(Matrix1);
    Matrix2 = es.eigenvalues();
    Matrix3 = es.eigenvectors();

    if (m==0)
    {
        return Matrix2;
    }
    if (m==1)
    {
        return Matrix3;
    }
}

MatrixXd Testing_Matrix_Even(int n, double r, int m)
{
    MatrixXd Matrix1(n,n);
    MatrixXd Matrix2;
    MatrixXd Matrix3;

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

    SelfAdjointEigenSolver<MatrixXd> es(Matrix1);
    Matrix2 = es.eigenvalues();
    Matrix3 = es.eigenvectors();

    if (m==0)
    {
        return Matrix2;
    }
    if (m==1)
    {
        return Matrix3;
    }
}

//3d로 되어 있는 것들 Xd로 바꿔야 함. 어떻게 바꾸지?

MatrixXd Testing_Hamiltonian_N(MatrixXd even, MatrixXd odd, double g)
{
    //수정 필요. 도대체 어디서 틀린거임 아오...
    // eigenvector는 각 케이스(even, odd)에 맞게 정확하게 나오고 있음. 근데 계산 결과가 이상함. 차이가 심함... 아니면 파이썬 N이 잘못된건가
    // MatrixXd로 받는 행렬은 return 조건을 꼭 m = 1로 둘 것.
    MatrixXd odd_eigenvec;
    MatrixXd even_eigenvec;

    odd_eigenvec = odd;
    even_eigenvec = even;

    //cout << odd_eigenvec << endl;
    //cout << even_eigenvec << endl;

    MatrixXd c;
    c = odd_eigenvec * even_eigenvec;

    //cout << c << endl;

    Matrix3d d = Matrix3d::Zero(3,3);

    d(0,1) = g * c(0,0);
    d(1,0) = g * c(0,0);
    d(1,2) = g * c(0,1);
    d(2,1) = g * c(0,1);

    return d;
}

MatrixXd Testing_Hamiltonian_exp(MatrixXd a, MatrixXd b, double temp)
{
    //g_0 
    MatrixXd Even = a;
    MatrixXd Odd = b;

    double zeroth = exp(Even(0));
    double first = exp(Odd(0));
    double second = exp(Even(1));

    MatrixXd Hamiltonian_exp;
    Hamiltonian_exp.setZero();

    Hamiltonian_exp(0,0) = temp * zeroth;
    Hamiltonian_exp(1,1) = temp * first;
    Hamiltonian_exp(2,2) = temp * second;

    return Hamiltonian_exp;
}

//time step (interval) 에 대한 propagater의 계산. Matrix3d로 되어 있는 것들 Xd로 바꿔야 함.

Matrix3d Testing_Hamiltonian_loc(MatrixXd a, MatrixXd b)
{
    //아 로컬해밀토니안 만든다고 했구나
    Matrix3d Hamiltonian;
    Hamiltonian.setZero();

    Hamiltonian(0,0) = a(0);
    Hamiltonian(1,1) = b(0);
    Hamiltonian(2,2) = a(1);

    return Hamiltonian;
}

MatrixXd round_propagater_ite_1(MatrixXd loc, MatrixXd sigma, int n, double tau)
{
    MatrixXd prop_iter_zero = MatrixXd::Identity(n,n);
    MatrixXd H_local = loc;
    MatrixXd Test_equation;

    Test_equation = -1 * H_local * prop_iter_zero - tau * sigma * prop_iter_zero;

    return Test_equation;
}

MatrixXd Testing_Sigma(MatrixXd N, MatrixXd prop_zero, double V)
{
    MatrixXd Test_N = N;
    MatrixXd Test_g_0 = prop_zero;
    double k_sum = V;

    MatrixXd Test_equation;
    Test_equation = 0.5 * k_sum * Test_N * Test_g_0 * Test_N;

    return Test_equation;
}




int main()
{
    MatrixXd prop_iter_one;
    MatrixXd prop_round;
    MatrixXd prop_Identity;
    double tau_interval;

    // 와 조립하는 것도 한세월... 너무 멍청하게 짠거 같은데?

    prop_Identity = MatrixXd::Identity(3,3);
    tau_interval = const_tau;
    prop_round = round_propagater_ite_1(Testing_Hamiltonian_loc(Testing_Matrix_Even(3,1,0),Testing_Matrix_Odd(3,1,0)),
                Testing_Sigma(Testing_Hamiltonian_N(Testing_Matrix_Even(3,1,1),Testing_Matrix_Odd(3,1,1),1),
                                Testing_Hamiltonian_exp(Testing_Matrix_Even(3,1,0),Testing_Matrix_Odd(3,1,0),const_tau),
                                Testing_Interact(Testing_coupling(0.2,1,10,10),const_linspace_test,Testing_green(const_linspace_test,10,const_tau),10)),3,const_tau);
    prop_iter_one = prop_Identity - const_tau * prop_round;

    cout << prop_iter_one << endl;

    return 0;
}