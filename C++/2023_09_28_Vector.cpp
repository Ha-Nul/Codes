#include <cmath>
#include <iostream>
#include <vector>
#include <Eigen/Dense>
using namespace std;
using namespace Eigen;

//linspace 함수는 인터넷에서 긁어옴...
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

vector<double> Testing_green(vector<double> x, int n, double tau) //vector<double> tau) 
{
    double T = 273; // You can adjust the value of T as needed
    vector<int> one_vec(n,1);
    vector<double> bose_dist(n);

    cout << "Testing_bose" << endl;

    for (int i = 0; i < n; i++)
    {
        bose_dist[i]=one_vec[i]/(exp(2 * x[i])-1);
        //result[i] = one_vec[i] / exp(x[i]);
        cout << bose_dist[i] << endl;
    }

    vector<double> Test_green(n);

    cout << "Test_Green\n" << endl;

    for (int j = 0; j < n; j++)
    {
        Test_green[j] = ((bose_dist[j] + 1)*exp(-1 * x[j] * tau) + (bose_dist[j])*exp(x[j] * tau));
        //Test_green[j] = ((bose_dist[j] + 1)*exp(-1 * x[j] * tau[j]) + (bose_dist[j])*exp(x[j] * tau[j]));
        cout << Test_green[j] << endl;
    }


    return Test_green;
}

vector<double> Testing_coupling(double v, double g, double W, int n)
{
    vector<double> v_array(n,v);
    vector<double> g_array(n,g);
    vector<double> W_array(n,W);
    vector<double> coupling_array(n);

    cout << "Testing_coupling\n" << endl;

    for (int i = 0; i < n ; i++)
    {
        coupling_array[i] = g_array[i] * sqrt(abs(const_linspace_test[i]) * v_array[i]/(1 + pow((abs(const_linspace_test[i]) * v_array[i]/W_array[i]),2)));
        cout << coupling_array[i] << endl;
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
    }

    cout << blank_matrix << endl;
    //역시 경향성은 비슷하게 나오는데 값이 다름. 왜 이런거지? Coupling 에서부터 뭐가 달라진 것 같긴 한데
    blank_factor = blank_matrix.sum();

    return blank_factor;
}

/*
Eigen이랑 Vector 헤더는 서로 Syntax가 맞지 않는걸로. 아 근데 그러면 g(1)계산은 어떻게 함? ㅎㅎ.......
Matrix<double, 10, 1> compatible_test(vector<double> a)
{
    Matrix<double, 10, 1> Lets_check;
    for (int i = 0; i < 10; i++)
    {
        Lets_check << a[i];
    }

    return Lets_check;
}
*/
int main()
{
    vector<double> check = Testing_coupling(1,0.2,10,10);
    //Chatgpt : C++에서는 파이썬처럼 배열을 한번에 출력하는 기능은 없다고 함. For loop을 사용해서 출력.

    double check2 = Testing_Interact(Testing_coupling(1,0.2,10,10),const_linspace_test,Testing_green(const_linspace_test,10,1),10);
    /*
    for (double value : check)//c++11 syntax인데 경고표시 뜸... 
    {
        cout << value << " ";
    }
    */

   cout << check2 << endl;


    return 0; 
}

