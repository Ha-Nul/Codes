#include<iostream>
#include<Eigen/Dense>
#include<vector>
#include<unsupported/Eigen/KroneckerProduct> // Eigen의 내장 크로네커 곱 함수 사용
#include<cmath>

using namespace Eigen;
using namespace std;

MatrixXd max_Entangled_Vec(int n){
    MatrixXd loc_SYS = MatrixXd::Zero(1,n);
    MatrixXd res_SYS = MatrixXd::Zero(1,n);

    MatrixXd store = MatrixXd::Zero(1,n*n);

    for (int i = 0 ; i < n ; i++){
        loc_SYS(0,i) = i;
        res_SYS(0,i) = i;
    }

    // 직접 구현한 kron 대신 Eigen의 내장 함수 kroneckerProduct 사용
    store = kroneckerProduct(res_SYS,loc_SYS).eval();
    
    return store; // 함수가 MatrixXd를 반환하도록 return 문 추가
}

int main()
{
    MatrixXd result = max_Entangled_Vec(3);
    cout << "max_Entangled_Vec(3) 결과:\n" << result << endl;
    return 0;
}
