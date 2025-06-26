#include<iostream>
#include<eigen3/Eigen/Dense>
#include<cmath>

using namespace Eigen;
using namespace std;

vector<double> value_list = {1,2,3,4,5,6,7,8,9};

MatrixXd Hankel(vector<double> list){
    vector<double> new_arr(18,0);

    for (int i = 0; i<9; i++){
        new_arr[i] = list[i];
        new_arr[i+9] = list[i] + list[i];
    }

    //test list element

    for (int i = 0; i < 18; i++){
        cout << new_arr[i] << endl;
    }

    MatrixXd blank = MatrixXd::Zero(9,9);

    for (int i = 0; i< 9; i++){
        for (int j = 0; j < 9; j++){
            blank(i,j) = new_arr[i+j];
        }
    }

    return blank;
}

MatrixXd Circ(vector<double> lis){
    vector<double> list(10,0);
    for (int i = 0; i < 9; i++){
        list[i] = lis[i];
    }

    for (int n = 0; n < 10; n++){
        cout << list[n] << endl;
    }

    cout << "##" << endl;

    MatrixXd blank = MatrixXd::Zero(9,9);

    for (int i = 0; i< 9; i++){
        for (int j = 0; j < 9; j++){
            blank(i,j) = list[j];
        }

        list[9] = list[0];

        for (int k = 0; k < 10; k++){
            list[k] = list[k+1];
        }

        // check sort 
        for (int n = 0; n < 10; n++){
            cout << list[n] << endl;
        }
        cout << "---" << endl;
    }

    return blank;

}

int main()
{
    MatrixXd Han = Circ(value_list);
    cout << Han;
    return 0;

}