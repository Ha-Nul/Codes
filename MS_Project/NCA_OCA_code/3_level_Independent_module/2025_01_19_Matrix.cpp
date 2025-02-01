#include <iostream>
#include <fstream>
#include <Eigen/Dense>
#include <vector>

using namespace std;
using namespace Eigen;

MatrixXd Make_N_Matrix(int lineNumber){
    ifstream file("/Users/e2_602_qma/Documents/GitHub/Codes/Codes/MS_Project/Matheiu/12_09_Mathieu/Matrix_data/M_N_gam_0to1.txt");
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

MatrixXd Make_Loc_Matrix(int lineNumber){
    ifstream file("/Users/e2_602_qma/Documents/GitHub/Codes/Codes/MS_Project/Matheiu/12_09_Mathieu/Matrix_data/M_H_loc_gam_0to1.txt");
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


int main() {
    cout << "Test program" << endl;
    //cout << "---------------------------------------------" << endl;
    //Make_N_Matrix(2);

    cout << "---------------------------------------------" << endl;
    cout << Make_N_Matrix(2) << endl;

    cout << "---------------------------------------------" << endl;
    cout << Make_Loc_Matrix(2) << endl;

}
