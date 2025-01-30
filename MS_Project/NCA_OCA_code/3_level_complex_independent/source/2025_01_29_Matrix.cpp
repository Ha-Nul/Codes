#include <iostream>
#include <fstream>
#include <Eigen/Dense>
#include <complex>
#include <vector>

using namespace std;
using namespace Eigen;

MatrixXcd Make_N_Matrix(int lineNumber) {
    ifstream Refile("/Users/e2_602_qma/Documents/GitHub/Codes/Codes/MS_Project/Matheiu/12_09_Mathieu/Matrix_data/realtest.txt");
    ifstream Imfile("/Users/e2_602_qma/Documents/GitHub/Codes/Codes/MS_Project/Matheiu/12_09_Mathieu/Matrix_data/imagtest.txt");
    
    if (!Refile.is_open() || !Imfile.is_open()) {
        cerr << "Error opening files!" << endl;
        return MatrixXcd::Zero(3, 3); // Or handle the error in a different way
    }

    string realLine, imagLine;
    MatrixXcd Nmatrix = MatrixXcd::Zero(3, 3);

    // lineNumber 행까지 읽기 (두 파일 모두 lineNumber 행까지 읽어야 함)
    for (int i = 0; i < lineNumber; ++i) {
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
