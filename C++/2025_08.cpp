#include <iostream>
#include<Eigen/Dense>

using namespace std;
using namespace Eigen;

int main() {

// Cave construction
    int cave_Size = 0;

    cin >> cave_Size;

    if (cave_Size == 0)
    {
        None;
    }
    else
    {
        MatrixXd cave = MatrixXd::Zero(cave_Size, cave_Size);

        for (int i = 0 ; i< cave_Size; i++){
            for (int j = 0; j < cave_Size; j++){
                cin >> cave[i][j];
            }
        }
    }

// Game start
    for (int i = 0 ; i < cave_Size ; i++){
        if (i==0){
            if cave[i+1][i] < cave[i+1][i]
            
        }
    }

    return 0;

}
