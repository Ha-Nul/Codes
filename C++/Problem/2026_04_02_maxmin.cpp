#include<iostream>
#include<vector>

using namespace std;

int main(){
    vector<double> arra(10);

    for (int i = 0; i < arra.size() ; i++){

        arra[i] = 0;

        cout << arra[i] << endl;
        cout << i << endl;
    }

    cout << "arra.size()";
}