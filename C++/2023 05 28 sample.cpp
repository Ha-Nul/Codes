#include <iostream>
using namespace std;

int main(){
    int num = 1;

    while(num){

        cout << "정수를 입력하십시오. (0을 입력하면 종료합니다)\n";
        cin >> num;
        cout << num << "이(가) 입력되었습니다\n";

    }

    cout << "반복이 끝났습니다 \n";


    return 0;
}