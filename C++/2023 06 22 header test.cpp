#include<iostream>
using namespace std;

int main()
{

int a = 1;
int b = 1;

while(a!=0 & b!=0)
{

    cin >> a >> b;

    if(a==0 & b==0)
        break;
    
    cout << a+b <<"\n";

}

return 0;

}