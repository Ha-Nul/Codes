#include <iostream>
using namespace std;

int main()
{

    int b = 0;

    cin >> b;

    int a[b];

    for(int i; i < b; i++)
    {
        a[i] = i+1;
        cout << i << a[i];
    }
}