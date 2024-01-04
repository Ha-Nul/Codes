#include<iostream>
using namespace std;

int main()
{
    string a;
    getline(cin,a);

    int b = 0;
 
    for(int i = 0; i < a.length(); i++)
    {
        if(a[i]==' ' && a.length()!=1 )
        {
            b = b+1;
        }
    }

    if(a[0]==' ' && a[a.length()-1]==' ')
    {
        if(a.length()==1)
        {
            cout << b;
        }
        else
        {
            cout << b-1;
        }
    }

    else if(a[0]==' ')
        {
            cout << b << endl;
        }

    else if(a[a.length()-1]==' ')
        {
            cout << b << endl;
        }
    
    else
        {
            cout << b+1 << endl;
        }

    return 0;

}