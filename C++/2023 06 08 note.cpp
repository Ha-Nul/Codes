#include <iostream>
using namespace std;

/*
int main()
{
int test[5] = {80, 60, 22, 50, 75};

for(int i=0; i<5; i++)
{
    cout << i+1 << "번째 사람의 점수는" << test[i] << "입니다 \n";
}

return 0;

}
*/

/*
int main()
{
    const int num = 5;
    int test [num];
    cout << num << "명의 점수를 입력하십시오.\n";

    for(int i=0; i<num; i++)
    {
        cin >> test[i];
    }    

    for(int j=0; j<num; j++)
    {
        cout << j+1 << "번째 사람의 점수는" << test[j] << "입니다\n";
    }

    return 0;
    
    
}

*/


/*
int main()
{
    const int num = 5;
    int test[num];

    cout << num << "명의 점수를 입력하십시오 \n";
    for(int i=0; i<num; i++)
    {
        cin >> test[i];
    }

    for(int s=0; s<num-1; s++)
    {
        for(int t=s+1; t<num; t++)
        {
            if(test[t] > test[s])
            {
                int tmp = test[s];
                test[t] = test[s];
                test[s] = tmp;
            }
        }
    }

    for(int j=0; j<num; j++)
    {
        cout << j+1 << "번째 사람의 점수는" << test[j] << "입니다. \n";
    }

    return 0;
}

*/

/*
int main()
{
    int test[][5] =
    {{80,60,22,50,75},{90,55,68,72,58}};

    for(int i=0; i<5; i++)
    {
        cout << i+1 << "번째 사람의 국어 점수는" << test[0][i] << "이고 수학 점수는" << test[1][i] << "입니다\n";
    }

    return 0;
}
*/

int main()
{
    int a,b;
    cin >> a >> b;

    int c[a];

    int d,e,f;

    for(int i; i<b; i++)
    {
        cin >> d >> e >> f;
        for(int j; d<=j<=e;j++)
        {
            c[j] = f;
        }

    }
    
    for(int k=0;k<b;k++)
    {
        cout << c[k];
    }
    return 0;
}
