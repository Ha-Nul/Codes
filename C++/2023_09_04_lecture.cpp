#include<iostream>
#include<string.h>
using namespace std;
//크로스플랫폼 개발시 Windows 가 포함되어 있는 환경이라면 가급적 long 을 사용하지 말고 (int , long long) 또는 C99 표준에서 추가된 stdint.h 에 포함된 (int32_t, int64_t) 자료형을 사용하는걸 권장v

int main()
{

    /*int a;
    int b;
    
    cin >> a >> b;

    int c[a];

    for(int i=0;i<b;i++)
    {
        
    }
    */

    /*
    char name[10];
    name[0] = 'J';
    name[1] = 'o';
    name[2] = 'h'; 
    name[3] = 'n';
    name[4] = '\0'; 맥이랑 윈도우랑 다르게 출력됨 ㅁㄴㅇㄹ

    

    cout << name << endl;
    */

   /*
   char name[] = "John";
   char name2[] = {'J','o','h','n'};

   cout << name[0] << endl;
   cout << name2 << endl; 

   cout << sizeof(name2) << endl; // 역시 윈도우랑 맥이랑 다르게 나옴. 자료형 문제인건가?/
   cout << sizeof(name) << endl; // end of string이 할당되어야 쓰레기값이 나오지 않음. 총 5바이트 할당됨.
    */
   /*
    char str[20] = "Hello";
    char str2[] = "World";

    cout << strlen(str) << endl;
    cout << sizeof(str) << endl; 
    strncat_s(str, str2, 4); // 이거 윈도우에서는 돌아가는데 맥에서는 안돌아감 왜 그런거지
    cout << str << endl;
        if (strcmp(str, "HelloWorld") == 0)
        cout << "OK" << endl; ㄴ
        else
        cout << "Fail" << endl; char str01[] = "10";
    char str02[] = "20";

    cout << atoi(str01) * atof(str02) << endl;
    */

    return 0;
}