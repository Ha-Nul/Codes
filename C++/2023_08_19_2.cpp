#include<iostream>
using namespace std;

union Year
{
    int ad;
    int dangi;
};

int main()
{
    Year myyear;

    cout << "서기를 입력하십시오. \n";
    cin >> myyear.ad;

    cout << "서기" << myyear.ad << "년입니다.\n";
    cout << "단기 또한" << myyear.dangi << "년 입니다.\n";

    cout << "단기를 입력하십시오";
    cin >> myyear.dangi;

    cout << "단기" << myyear.dangi << "년입니다.\n";
    cout << "서기 또한" << myyear.ad << "년입니다.\n";

    return 0;
}