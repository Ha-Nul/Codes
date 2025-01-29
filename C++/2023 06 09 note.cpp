#include <iostream>
#include <complex>
#include <vector>

using namespace std;

int main() {

    vector<double> arr(2,0);

    cout << " ** First output : "; 
    cin >> arr[0];

    cout << " ** Second output : ";
    cin >> arr[1];

    std::complex<double> z1(arr[0], arr[1]);  // 3 + 4i
    std::complex<double> z2(arr[1], arr[0]);  // 1 + 2i

    // 덧셈
    std::complex<double> sum = z1 + z2;
    cout << "z1 + z2 = " << sum << endl;

    // 곱셈
    std::complex<double> product = z1 * z2;
    cout << "z1 * z2 = " << product << endl;

    // 나눗셈
    std::complex<double> quotient = z1 / z2;
    cout << "z1 / z2 = " << quotient << endl;

    // 절대값 (magnitude)
    double mag = abs(z1);
    cout << "|z1| = " << mag << endl;

    return 0;
}
