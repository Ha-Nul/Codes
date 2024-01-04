#include<iostream>
#include<vector>
#include<firstheader.hpp>


using namespace std;

template<typename T>
std::vector<T>
conv(std::vector<T> const &f, std::vector<T> const &g) {
  int const nf = f.size();
  int const ng = g.size();
  int const n  = nf + ng - 1;
  std::vector<T> out(n, T());
  for(auto i(0); i < n; ++i) {
    int const jmn = (i >= ng - 1)? i - (ng - 1) : 0;
    int const jmx = (i <  nf - 1)? i            : nf - 1;
    cout << "j min : " << jmn << endl;
    cout << "j max : " << jmx << endl;
    for(auto j(jmn); j <= jmx; ++j) {
      out[i] += (f[j] * g[i - j]);
      cout << "convolu" << out[i] << endl;
    }
  }
  return out; 
}

vector<double> linspace(const double &min,const double &max, int n)
        {
            vector<double> result;
            // vector iterator
            int iterator = 0;

            for (int i = 0; i <= n-2; i++)	
            {
                double temp = min + i*(max-min)/(floor((double)n) - 1);
                result.insert(result.begin() + iterator, temp);
                iterator += 1;
            }

            //iterator += 1;

            result.insert(result.begin() + iterator, max);
            return result;
        }


int main()
{
    int j = 0;
    vector<double> a = linspace(0,1,100);
    vector<double> b = linspace(0,1,100);

    auto v1 = conv(a,b);

    for(auto i: v1){
        cout << j++ << "\t" << i << " " << endl;
    }
}