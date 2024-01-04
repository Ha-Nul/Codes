    #include <iostream>
    #include <vector>
    #include <Eigen/Dense>
    #include <firstheader.hpp>

    using namespace std;
    using namespace Eigen;
    
    vector<MatrixXd> convolve(const vector<MatrixXd>& Signal,
              const vector<MatrixXd>& Kernel, int n, int i)
    {
        size_t SignalLen = i;
        size_t KernelLen = Kernel.size();
        size_t ResultLen = SignalLen + KernelLen - 1;

        vector<MatrixXd> Result(ResultLen,MatrixXd::Zero(n,n));

        for (size_t n = 0; n < ResultLen; ++n)
        {
            size_t kmin = (n >= KernelLen - 1) ? n - (KernelLen - 1) : 0;
            size_t kmax = (n < SignalLen - 1) ? n : SignalLen - 1;

            for (size_t k = kmin; k <= kmax; ++k)
            {
                Result[n] += Signal[k] * Kernel[n - k];
            }
        }

        return Result;
    }

    MatrixXd test(int n)
    {
        cout << "This is sigma index test with " << n << endl;
        MatrixXd bucket = MatrixXd::Zero(3,3);
        vector<MatrixXd> one(n+1,MatrixXd::Ones(3,3));

        for (int i=0; i<=n; i++)
        {

            bucket += 0.5 * (one[i] + one[i+1]);
            cout << "sigma bucket " << i << "th" << endl;
            cout << bucket << endl;
            if (i+1 == n)
            {
                break;
            }
        }
    
        return bucket;
    }

    vector<MatrixXd> Propagator(int k)
    {
        vector<MatrixXd> proparray(k,MatrixXd::Zero(3,3));
        MatrixXd Iden = MatrixXd::Identity(3,3);

        proparray[0] = Iden;

        MatrixXd bucket1 = MatrixXd::Zero(3,3);
        MatrixXd bucket2 = MatrixXd::Zero(3,3);
        MatrixXd bucket = MatrixXd::Zero(3,3);
        
        for(int i = 1; i < k; i++)
        { 
            cout << "this is prop test" << endl;
            if (i>1)
            {
                bucket1 = test(i);
                proparray[i] = proparray[i-1] + bucket1;
                cout << "bucket change" << endl;
                bucket2 = test(i+1);
                proparray[i] = proparray[i-1] + 0.5 * (bucket1 + bucket2);
            }
                proparray[i] = proparray[i-1] + bucket1;

        }

        return proparray;
    }
/*
MatrixXd round_propagater_ite(const MatrixXd &loc, const vector<MatrixXd> &sigma, const vector<MatrixXd> &ite, int n)
{
    //assume that array size of prop and sigma is same as k

    vector<MatrixXd> convol_array = convolve(sigma,ite,3);

    MatrixXd sigsum = MatrixXd::Zero(3,3);

    for (int j=0;j<convol_array.size();j++)
    {
        sigsum = sigsum + convol_array[j];
    }

    MatrixXd Bucket;
    Bucket = -loc * ite[1] + sigsum;
    
    //cout << -loc * ite << endl;
    return Bucket;
}
*/

int main()
{
    Testing test;

    vector<double> test_t = test.grid;

    for (int i = 0; i<100 ; i++)
    {
        cout << test_t[i] << endl;
    }
}


