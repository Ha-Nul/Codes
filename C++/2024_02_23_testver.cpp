#include<iostream>
#include <fstream>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <vector>
#include <cmath>
#include <iomanip>
#include <string>

using namespace std;
using namespace Eigen;

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

vector<double> grid = linspace(0,1,20);
double dtau = grid[1]-grid[0];

vector<double> testf()
{
    vector<double> value(grid.size(),0);
    
    for (int i=0;i<grid.size();i++)
    {
        value[i] = 0.6 * grid[i] + 2;
    } 
    
    return value;
}

vector<double> testg()
{
    return grid;
}


double weight(const int& n, const int& m){
	if(n==0) return 0.0;
	else if(m==0 or m==n) return 0.5;
	else return 1.0;
};

vector<double> updateGt_volterra(){

    vector<double> Gt(grid.size(),0);
    vector<double> dGt(grid.size(),0);
    vector<double> St = testf();

    Gt[0] = 1;
    dGt[0] = 1;
    double a = 0;
    vector<double> c(grid.size(),0);
    vector<double> d(grid.size(),0);

    // Green's function update
	for(int n=0; n<grid.size()-1; n++){
        double b = 0;
		//if(n==0) cout<<Gt[n]<<endl;
		Gt[n+1] = Gt[n];
		//if(n==0) cout<<Gt[n+1]<<endl;

		//reconstruct Roundpropagator
		Gt[n+1] -= 0.5*dtau*( 2.0*Gt[n] + dtau*dGt[n] );
		//if(n==0) cout<<Gt[n+1]<<endl;
		for(int p=0; p<n+1; p++){
            a = 0.5*pow(dtau,2)
				*( weight(n,p)*St[n-p] + weight(n+1,p)*St[n+1-p])*Gt[p];
			Gt[n+1] += a;
            b += a;
            c[n+1] += weight(n,p)*St[n-p] * Gt[p];
            d[n+1] += weight(n+1,p)*St[n+1-p]*Gt[p];
            //cout << "this is forloop : " << "\t" << n <<"\t" << a << endl; //weight(n,p) << "\t" << weight(n+1,p) << endl;
		}
         //weight(n,p) << "\t" << weight(n+1,p) << endl;
		
        //cout << "this is last term : " << Gt[n] + dtau*dGt[n] << endl;
		Gt[n+1] += 0.25*pow(dtau,2)*St[0]*( Gt[n] + dtau*dGt[n] );
        d[n+1] += 0.5 * St[0]*( Gt[n] + dtau*dGt[n] );

        cout << "this is forloop : " << n+1 <<"\t" << b << "\t" << dtau * c[n+1] << "\t"<< dtau * d[n+1] << endl;

		//if(n==0) cout<<Gt[n+1]<<endl;

        //round Green's function update
		dGt[n+1] = -Gt[n+1];
		for(int p=0; p<n+2; p++){
			dGt[n+1] += dtau*weight(n+1,p)*St[n+1-p]*Gt[p];
            //cout << "this is Gtp : " << "\t" << Gt[p] << endl;
        }
        //cout << "this is dGt : " <<  dGt[n+1] << endl;
	}

    return Gt;
};

double round_propagator_ite(const vector<double> &sigma, const vector<double> &ite, int n, int boolean)
{   
    double sigsum = 0;

    if (n == 1)
    {
        sigsum = 0.5 * dtau * (sigma[1]*ite[0] + sigma[0]*ite[1]);
    }
    else if (n > 1){
        for (int i = 0 ; i < n ; i++)
        {
            sigsum += 0.5 * dtau * (sigma[n-(i)] * ite[i] + sigma[n-(i+1)] * ite[i+1]);

            if (i+1 == n)
            {
                //cout << "this is last ite : " << ite[i+1] << endl;
                break;
            }

        }
    }

    double Bucket = 0;
    if (boolean == 0)
    {
        Bucket = -ite[n] + sigsum;
    }    
    else if (boolean == 1)
    {
        Bucket = sigsum;
    }

    return Bucket;
}



vector<double> Propagator()
{
    vector<double> P_arr(grid.size(),0);
    vector<double> sig = testf();
    vector<double> sigP(grid.size(),0);

    double sig_form;
    double sig_late;
    double a = 0;

    P_arr[0] = 1;
    sigP[0] = 1;

    for (int n=1; n<grid.size();n++)
    {
        if (n == 1)
        {
            P_arr[1] = P_arr[0];
            a = 0.5 * dtau * ( 0.5 * dtau * (sig[1]*P_arr[0] + sig[0]*(P_arr[0] + dtau * 1)) );
            //cout << "a value" << a << endl;
            P_arr[1] = P_arr[0] - 0.5 * dtau * (2 * P_arr[0] + dtau * 1) + a;
            sigP[1] = P_arr[1];
    
            //cout << "this is P_arr : " << "\t" << n << "\t" << P_arr[n] << endl;
        }
        else if ( n > 1)
        {
            sig_form = round_propagator_ite(sig,P_arr,n-1,0);
            //cout << "this is sigform : " << sig_form << endl;
            sigP[n] = P_arr[n-1] + dtau * sig_form;
            //cout << "--------------dummy---------------" << endl;
            sig_late = 0.5 * dtau * (round_propagator_ite(sig,P_arr,n-1,1) + round_propagator_ite(sig,sigP,n,1));
            cout << "this is sig_late : " << n << "\t" << sig_late << "\t" << round_propagator_ite(sig,P_arr,n-1,1) << "\t" << round_propagator_ite(sig,sigP,n,1) << endl;

            P_arr[n] = P_arr[n-1] - 0.5 * dtau * (2 * P_arr[n-1] + dtau * sig_form) + sig_late;
            //cout << "this is P_arr : " << "\t" << n << "\t" << P_arr[n] << endl;

            //cout << "this is sig_late : " << "\t" << n << "\t" << sig_late << endl;
            //P_arr[n] += 0.25 * pow(dtau,2) * sig[0] * (P_arr[n-1] + dtau * sig_form);
        }
    }

    return P_arr;
}


/*
double loop(int n)
{   
    if (n == 1)
    {
        cout << "index :" << n << "," << n-1 << endl;
    }
    else if (n > 1){
        for (int i = 0 ; i < n ; i++)
        {
            cout << "index :" << n << "," << i << endl;

            if (i+1 == n)
            {
                break;
            }

        }
    }
    return 0;

};



double weightint(int n)
{
    double weigh_int = 0;

    for (int p=0; p<n+2; p++)
    {
        weigh_int += inter*weight(n+1,p)*testf()[n+1-p]*testg()[p];
        cout << "index :" << n+1 << "," << p << endl;
    }

    return weigh_int;
}

double trapint(int n)
{
    double bb = 0;

    if (n == 1)
    {
        bb += testf()[0]*testg()[1] + testf()[1]*testg()[0];
        cout << "index :" << n << "," << n-1 << endl;
    }
    else if (n > 1){
        for (int i = 0 ; i < n ; i++)
        {
            cout << "index :" << n << "," << i << "," << i+1 << endl;
            bb += 0.5 * inter * (testf()[n-(i)] * testg()[i] + testf()[n-(i+1)] * testg()[i+1]);
            
            if (i+1 == n)
            {
                //break;
            }

        }
    }
    //return bb;
}

double indexcount1(int k)
{
    for (int n = 0; n<k ;n++)
    {
        for (int p=0; p<n+2;p++)
        {
            cout << "index :" << n+1 << "," << p << endl;
        }
    }
};

double indexcount2(int k)
{
    for (int i = 0 ; i<k ; i++)
    {
        trapint(i);
    }
};

*/

int main()
{
    /*
    for (int i = 0 ; i<20;  i++)
    {
        for (int j = 0; j < i+1 ; j++ )
        {
            cout << "index :" << i << "," << j << "\t" << "value : " <<weight(i,j) << endl;
        }
    }

    cout << "-------------------index compare------------------" << endl;

    loop(19);
    */


    
    //cout << weightint(8) << endl;
    //cout << trapint(9) << endl;
    
   
    //indexcount1(9);
    //indexcount2(10);

    vector<double> a = updateGt_volterra();
    vector<double> b = Propagator();

    for (int i = 0; i < grid.size(); i++)
    {
        cout << a[i] << endl;
    }

    cout << "------------------------------------" << endl;

    for (int i = 0; i < grid.size(); i++)
    {
        cout << b[i] << endl;
    }
    
    return 0;
}



