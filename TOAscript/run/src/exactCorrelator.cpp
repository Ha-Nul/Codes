#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <Eigen/Dense>
#include "json.hpp"
//#include "dyson.hpp"
using namespace std;

int main(int argc, char* argv[]){
	cout<<"***********************"<<endl;
	cout<<"** Exact calculation **"<<endl;
	cout<<"***********************"<<endl;
	nlohmann::json input;

	ifstream para(argv[1]);
	para >> input;
	para.close();

	int Ntau = input["Ntau"];
    int nB_level = input["nB_level"];
	double beta = input["beta"];
    double gCoupling = input["gCoupling"];
    double w_boson = input["w_boson"];

	cout<<"** parameters **"<<endl;
	cout<<"   Ntau = "<<Ntau<<endl;
	cout<<"   nB_level = "<<nB_level<<endl;
	cout<<"   beta = "<<beta<<endl;
	cout<<"   gCoupling = "<<gCoupling<<endl;
	cout<<"****************"<<endl;

    const int NHilbertLoc = 3;

    Eigen::MatrixXcd Hloc(NHilbertLoc,NHilbertLoc);

    Hloc << -3.7841100000e-01 , 0.0000000000e+00 ,0.0000000000e+00 
        , 0.0000000000e+00 , 9.1806000000e-01 , 0.0000000000e+00
        , 0.0000000000e+00 , 0.0000000000e+00 , 1.2940100000e+00;

    cout<<Hloc<<endl;

    Eigen::MatrixXcd Ntilde(NHilbertLoc,NHilbertLoc);
    Ntilde << 0 , -0.482331 , 0
        , 0.482331 , 0 , -0.906537
        , 0 , 0.906537 , 0;

    cout<<Ntilde<<endl;

	int dim = NHilbertLoc * nB_level;

    Eigen::MatrixXcd H = Eigen::MatrixXd::Zero(dim,dim);
    for(int i=0; i<nB_level; i++) for(int j=0; j<nB_level; j++){
        if(i==j){
            H.block(i*NHilbertLoc,j*NHilbertLoc,NHilbertLoc,NHilbertLoc) = Hloc + i*w_boson*Eigen::MatrixXd::Identity(NHilbertLoc,NHilbertLoc);
        }
        if(j==i+1 && j<nB_level){
            H.block(i*NHilbertLoc,j*NHilbertLoc,NHilbertLoc,NHilbertLoc) = sqrt(i+1)*gCoupling*complex<double>(0.0,-1.0)*Ntilde;
        }
        if(j==i-1 && i<nB_level){
            H.block(i*NHilbertLoc,j*NHilbertLoc,NHilbertLoc,NHilbertLoc) = sqrt(j+1)*gCoupling*complex<double>(0.0,1.0)*Ntilde;
        }
    }

    cout<<H<<endl;

	/*Eigen::MatrixXd aa = Eigen::MatrixXd::Zero(dim,dim);
	Eigen::MatrixXd ac = Eigen::MatrixXd::Zero(dim,dim);
	Eigen::MatrixXd n_boson = Eigen::MatrixXd::Zero(dim,dim);
	Eigen::MatrixXd Sz = Eigen::MatrixXd::Zero(dim,dim);
	Eigen::MatrixXd Sx = Eigen::MatrixXd::Zero(dim,dim);
	for(int sz=0; sz<NHilbertLoc; sz++){
		for(int n=0; n<nB_level; n++){
			int i=n*NHilbertLoc+sz;
			for(int sz1=0; sz1<NHilbertLoc; sz1++){
				for(int n1=0; n1<nB_level; n1++){
					int i1=n1*NHilbertLoc+sz1;
					if(n1==n+1 && sz1==sz) ac(i1,i)=sqrt(n+1);
					if(n1==n-1 && sz1==sz) aa(i1,i)=sqrt(n);
					if(n1==n && sz1==sz) n_boson(i1,i)=n;
					if(n1==n && sz1==sz) Sz(i1,i)=(sz1==1 ? 1.0 : -1.0);
					if(n1==n && sz1!=sz) Sx(i1,i)=1.0;
				}
			}
		}
	}
	Eigen::MatrixXd H = epsilon*Sz+B_field*Sx+g_coupling*(aa+ac)*Sz+w_boson*n_boson;
	//cout<<"* Hamiltonian *"<<endl<<H<<endl;

	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(H);

	Eigen::MatrixXd evec = solver.eigenvectors();
	Eigen::VectorXd eval = solver.eigenvalues();

	//cout<<"eigenvalue: "<<endl<<eval.transpose()<<endl;

	//hdiag.setZero();
	//for(int i=0;i<dim;i++) hdiag(i,i)=eval(i);
	double Z = 0.0;
	double E0 = eval(0);
	for(int i=0;i<dim;i++) Z += exp(-beta*(eval(i)-E0));

	double tau_tmp;
	double dtau = beta/Ntau;

	Eigen::MatrixXd diag1;
	Eigen::MatrixXd diag2;

	stringstream ss;
	ss<<"StS0_b"<<beta<<"ep"<<epsilon<<"B"<<B_field<<"g"<<g_coupling<<"w"<<w_boson<<"_nB"<<nB_level<<".exact";
	ofstream ofstr(ss.str().c_str());
	ofstr<<scientific;
	vector<double> SzSz(Ntau+1);
	for(int m=0; m<Ntau+1; m++){
	    tau_tmp = m*dtau;
	    diag1 = Eigen::MatrixXd::Zero(dim,dim);
	    for(int i=0; i<dim; i++) diag1(i,i) = exp(-tau_tmp*(eval(i)-E0));
	    diag2 = Eigen::MatrixXd::Zero(dim,dim);
	    for(int i=0; i<dim; i++) diag2(i,i) = exp(-(beta-tau_tmp)*(eval(i)-E0));
	    double SzSz_tmp = (evec*diag2*evec.adjoint()*Sz*evec*diag1*evec.adjoint()*Sz).trace();
	    SzSz[m] = SzSz_tmp/Z;

	    ofstr<<left<<setw(30)<<tau_tmp
		    <<right<<setw(30)<<SzSz[m]<<endl;
	}
	ofstr.close();*/

	return 0;
};

