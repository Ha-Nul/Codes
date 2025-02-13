//#include <mpi.h>
#include <cmath>
#include <iomanip>
#include <string>
#include <sstream>
#include "dyson.hpp"
//#include "parameters.hpp"
//#include "diagram.hpp"
#include "green.hpp"
using namespace std;
//using namespace parameters;

dyson::dyson(const int& Nflavor_i, const double& beta_i, const int& Ntau_i)
	: Nflavor(Nflavor_i), Ntau(Ntau_i), beta(beta_i), dtau(beta_i/Ntau_i){
	ZeroM = Eigen::MatrixXd::Zero(Nflavor,Nflavor);
	ZeroV = Eigen::MatrixXd::Zero(2*Nflavor,1);
	IdentityM = Eigen::MatrixXd::Identity(Nflavor,Nflavor);

	/*Eigen::Matrix2d Sz_tmp;
	Sz_tmp << 1.0, 0.0,
	       0.0, -1.0;
	Sz = Sz_tmp;

	Eigen::Matrix2d Hatom_tmp;
	Hatom_tmp << ep_i, B_i,
	  B_i, -ep_i;
	Hatom = Hatom_tmp;

	Zero2x2 = Eigen::MatrixXd::Zero(2,2);
	Zero4x1 = Eigen::MatrixXd::Zero(4,1);*/

	/*
	Tt_amp_up.resize(Ntau+1);
	Tt_amp_dn.resize(Ntau+1);
	Tt_amp.resize(Ntau+1);
	for(int n=0; n<Ntau+1; n++)
		Tt_amp[n].resize(n+1);
	*/
};
void dyson::readWt(ifstream& ifstr){
	string line;
	stringstream ss;
	string sharp, index;

	int Nt_tmp = 0;
	while( getline(ifstr,line) ){
		Nt_tmp++;
	}
	double dt_tmp = beta/(Nt_tmp-1);

	ifstr.clear();
	ifstr.seekg(0);

	double tau;
	vector<double> Wt_tmp(Nt_tmp);
	for(int i=0; i<Nt_tmp; i++){
		ifstr >> tau >> Wt_tmp[i];
	}

	int it_tmp;
	Wt.resize(Ntau+1);
	Wt[0] = Wt_tmp[0];
	Wt[Ntau] = Wt_tmp[Nt_tmp-1];
	for(int i=1; i<Ntau; i++){
		tau = i*dtau;
		it_tmp = static_cast<int>(tau/dt_tmp);
		Wt[i] = (
				((it_tmp+1)*dt_tmp-tau)*Wt_tmp[it_tmp]
				+ (tau - it_tmp*dt_tmp)*Wt_tmp[it_tmp+1]
			) / dt_tmp;
	}
};
void dyson::readHloc(ifstream& ifstr){
	Hloc.resize(Nflavor,Nflavor);

	for(int i=0; i<Nflavor; i++) for(int j=0; j<Nflavor; j++){
		ifstr >> Hloc(i,j);
	}
};
void dyson::readVOp(ifstream& ifstr){
	VOp.resize(Nflavor,Nflavor);

	for(int i=0; i<Nflavor; i++) for(int j=0; j<Nflavor; j++){
		ifstr >> VOp(i,j);
	}
};

//int dyson::gen4index(const int& a, const int& b){ return a*Nflavor+b; };
//vector<int> dyson::gen22index(const int& ab){
//	vector<int> twoIndices(2);
//	twoIndices[0] = ab/Nflavor;
//	twoIndices[1] = ab%Nflavor;
//	return twoIndices;
//};
//Eigen::MatrixXd dyson::transformTo2x2(const Eigen::MatrixXd& V4x1){
//	Eigen::Matrix2d	Mtmp;
//	Mtmp(0,0) = V4x1(0);
//	Mtmp(1,0) = V4x1(1);
//	Mtmp(0,1) = V4x1(2);
//	Mtmp(1,1) = V4x1(3);
//
//	return Mtmp;
//};
//Eigen::MatrixXd dyson::transformTo4x1(const Eigen::MatrixXd& M2x2){
//	Eigen::Vector4d	Vtmp;
//
//	Vtmp(0) = M2x2(0,0);
//	Vtmp(1) = M2x2(1,0);
//	Vtmp(2) = M2x2(0,1);
//	Vtmp(3) = M2x2(1,1);
//
//	return Vtmp;
//};
//Eigen::Matrix4d dyson::tensorProd(const Eigen::MatrixXd& M2x2_1, const Eigen::MatrixXd& M2x2_2){
//	Eigen::Matrix4d	Mtmp;
//
//	Mtmp.topLeftCorner(2,2) = M2x2_1(0,0)*M2x2_2;
//	Mtmp.topRightCorner(2,2) = M2x2_1(1,0)*M2x2_2;
//	Mtmp.bottomLeftCorner(2,2) = M2x2_1(0,1)*M2x2_2;
//	Mtmp.bottomRightCorner(2,2) = M2x2_1(1,1)*M2x2_2;
//
//	return Mtmp;
//};
//void dyson::initWt(const double& gCouple, const double& wBoson){
//	double tau;
//
//	mWt.resize(Ntau+1);
//	for(int i=0; i<Ntau+1; i++){
//		tau = i*dtau;
//		mWt[i] = pow(gCouple,2)*cosh( (tau-0.5*beta)*wBoson ) / sinh( 0.5*beta*wBoson );
//		//mWt[i] = -pow(gCouple,2)*cosh( (tau-0.5*beta)*wBoson ) / sinh( 0.5*beta*wBoson );
//	}
//};
//int dyson::getSimpsonFactor(const int& ix, const int& Nx){
//	if(ix==0 || ix==Nx) return 1;
//	else if(ix%2==1) return 4;
//	else return 2;
//};
//void dyson::initWtSingleModeWaveguide(const double& g0, const double& Omega11, const double& omega_c, const int& Nkx){
//	double preFactor = g0*pow(Omega11,4);
//	double kc = sqrt(pow(omega_c/Omega11,2)-1);
//	double betaTilde = beta*Omega11;
//	double dtauTilde = (this->dtau)*Omega11;
//	double NkxRegular;
//	if(Nkx%2==1) NkxRegular = Nkx+1;
//	else NkxRegular = Nkx;
//	double dkx = kc/NkxRegular;
//
//	double kx;
//	double simpsonFactor, intFactor, omega11;
//	mWt.resize(Ntau+1);
//	for(int it=0; it<Ntau+1; it++){
//		mWt[it] = 0.0;
//	}
//	for(int ikx=0; ikx<NkxRegular+1; ikx++){
//		kx = ikx*dkx;
//		omega11 = sqrt(kx*kx+1.0);
//		//cout<<"for ikx="<<ikx<<", omega11 = "<<omega11<<endl;
//		intFactor = 1.0/omega11/sinh(0.5*betaTilde*omega11);
//		simpsonFactor = getSimpsonFactor(ikx,Nkx);
//		for(int it=0; it<Ntau+1; it++){
//			mWt[it] += simpsonFactor*intFactor*cosh(omega11*(it*dtauTilde - 0.5*betaTilde));
//		}
//	}
//	for(int it=0; it<Ntau+1; it++){
//		mWt[it] *= preFactor*2.0*dkx/3.0;
//	}
//};
void dyson::initGt(){
	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver(Hloc);

	Eigen::VectorXd eval = eigensolver.eigenvalues();
	Eigen::MatrixXd evec = eigensolver.eigenvectors();

	lambda = 0.0;
	Eigen::MatrixXd expFactor = ZeroM;
	for(size_t i=0; i<eval.size(); i++) expFactor(i,i) = exp(-dtau*(eval(i)-lambda));
	Eigen::MatrixXd propagator = IdentityM;

	//G0t.resize(Ntau+1);
	Gt.resize(Ntau+1);
	for(int ti=0; ti<Ntau+1; ti++){
		//G0t[ti] = evec*propagator*evec.adjoint();
		//Gt[ti] = G0t[ti];
		Gt[ti] = evec*propagator*evec.adjoint();
		propagator = propagator*expFactor;
	}

	/*double tau;
	double G0t_tmp;

	lambda = 0.0;
	G0t.resize(Ntau+1);
	Gt.resize(Ntau+1);
	for(int i=0; i<Ntau+1; i++){
		tau = i*dtau;
		G0t_tmp = exp(-tau*(epsilon - lambda));
		G0t[i] = G0t_tmp;
		Gt[i] = G0t_tmp;
	}*/
};
void dyson::readGt(ifstream& ifstr){
	string line;
	stringstream ss;
	string sharp, index;

	getline(ifstr,line);
	ss<<line;
	ss>>sharp>>index>>lambda;
	cout<<"lambda = "<<lambda<<endl;

	int Nt_tmp = 0;
	while( getline(ifstr,line) ){
		Nt_tmp++;
	}
	double dt_tmp = beta/(Nt_tmp-1);

	ifstr.clear();
	ifstr.seekg(0);

	double tau;
	vector<Eigen::MatrixXd> Gt_tmp(Nt_tmp);
	getline(ifstr,line);
	for(int i=0; i<Nt_tmp; i++){
		ifstr>>tau;
        Gt_tmp[i] = ZeroM;
        for(int iflavor=0; iflavor<Nflavor; iflavor++)
        for(int jflavor=0; jflavor<Nflavor; jflavor++)
            ifstr>>Gt_tmp[i](iflavor,jflavor);
	}

	int it_tmp;
	Gt.resize(Ntau+1);
	Gt[0] = Gt_tmp[0];
	Gt[Ntau] = Gt_tmp[Nt_tmp-1];
	for(int i=1; i<Ntau; i++){
		tau = i*dtau;
		it_tmp = static_cast<int>(tau/dt_tmp);
		Gt[i] = (
				((it_tmp+1)*dt_tmp-tau)*Gt_tmp[it_tmp]
				+ (tau - it_tmp*dt_tmp)*Gt_tmp[it_tmp+1]
			) / dt_tmp;
		//Gt[i] = G0t_tmp;
	}

	//string sharp, index;
	//ifstr>>sharp>>index>>lambda;
	//double tau;
	//Gt.resize(Ntau+1);
	//for(int i=0; i<Ntau+1; i++){
	//	ifstr>>tau>>Gt[i];
	//	//Gt[i] = G0t_tmp;
	//}
};
////void dyson::setGtGlobal(){
////	parameters::Ntime = Ntau;
////	dt = dtau;
////	//cout<<"dt = "<<dt<<endl<<flush;
////
////	Eigen::MatrixXd Gt_tmp(Ntime+1,1);
////	for(int i=0; i<Ntime+1; i++)
////		Gt_tmp(i,0) = Gt[i];
////
////	parameters::Gt[0].setOtx(Gt_tmp);
////};
////
/////*
////void dyson::initGtexact(const double& epsilon, const double& gCouple, const double wBoson, const int& nBcut){
////	int Nh = nBcut;
////
////	Eigen::MatrixXd H = Eigen::MatrixXd::Zero(Nh,Nh);
////	H(0,0) = epsilon;
////	for(int i=1; i<Nh; i++){
////		H(i,i) = epsilon + i*wBoson;
////		H(i-1,i) = gCouple;
////		H(i,i-1) = gCouple;
////	}
////
////	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver(H);
////
////	Eigen::MatrixXd Evec = eigensolver.eigenvectors();
////	Eigen::VectorXd Eval = eigensolver.eigenvalues();
////
////	Eigen::VectorXd exp_mdtH(Nh);
////	Eigen::VectorXd exp_mtH(Nh);
////
////	double Z = 0.0;
////	for(int i=0; i<Eval.size(); i++){
////		Z += exp(-beta*Eval(i));
////		exp_mdtH(i) = exp( -dtau*Eval(i) );
////		exp_mtH(i) = 1.0;
////	}
////
////	double lambda = -log(Z)/beta;
////	double exp_dtlambda = exp( dtau*lambda );
////	double exp_tlambda = 1.0;
////
////	double tau;
////
////	for(int i=0; i<Ntau+1; i++){
////		tau = i*dtau;
////		Gt[i] = exp_mtH.sum()*exp_dtlambda/static_cast<double>(Nh);
////		exp_mtH = exp_mtH.cwiseProduct(exp_mdtH).eval();
////	}
////
////};
////*/
//
//void dyson::initQtReducible(){
//	Qt_reducible.resize(Ntau+1);
//	for(int it=0; it<Ntau+1; it++){
//		Qt_reducible[it].resize(it+1);
//		for(int it1=0; it1<it+1; it1++){
//			Qt_reducible[it][it1].resize(it1+1);
//			for(int it2=0; it2<it1+1; it2++){
//				Qt_reducible[it][it1][it2] = Eigen::MatrixXd::Zero(4,4);
//			}
//		}
//	}
//};
//void dyson::initQt2ndOrder(){
//	Qt.resize(Ntau+1);
//	for(int it=0; it<Ntau+1; it++){
//		Qt[it].resize(it+1);
//		for(int it2=0; it2<it+1; it2++){
//			Qt[it][it2].resize(it2+1);
//			for(int it1=0; it1<it2+1; it1++){
//				Qt[it][it2][it1] = mWt[it2]*mWt[it-it1]*tensorProd(Sz*Gt[it1]*Sz,Sz*Gt[it-it2]*Sz);
//			}
//		}
//	}
//
//	/*Eigen::MatrixXd Qt_n_tmp;
//
//	for(int n=0; n<Ntau+1; n++){
//		Qt_n_tmp.resize(n+1,n+1);
//		for(int i=0; i<n+1; i++) for(int j=i; j<n+1; j++){
//			Qt_n_tmp(j,i) = Gt[i]*Gt[n-j]*mWt[j]*mWt[n-i];
//		}
//		Qt.setQt(n,Qt_n_tmp);
//	}*/
//};
//void dyson::initQt3rdOrder(){
//	Qt.resize(Ntau+1);
//
//	Eigen::Matrix4d Qt_tmp;
//	for(int n=0; n<Ntau+1; n++){
//		Qt[n].resize(n+1);
//		for(int m2=0; m2<n+1; m2++){
//			Qt[n][m2].resize(m2+1);
//			for(int m1=0; m1<m2+1; m1++){
//				// Q2 //
//				Qt[n][m2][m1] = mWt[m2]*mWt[n-m1]*tensorProd(Sz*Gt[m1]*Sz,Sz*Gt[n-m2]*Sz);
//
//				// Q3 //
//				Qt_tmp = Eigen::MatrixXd::Zero(4,4);
//				for(int p=0; p<m1+1; p++) for(int q=m2; q<n+1; q++){
//					Qt_tmp += weight(m1,p)*weight(n-m2,q-m2)
//						* (mWt[m2]*mWt[n-p]*mWt[q-m1] + mWt[q]*mWt[m2-p]*mWt[n-m1] + mWt[m2]*mWt[q-p]*mWt[n-m1])
//						* tensorProd(Sz*Gt[m1-p]*Sz*Gt[p]*Sz,Sz*Gt[n-q]*Sz*Gt[q-m2]*Sz);
//				}
//
//				for(int p=m2; p<n+1; p++) for(int q=p; q<n+1; q++){
//					Qt_tmp += weight(n-m2,p-m2)*weight(n-p,q-p)
//						* (mWt[q]*mWt[p-m1]*mWt[n-m2] + mWt[p]*mWt[q-m1]*mWt[n-m2] + mWt[p]*mWt[q-m2]*mWt[n-m1] + mWt[m2]*mWt[q-m1]*mWt[n-p])
//						* tensorProd(Sz*Gt[m1]*Sz,Sz*Gt[n-q]*Sz*Gt[q-p]*Sz*Gt[p-m2]*Sz);
//				}
//
//				for(int p=0; p<m1+1; p++) for(int q=p; q<m1+1; q++){
//					Qt_tmp += weight(m1,p)*weight(m1-p,q-p)
//						* (mWt[m1]*mWt[n-p]*mWt[m2-q] + mWt[m1]*mWt[m2-p]*mWt[n-q] + mWt[m2]*mWt[m1-p]*mWt[n-q] + mWt[q]*mWt[m2-p]*mWt[n-m1])
//						* tensorProd(Sz*Gt[m1-q]*Sz*Gt[q-p]*Sz*Gt[p]*Sz,Sz*Gt[n-m2]*Sz);
//				}
//				Qt[n][m2][m1] += pow(dtau,2)*Qt_tmp;
//			}
//		}
//	}
//};
//
////void dyson::diagMC_Qt(const int& Nordermax){
////	//cout<<"** diagMC_Qt **"<<endl<<flush;
////	const int root = 0;
////	const int rank = MPI::COMM_WORLD.Get_rank();
////
////	//Qt.resize(Ntau+1);
////
////	//cout<<"Initialize"<<endl<<flush;
////	diagram DiagMC;
////	DiagMC.initDiagram();
////
////	//cout<<"Sampling"<<endl<<flush;
////	vertex_4pt Qt_tmp, Qt_acc, Qt_sqacc;
////	Qt_acc.initZero();
////	Qt_sqacc.initZero();
////	for(int i=0; i<Nerr; i++){
////		DiagMC.MCsampling();
////		DiagMC.partialSum();
////
////		Qt_tmp = DiagMC.getQt();
////		Qt_acc += Qt_tmp;
////		Qt_sqacc += Qt_tmp*Qt_tmp;
////		//for(int n=0; n<Ntau+1; n++){
////		//	Qt[n] = DiagMC.getQt(n);
////		//}
////	}
////
////	double *source, *result, *output1, *output2;
////
////	int Ndata = 0;
////	for(int it=0; it<Ntau+1; it++)
////		Ndata += (it+1)*(it+1);
////
////	source = new double[2*Ndata];
////	result = new double[2*Ndata];
////	output1 = new double[Ndata];
////	output2 = new double[Ndata];
////
////	//cout<<"MPI Pack"<<endl<<flush;
////	int position = 0;
////	for(int it=0; it<Ntau+1; it++){
////		MPI::DOUBLE.Pack(Qt_acc.getQt(it).data(),(it+1)*(it+1),source,16*Ndata,position,MPI::COMM_WORLD);
////		MPI::DOUBLE.Pack(Qt_sqacc.getQt(it).data(),(it+1)*(it+1),source,16*Ndata,position,MPI::COMM_WORLD);
////	}
////
////	//cout<<"MPI Reduce and Bcast"<<endl<<flush;
////	MPI::COMM_WORLD.Reduce(source,result,2*Ndata,MPI::DOUBLE,MPI::SUM,root);
////	MPI::COMM_WORLD.Bcast(result,2*Ndata,MPI::DOUBLE,root);
////
////	//cout<<"MPI Unpack"<<endl<<flush;
////	position = 0;
////	for(int it=0; it<Ntau+1; it++){
////		MPI::DOUBLE.Unpack(result,16*Ndata,output1,(it+1)*(it+1),position,MPI::COMM_WORLD);
////		MPI::DOUBLE.Unpack(result,16*Ndata,output2,(it+1)*(it+1),position,MPI::COMM_WORLD);
////
////		for(int j=0; j<(it+1)*(it+1); j++){
////			output1[j] /= NerrTot;
////			output2[j] = sqrt(abs(output2[j]/NerrTot - output1[j]*output1[j])/NerrTot);
////		}
////		Qt.setQt(it,output1);
////		Qt_err.setQt(it,output2);
////	}
////
////	delete[] source;
////	delete[] result;
////	delete[] output1;
////	delete[] output2;
////};
////
double dyson::weight(const int& n, const int& m){
	if(n==0) return 0.0;
	else if(m==0 or m==n) return 0.5;
	else return 1.0;
};
//void dyson::initTt(){
//	Tt.resize(Ntau+1);
//	for(int i=0; i<Tt.size(); i++){
//		Tt[i].resize(i+1);
//	}
//	//Tt = Eigen::MatrixXd::Zero(Ntau+1,Ntau+1);
//	Tt[0][0] = mWt[0]*Sz;
//};
//void dyson::initTt_amp(){
//	Tt_amp.resize(Ntau+1);
//	Tt_amp_up.resize(Ntau+1);
//	Tt_amp_dn.resize(Ntau+1);
//	Tt_amp_W.resize(Ntau+1);
//	for(int i=0; i<Tt.size(); i++){
//		Tt_amp[i].resize(i+1);
//		Tt_amp_up[i].resize(i+1);
//		Tt_amp_dn[i].resize(i+1);
//		Tt_amp_W[i].resize(i+1);;
//	}
//	Tt_amp[0][0] = Zero2x2;
// 	for(int n=0; n<Ntau+1; n++){
//		Tt_amp[n][0] = Zero2x2;
//		Tt_amp[n][n] = Zero2x2;
//		Tt_amp_up[n][n] = Zero2x2;
//		Tt_amp_dn[n][0] = Zero2x2;
//		Tt_amp_W[n][0] = Zero2x2;
//	}
//};
//double dyson::calculateBulkM_TCA(const int& n){
//	// T5 and T7 contribution
//	//return 0.25*dtau*dtau*mWt[n]*(1.0 + 0.25*dtau*dtau);
//	return 0.25*dtau*dtau*mWt[n];
//};
//double dyson::calculateBulkM_Q2(const int& n){
//	//return 0.25*dtau*dtau*mWt[n];
//	return pow(0.25*dtau*dtau*mWt[n],2);
//};
Eigen::MatrixXd dyson::calculateBulkB_OCA(const int& n, const int& m){
	return Wt[n]*VOp*Gt[n-m]*VOp*Gt[m]*VOp;;
}
////double dyson::calculateBulkB_Ladder(const int& n, const int& m){
////	double Bacc = calculateBulkB_OCA(n,m);
////	double Btmp;
////
////	// B5 //
////	Btmp = 0.0;
////	for(int p=1; p<m+1; p++){
////		Btmp += 0.5*weight(m,p)*Gt[p]*Gt[0]*Tt(n-p,m-p);
////	}
////	for(int q=m; q<n; q++){
////		Btmp += 0.5*weight(n-m,q-m)*Gt[0]*Gt[n-q]*Tt(q,m);
////	}
////	for(int p=1; p<m+1; p++) for(int q=m; q<n; q++){
////		Btmp += weight(m,p)*weight(n-m,q-m)*Gt[p]*Gt[n-q]*Tt(q-p,m-p);
////	}
////	Bacc += pow(dtau,2)*mWt[n]*Btmp;
////	//cout<<"after B5: Bacc = "<<Bacc<<endl;
////
////	return Bacc;
////}
////
//void dyson::calculateBulkTtAmputatedGUp(const int& m){
//	for(int p=0; p<m; p++){
//		Tt_amp_up[m][p] = Zero2x2;
//		for(int q=p; q<m+1; q++){
//			Tt_amp_up[m][p] += weight(m-p,q-p)
//				* Gt[m-q]*Tt[q][p];
//		}
//	}
//};
//void dyson::calculateBulkTtAmputatedGDn(const int& nm){
//	for(int q=1; q<nm+1; q++){
//		Tt_amp_dn[nm][q] = Zero2x2;
//		for(int p=0; p<q+1; p++){
//			Tt_amp_dn[nm][q] += weight(q,p)
//				* Tt[nm-p][q-p]*Gt[p];
//		}
//	}
//};
//void dyson::calculateBulkTtAmputatedW(const int& m){
//	if(m!=0){
//		for(int q=m+1; q<Ntau+1; q++){
//			Tt_amp_W[q][m] = Zero2x2;
//			for(int p=0; p<m+1; p++){
//				Tt_amp_W[q][m] += weight(m,p)
//					* mWt[q-p]*Tt[m][p];
//			}
//		}
//	}
//};
//void dyson::calculateBulkQtReducible(const int& n){
//	for(int m1=1; m1<n; m1++) for(int m2=m1; m2<n; m2++){
//		Qt_reducible[n][m2][m1] = Eigen::MatrixXd::Zero(4,4);
//		for(int p=m2; p<n+1; p++){
//			Qt_reducible[n][m2][m1] += weight(n-m2,p-m2)*tensorProd(Tt_amp_W[p][m1],Tt[n-m2][p-m2]);
//		}
//	}
//	/*for(int nm=1; nm<n; nm++){
//		Qt_reducible[n][n-nm][n] = Eigen::MatrixXd::Zero(4,4);
//		for(int p=0; p<nm+1; p++){
//			Qt_reducible[n][n-nm][n] += weight(nm,p)*tensorProd(Tt[nm][p],Tt_amp_W[n-nm+p][n]);
//		}
//	}*/
//};
//void dyson::calculateBulkTtAmputated(const int& n, const int& m){
//	Tt_amp[n][m] = Eigen::MatrixXd::Zero(2,2);
//	/*for(int q=m; q<n; q++){
//		Tt_amp[n][m] += 0.5*weight(n-m,q-m)
//			* Gt[n-q]*Tt[q][m];
//	}
//	for(int p=1; p<m+1; p++){
//		Tt_amp[n][m] += 0.5*weight(m,p)
//			* Tt[n-p][m-p]*Gt[p];
//	}*/
//	for(int p=0; p<m+1; p++) for(int q=m; q<n+1; q++){
//		Tt_amp[n][m] += weight(m,p)*weight(n-m,q-m)
//			* Gt[n-q]*Tt[q-p][m-p]*Gt[p];
//	}
//};
//Eigen::Matrix2d dyson::calculateBulkB_TCA(const int& n, const int& m){
//	// B1 //
//	Eigen::Matrix2d Bacc = calculateBulkB_OCA(n,m);
//	//if(n==2 && m==1) cout<<"after B1: Bacc = "<<Bacc<<endl;
//
//	Eigen::Matrix2d Btmp;
//
//	// B2 //
//	Btmp = Zero2x2;
//	for(int p=m+1; p<n+1; p++)
//		Btmp += weight(n-m,p-m)*mWt[p]*Tt_amp_dn[n-m][p-m]*Sz*Gt[m]*Sz;
//	Bacc += pow(dtau,2)*Btmp;
//	/*Btmp = Zero2x2;
//	for(int q=m; q<n+1; q++) for(int p=m; p<q+1; p++){
//		//cout<<"("<<Ntau<<","<<n-p<<","<<q-p<<")";
//		Btmp += weight(n-m,q-m)*weight(q-m,p-m)
//			* mWt[q]
//			* Tt[n-p][q-p]*Gt[p-m];
//	}
//	Bacc += pow(dtau,2)*Btmp*Sz*Gt[m]*Sz;*/
//	//if(n==2 && m==1) cout<<"after B2: Bacc = "<<Bacc<<endl;
//
//	// B3 //
//	Btmp = Zero2x2;
//	for(int p=0; p<m+1; p++)
//		Btmp += weight(m,p)*mWt[n-p]*Sz*Gt[n-m]*Sz*Tt_amp_up[m][p];
//	Bacc += pow(dtau,2)*Btmp;
//	/*Btmp = Zero2x2;
//	for(int p=0; p<m+1; p++) for(int q=p; q<m+1; q++) {
//		Btmp += weight(m,p)*weight(m-p,q-p)
//			* mWt[n-p]
//			* Gt[m-q]*Tt[q][p];
//	}
//	Bacc += pow(dtau,2)*Sz*Gt[n-m]*Sz*Btmp;*/
//	//if(n==2 && m==1) cout<<"after B3: Bacc = "<<Bacc<<endl;
//
//	// B4 //
//	Btmp = Zero2x2;
//	for(int p=0; p<m+1; p++) for(int q=m; q<n+1; q++)
//		Btmp += weight(m,p)*weight(n-m,q-m)*mWt[q-p]*Tt_amp_dn[n-m][q-m]*Sz*Tt_amp_up[m][p];
//	Bacc += pow(dtau,4)*Btmp;
//	/*Btmp = Zero2x2;
//	for(int p=0; p<m+1; p++) for(int q=p; q<m+1; q++) {
//		for(int s=m; s<n+1; s++) for(int r=m; r<s+1; r++) {
//		       Btmp += weight(m,p)*weight(m-p,q-p)*weight(n-m,s-m)*weight(s-m,r-m)
//			       * mWt[s-p]
//			       * Tt[n-r][s-r]*Gt[r-m]*Sz*Gt[m-q]*Tt[q][p];
//		}
//	}
//	Bacc += pow(dtau,4)*Btmp;*/
//	//if(n==2 && m==1) cout<<"after B4: Bacc = "<<Bacc<<endl;
//
//	// B5 //
//	Btmp = Zero2x2;
//	for(int p=1; p<m+1; p++){
//		Btmp += 0.5*weight(m,p)*Sz*Tt[n-p][m-p]*Gt[p]*Sz;
//	}
//	for(int q=m; q<n; q++){
//		Btmp += 0.5*weight(n-m,q-m)*Sz*Gt[n-q]*Tt[q][m]*Sz;
//	}
//	for(int p=1; p<m+1; p++) for(int q=m; q<n; q++){
//		Btmp += weight(m,p)*weight(n-m,q-m)*Sz*Gt[n-q]*Tt[q-p][m-p]*Gt[p]*Sz;
//	}
//	Bacc += pow(dtau,2)*mWt[n]*Btmp;
//	//if(n==2 && m==1) cout<<"after B5: Bacc = "<<Bacc<<endl;
//
//	// B6 //
//	Btmp = Zero2x2;
//	for(int q=m; q<n+1; q++) for(int p=m; p<q+1; p++){
//		if(!(p==n)){
//			Btmp += weight(n-m,q-m)*weight(q-m,p-m)
//				* mWt[q]
//				* Tt[n-p][q-p]*Tt_amp[p][m]*Sz;
//		}
//	}
//	Bacc += pow(dtau,4)*Btmp;
//	//if(n==2 && m==1) cout<<"after B6: Bacc = "<<Bacc<<endl;
//
//	// B7 //
//	Btmp = Zero2x2;
//	for(int p=0; p<m+1; p++) for(int q=p; q<m+1; q++){
//		if(!(q==0)){
//			Btmp += weight(m,p)*weight(m-p,q-p)
//				* mWt[n-p]
//				* Sz*Tt_amp[n-q][m-q]*Tt[q][p];
//		}
//	}
//	Bacc += pow(dtau,4)*Btmp;
//	//if(n==2 && m==1) cout<<"after B7: Bacc = "<<Bacc<<endl;
//
//	// B8 //
//	Eigen::MatrixXd Btmp_4x1 = Zero4x1;
//	for(int p=1; p<m+1; p++) for(int q=m; q<n; q++){
//		Btmp_4x1 += weight(m,p)*weight(n-m,q-m)
//			* Qt_reducible[n][q][p]*transformTo4x1(Tt_amp[q-p][m-p]);
//	}
//	Bacc += pow(dtau,6)*transformTo2x2(Btmp_4x1);
//
//	/*
//	Btmp = Zero2x2;
//	for(int r=0; r<m+1; r++) for(int q=0; q<r+1; q++) for(int p=0; p<q+1; p++){
//		for(int s=m; s<n+1; s++) for(int t=s; t<n+1; t++) for(int u=t; u<n+1; u++){
//			Btmp += weight(q,p)*weight(r,q)*weight(m,r)
//				* weight(n-m,s-m)*weight(n-s,t-s)*weight(n-t,u-t)
//				* mWt[u-p]
//				* Tt[n-t][u-t]*Gt[t-s]*Tt[s-r][m-r]*Gt[r-q]*Tt[q][p];
//		}
//	}
//	Bacc += pow(dtau,6)*Btmp;
//	//if(n==2 && m==1) cout<<"after B8: Bacc = "<<Bacc<<endl;
//	*/
//
//	return Bacc;
//}
//Eigen::Matrix2d dyson::calculateBulkB(const int& n, const int& m){
//	Eigen::MatrixXd Btmp;
//
//	Eigen::MatrixXd Bacc = transformTo4x1(calculateBulkB_TCA(n,m));
//
//	// B9 //
//	Btmp = Zero4x1;
//	for(int p=0; p<m+1; p++) for(int q=m; q<n+1; q++){
//		Btmp += weight(m,p)*weight(n-m,q-m)
//			*Qt[n][q][p]*transformTo4x1(Gt[q-m]*Sz*Gt[m-p]);
//	}
//	Bacc += pow(dtau,2)*Btmp;
//	//cout<<"after B9: Bacc = "<<Bacc<<endl;
//
//	// B10 //
//	Btmp = Zero4x1;
//	//for(int r=m; r<n; r++) for(int s=r; s<n+1; s++){
//	//	Btmp += 0.25*weight(n-m,r-m)*weight(n-r,s-r)
//	//		* Qt[n][s][0]*transformTo4x1(Gt[s-r]*Tt[r][m]*Gt[0]);
//	//}
//	//for(int q=1; q<m+1; q++) for(int p=0; p<q+1; p++){
//	//	Btmp += 0.25*weight(q,p)*weight(m,q)
//	//		* Qt[n][n][p]*transformTo4x1(Gt[0]*Tt[n-q][m-q]*Gt[q-p]);
//	//}
//	for(int p=1; p<m+1; p++){
//		Btmp += 0.125*weight(m,p)
//			* Qt[n][n][0]*transformTo4x1(Tt_amp_up[n-p][m-p]*Gt[p]);
//		Btmp += 0.5*weight(m,p)
//			* Qt[n][n][p]*transformTo4x1(Tt_amp[n-p][m-p]);
//	}
//	for(int q=m; q<n; q++){
//		Btmp += 0.125*weight(n-m,n-q)
//			* Qt[n][n][0]*transformTo4x1(Gt[n-q]*Tt_amp_dn[q][m]);
//		Btmp += 0.5*weight(n-m,q-m)
//			* Qt[n][q][0]*transformTo4x1(Tt_amp[q][m]);
//	}
//	for(int p=1; p<m+1; p++) for(int q=m; q<n; q++) {
//		Btmp += weight(m,p)*weight(n-m,q-m)
//			* Qt[n][q][p]*transformTo4x1(Tt_amp[q-p][m-p]);
//	}
//	Bacc += pow(dtau,4)*Btmp;
//	//cout<<"after B10: Bacc = "<<Bacc<<endl;
//
//	return transformTo2x2(Bacc);
//};
vector<Eigen::MatrixXd> dyson::calculateBoundaryB_OCA(const int& n){
	//cout<<"* calculateBoundaryB *: n="<<n<<endl;
	vector<Eigen::MatrixXd> Bacc(2);
	Bacc[0] = Wt[n]*VOp*Gt[n]*VOp*VOp;
	Bacc[1] = Wt[n]*VOp*VOp*Gt[n]*VOp;
	//cout<<"  T1 contribution: "<<endl<<Bacc<<endl;

	return Bacc;
};
//vector<Eigen::Matrix2d> dyson::calculateBoundaryB(const int& n){
//	//cout<<"* calculateBoundaryB *: n="<<n<<endl;
//	//cout<<"  T1 contribution: "<<endl<<Bacc<<endl;
//	vector<Eigen::Matrix2d> Bacc = calculateBoundaryB_OCA(n);
//
//	//cout<<"  T2 contribution: "<<endl;
//	Eigen::Matrix2d Btmp = Eigen::MatrixXd::Zero(2,2);
//	for(int q=1; q<n; q++){
//		Btmp += 0.5*mWt[q]*Tt[n][q];
//	}
//	for(int p=1; p<n+1; p++) for(int q=p; q<n+1; q++){
//		Btmp += weight(n,p)*weight(n-p,q-p)*mWt[q]*Tt[n-p][q-p]*Gt[p];
//	}
//	Bacc[0] += pow(dtau,2)*Btmp;
//	//cout<<"  T2 contribution: "<<pow(dtau,2)*Btmp<<endl;
//
//	//cout<<"  T3 contribution: "<<pow(dtau,2)*Btmp<<endl;
//	Btmp = Eigen::MatrixXd::Zero(2,2);
//	for(int p=1; p<n; p++)
//		Btmp += 0.5*mWt[n-p]*Tt[n][p];
//	for(int q=0; q<n; q++) for(int p=0; p<q+1; p++){
//		Btmp += weight(q,p)*weight(n,q)*mWt[n-p]*Gt[n-q]*Tt[q][p];
//	}
//	Bacc[1] += pow(dtau,2)*Btmp;
//	//cout<<"  T3 contribution: "<<pow(dtau,2)*Btmp<<endl;
//
//	//cout<<"**********************"<<endl;
//
//	return Bacc;
//};
//Eigen::Matrix2d dyson::calculateBoundaryM(const int& n){
//	Eigen::Matrix2d Mtmp;
//	Mtmp(0,0) = mWt[0];
//	Mtmp(0,1) = mWt[n];
//	Mtmp(1,0) = mWt[n];
//	Mtmp(1,0) = mWt[0];
//
//	//Mtmp(0,0) = 0.0;
//	//Mtmp(0,1) = mWt[n];
//	//Mtmp(1,0) = mWt[n];
//	//Mtmp(1,0) = mWt[0];
//
//	//Mtmp(1,1) = 0.0;
//
//	return 0.25*pow(dtau,2)*Mtmp;
//};
//void dyson::calculateBoundaryTtAmputated(const int& n){
//	// Bulk Tt_amp //
//	/*
//	for(int p=0; p<n; p++){
//		Tt_amp_up[n][p] = Eigen::MatrixXd::Zero(2,2);
//		for(int q=p; q<n+1; q++){
//			Tt_amp_up[n][p] += weight(n-p,q-p)
//				* Gt[n-q]*Tt[q][p];
//		}
//	}
//
//	int nm = n-m;
//	for(int q=0; q<nm+1; q++){
//		Tt_amp_dn[nm][q] = Eigen::MatrixXd::Zero(2,2);
//		for(int p=0; p<q+1; p++){
//			Tt_amp_dn[nm][q] += weight(q,p)
//				* Tt[nm-p][q-p]*Gt[p];
//		}
//	}
//	*/
//
//	Tt_amp[n][0] = Eigen::MatrixXd::Zero(2,2);
//	Tt_amp[n][n] = Eigen::MatrixXd::Zero(2,2);
//};
//
////void dyson::TupdateBruteForce_OCA(){
////	//Tupdate_OCA();
////	//Eigen::MatrixXd Tt_new;
////
////	Tt = Eigen::MatrixXd::Zero(Ntau+1,Ntau+1);
////
////	for(int m=0; m<Ntau+1; m++) for(int n=m; n<Ntau+1; n++){
////		Tt(n,m) = mWt[n]*Gt[m]*Gt[n-m];
////		//if(n==1) cout<<T(n,m)<<" ";
////	}
////	//cout<<endl;
////	//Bacc += (mWt[n])*Gt[m]*Gt[n-m];
////	//Bacc << mWt[n]*Gt[n] , mWt[n]*Gt[n];
////
////	/*int Niter = 100;
////	for(int i=0; i<Niter; i++){
////		Tt_new = Eigen::MatrixXd::Zero(Ntau+1,Ntau+1);
////		for(int m=0; m<Ntau+1; m++) for(int n=m; n<Ntau+1; n++){
////			Tt_new(n,m) = mWt[n]*Gt[m]*Gt[n-m];
////			for(int p=0; p<m+1; p++) for(int q=m; q<n+1; q++){
////				Tt_new(n,m) += dtau*dtau*weight(m,p)*weight(n-m,q-m)
////						*mWt[n]*Gt[p]*Gt[n-q]*Tt(q-p,m-p);
////			}
////		}
////		Tt = Tt_new;
////	}*/
////};
////void dyson::TupdateBruteForce_Ladder(){
////	Tupdate_OCA();
////	Eigen::MatrixXd Tt_new;
////
////	int Niter = 100;
////	for(int i=0; i<Niter; i++){
////		Tt_new = Eigen::MatrixXd::Zero(Ntau+1,Ntau+1);
////		for(int m=0; m<Ntau+1; m++) for(int n=m; n<Ntau+1; n++){
////			Tt_new(n,m) = mWt[n]*Gt[m]*Gt[n-m];
////			for(int p=0; p<m+1; p++) for(int q=m; q<n+1; q++){
////				Tt_new(n,m) += dtau*dtau*weight(m,p)*weight(n-m,q-m)
////						*mWt[n]*Gt[p]*Gt[n-q]*Tt(q-p,m-p);
////			}
////		}
////		Tt = Tt_new;
////	}
////};
////void dyson::TupdateBruteForce_TCA(){
////	Tupdate_OCA();
////	Eigen::MatrixXd Tt_new;
////
////	int Niter = 300;
////	for(int i=0; i<Niter; i++){
////		Tt_new = Eigen::MatrixXd::Zero(Ntau+1,Ntau+1);
////		for(int m=0; m<Ntau+1; m++) for(int n=m; n<Ntau+1; n++){
////			// T1 //
////			Tt_new(n,m) = mWt[n]*Gt[m]*Gt[n-m];
////
////			// T2 //
////			for(int p=m; p<n+1; p++) for(int q=p; q<n+1; q++){
////				Tt_new(n,m) += pow(dtau,2)*weight(n-m,p-m)*weight(n-p,q-p)
////						*mWt[q]*Gt[m]*Gt[p-m]*Tt(n-p,q-p);
////			}
////
////			// T3 //
////			for(int p=0; p<m+1; p++) for(int q=p; q<m+1; q++){
////				Tt_new(n,m) += pow(dtau,2)*weight(m,p)*weight(m-p,q-p)
////						*mWt[n-p]*Gt[m-q]*Gt[n-m]*Tt(q,p);
////			}
////
////			// T4 //
////			for(int p=0; p<m+1; p++) for(int q=p; q<m+1; q++)
////			for(int r=m; r<n+1; r++) for(int s=r; s<n+1; s++){
////				Tt_new(n,m) += pow(dtau,4)*weight(m,p)*weight(m-p,q-p)*weight(n-m,r-m)*weight(n-r,s-r)
////						*mWt[s-p]*Gt[m-q]*Gt[r-m]*Tt(q,p)*Tt(n-r,s-r);
////			}
////
////			// T5 //
////			for(int p=0; p<m+1; p++) for(int q=m; q<n+1; q++){
////				Tt_new(n,m) += pow(dtau,2)*weight(m,p)*weight(n-m,q-m)
////						*mWt[n]*Gt[p]*Gt[n-q]*Tt(q-p,m-p);
////			}
////
////			// T6 //
////			for(int p=0; p<m+1; p++) for(int q=m; q<n+1; q++)
////			for(int r=q; r<n+1; r++) for(int s=r; s<n+1; s++){
////				Tt_new(n,m) += pow(dtau,4)*weight(m,p)*weight(n-m,q-m)*weight(n-q,r-q)*weight(n-r,s-r)
////						*mWt[s]*Gt[p]*Gt[r-q]*Tt(q-p,m-p)*Tt(n-r,s-r);
////			}
////
////			// T7 //
////			for(int p=0; p<m+1; p++) for(int q=p; q<m+1; q++)
////			for(int r=q; r<m+1; r++) for(int s=m; s<n+1; s++){
////				Tt_new(n,m) += pow(dtau,4)*weight(m,p)*weight(m-p,q-p)*weight(m-q,r-q)*weight(n-m,s-m)
////						*mWt[n-p]*Gt[r-q]*Gt[n-s]*Tt(q,p)*Tt(s-r,m-r);
////			}
////
////			// T8 //
////			for(int p=0; p<m+1; p++) for(int q=p; q<m+1; q++) for(int r=q; r<m+1; r++) 
////			for(int s=m; s<n+1; s++) for(int t=s; t<n+1; t++) for(int u=t; u<n+1; u++){
////				Tt_new(n,m) += pow(dtau,6)*weight(m,p)*weight(m-p,q-p)*weight(m-q,r-q)
////						*weight(n-m,s-m)*weight(n-s,t-s)*weight(n-t,u-t)
////						*mWt[u-p]*Gt[r-q]*Gt[t-s]*Tt(q,p)*Tt(s-r,m-r)*Tt(n-t,u-t);
////			}
////		}
////		Tt = Tt_new;
////	}
////};
void dyson::Tupdate_TOA(){
	Tupdate_OCA();

	std::vector<std::vector<Eigen::MatrixXd> > Tt_new(Ntau+1);
    for(int n=0; n<Ntau+1; n++) Tt_new[n].resize(n+1);

    Tt_new[0][0] = Tt[0][0];
	for(int n=1; n<Ntau+1; n++){
        for(int m=0; m<n+1; m++){
            // T1 //
            Tt_new[n][m] = Tt[n][m];

            // T2 //
            for(int p=m; p<n+1; p++) for(int q=p; q<n+1; q++){
                Tt_new[n][m] += pow(dtau,2)*weight(n-m,p-m)*weight(n-p,q-p)
                        *Wt[q]*Tt[n-p][q-p]*Gt[p-m]*VOp*Gt[m]*VOp;
            }

            // T3 //
            for(int p=0; p<m+1; p++) for(int q=p; q<m+1; q++){
                Tt_new[n][m] += pow(dtau,2)*weight(m,p)*weight(m-p,q-p)
                        *Wt[n-p]*VOp*Gt[n-m]*VOp*Gt[m-q]*Tt[q][p];
            }

            // T5 //
            for(int p=0; p<m+1; p++) for(int q=m; q<n+1; q++){
                Tt_new[n][m] += pow(dtau,2)*weight(m,p)*weight(n-m,q-m)
                        *Wt[n]*VOp*Gt[n-q]*Tt[q-p][m-p]*Gt[p]*VOp;
            }

            // T9 //
            for(int p=0; p<m+1; p++) for(int q=m; q<n+1; q++){
                Tt_new[n][m] += pow(dtau,2)*weight(m,p)*weight(n-m,q-m)
                    *Wt[q]*Wt[n-p]*VOp*Gt[n-q]*VOp*Gt[q-m]*VOp*Gt[m-p]*VOp*Gt[p]*VOp;
            }
        }
	}
	Tt = Tt_new;
};
void dyson::Tupdate_TOA(const Eigen::MatrixXd& X){
	Tupdate_OCA();

	std::vector<std::vector<Eigen::MatrixXd> > Tt_new(Ntau+1);
    for(int n=0; n<Ntau+1; n++) Tt_new[n].resize(n+1);
    Tt_new[0][0] = Wt[0]*VOp*X*VOp;
	for(int n=1; n<Ntau+1; n++) for(int m=0; m<n+1; m++){
		// T1 //
		Tt_new[n][m] = Wt[n]*VOp*Gt[n-m]*X*Gt[m]*VOp;

		// T2 //
		for(int p=m; p<n+1; p++) for(int q=p; q<n+1; q++){
			Tt_new[n][m] += pow(dtau,2)*weight(n-m,p-m)*weight(n-p,q-p)
					*Wt[q]*Tt[n-p][q-p]*Gt[p-m]*X*Gt[m]*VOp;
		}

		// T3 //
		for(int p=0; p<m+1; p++) for(int q=p; q<m+1; q++){
			Tt_new[n][m] += pow(dtau,2)*weight(m,p)*weight(m-p,q-p)
					*Wt[n-p]*VOp*Gt[n-m]*X*Gt[m-q]*Tt[q][p];
		}

		// T5 + T9 //
		for(int p=0; p<m+1; p++) for(int q=m; q<n+1; q++){
			Tt_new[n][m] += pow(dtau,2)*weight(m,p)*weight(n-m,q-m)
                *(Wt[n]*Wt[q-p] + Wt[q]*Wt[n-p])
                *VOp*Gt[n-q]*VOp*Gt[q-m]*X*Gt[m-p]*VOp*Gt[p]*VOp;
		}

		// T9 //
		//for(int p=0; p<m+1; p++) for(int q=m; q<n+1; q++){
		//	Tt_new[n][m] += pow(dtau,2)*weight(m,p)*weight(n-m,q-m)
		//		*Wt[q]*Wt[n-p]*VOp*Gt[n-q]*VOp*Gt[q-m]*X*Gt[m-p]*VOp*Gt[p]*VOp;
		//}
	}
	Tt = Tt_new;
};

void dyson::Tupdate_OCA(){
	Tt.resize(Ntau+1);
	for(int i=0; i<Tt.size(); i++){
		Tt[i].resize(i+1);
	}

	//Tt[0][0] = mWt[0]*Eigen::MatrixXd::Identity(2,2);
	Tt[0][0] = Wt[0]*VOp*VOp*VOp;

	vector<Eigen::MatrixXd> Tboundary;
	for(int n=1; n<Ntau+1; n++){
		for(int m=1; m<n; m++){
			Tt[n][m] = calculateBulkB_OCA(n,m);
		}
		Tboundary = calculateBoundaryB_OCA(n);
		Tt[n][0] = Tboundary[0];
		Tt[n][n] = Tboundary[1];
		//if(n==1) cout<<Tboundary(0)<<" "<<Tboundary(1)<<endl;
	}
};

void dyson::Tupdate_OCA(const Eigen::MatrixXd& X){
	Tt.resize(Ntau+1);
	for(int i=0; i<Tt.size(); i++){
		Tt[i].resize(i+1);
	}
	for(int n=0; n<Ntau+1; n++) for(int m=0; m<n+1; m++){
        Tt[n][m] = Wt[n]*VOp*Gt[n-m]*X*Gt[m]*VOp;
	}
};

////void dyson::Tupdate_Ladder(){
////	Tt = Eigen::MatrixXd::Zero(Ntau+1,Ntau+1);
////	Tt(0,0) = mWt[0];
////
////	Eigen::MatrixXd Tboundary;
////	for(int n=1; n<Ntau+1; n++){
////		for(int m=1; m<n; m++){
////			Tt(n,m) = calculateBulkB_Ladder(n,m)/(1.0 - calculateBulkM_TCA(n));
////		}
////		Tboundary = calculateBoundaryB_OCA(n);
////		Tt(n,0) = Tboundary(0);
////		Tt(n,n) = Tboundary(1);
////	}
////};
//void dyson::Tupdate_TCA0(){
//	Tupdate_OCA();
//
//	for(int n=0; n<Ntau+1; n++){
//		for(int m=0; m<n+1; m++){
//			for(int p=0; p<m+1; p++) for(int q=m; q<n+1; q++){
//				Tt[n][m] += pow(dtau,2)
//					* weight(m,p)*weight(n-m,q-m)
//					* mWt[n]*mWt[q-p]
//					* Sz*Gt[n-q]*Sz*Gt[q-m]*Sz*Gt[m-p]*Sz*Gt[p]*Sz;
//			}
//			for(int q=m; q<n+1; q++) for(int p=m; p<q+1; p++){
//				Tt[n][m] += pow(dtau,2)
//					* weight(n-m,q-m)*weight(q-m,p-m)
//					* mWt[q]*mWt[n-p]
//					* Sz*Gt[n-q]*Sz*Gt[q-p]*Sz*Gt[p-m]*Sz*Gt[m]*Sz;
//			}
//			for(int p=0; p<m+1; p++) for(int q=p; q<m+1; q++){
//				Tt[n][m] += pow(dtau,2)
//					* weight(m,p)*weight(m-p,q-p)
//					* mWt[q]*mWt[n-p]
//					* Sz*Gt[n-m]*Sz*Gt[m-q]*Sz*Gt[q-p]*Sz*Gt[p]*Sz;
//			}
//		}
//	}
//};
//void dyson::Tupdate_Q2_0(){
//	Tupdate_OCA();
//
//	for(int n=0; n<Ntau+1; n++){
//		for(int m=0; m<n+1; m++){
//			for(int p=0; p<m+1; p++) for(int q=m; q<n+1; q++){
//				Tt[n][m] += pow(dtau,2)
//					* weight(m,p)*weight(n-m,q-m)
//					* (mWt[n]*mWt[q-p] + mWt[q]*mWt[n-p])
//					* Sz*Gt[n-q]*Sz*Gt[q-m]*Sz*Gt[m-p]*Sz*Gt[p]*Sz;
//			}
//			for(int q=m; q<n+1; q++) for(int p=m; p<q+1; p++){
//				Tt[n][m] += pow(dtau,2)
//					* weight(n-m,q-m)*weight(q-m,p-m)
//					* mWt[q]*mWt[n-p]
//					* Sz*Gt[n-q]*Sz*Gt[q-p]*Sz*Gt[p-m]*Sz*Gt[m]*Sz;
//			}
//			for(int p=0; p<m+1; p++) for(int q=p; q<m+1; q++){
//				Tt[n][m] += pow(dtau,2)
//					* weight(m,p)*weight(m-p,q-p)
//					* mWt[q]*mWt[n-p]
//					* Sz*Gt[n-m]*Sz*Gt[m-q]*Sz*Gt[q-p]*Sz*Gt[p]*Sz;
//			}
//		}
//	}
//};
//
//void dyson::Tupdate_TCA(){
//	initTt();
//	initTt_amp();
//	initQtReducible();
//
//	Eigen::Matrix2d BulkB_tmp;
//	vector<Eigen::Matrix2d> BoundaryB_tmp;
//	Eigen::Matrix2d Kernel;
//	for(int n=1; n<Ntau+1; n++){
//		calculateBulkQtReducible(n);
//
//		/*if(n==2){
//			cout<<endl;
//			for(int m2=0; m2<n+1; m2++){
//				for(int m1=0; m1<m2+1; m1++){
//					cout<<"("<<n<<","<<m2<<","<<m1<<")"<<endl;
//					cout<<Qt_reducible[n][m2][m1]<<endl;
//				}
//				cout<<endl;
//			}
//			cout<<endl;
//		}*/
//
//		for(int m=1; m<n; m++){
//			//Tt[n][m] = calculateBulkB_TCA(n,m)/(1.0 - calculateBulkM_TCA(n));
//			BulkB_tmp = calculateBulkB_TCA(n,m);
//
//			Tt[n][m](0,0) = BulkB_tmp(0,0)/(1.0 - calculateBulkM_TCA(n));
//			Tt[n][m](0,1) = BulkB_tmp(0,1)/(1.0 + calculateBulkM_TCA(n));
//			Tt[n][m](1,0) = BulkB_tmp(1,0)/(1.0 + calculateBulkM_TCA(n));
//			Tt[n][m](1,1) = BulkB_tmp(1,1)/(1.0 - calculateBulkM_TCA(n));
//
//			calculateBulkTtAmputated(n,m);
//		}
//		BoundaryB_tmp = calculateBoundaryB(n);
//		Kernel = ( Eigen::MatrixXd::Identity(2,2) - calculateBoundaryM(n) ).inverse();
//
//		Tt[n][0] = Kernel(0,0)*BoundaryB_tmp[0] + Kernel(0,1)*BoundaryB_tmp[1];
//		Tt[n][n] = Kernel(1,0)*BoundaryB_tmp[0] + Kernel(1,1)*BoundaryB_tmp[1];
//
//		calculateBulkTtAmputatedGUp(n);
//		calculateBulkTtAmputatedGDn(n);
//		calculateBulkTtAmputatedW(n);
//
//		//calculateBoundaryTtAmputated(n);
//		//calculateBoundaryTtAmputated(n);
//	}
//};
//void dyson::Tupdate(){
//	/*Tt.resize(Ntau+1);
//	for(int i=0; i<Tt.size(); i++){
//		Tt[i].resize(i+1);
//	}
//	Tt[0][0] = mWt[0]*Sz;*/
//	initTt();
//	initTt_amp();
//	initQtReducible();
//
//	double m5, m10;
//	Eigen::Matrix2d BulkB_tmp;
//	vector<Eigen::Matrix2d> BoundaryB_tmp;
//	Eigen::Matrix2d Kernel;
//	for(int n=1; n<Ntau+1; n++){
//		calculateBulkQtReducible(n);
//
//		m5 = calculateBulkM_TCA(n);
//		m10 = calculateBulkM_Q2(n);
//		for(int m=1; m<n; m++){
//			//cout<<"("<<n<<","<<m<<")"<<flush;
//			BulkB_tmp = calculateBulkB(n,m);
//
//			Tt[n][m](0,0) = BulkB_tmp(0,0)/(1.0 - m5 - m10);
//			Tt[n][m](0,1) = BulkB_tmp(0,1)/(1.0 + m5 - m10);
//			Tt[n][m](1,0) = BulkB_tmp(1,0)/(1.0 + m5 - m10);
//			Tt[n][m](1,1) = BulkB_tmp(1,1)/(1.0 - m5 - m10);
//
//			calculateBulkTtAmputated(n,m);
//		}
//		BoundaryB_tmp = calculateBoundaryB(n);
//		Kernel = ( Eigen::MatrixXd::Identity(2,2) - calculateBoundaryM(n) ).inverse();
//
//		Tt[n][0] = Kernel(0,0)*BoundaryB_tmp[0] + Kernel(0,1)*BoundaryB_tmp[1];
//		Tt[n][n] = Kernel(1,0)*BoundaryB_tmp[0] + Kernel(1,1)*BoundaryB_tmp[1];
//
//		calculateBulkTtAmputatedGUp(n);
//		calculateBulkTtAmputatedGDn(n);
//		calculateBulkTtAmputatedW(n);
//	}
//};
////void dyson::TupdateDebug(){
////	Tt = Eigen::MatrixXd::Zero(Ntau+1,Ntau+1);
////	Tt(0,0) = mWt[0];
////	cout<<"Tt(0,0) = "<<Tt(0,0)<<endl;
////
////	Eigen::MatrixXd Tboundary;
////	for(int n=1; n<2; n++){
////		for(int m=1; m<n; m++){
////			cout<<"("<<n<<","<<m<<")";
////			cout<<", BuldB: "<<calculateBulkB(n,m)<<endl;
////			cout<<", BuldM: "<<calculateBulkM(n)<<endl;
////			Tt(n,m) = calculateBulkB(n,m)/(1.0 - calculateBulkM(n));
////			cout<<"Tt = "<<Tt(n,m)<<endl;
////		}
////		cout<<"[("<<n<<","<<0<<"),("<<n<<","<<n<<")] "<<endl;
////		cout<<"BoundaryB: "<<endl
////			<<calculateBoundaryB(n)<<endl;
////		cout<<"BoundaryM: "<<endl
////			<<calculateBoundaryM(n)<<endl;
////		cout<<endl;
////
////		cout<<"(I-BoundaryM).inverse() = "<<endl
////			<<( Eigen::MatrixXd::Identity(2,2) - calculateBoundaryM(n) ).inverse()<<endl;
////
////		Tboundary = ( Eigen::MatrixXd::Identity(2,2) - calculateBoundaryM(n) ).inverse()*calculateBoundaryB(n);
////		Tt(n,0) = Tboundary(0);
////		cout<<"Tt = "<<setprecision(15)<<Tt(n,0)<<endl;
////		Tt(n,n) = Tboundary(1);
////		cout<<"Tt = "<<setprecision(15)<<Tt(n,n)<<endl;
////	}
////};
////
////void dyson::setSt(const double& s_value){
////	St.resize(Ntau+1);
////
////	for(int n=0; n<Ntau+1; n++){
////		St[n] = 0.0;
////	}
////};
////
void dyson::updateSt_NCA(){
	St.resize(Ntau+1);

	for(int n=0; n<Ntau+1; n++){
		St[n] = Wt[n]*VOp*Gt[n]*VOp;
	}
};

////void dyson::updateSt_OCA(){
////	St.resize(Ntau+1);
////
////	double St_tmp;
////	for(int n=0; n<Ntau+1; n++){
////		St_tmp = 0.0;
////		for(int p=0; p<n+1; p++) for(int q=p; q<n+1; q++){
////			St_tmp += weight(n,p)*weight(n-p,q-p)*Gt[p]*Gt[q-p]*Gt[n-q]*mWt[q]*mWt[n-p];
////		}
////		St[n] = mWt[n]*Gt[n] + pow(dtau,2)*St_tmp;
////	}
////};
void dyson::updateSt(){
	St.resize(Ntau+1);

	Eigen::MatrixXd St_tmp;
	for(int n=0; n<Ntau+1; n++){
		St_tmp = ZeroM;
		for(int p=0; p<n+1; p++) for(int q=p; q<n+1; q++){
			St_tmp += weight(n,p)*weight(n-p,q-p)*Wt[n-p]*VOp*Gt[n-q]*Tt[q][p];
		}
		St[n] = Wt[n]*VOp*Gt[n]*VOp + pow(dtau,2)*St_tmp;
	}
};
////void dyson::updateGt_dyson(const double& epsilon, const int& Niter){
////	//cout<<"dyson equation with lambda = "<<lambda<<endl;
////	//vector<double> Gt0(Ntau+1);
////	Gt_prev.resize(Ntau+1);
////	for(int n=0; n<Ntau+1; n++){
////		Gt_prev[n] = Gt[n];
////	}
////
////	for(int n=0; n<Ntau+1; n++){
////		//Gt0[n] = exp(-epsilon*n*dtau);
////		Gt[n] = G0t[n];
////	}
////
////	vector<double> Gt_new(Ntau+1);
////	for(int i=0; i<Niter; i++){
////		for(int n=0; n<Ntau+1; n++){
////			Gt_new[n] = 0.0;
////			for(int p=0; p<n+1; p++) for(int q=p; q<n+1; q++){
////				Gt_new[n] += weight(n,p)*weight(n-p,q-p)*G0t[n-q]*St[q-p]*Gt[p];
////			}
////
////			Gt_new[n] *= pow(dtau,2);
////			Gt_new[n] += G0t[n];
////		}
////
////		for(int n=0; n<Ntau+1; n++)
////			Gt[n] = Gt_new[n];
////	}
////};
////
void dyson::updateGt_volterra(){
	//cout<<"volterra equation with lambda = "<<lambda<<endl;
	Gt_prev.resize(Ntau+1);
	for(int n=0; n<Ntau+1; n++){
		Gt_prev[n] = Gt[n];
	}

	Eigen::MatrixXd Katom = Hloc - lambda*IdentityM;
	vector<Eigen::MatrixXd> dGt(Ntau+1);
	Gt[0] = IdentityM;
	dGt[0] = -Katom*Gt[0];

	for(int n=0; n<Ntau; n++){
		//if(n==0) cout<<Gt[n]<<endl;
		Gt[n+1] = Gt[n];
		//if(n==0) cout<<Gt[n+1]<<endl;
		Gt[n+1] -= 0.5*Katom*dtau*( 2.0*Gt[n] + dtau*dGt[n] );
		//if(n==0) cout<<Gt[n+1]<<endl;
		for(int p=0; p<n+1; p++){
			Gt[n+1] += 0.5*pow(dtau,2)
				*( weight(n,p)*St[n-p] + weight(n+1,p)*St[n+1-p])*Gt[p];
			//if(n==0) cout<<Gt[n+1]<<endl;
		}
		Gt[n+1] += 0.25*pow(dtau,2)*St[0]*( Gt[n] + dtau*dGt[n] );
		//if(n==0) cout<<Gt[n+1]<<endl;

		dGt[n+1] = -Katom*Gt[n+1];
		for(int p=0; p<n+2; p++)
			dGt[n+1] += dtau*weight(n+1,p)*St[n+1-p]*Gt[p];
	}
};
void dyson::normalizeGt(){
	lambda_prev = lambda;
	double dlambda = -log(abs(Gt[Ntau].trace()))/beta;
	double expDtauDLambda = exp(dtau*dlambda);

	double factor = 1.0;
	for(int n=0; n<Ntau+1; n++){
		//G0t[n] *= factor;
		Gt[n] *= factor;
		factor *= expDtauDLambda;
	}

	lambda += dlambda;
	//cout<<"* normalization *"<<endl;
	//cout<<"  lambda = "<<lambda<<endl;
	//cout<<"  lambda_prev = "<<lambda_prev<<endl;
};
//
//void dyson::mixingGt(){
//	//lambda = 0.5*(lambda + lambda_prev);
//	lambda = lambda + lambda_prev;
//
//	for(int n=0; n<Ntau+1; n++){
//		Gt[n] = 0.5*(Gt[n] + Gt_prev[n]);
//	}
//};
//
void dyson::calculateStS0(const Eigen::MatrixXd& X){
	//cout<<"*  calculateStS0 **"<<endl;
	double StS0_tmp;
	StS0.resize(Ntau+1);
	for(int i=0; i<Ntau+1; i++){
		//cout<<i<<" "<<flush;
		StS0[i] = (Gt[Ntau-i]*X*Gt[i]*X).trace();
		StS0_tmp = 0.0;
		for(int p=0; p<i+1; p++) for(int q=i; q<Ntau+1; q++)
			StS0_tmp += weight(i,p)*weight(Ntau-i,q-i)*(Gt[Ntau-q]*Tt[q-p][i-p]*Gt[p]*X).trace();
		StS0[i] += pow(dtau,2)*StS0_tmp;
	}
	//cout<<endl;
};

void dyson::calculateStS0_NCA(const Eigen::MatrixXd& X){
	StS0.resize(Ntau+1);
	for(int i=0; i<Ntau+1; i++){
		StS0[i] = (Gt[Ntau-i]*X*Gt[i]*X).trace();
		//StS0[i] = (Gt[Ntau-i]*Gt[i]).trace();
	}
};

////void dyson::calculateStS0_OCA(){
////	//cout<<"*  calculateStS0 **"<<endl;
////	double StS0_tmp;
////	StS0.resize(Ntau+1);
////	for(int i=0; i<Ntau+1; i++){
////		//cout<<i<<" "<<flush;
////		StS0[i] = Gt[i]*Gt[Ntau-i];
////		StS0_tmp = 0.0;
////		for(int p=0; p<i+1; p++) for(int q=i; q<Ntau+1; q++)
////			StS0_tmp += weight(i,p)*weight(Ntau-i,q-i)*Gt[p]*Gt[i-p]*Gt[q-i]*Gt[Ntau-q]*mWt[q-p];
////		StS0[i] += pow(dtau,2)*StS0_tmp;
////	}
////	//cout<<endl;
////};
////
void dyson::printWt(ostream& ostr){
	double tau;
	ostr<<scientific;
	for(int i=0; i<Ntau+1; i++){
		tau = i*dtau;
		ostr<<setw(30)<<left<<setprecision(10)<<tau
			<<setw(30)<<right<<setprecision(10)<<Wt[i]<<endl;
	}
};
void dyson::printHloc(std::ostream& ostr){
	for(int i=0; i<Nflavor; i++){
		for(int j=0; j<Nflavor; j++) ostr<<setw(30)<<right<<Hloc(i,j);
		ostr<<endl;
	}
};
void dyson::printVOp(std::ostream& ostr){
	for(int i=0; i<Nflavor; i++){
		for(int j=0; j<Nflavor; j++) ostr<<setw(30)<<right<<VOp(i,j);
		ostr<<endl;
	}
};

void dyson::printGt(ostream& ostr){
	double tau;
	ostr<<scientific;
	ostr<<setw(10)<<left<<"# lambda: "
		<<setw(30)<<left<<setprecision(10)<<lambda<<endl;
	for(int i=0; i<Ntau+1; i++){
		tau = i*dtau;
		ostr<<setw(30)<<left<<setprecision(10)<<tau;
		for(int a=0; a<Nflavor; a++) for(int b=0; b<Nflavor; b++) 
			ostr<<setw(30)<<right<<setprecision(10)<<Gt[i](a,b);
		ostr<<endl;
	}
};
//void dyson::printQt(ostream& ostr, const int& n, const int& iflavor, const int& jflavor){
//	double tau;
//	ostr<<scientific;
//	for(int i=0; i<n+1; i++){
//		for(int j=0; j<i+1; j++){
//			ostr<<setw(30)<<right<<setprecision(15)<<Qt[n][i][j](iflavor,jflavor);
//		}
//		ostr<<endl;
//	}
//};
////void dyson::printQt_err(ostream& ostr, const int& n){
////	double tau;
////	ostr<<scientific;
////	for(int i=0; i<n+1; i++){
////		for(int j=0; j<n+1; j++){
////			ostr<<setw(30)<<right<<setprecision(15)<<Qt_err.getQt(n,i,j);
////		}
////		ostr<<endl;
////	}
////};
//void dyson::printTt(ostream& ostr, const int& iflavor, const int& jflavor){
//	double tau;
//	ostr<<scientific;
//	for(int n=0; n<Ntau+1; n++){
//		for(int m=0; m<n+1; m++){
//			ostr<<setw(15)<<right<<setprecision(7)<<Tt[n][m](iflavor,jflavor);
//		}
//		ostr<<endl;
//	}
//};
void dyson::printSt(ostream& ostr){
	double tau;
	ostr<<scientific;
	for(int i=0; i<Ntau+1; i++){
		tau = i*dtau;
		ostr<<setw(30)<<left<<setprecision(10)<<tau;
		for(int a=0; a<Nflavor; a++) for(int b=0; b<Nflavor; b++) 
			ostr<<setw(30)<<right<<setprecision(10)<<St[i](a,b);
		ostr<<endl;
	}
};
void dyson::printStS0(ostream& ostr){
	double tau;
	ostr<<scientific;
	for(int i=0; i<Ntau+1; i++){
		tau = i*dtau;
		ostr<<setw(30)<<left<<setprecision(10)<<tau
			<<setw(30)<<right<<setprecision(10)<<StS0[i]<<endl;
	}
};
void dyson::printG05beta(std::ostream& ostr, const int& iter){
	ostr<<scientific;
	ostr<<setw(30)<<left<<setprecision(10)<<iter;
	for(int a=0; a<Nflavor; a++) for(int b=0; b<Nflavor; b++){
        ostr<<setw(30)<<right<<setprecision(10)<<Gt[Ntau/2](a,b);
    }
    ostr<<endl;
};
