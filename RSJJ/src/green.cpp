#include <cmath>
#include <vector>
#include <iostream>
#include <iomanip>
#include <sstream>
#include "green.hpp"
#include "parameters.hpp"
using namespace std;
using namespace parameters;


////////////////////////////
// interaction_t function //
////////////////////////////

const double interaction_t::eps = 1.0e-15;

interaction_t::interaction_t(){};
void interaction_t::init(const std::vector<double>& mWt_i){
	NtW = mWt_i.size();
	mWt.resize( NtW );
	dt = parameters::beta/(NtW-1);
	for(int i=0; i<mWt.size(); i++){
		mWt[i] = mWt_i[i];
	}
};
/*void interaction_t::initWt_spinBoson(){
	//cout<<"* initW2 module *"<<endl;
	
	Ot.resize(Ntime+1);

	for(int ti=0; ti<Ntime+1; ti++){
		//cout<<ti<<" ";
		Ot[ti] = -g_coupling*g_coupling*cosh( w_boson*(ti*parameters::dt-0.5*beta) ) / sinh( 0.5*beta*w_boson );
	}

	//cout<<"Otk(0,0) = "<<Otk(0,0)<<endl<<flush;
	//cout<<"Otx(0,0) = "<<Otx(0,0)<<endl<<flush;

	//cout<<"* initW2 module *"<<endl;
};*/
void interaction_t::initWt_spinBoson(const int& NtW_i, const double& gCouple, const double& wBoson){
	NtW = NtW_i;
	dt = parameters::beta/NtW;

	double tau;
	mWt.resize(NtW+1);
	for(int i=0; i<NtW+1; i++){
		tau = i*dt;
		mWt[i] = pow(gCouple,2)*cosh( (tau-0.5*parameters::beta)*wBoson ) / sinh( 0.5*parameters::beta*wBoson );
	}
};
void interaction_t::initWtSingleModeWaveguide(const int& NtW_i, const double& g0, const double& Omega11, const double& omega_c, const int& Nkx){
	NtW = NtW_i;
	dt = parameters::beta/NtW;

	double preFactor = g0*pow(Omega11,4);
	double kc = sqrt(pow(omega_c/Omega11,2)-1.0);
	double betaTilde = parameters::beta*Omega11;
	double dtauTilde = (this->dt)*Omega11;
	double NkxRegular;
	if(Nkx%2==1) NkxRegular = Nkx+1;
	else NkxRegular = Nkx;
	double dkx = kc/NkxRegular;

	double kx;
	double simpsonFactor, intFactor, omega11;
	mWt.resize(Ntau+1);
	for(int it=0; it<Ntau+1; it++){
		mWt[it] = 0.0;
	}
	for(int ikx=0; ikx<NkxRegular+1; ikx++){
		kx = ikx*dkx;
		omega11 = sqrt(kx*kx+1.0);
		intFactor = 1.0/omega11/sinh(0.5*betaTilde*omega11);
		simpsonFactor = getSimpsonFactor(ikx,Nkx);
		for(int it=0; it<Ntau+1; it++){
			mWt[it] += simpsonFactor*intFactor*cosh(omega11*(it*dtauTilde - 0.5*betaTilde));
		}
	}
	for(int it=0; it<Ntau+1; it++){
		mWt[it] *= preFactor*2.0*dkx/3.0;
	}
};

double interaction_t::getWt(const coordinate& xi, const coordinate& xj) const{
	////const double eps =1.e-10;
	//const double eps =1.e-15;
	double tji = xj.time - xi.time;
	int it = static_cast<int>(tji/dt);
	if(it==0){ return mWt[0]; }
	else if(it==NtW-1){ return mWt[it]; }
	else return (((it+1)*dt-tji)*mWt[it] + (tji-it*dt)*mWt[it+1])/dt;
}
void interaction_t::print(ostream& ostr){
	//for(int i=0; i<mWt.size(); i++)
	//	ostr<<mWt[i]<<endl;
	double tau;
	ostr<<scientific;
	for(int i=0; i<NtW+1; i++){
		tau = i*dt;
		ostr<<setw(30)<<left<<setprecision(10)<<tau
			<<setw(30)<<right<<setprecision(10)<<-mWt[i]<<endl;
	}
}


////////////////////////////////////
// pseudoparticle Grren  function //
////////////////////////////////////

const double pseudogreen_t::eps = 1.0e-15;

pseudogreen_t::pseudogreen_t(){
	Eigen::Matrix2d Sz_tmp;
	Sz_tmp << 1.0, 0.0,
	       0.0, -1.0;

	Sz = Sz_tmp;
};
void pseudogreen_t::initOp( const Eigen::MatrixXd& Op_i){
	Sz = Op_i;
};
void pseudogreen_t::init(const std::vector<Eigen::Matrix2d>& Gt_i){
	NtG = Gt_i.size();
	dt = parameters::beta/(NtG-1);

	Gt.resize( NtG );
	for(int i=0; i<Gt.size(); i++){
		Gt[i] = Gt_i[i];
	}
};
//void pseudogreen_t::initGt_atomic(){
//	Eigen::Matrix2d H;
//	H << epsilon, B_field,
//	  B_field, -epsilon;
//
//	Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> eigensolver(H);
//
//	Eigen::Vector2d eval = eigensolver.eigenvalues();
//	Eigen::Matrix2d evec = eigensolver.eigenvectors();
//
//	Eigen::Matrix2d expFactor = Eigen::MatrixXd::Zero(2,2);
//	expFactor(0,0) = exp(-dt*eval(0));
//	expFactor(1,1) = exp(-dt*eval(1));
//	Eigen::Matrix2d propagator = Eigen::MatrixXd::Identity(2,2);
//
//	Ot.resize(Ntime+1);
//	for(int ti=0; ti<Ntime+1; ti++){
//		Ot[ti] = evec*propagator*evec.adjoint();
//		propagator = propagator*expFactor;
//	}
//};
Eigen::Matrix2d pseudogreen_t::getGtSz(const coordinate& xi, const coordinate& xj) const {
	////const double eps =1.e-10;
	//const double eps =1.e-15;
	double tji = xj.time - xi.time;
	int it = static_cast<int>(tji/dt);
	//if(it==0){ return Gt[0]*Sz; }
	//else if(it==NtG-1){ return Gt[it]*Sz; }
	//else return ((((it+1)*dt-tji)*Gt[it] + (tji-it*dt)*Gt[it+1])/dt)*Sz;
	return ((((it+1)*dt-tji)*Gt[it] + (tji-it*dt)*Gt[it+1])/dt)*Sz;

	//double tji = xj.time - xi.time;
	//if(tji<eps){ tji += beta; }
	//int it = static_cast<int>(tji/dt);

	//if(it>Ntau-1){ return Sz*Gt[it]*Sz; }
	//else return Sz*((((it+1)*dt-tji)*Gt[it] + (tji-it*dt)*Gt[it+1])/dt)*Sz;
}
Eigen::Matrix2d pseudogreen_t::getSzGtSz(const coordinate& xi, const coordinate& xj) const {
	////const double eps =1.e-10;
	//const double eps =1.e-15;
	double tji = xj.time - xi.time;
	int it = static_cast<int>(tji/dt);
	//if(it==0){ return Sz*Gt[0]*Sz; }
	//else if(it==NtG-1){ return Sz*Gt[it]*Sz; }
	//else return Sz*( (((it+1)*dt-tji)*Gt[it] + (tji-it*dt)*Gt[it+1])/dt )*Sz;
	return Sz*( (((it+1)*dt-tji)*Gt[it] + (tji-it*dt)*Gt[it+1])/dt )*Sz;

	//double tji = xj.time - xi.time;
	//if(tji<eps){ tji += beta; }
	//int it = static_cast<int>(tji/dt);

	//if(it>Ntau-1){ return Sz*Gt[it]*Sz; }
	//else return Sz*((((it+1)*dt-tji)*Gt[it] + (tji-it*dt)*Gt[it+1])/dt)*Sz;
}
void pseudogreen_t::print(ostream& ostr){
	double tau;
	ostr<<scientific;
	for(int i=0; i<NtG+1; i++){
		tau = i*dt;
		ostr<<setw(30)<<left<<setprecision(10)<<tau
			<<setw(30)<<right<<setprecision(10)<<Gt[i](0,0)
			<<setw(30)<<right<<setprecision(10)<<Gt[i](0,1)
			<<setw(30)<<right<<setprecision(10)<<Gt[i](1,0)
			<<setw(30)<<right<<setprecision(10)<<Gt[i](1,1)<<endl;
	}
}

////////////////////
// 4-point vertex //
////////////////////

vertex_4pt::vertex_4pt(){};
void vertex_4pt::initZero(const int& Nt_i){
	Qt.resize(Nt_i);
	for(int n=0; n<Nt_i; n++){
		Qt[n].resize(n+1);
		for(int m2=0; m2<n+1; m2++){
			Qt[n][m2].resize(m2+1);
			for(int m1=0; m1<m2+1; m1++){
				Qt[n][m2][m1] = Eigen::MatrixXd::Zero(4,4);
				//cout<<"("<<n<<","<<m2<<","<<m1<<")";
			}
		}
	}
	//cout<<endl;
};
void vertex_4pt::accumulate(const vector<int>& iflavor_i, const vector<int>& it_i, const double& value_i){
	//for(int i=0; i<3; i++)
	//	cout<<it_i[i]<<" ";
	//cout<<endl;

	Qt[it_i[0]][it_i[1]][it_i[2]](iflavor_i[0],iflavor_i[1]) += value_i;

	//cout<<endl<<Qt[it_i[0]][it_i[1]][it_i[2]]<<endl;
};
void vertex_4pt::setQt(const int& it, const int& it2, const int& it1, const Eigen::MatrixXd& Qt_tt2t1){
	Qt[it][it2][it1] = Qt_tt2t1;
};
void vertex_4pt::setQt(const int& it, const int& it2, const int& it1, const double* Qt_tt2t1){
	Qt[it][it2][it1] = Eigen::Map<const Eigen::MatrixXd>(Qt_tt2t1,4,4);
};
Eigen::Matrix4d vertex_4pt::getQt(const int& it, const int& it2, const int& it1) const{
	return Qt[it][it2][it1];
};
double vertex_4pt::getQt(const int& it, const int& it2, const int& it1, const int& iout, const int& iin) const{
	return Qt[it][it2][it1](iout,iin);
};
void vertex_4pt::print(std::ostream& ostr,const int& it, const int& iout, const int& iin){
	ostr<<scientific;
	//ostr<<"# end time point: "<<it*dtout<<" #"<<endl;
	for(int it2=0; it2<it+1; it2++){
		for(int it1=0; it1<it2+1; it1++){
			ostr<<setw(30)<<left<<Qt[it][it2][it1](iout,iin);
		}
		ostr<<endl;
	}
};
void vertex_4pt::operator/=(const double& denominator_i){
	for(int it=0; it<Qt.size(); it++){
		for(int it2=0; it2<it+1; it2++){
			for(int it1=0; it1<it2+1; it1++)
				Qt[it][it2][it1] /= denominator_i;
		}
	}
};
void vertex_4pt::operator+=(const vertex_4pt& Qt_i){
	for(int it=0; it<Qt.size(); it++){
		for(int it2=0; it2<it+1; it2++){
			for(int it1=0; it1<it2+1; it1++){
				Qt[it][it2][it1] += Qt_i.getQt(it,it2,it1);
			}
		}
	}
};
vertex_4pt vertex_4pt::operator*(const vertex_4pt& Qt_i){
	vertex_4pt Q_prod;
	Q_prod.initZero(Qt.size());

	for(int it=0; it<Qt.size(); it++){
		for(int it2=0; it2<it+1; it2++){
			for(int it1=0; it1<it2+1; it1++){
				Q_prod.setQt(it,it2,it1,Qt[it][it2][it1].cwiseProduct(Qt_i.Qt[it][it2][it1]));
			}
		}
	}

	return Q_prod;
};


/*
/////////////////////////////////////
// general imaginary-time function //
/////////////////////////////////////

const double function_t::eps = 1.0e-15;

function_t::function_t(){};
void function_t::initZero(){
	Otk.resize(Ntime,Nx);
	Otx.resize(Ntime,Nx);
	for(int i=0; i<Ntime; i++){
		for(int j=0; j<Nx; j++){
			Otk(i,j) = 0.;
			Otx(i,j) = 0.;
		}
	}
	Otx.conservativeResize(Ntime+1,Nx);
	for(int xi=0; xi<Nx; xi++) Otx(Ntime,xi) = Otx(0,xi);
};
Eigen::MatrixXd function_t::getOtx(){ return Otx; };
Eigen::MatrixXd function_t::getOtk(){ return Otk; };
void function_t::setOtx(const Eigen::MatrixXd& Otx_i){
	Otx = Otx_i; 
	Otx.conservativeResize(Ntime+1,Nx);
	for(int xi=0; xi<Nx; xi++) Otx(Ntime,xi) = Otx(0,xi);
};
void function_t::setOtk(const Eigen::MatrixXd& Otk_i){ Otk = Otk_i; };
void function_t::operator+=(const function_t& ft_i){
	Otk += ft_i.Otk;
	Otx += ft_i.Otx;
};
void function_t::printX(ostream& ostr){ ostr<<Otx; };


//////////////////////
// green_t function //
//////////////////////

green_t::green_t(){};
void green_t::initGt0(const lattice& Lattice_i){
	vector<double> ek = Lattice_i.getenergy();
	int Ngrid = ek.size();
	Eigen::MatrixXd ker = Lattice_i.getkernel();

	Otk.resize(Ntime,Ngrid);
	Otx.resize(Ntime,Ngrid);

	double tau, ektmp;
	for(int ti=0; ti<Ntime; ti++){
		tau = ti*parameters::dt;

		for(int ki=0; ki<Ngrid; ki++){
			ektmp = ek[ki];
			Otk(ti,ki) = -exp(-tau*ektmp)/(1.+exp(-beta*ektmp));
		}
	}

	Otx.transpose() = ker*Otk.transpose();

	Otx.conservativeResize(Ntime+1,Ngrid);
	Otx(Ntime,0) =  -1. - Otx(0,0);
	for(int xi=1; xi<Ngrid; xi++) Otx(Ntime,xi) =  -Otx(0,xi);
};
void green_t::initGt_spinBoson(){
	//vector<double> ek = Lattice_i.getenergy();
	//int Ngrid = ek.size();
	//Eigen::MatrixXd ker = Lattice_i.getkernel();

	Otx.resize(Ntime+1,Nx);
	//Otx.resize(Ntime,Ngrid);

	for(int xi=0; xi<Nx; xi++){
		for(int ti=0; ti<Ntime+1; ti++){
			//tau = ti*parameters::dt;
			Otx(ti,xi) = -exp(+beta*epsilon)/( exp(beta*epsilon)+exp(-beta*epsilon) );
		}
	}

	//Otx.transpose() = ker*Otk.transpose();
	//Otx.conservativeResize(Ntime+1,Ngrid);
	//Otx(Ntime,0) =  -1. - Otx(0,0);
	//for(int xi=1; xi<Ngrid; xi++) Otx(Ntime,xi) =  -Otx(0,xi);
};

void green_t::initGtSpinlessAtom(const lattice& Lattice_i){
	vector<double> ek = Lattice_i.getenergy();
	int Ngrid = ek.size();
	Eigen::MatrixXd ker = Lattice_i.getkernel();

	Otk.resize(Ntime,Ngrid);
	Otx.resize(Ntime,Ngrid);

	double tau, mu_tmp;
	for(int ti=0; ti<Ntime; ti++){
		tau = ti*parameters::dt;
		for(int ki=0; ki<Ngrid; ki++)
			Otk(ti,ki) = -exp(-tau*ek[ki])/(1.0 + exp(-beta*ek[ki]));
	}

	Otx.transpose() = ker*Otk.transpose();

	Otx.conservativeResize(Ntime+1,Ngrid);
	Otx(Ntime,0) =  -1. - Otx(0,0);
	for(int xi=1; xi<Ngrid; xi++) Otx(Ntime,xi) =  -Otx(0,xi);
};
void green_t::initGtHubbardAtom(const lattice& Lattice_i){
	vector<double> ek = Lattice_i.getenergy();
	int Ngrid = ek.size();
	Eigen::MatrixXd ker = Lattice_i.getkernel();

	Otk.resize(Ntime,Ngrid);
	Otx.resize(Ntime,Ngrid);

	double tau, ektmp;
	double Uover2 = U/2.;
	for(int ti=0; ti<Ntime; ti++){
		tau = ti*parameters::dt;
		for(int ki=0; ki<Ngrid; ki++)
			Otk(ti,ki) = -(exp((beta-tau)*Uover2)+exp(tau*Uover2))/2./(exp(beta*Uover2)+1.);
	}

	Otx.transpose() = ker*Otk.transpose();

	Otx.conservativeResize(Ntime+1,Ngrid);
	Otx(Ntime,0) =  -1. - Otx(0,0);
	for(int xi=1; xi<Ngrid; xi++) Otx(Ntime,xi) =  -Otx(0,xi);
};

double green_t::getXvalue(const coordinate& xi, const coordinate& xj){
	//const double eps =1.e-10;
	//const double eps =1.e-15;
	double Tsign = 1;
	double tji = xj.time - xi.time;
	if(tji<eps){
		tji += beta;
		Tsign = -1;
	}
	int it = static_cast<int>(tji/dt);
	int ix = Lattice->getXindex(xj.space - xi.space);

	if(it>Ntime-1){ return Tsign*Otx(it,ix); }
	else return Tsign*(((it+1)*dt-tji)*Otx(it,ix) + (tji-it*dt)*Otx(it+1,ix))/dt;
}

////////////////////////////
// interaction_t function //
////////////////////////////

interaction_t::interaction_t(){};
void interaction_t::set_const(const double& W_i){
	Otx = Eigen::MatrixXd::Constant(Ntime+1,Nx,W_i);
}
void interaction_t::initWt_spinBoson(){
	//cout<<"* initW2 module *"<<endl;
	
	Otx.resize(Ntime+1,Nx);

	for(int xi=0; xi<Nx; xi++){
		for(int ti=0; ti<Ntime+1; ti++){
			//cout<<ti<<" ";
			Otx(ti,xi) = -g_coupling*g_coupling*cosh( w_boson*(ti*parameters::dt-0.5*beta) ) / sinh( 0.5*beta*w_boson );
		}
	}

	//cout<<"Otk(0,0) = "<<Otk(0,0)<<endl<<flush;
	//cout<<"Otx(0,0) = "<<Otx(0,0)<<endl<<flush;

	//cout<<"* initW2 module *"<<endl;
};
void interaction_t::initWt0(){
	Otx.resize(Ntime,Nx); 
	Otk.resize(Ntime,Nx); 
	for(int ti=0; ti<Ntime; ti++) for(int xki=0; xki<Nx; xki++){
		Otx(ti,xki) = 0.;
		Otk(ti,xki) = U;
	}
	Otx(0,0) = U;
	Otx.conservativeResize(Ntime+1,Nx);
	for(int xi=0; xi<Nx; xi++) Otx(Ntime,xi) = Otx(0,xi);
};
void interaction_t::initWt1(const lattice& Lattice_i){
	//cout<<"* initW2 module *"<<endl;
	const double eps = 1.e-10;
	vector<Eigen::VectorXd> kgrid = Lattice_i.getkgrid();
	Eigen::MatrixXd ker = Lattice_i.getkernel();

	double tau;
	double ep_kpq, ep_k, exptmp;;

	Otk.resize(Ntime,Nx); 
	Otx.resize(Ntime,Nx);

	for(int ti=0; ti<Ntime; ti++){
		//cout<<ti<<" ";
		tau = ti*parameters::dt;

		for(int qi=0; qi<Nx; qi++){
			Otk(ti,qi) = 0.;
			for(int ki=0; ki<Nx; ki++){
				ep_kpq = Lattice_i.dispersion(kgrid[ki]+kgrid[qi]);
				ep_k = Lattice_i.dispersion(kgrid[ki]);

				if(abs(ep_kpq-ep_k)>eps){
					Otk(ti,qi) -= (nF(ep_kpq)-nF(ep_k)) 
						* exp(-tau*(ep_kpq-ep_k)) / (exp(-beta*(ep_kpq-ep_k)) - 1.);
				}
				else{
					exptmp = exp(beta*ep_k);
					Otk(ti,qi) -= exptmp/(exptmp+1.)/(exptmp+1.);
					//if(ti==0 && qi==0) cout<<exptmp<<endl;
				}
			}
			Otk(ti,qi) *= U*U/Nx;
		}
	}

	Otx.transpose() = ker*Otk.transpose();

	Otx.conservativeResize(Ntime+1,Nx);
	for(int xi=0; xi<Nx; xi++) Otx(Ntime,xi) = Otx(0,xi);

	//cout<<"Otk(0,0) = "<<Otk(0,0)<<endl<<flush;
	//cout<<"Otx(0,0) = "<<Otx(0,0)<<endl<<flush;

	//cout<<"* initW2 module *"<<endl;
};
void interaction_t::initWt2HubbardAtom(){
	//cout<<"* initW2 module *"<<endl;
	Otk.resize(Ntime,Nx); 
	Otx.resize(Ntime,Nx);

	double tau;
	double factor = U*U/2.*tanh(beta*U/4.);
	double constant = U*U/8./cosh(beta*U/4.)/cosh(beta*U/4.);
	for(int xi=0; xi<Nx; xi++){
		for(int ti=0; ti<Ntime; ti++){
			//cout<<ti<<" ";
			tau = ti*parameters::dt;
			Otk(ti,xi) = -(factor*(exp(-tau*U)/2. + cosh(tau*U)/(exp(beta*U)-1.)) + constant);
			Otx(ti,xi) = Otk(ti,xi);
		}
	}

	Otx.conservativeResize(Ntime+1,Nx);
	for(int xi=0; xi<Nx; xi++) Otx(Ntime,xi) = Otx(0,xi);

	//cout<<"Otk(0,0) = "<<Otk(0,0)<<endl<<flush;
	//cout<<"Otx(0,0) = "<<Otx(0,0)<<endl<<flush;

	//cout<<"* initW2 module *"<<endl;
};
void interaction_t::initWt_upup_exactHubbardAtom(){
	//cout<<"* initW2 module *"<<endl;
	Otk.resize(Ntime,Nx); 
	Otx.resize(Ntime,Nx);

	double chi_upup_q0 = beta/4.;
	double W_upup_q0 = -U*U*chi_upup_q0;

	for(int xi=0; xi<Nx; xi++){
		for(int ti=0; ti<Ntime; ti++){
			//cout<<ti<<" ";
			Otk(ti,xi) = W_upup_q0/beta;
			Otx(ti,xi) = Otk(ti,xi);
		}
	}

	Otx.conservativeResize(Ntime+1,Nx);
	for(int xi=0; xi<Nx; xi++) Otx(Ntime,xi) = Otx(0,xi);

	//cout<<"Otk(0,0) = "<<Otk(0,0)<<endl<<flush;
	//cout<<"Otx(0,0) = "<<Otx(0,0)<<endl<<flush;

	//cout<<"* initW2 module *"<<endl;
};
void interaction_t::initWt_updn_exactHubbardAtom(){
	//cout<<"* initW2 module *"<<endl;
	Otk.resize(Ntime,Nx); 
	Otx.resize(Ntime,Nx);

	double chi_updn_q0 = beta*(0.5/(1.+exp(beta*U/2.))-0.25);
	double W_updn_q0 = -U*U*chi_updn_q0;

	for(int xi=0; xi<Nx; xi++){
		for(int ti=0; ti<Ntime; ti++){
			//cout<<ti<<" ";
			Otk(ti,xi) = W_updn_q0/beta;
			Otx(ti,xi) = Otk(ti,xi);
		}
	}

	Otx.conservativeResize(Ntime+1,Nx);
	for(int xi=0; xi<Nx; xi++) Otx(Ntime,xi) = Otx(0,xi);

	//cout<<"Otk(0,0) = "<<Otk(0,0)<<endl<<flush;
	//cout<<"Otx(0,0) = "<<Otx(0,0)<<endl<<flush;

	//cout<<"* initW2 module *"<<endl;
};*/

/*void interaction_t::initWt_upup_rpaHubbardAtom(){
	//cout<<"* initW2 module *"<<endl;
	polarization_w P0;
	P0.init0orderHubbardAtom();

	Otk.resize(Ntime,Nx); 
	Otx.resize(Ntime,Nx);

	double chi_upup_q0 = beta/4.;
	double W_upup_q0 = -U*U*chi_upup_q0;

	for(int xi=0; xi<Nx; xi++){
		for(int ti=0; ti<Ntime; ti++){
			//cout<<ti<<" ";
			Otk(ti,xi) = W_upup_q0/beta;
			Otx(ti,xi) = Otk(ti,xi);
		}
	}

	Otx.conservativeResize(Ntime+1,Nx);
	for(int xi=0; xi<Nx; xi++) Otx(Ntime,xi) = Otx(0,xi);

	//cout<<"Otk(0,0) = "<<Otk(0,0)<<endl<<flush;
	//cout<<"Otx(0,0) = "<<Otx(0,0)<<endl<<flush;

	//cout<<"* initW2 module *"<<endl;
};
void interaction_t::initWt_updn_rpaHubbardAtom(){
	//cout<<"* initW2 module *"<<endl;
	Otk.resize(Ntime,Nx); 
	Otx.resize(Ntime,Nx);

	double chi_updn_q0 = beta*(0.5/(1.+exp(beta*U/2.))-0.25);
	double W_updn_q0 = -U*U*chi_updn_q0;

	for(int xi=0; xi<Nx; xi++){
		for(int ti=0; ti<Ntime; ti++){
			//cout<<ti<<" ";
			Otk(ti,xi) = W_updn_q0/beta;
			Otx(ti,xi) = Otk(ti,xi);
		}
	}

	Otx.conservativeResize(Ntime+1,Nx);
	for(int xi=0; xi<Nx; xi++) Otx(Ntime,xi) = Otx(0,xi);

	//cout<<"Otk(0,0) = "<<Otk(0,0)<<endl<<flush;
	//cout<<"Otx(0,0) = "<<Otx(0,0)<<endl<<flush;

	//cout<<"* initW2 module *"<<endl;
};

double interaction_t::getXvalue(const coordinate& xi, const coordinate& xj){
	////const double eps =1.e-10;
	//const double eps =1.e-15;
	double tji = xj.time - xi.time;
	if(tji<eps){ tji += beta; }
	int it = static_cast<int>(tji/dt);
	int ix = Lattice->getXindex(xj.space - xi.space);

	if(it>Ntime-1){ return Otx(it,ix); }
	else return (((it+1)*dt-tji)*Otx(it,ix) + (tji-it*dt)*Otx(it+1,ix))/dt;
}


////////////////////
// 4-point vertex //
////////////////////

vertex_4pt::vertex_4pt(){
	Qt.resize(Ntime+1);
};
vertex_4pt::~vertex_4pt(){
	delete[] Qt;
};
void vertex_4pt::initZero(){
	for(int it=0; it<Qt.size(); it++){
		Qt[it] = Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic>::Zero(it+1,it+1);
	}
};
void vertex_4pt::accumulate(const int& it, const int& it2, const int& it1, const double& value_i){
	Qt[it](it2,it1) += value_i;
};
void vertex_4pt::setQt(const std::vector<Eigen::MatrixXd>& Qt_i){
	for(int i=0; i<Qt.size(); i++)
		Qt[i] = Qt_i[i];
};
void vertex_4pt::setQt(const int& n, const Eigen::MatrixXd& Qt_n){
	Qt[n] = Qt_n;
};
void vertex_4pt::setQt(const int& n, const double* Qt_i){
	Qt[n] = Eigen::Map<const Eigen::MatrixXd>(Qt_i,n+1,n+1);
};
Eigen::MatrixXd vertex_4pt::getQt(const int& n) const{
	return Qt[n];
};
double vertex_4pt::getQt(const int& n, const int& i, const int& j) const{
	return Qt[n](i,j);
};
void vertex_4pt::print(std::ostream& ostr,const int& it){
	ostr<<"# end time point: "<<it*dtout<<" #"<<endl;
	ostr<<Qt[it]<<endl;
};
void vertex_4pt::operator=(vertex_4pt Q_i){
	for(int it=0; it<Nout+1; it++){
		Qt[it]  = Q_i.Qt[it];
	}
};
void vertex_4pt::operator/=(const double& denominator_i){
	for(int it=0; it<Qt.size(); it++){
		Qt[it] /= denominator_i;
	}
};
vertex_4pt vertex_4pt::operator+(const vertex_4pt& vertex_i){
	vertex_4pt vertex_sum;
	for(int i=0; i<Qt.size(); i++){
		vertex_sum.Qt[i] = Qt[i] + vertex_i.Qt[i];
	}
	return vertex_sum;
};
void vertex_4pt::operator+=(const vertex_4pt& v_i){
	for(int i=0; i<Qt.size(); i++){
		Qt[i] += v_i.getQt(i);
	}
};
vertex_4pt vertex_4pt::operator*(const vertex_4pt& Qt_i){
	vertex_4pt Q_prod;

	for(int it=0; it<Qt.size(); it++)
		Q_prod.setQt(it, Qt[it].cwiseProduct(Qt_i.Qt[it]));

	return Q_prod;
};

*/

/////////////////////////////////////
// general imaginary-freq function //
/////////////////////////////////////

/*function_w::function_w(){};
void function_w::init(){
	Owx.resize(Nw,Nx); 
	Owk.resize(Nw,Nx); 
	for(int i=0; i<Nw; i++) for(int j=0; j<Nx; j++){
		Owx(i,j) = 0.;
		Owk(i,j) = 0.;
	}
};
void function_w::setZero(){
	for(int i=0; i<Owx.rows(); i++) for(int j=0; j<Owx.cols(); j++) Owx(i,j) = 0.;
	for(int i=0; i<Owk.rows(); i++) for(int j=0; j<Owk.cols(); j++) Owk(i,j) = 0.;
};
Eigen::Matrix<std::complex<double>,Eigen::Dynamic,Eigen::Dynamic> function_w::get(){ return Owx; };
Eigen::MatrixXcd function_w::getOwx(){ return Owx; };
Eigen::MatrixXcd function_w::getOwk(){ return Owk; };
complex<double> function_w::getOwx(const int& n_i, const int& x_i){ return Owx(n_i,x_i); };
complex<double> function_w::getOwk(const int& n_i, const int& k_i){ return Owk(n_i,k_i); };
void* function_w::getbufferOwx(){ return &Owx(0,0); };
void* function_w::getbufferOwk(){ return &Owk(0,0); };
void function_w::setOwx(const complex<double>* Owx_i){ Owx = Eigen::Map<const Eigen::MatrixXcd>(Owx_i,Nw,Nx); };
void function_w::setOwk(const complex<double>* Owk_i){ Owk = Eigen::Map<const Eigen::MatrixXcd>(Owk_i,Nw,Nx); };
void function_w::setOwx(const Eigen::MatrixXcd& Owx_i){ Owx = Owx_i; };
void function_w::setOwk(const Eigen::MatrixXcd& Owk_i){ Owk = Owk_i; };
void function_w::setOwx(const int n, const int x, const complex<double> Owx_i){ Owx(n,x) = Owx_i; };
void function_w::setOwk(const int n, const int k, const complex<double> Owk_i){ Owk(n,k) = Owk_i; };
void function_w::accumulateOwx(const int& ix, const Eigen::VectorXcd& Owx_i){ Owx.col(ix) += Owx_i; };
void function_w::operator+= (const function_w& Ow_i){
	Owx += Ow_i.Owx; 
	Owk += Ow_i.Owk; 
};
void function_w::operator-= (const function_w& Ow_i){
	Owx -= Ow_i.Owx; 
	Owk -= Ow_i.Owk; 
}
void function_w::operator/= (const double& denominator_i){ 
	Owx /= denominator_i; 
	Owk /= denominator_i; 
};
function_w function_w::operator* (const double& scalar_i){
	function_w tmp;
	tmp.Owx = scalar_i*Owx;
	tmp.Owk = scalar_i*Owk;
       	
	return tmp;
};
function_w function_w::operator* (const function_w& Ow_i){
	function_w tmp;
	tmp.setOwk( (Owk.array()*Ow_i.Owk.array()).matrix() );
       	return tmp; 
};
void function_w::operator*= (const function_w& Ow_i){
	setOwk( (Owk.array()*Ow_i.Owk.array()).matrix() );
};
void function_w::fourier_xtok(){
	Owk.resize(Nw,Nx);

	Eigen::MatrixXd ker = Lattice.getkernel();
	for(int n=0; n<Nw; n++) Owk.row(n) = Nx*Owx.row(n)*ker;
};
void function_w::fourier_ktox(){
	Owx.resize(Nw,Nx);

	Eigen::MatrixXd ker = Lattice.getkernel();
	Owx.transpose() = ker*Owk.transpose();
};
void function_w::print_x(std::ostream& ostr){
	ostr<<Owx<<std::endl;
};
void function_w::print_k(std::ostream& ostr){
	ostr<<Owk<<std::endl;
};
void function_w::printOwk(std::ostream& ostr){
	stringstream ss;
	ostr<<setw(25)<<left<<"# w_n";
	for(int ik=0; ik<Nx; ik++){
		ss.str(""); ss<<"("<<Lattice.getkgrid(ik)(0)<<","<<Lattice.getkgrid(ik)(1)<<")";
		ostr<<setw(50)<<right<<ss.str();
	}
	ostr<<endl;

	ostr<<setprecision(15)<<scientific;
	for(int n=0; n<Nw; n++){
		ostr<<setw(25)<<left<<w(n);
		for(int ik=0; ik<Nx; ik++){
			ostr<<setw(25)<<right<<Owk(n,ik).real()
				<<setw(25)<<right<<Owk(n,ik).imag();
		}
		ostr<<endl;
	}

	//for(int n=0; n<Owk.rows(); n++){
	//	ostr<<setw(15)<<left<<w(n);
	//	for(int ix=0; ix<Owk.cols(); ix++){
	//		ostr<<setw(15)<<right<<Owk(n,ix).real()
	//			<<setw(15)<<right<<Owk(n,ix).imag();
	//	}
	//	ostr<<endl;
	//}
};
void function_w::printOwx(std::ostream& ostr){
	stringstream ss;
	ostr<<setw(25)<<left<<"# w_n";
	for(int ik=0; ik<Nx; ik++){
		ss.str(""); ss<<"("<<Lattice.getkgrid(ik)(0)<<","<<Lattice.getkgrid(ik)(1)<<")";
		ostr<<setw(50)<<right<<ss.str();
	}
	ostr<<endl;

	ostr<<setprecision(15)<<scientific;
	for(int n=0; n<Nw; n++){
		ostr<<setw(25)<<left<<w(n);
		for(int ix=0; ix<Nx; ix++){
			ostr<<setw(25)<<right<<Owx(n,ix).real()
				<<setw(25)<<right<<Owx(n,ix).imag();
		}
		ostr<<endl;
	}

	//for(int n=0; n<Owk.rows(); n++){
	//	ostr<<setw(15)<<left<<w(n);
	//	for(int ix=0; ix<Owk.cols(); ix++){
	//		ostr<<setw(15)<<right<<Owk(n,ix).real()
	//			<<setw(15)<<right<<Owk(n,ix).imag();
	//	}
	//	ostr<<endl;
	//}
};


/////////////
// green_w //
/////////////

green_w::green_w(){};
void green_w::initGw0(const lattice& Lattice_i){
	vector<double> ek = Lattice_i.getenergy();
	int Ngrid = ek.size();
	Eigen::MatrixXd ker = Lattice_i.getkernel();

	Owk.resize(Nw,Nx);
	Owx.resize(Nw,Nx);

	double wn;
	for(int n=0; n<Nw; n++){
		wn = w(n);
		for(int ki=0; ki<Ngrid; ki++)
			Owk(n,ki) = 1./complex<double>(-ek[ki],w(n));
	}
	Owx.transpose() = ker*Owk.transpose();
};
void green_w::initGwHubbardAtomexact(){
	Owk.resize(Nw,1);
	Owx.resize(Nw,1);

	double wn;
	for(int n=0; n<Nw; n++){
		wn = w(n);
		Owk(n,0) = complex<double>(0.,-4.*wn)/(4.*wn*wn+U*U);
	}
	Owx.transpose() = Owk.transpose();
};


///////////////////
// interaction_w //
///////////////////

interaction_w::interaction_w(){};
void interaction_w::initW1(){
	Owx.resize(Nw,Nx); 
	Owk.resize(Nw,Nx);
	for(int n=0; n<Nw; n++){
		for(int ki=0; ki<Nx; ki++){
			Owk(n,ki) = U;
			Owx(n,ki) = 0.;
		}
		Owx(n,0) = U;
	}
};
void interaction_w::initW2(){
	Owx.resize(Nw,Nx); 
	Owk.resize(Nw,Nx);
	for(int ki=0; ki<Nx; ki++){
		for(int n=0; n<Nw; n++){
			Owk(n,ki) = -0.5*tanh(beta*U/4.)*U*U*U/(wB(n)*wB(n)+U*U);
			Owx(n,ki) = 0.;
		}
	}
	for(int ki=0; ki<Nx; ki++){
		Owk(0,ki) -= U*U*beta/8./cosh(beta*U/4.)/cosh(beta*U/4.);
	}

};

//////////////////
// selfenergy_w //
//////////////////

selfenergy_w::selfenergy_w(){};
void selfenergy_w::init2order(const lattice& Lattice_i){
	cout<<"* 2nd order analytic calculation *"<<endl;
	const double eps = 1.e-10;
	vector<Eigen::VectorXd> kgrid = Lattice_i.getkgrid();

	double ep_kpq, ep_kppq, ep_kp, exptmp;

	Owk.resize(Nw,Nx); 
	for(int n=0; n<Nw; n++) for(int ki=0; ki<Nx; ki++){
		Owk(n,ki) = 0.;
		for(int qi=0; qi<Nx; qi++) for(int kj=0; kj<Nx; kj++){
			ep_kpq = Lattice_i.dispersion(kgrid[ki]+kgrid[qi]);
			ep_kppq = Lattice_i.dispersion(kgrid[kj]+kgrid[qi]);
			ep_kp = Lattice_i.dispersion(kgrid[kj]);

			if(abs(ep_kppq-ep_kp)>eps){
				Owk(n,ki) -= (nF(ep_kpq)+nB(ep_kppq-ep_kp))
					  *  (nF(ep_kppq)-nF(ep_kp))
					  /  ( complex<double>(0.,w(n)) - (ep_kpq-ep_kppq+ep_kp));
			}
			else{
				exptmp = exp(beta*ep_kp);
				Owk(n,ki) += exptmp/(exptmp+1)/(exptmp+1)
					  /  (complex<double>(0.,w(n)) - ep_kpq);
			}
		}
		Owk(n,ki) *= U*U/Nx/Nx;
	}
	cout<<"**********************************"<<endl;
};*/
/*void sigma::set(const std::list<line>::iterator& measuringline_i){
	//cout<<measuringline_i->getxj()<<endl<<flush;
	coordinate xi( measuringline_i->getxj() );
	coordinate xj( measuringline_i->getxi() );

	double tji = xj.time - xi.time;
	Eigen::VectorXi xji = xj.space - xi.space;
	int Lgrid = Lattice.getLgrid();
	int dim = Lattice.getdim();

	int ix = 0;
	for(int i=0; i<xji.size(); i++){
		xji(i) = (xji(i) + Lgrid)%Lgrid;
		ix += xji(i)*pow(Lgrid,dim-i-1);
	}

	double wt;
	for(int n=0; n<Swx.rows(); n++){
		for(int jx=0; jx<Swx.cols(); jx++) Swx(n,jx) = 0.;
		wt = w(n)*tji;
		Swx(n,ix) = complex<double>(cos(wt),sin(wt));
	}
};
*/

////////////////////
// polarization_w //
////////////////////

/*polarization_w::polarization_w(){};
void polarization_w::init0orderHubbardAtom(){
	Owk.resize(Nw,1);
	Owx.resize(Nw,1);
	Owk(0,0) = beta/4.;
	Owx(0,0) = beta/4.;
	for(int n=1; n<Nw; n++){
		Owk(n,0) = 0.;
		Owx(n,0) = 0.;
	}
};*/
/*void polarization_w::init0orderHubbardAtom(){
	Owk.resize(Nw,1);
	Owx.resize(Nw,1);
	double factor = U*tanh(beta*U/4.)/2.;
	for(int n=0; n<Nw; n++){
		Owk(n,0) = factor/(wB(n)*wB(n)+U*U);
		Owx(n,0) = factor/(wB(n)*wB(n)+U*U);
	}
	Owk(0,0) += beta/cosh(beta*U/4.)/cosh(beta*U/4.)/8.;
	Owx(0,0) += beta/cosh(beta*U/4.)/cosh(beta*U/4.)/8.;
};*/
