#include "controller.hpp"
#include <mpi.h>
#include <complex>
#include <fftw3.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include "parameters.hpp"
#include "diagram.hpp"
using namespace std;
using namespace parameters;

/////////////////////////////
// miscellaneous functions //
/////////////////////////////

bool all_of_smaller_than(const vector<int>& v, const int& m){
	for(auto& it : v)
		if(it>m) return false;
	return true;
};
int factorial(const int& i){
	if(i>1) return i*factorial(i-1);
	else if(i==1) return 1;
	else exit(EXIT_FAILURE);
};
int count_combination(const vector<int>& v){
	int num = v.size();
	vector<int> den;
	int tmp;
	tmp = 1;
	for(int i=1; i<v.size(); i++){
		if(v[i-1]==v[i]){
			tmp++;
			if(i==v.size()-1) den.push_back(tmp);
		}
		else{
			den.push_back(tmp);
			tmp = 1;
		}
	}
	int out = factorial(num);
	for(auto& it : den) out /= factorial(it);

	return out;
};


///////////////////////////////
// controller implementation //
///////////////////////////////

controller::controller(){
	//cout<<"*** controller constructor ***"<<endl<<flush;
	for(int sp=0; sp<2; sp++){
		normalization[sp].resize(Nordermax+1); 
		dnormalization[sp].resize(Nordermax+1);
	}
	Ploc.resize(Nordermax); dPloc.resize(Nordermax);
	Sw = new selfenergy_w**[Nordermax];	dSw = new selfenergy_w**[Nordermax];
	Pw = new polarization_w**[Nordermax];	dPw = new polarization_w**[Nordermax];
	sw = new selfenergy_w**[Nordermax];	dsw = new selfenergy_w**[Nordermax];
	pw = new polarization_w**[Nordermax];	dpw = new polarization_w**[Nordermax];

	for(int iorder=0; iorder<Nordermax; iorder++){
		Sw[iorder] = new selfenergy_w*[Nflavor];
		dSw[iorder] = new selfenergy_w*[Nflavor];
		sw[iorder] = new selfenergy_w*[Nflavor];
		dsw[iorder] = new selfenergy_w*[Nflavor];
		for(int iflavor=0; iflavor<Nflavor; iflavor++){
			Sw[iorder][iflavor] = new selfenergy_w[Nflavor];
			dSw[iorder][iflavor] = new selfenergy_w[Nflavor];
			sw[iorder][iflavor] = new selfenergy_w[Nflavor];
			dsw[iorder][iflavor] = new selfenergy_w[Nflavor];
		}

		Pw[iorder] = new polarization_w*[Nflavor];
		dPw[iorder] = new polarization_w*[Nflavor];
		pw[iorder] = new polarization_w*[Nflavor];
		dpw[iorder] = new polarization_w*[Nflavor];
		for(int iflavor=0; iflavor<Nflavor; iflavor++){
			Pw[iorder][iflavor] = new polarization_w[Nflavor];
			dPw[iorder][iflavor] = new polarization_w[Nflavor];
			pw[iorder][iflavor] = new polarization_w[Nflavor];
			dpw[iorder][iflavor] = new polarization_w[Nflavor];
		}
	}

	Uw.initW1();
	//Uw.initW1(); W2w.initW2();
	//ofstream ofstr("W2w.dat"); W2w.printOwk(ofstr); ofstr.close();

	SwSum = new selfenergy_w*[Nflavor];	dSwSum = new selfenergy_w*[Nflavor];
	PwSum = new polarization_w*[Nflavor];	dPwSum = new polarization_w*[Nflavor];
	for(int iflavor=0; iflavor<Nflavor; iflavor++){
		SwSum[iflavor] = new selfenergy_w[Nflavor];
		SwSum[iflavor][iflavor].init();
		SwSum[iflavor][(iflavor^1)].init();

		dSwSum[iflavor] = new selfenergy_w[Nflavor];
		dSwSum[iflavor][iflavor].init();
		dSwSum[iflavor][(iflavor^1)].init();

		PwSum[iflavor] = new polarization_w[Nflavor];
		PwSum[iflavor][iflavor].init0orderHubbardAtom();
		PwSum[iflavor][(iflavor^1)].init();

		dPwSum[iflavor] = new polarization_w[Nflavor];
		dPwSum[iflavor][iflavor].init();
		dPwSum[iflavor][(iflavor^1)].init();
	}

	WwSum = new interaction_w*[Nflavor];
	WwSum_prev = new interaction_w*[Nflavor];
	for(int iflavor=0; iflavor<Nflavor; iflavor++){
		WwSum[iflavor] = new interaction_w[Nflavor];
		WwSum_prev[iflavor] = new interaction_w[Nflavor];
		for(int jflavor=0; jflavor<Nflavor; jflavor++){
			WwSum[iflavor][jflavor].init();
			WwSum_prev[iflavor][jflavor].init();
		}
	}

	/*dysonequation_W();
	WwSum = WwSum - Uw;

	for(int iflavor=0; iflavor<Nflavor; iflavor++) for(int jflavor=0; jflavor<Nflavor; jflavor++){
		if(iflavor==jflavor){
			Wt[0][iflavor][jflavor].InverseTimeFourier(WwSum[iflavor][jflavor].getOwk());
			Gt[0][iflavor][jflavor].initGtHubbardAtom(Lattice);
		}
		else{
			Wt[0][iflavor][jflavor].InverseTimeFourier(WwSum[iflavor][jflavor].getOwk());
			Gt[0][iflavor][jflavor].initZero();
		}
	}*/

	//cout<<"iterative scheme"<<endl;
	dysonequation_W();

	//interaction_t W2t;
	//W2t.initWt2HubbardAtom();
	for(int iflavor=0; iflavor<Nflavor; iflavor++) for(int jflavor=0; jflavor<Nflavor; jflavor++){
		WwSum_prev[iflavor][jflavor] = WwSum[iflavor][jflavor];

		if(iflavor==jflavor){
			//WwSum[iflavor][iflavor] -= W2w;
			Gt[0][iflavor][iflavor].initGt0(Lattice);
		}
		else{
			WwSum[iflavor][jflavor] -= Uw;
			Gt[0][iflavor][jflavor].initZero();
		}

		WwSum[iflavor][jflavor].fourier_ktox();
		Wt[0][iflavor][jflavor].initZero();
		Wt[0][iflavor][jflavor].setOtx(InverseTimeFourier(WwSum[iflavor][jflavor].getOwx()));

		//if(iflavor==jflavor) Wt[0][iflavor][jflavor] += W2t;
	}

	/*selfenergy_w Sw1;
	Sw1.setOwx( construct_Sw1(Gt[0][0][0].getOtx(), Wt[0][0][0].getOtx()) );

	for(int iflavor=0; iflavor<Nflavor; iflavor++) for(int jflavor=0; jflavor<Nflavor; jflavor++){
		WwSum_prev[iflavor][jflavor] = WwSum[iflavor][jflavor];

		if(iflavor==jflavor){
			Wt[0][iflavor][jflavor].initWt_upup_exactHubbardAtom();
			Gt[0][iflavor][iflavor].initGt0(Lattice);
			//Gt[0][iflavor][iflavor].initGtHubbardAtom(Lattice);
		}
		else{
			Wt[0][iflavor][jflavor].initWt_updn_exactHubbardAtom();
			Gt[0][iflavor][jflavor].initZero();
		}
	}*/


	/*normalization.resize(Nordermax);
	dnormalization.resize(Nordermax);
	normalization[0] = U*n0;
	dnormalization[0] = 0.0;

	init_order_combination();

	Sw = new selfenergy_w**[Nordermax];
	dSw = new selfenergy_w**[Nordermax];

	Gw = new green_w**[Nordermax];
	for(int iorder=0; iorder<Nordermax; iorder++){
		Sw[iorder] = new selfenergy_w*[Nflavor];
		dSw[iorder] = new selfenergy_w*[Nflavor];

		Gw[iorder] = new green_w*[Nflavor];
		for(int iflavor=0; iflavor<Nflavor; iflavor++){
			Sw[iorder][iflavor] = new selfenergy_w[Nflavor];
			dSw[iorder][iflavor] = new selfenergy_w[Nflavor];

			Gw[iorder][iflavor] = new green_w[Nflavor];
		}
	}

	for(int iflavor=0; iflavor<Nflavor; iflavor++)
		Gw[0][iflavor][iflavor].initGw0(Lattice);*/

	printinput(global.getdir());	

	/*ofstream ofstr; stringstream ss;
	ss.str(""); ss<<global.getdir()<<"/Swx_o1.exact";
	ofstr.open(ss.str().c_str());
	Sw1.printOwx(ofstr); 
	ofstr.close();*/

	//cout<<"*** controller constructor completed ***"<<endl<<flush;
};
void controller::init_order_combination(){
	order_combination.resize(Nordermax);

	vector<int> head(1);
	vector<int> body;

	for(int n=2; n<Nordermax+1; n++){
		head[0] = n;
		order_combination[n-1].push_back(head);
		for(int l=1; l<n-1; l++){
			head[0] = n-l;
			for(int m=0; m<order_combination[l-1].size(); m++){
				if( all_of_smaller_than(order_combination[l-1][m],n-l) ){
					body.resize(order_combination[l-1][m].size()+1);
					body[0] = head[0];
					for(int i=1; i<body.size(); i++)
						body[i] = order_combination[l-1][m][i-1];
					order_combination[n-1].push_back(body);
				}
			}
		}

		//for(auto& it0 : order_combination[n-1]){
		//	for(auto& it1 : it0) cout<<it1;
		//	cout<<"("<<count_combination(it0)<<") ";
		//}
		//cout<<endl;
	}
};
void controller::debug(){
	printinput(global.getdir());

	diagram Diag;
	Diag.debug();
	//Diag.normalize_S(0,0,U*n0);
	//Diag.normalize_P(0,0,U*n0);

	//Diag.printoutput_S(global.getdir());
	//Diag.printoutput(global.getdir());
};
void controller::construct_W2(){
	WwSum_prev = new interaction_w*[Nflavor];
	for(int iflavor=0; iflavor<Nflavor; iflavor++){
		WwSum_prev[iflavor] = new interaction_w[Nflavor];
		WwSum_prev[iflavor][iflavor].initW2();
		WwSum_prev[iflavor][(iflavor^1)].init();
	}

	/*Eigen::MatrixXcd I = Eigen::MatrixXcd::Identity(Nflavor,Nflavor);
	Eigen::MatrixXcd W0 = Eigen::MatrixXcd::Zero(Nflavor,Nflavor);
	Eigen::MatrixXcd W; W0(0,1) = -U; W0(1,0) = -U;
	Eigen::MatrixXcd Pi(Nflavor,Nflavor);
	for(int n=0; n<Nw; n++) for(int ik=0; ik<Nx; ik++) {
		for(int iflavor=0; iflavor<Nflavor; iflavor++) for(int jflavor=0; jflavor<Nflavor; jflavor++)
			Pi(iflavor,jflavor) = PwSum[iflavor][jflavor].getOwk(n,ik);

		//W = (I - W0*Pi).inverse()*W0;
		W = W0 + W0*Pi*W0;

		for(int iflavor=0; iflavor<Nflavor; iflavor++) for(int jflavor=0; jflavor<Nflavor; jflavor++)
			WwSum_prev[iflavor][jflavor].setOwk(n,ik,W(iflavor,jflavor));
	}*/
};
Eigen::MatrixXcd
controller::construct_Sw1(const Eigen::MatrixXd& Gt, const Eigen::MatrixXd& Wt){
	//cout<<"* construct Sw1 *"<<endl<<flush;
	Eigen::MatrixXcd Swx1(Ntime,Nx);

	//cout<<"dt = "<<dt<<endl<<flush;
	//cout<<"Nx = "<<Nx<<endl;
	//cout<<"Nw = "<<Nw<<endl;
	//cout<<"Gt: ("<<Gt.rows()<<","<<Gt.cols()<<")"<<endl<<flush;
	//cout<<"Wt: ("<<Wt.rows()<<","<<Wt.cols()<<")"<<endl<<flush;

	fftw_complex* in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*(Ntime));
	fftw_complex* out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*(Ntime));

	fftw_plan p = fftw_plan_dft_1d(Ntime,in,out,FFTW_BACKWARD,FFTW_ESTIMATE);

	for(int x=0; x<Nx; x++){
		for(int i=0; i<Ntime; i++){
			in[i][0] = -Gt(i,x)*Wt(i,x)*cos(M_PI*i/Ntime);
			in[i][1] = -Gt(i,x)*Wt(i,x)*sin(M_PI*i/Ntime);
			//cout<<-Gt(i,x)*Wt(i,x)*cos(M_PI*i*dt)<<"    "<<-Gt(i,x)*Wt(i,x)*sin(M_PI*i*dt)<<endl<<flush;
		}

		fftw_execute(p);

		for(int i=0; i<Ntime; i++) Swx1(i,x) = dt*complex<double>(out[i][0],out[i][1]);
	}


	fftw_destroy_plan(p);
	fftw_free(in); fftw_free(out);

	return Swx1;
};
void controller::dysonequation_S(){
	//Eigen::MatrixXcd Gwtmp(Nw,Nx);
	//complex<double> G0wtmp, Swtmp;
	Eigen::MatrixXcd Gwtildetmp(Nw,Nx);
	complex<double> G0wtmp, Swtmp;
	for(int n=0; n<Nw; n++) for(int ik=0; ik<Nx; ik++){
		G0wtmp = 0.0; Swtmp = 0.0;
		// spin symmetrization
		for(int iflavor=0; iflavor<Nflavor; iflavor++){
			G0wtmp += Gw[0][iflavor][iflavor].getOwk(n,ik);
			for(int iorder=0; iorder<Nordermax; iorder++){
				Swtmp += Sw[iorder][iflavor][iflavor].getOwx(n,ik);
			}
		}
		G0wtmp /= 2.0; Swtmp /= 2.0;

		//Gwtmp(n,ik) = G0wtmp/(1.0 - G0wtmp*Swtmp);
		Gwtildetmp(n,ik) = G0wtmp/(1.0 - G0wtmp*Swtmp);
	}
	/*for(int iflavor=0; iflavor<Nflavor; iflavor++){
		Gw[Nordermax-1][iflavor][iflavor].setOwx(Gwtmp);
		Gw[Nordermax-1][iflavor][iflavor].setOwk(Gwtmp);
	}*/

	for(int iflavor=0; iflavor<Nflavor; iflavor++){
		GwSum[iflavor][iflavor].setOwk(Gwtildetmp);
	}

	/*green_w Gwtilde;
	Gwtilde.setOwx(Gwtildetmp);
	Gwtilde.setOwk(Gwtildetmp);
	Gw_memory.push_back(Gwtilde);*/
};
void controller::dysonequation_W(){
	Eigen::MatrixXcd I = Eigen::MatrixXcd::Identity(Nflavor,Nflavor);
	Eigen::MatrixXcd W0 = Eigen::MatrixXcd::Zero(Nflavor,Nflavor);
	Eigen::MatrixXcd W; W0(0,1) = U; W0(1,0) = U;
	Eigen::MatrixXcd Pi(Nflavor,Nflavor);
	for(int n=0; n<Nw; n++) for(int ik=0; ik<Nx; ik++) {
		for(int iflavor=0; iflavor<Nflavor; iflavor++) for(int jflavor=0; jflavor<Nflavor; jflavor++)
			Pi(iflavor,jflavor) = PwSum[iflavor][jflavor].getOwk(n,ik);

		W = (I + W0*Pi).inverse()*W0;

		//check for limit case
		//W = W0 - W0*Pi*W0; 

		for(int iflavor=0; iflavor<Nflavor; iflavor++) for(int jflavor=0; jflavor<Nflavor; jflavor++)
			WwSum[iflavor][jflavor].setOwk(n,ik,W(iflavor,jflavor));
	}
};
Eigen::MatrixXd controller::BosonicInverseTimeFourier(const Eigen::MatrixXcd& Owk){
        Eigen::MatrixXd Otk(Ntime,Nx);

        fftw_complex* in; double* out;
        in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*(Ntime/2+1));
        out = (double*) fftw_malloc(sizeof(double)*Ntime);

        //double factor = 0.5*U*U*U*tanh(beta*U/4.);
        fftw_plan p = fftw_plan_dft_c2r_1d(Ntime,in,out,FFTW_ESTIMATE);
        for(int k=0; k<Nx; k++){
                for(int i=0; i<Nw; i++){ in[i][0] = Owk(i,k).real();  in[i][1] = -Owk(i,k).imag();}
                for(int i=Nw; i<(Ntime/2+1); i++){ in[i][0] = 0.;  in[i][1] = 0.; }
                //for(int i=Nw; i<(Ntime/2+1); i++){ in[i][0] = factor/(wB(i)*wB(i)+U*U);  in[i][1] = 0.; }

                fftw_execute(p);

                for(int i=0; i<Ntime; i++) Otk(i,k) = out[i]/beta;
        }

        fftw_destroy_plan(p);
        fftw_free(in); fftw_free(out);

        return Otk;

        //int Nfreq = freq.size();
        //int Ntau = 2*Nfreq;
        //WrapVecDoub data(2*Ntau);
        //std::vector<complex<double> > time(Ntau);
        //for(int i=0; i<Nfreq; i++){
        //        data[i] = Freq[i];
        //        data[-i-1] = conj(Freq[i]);
        //}

        //four1(data,-1);

        //complex<double> temp    = 1.;
        //complex<double> exp(cos(M_PI/Ntau),-sin(M_PI/Ntau));
        //for(int i=0; i<Ntau; i++){
        //        Time[i] = complex<double>((temp*data[i]).real()/beta,0.);
        //        temp    *= exp;
        //}

        //return Time;
};

Eigen::MatrixXd controller::InverseTimeFourier(const Eigen::MatrixXcd& Owk){
	Eigen::MatrixXd Otk(Ntime,Nx);

	fftw_complex* in; double* out;
	in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*(Ntime/2+1));
	out = (double*) fftw_malloc(sizeof(double)*Ntime);

	//double factor = 0.5*U*U*U*tanh(beta*U/4.);
	fftw_plan p = fftw_plan_dft_c2r_1d(Ntime,in,out,FFTW_ESTIMATE);
	for(int k=0; k<Nx; k++){
		for(int i=0; i<Nw; i++){ in[i][0] = Owk(i,k).real();  in[i][1] = -Owk(i,k).imag();}
		for(int i=Nw; i<(Ntime/2+1); i++){ in[i][0] = 0.;  in[i][1] = 0.; }
		//for(int i=Nw; i<(Ntime/2+1); i++){ in[i][0] = factor/(wB(i)*wB(i)+U*U);  in[i][1] = 0.; }
	
		fftw_execute(p);

		for(int i=0; i<Ntime; i++) Otk(i,k) = out[i]/beta;
	}


	fftw_destroy_plan(p);
	fftw_free(in); fftw_free(out);
	
	return Otk;

        //int Nfreq = freq.size();         
        //int Ntau = 2*Nfreq;             
        //WrapVecDoub data(2*Ntau);
        //std::vector<complex<double> > time(Ntau);
        //for(int i=0; i<Nfreq; i++){            
        //        data[i] = Freq[i];
        //        data[-i-1] = conj(Freq[i]);
        //}

        //four1(data,-1);

        //complex<double> temp    = 1.;  
        //complex<double> exp(cos(M_PI/Ntau),-sin(M_PI/Ntau));
        //for(int i=0; i<Ntau; i++){             
        //        Time[i] = complex<double>((temp*data[i]).real()/beta,0.);
        //        temp    *= exp;
        //}

        //return Time;
};
void controller::iterate_W(const double& alpha_i){
	const int root = 0;
	const int rank = MPI::COMM_WORLD.Get_rank();
	if(rank==root) cout<<"iterate_W module"<<endl<<flush;

	Eigen::MatrixXcd I = Eigen::MatrixXcd::Identity(Nflavor,Nflavor);
	Eigen::MatrixXcd W0 = Eigen::MatrixXcd::Zero(Nflavor,Nflavor);
	W0(0,1) = U; W0(1,0) = U;
	Eigen::MatrixXcd W(Nflavor,Nflavor);
	Eigen::MatrixXcd Pi(Nflavor,Nflavor);
	for(int n=0; n<Nw; n++) for(int ik=0; ik<Nx; ik++) {
		//cout<<"("<<n<<","<<ik<<")"<<flush;
		for(int iflavor=0; iflavor<Nflavor; iflavor++) for(int jflavor=0; jflavor<Nflavor; jflavor++){
			Pi(iflavor,jflavor) = PwSum[iflavor][jflavor].getOwk(n,ik);
			W(iflavor,jflavor) = WwSum_prev[iflavor][jflavor].getOwk(n,ik);
		}

		if(rank==root && (n==0 && ik==0)) cout<<alpha_i*((I + W0*Pi)*W - W0).eval()<<endl;
		W -= alpha_i*((I + W0*Pi)*W - W0).eval();

		for(int iflavor=0; iflavor<Nflavor; iflavor++) for(int jflavor=0; jflavor<Nflavor; jflavor++)
			WwSum[iflavor][jflavor].setOwk(n,ik,W(iflavor,jflavor));
	}
};
//void controller::flavorsymmetrize(){ };
void controller::bolditeration(){
	const int root = 0;
	const int rank = MPI::COMM_WORLD.Get_rank();
	if(rank==root) cout<<"*** bold iteration ***"<<endl<<flush;

	//const double t = 0.5;

	//double alpha;
	//double alpha_numerator;
	//double alpha_denomination;

	if(rank==root) printinput(global.getdir());	

	//alpha_denomination = 0.;
	for(int iiteration=1; iiteration<Niteration+1; iiteration++){
		if(rank==root) cout<<" ** "<<iiteration<<"th iteration **"<<endl<<flush;
		if(rank==root) cout<<"  * sampling *"<<endl<<flush;
		MPI::COMM_WORLD.Barrier();
		sampling();
		//flavorsymmetrize();

		/*if(rank==root) cout<<"  * stepsize (alpha) adjust *"<<endl<<flush;
		alpha_numerator = pow(static_cast<double>(iiteration),t);
		alpha_denomination += alpha_numerator;
		alpha = alpha_numerator/alpha_denomination;
		if(rank==root) cout<<"    alpha : "<<alpha<<endl;

		if(rank==root) cout<<"  * iterate Ww *"<<endl<<flush;
		iterate_W(alpha);

		if(rank==root) cout<<"  * fourier transform into Wt *"<<endl<<flush;
		interaction_t W2t;
		W2t.initWt2HubbardAtom();
		for(int iflavor=0; iflavor<Nflavor; iflavor++) for(int jflavor=0; jflavor<Nflavor; jflavor++){
			WwSum_prev[iflavor][jflavor] = WwSum[iflavor][jflavor];

			if(iflavor==jflavor){ WwSum[iflavor][iflavor] -= W2w; }
			else{ WwSum[iflavor][jflavor] -= Uw; }

			WwSum[iflavor][jflavor].fourier_ktox();
			Wt[0][iflavor][jflavor].initZero();
			Wt[0][iflavor][jflavor].setOtx(InverseTimeFourier(WwSum[iflavor][jflavor].getOwx()));

			if(iflavor==jflavor) Wt[0][iflavor][jflavor] += W2t;
		}*/

		//printiterativeoutput(global.getdir(),iiteration);
		//printoutput(global.getdir());
		if(rank==root) cout<<endl<<endl;
	}
	if(rank==root) printoutput(global.getdir());
};
void controller::bolditeration_S(){
	const int root = 0;
	const int rank = MPI::COMM_WORLD.Get_rank();
	if(rank==root) cout<<"*** bold iteration ***"<<endl<<flush;

	const double t = 0.5;

	double alpha;
	double alpha_numerator;
	double alpha_denomination;

	if(rank==root) printinput(global.getdir());	

	alpha_denomination = 0.;
	for(int iiteration=1; iiteration<Niteration+1; iiteration++){
		if(rank==root) cout<<" ** "<<iiteration<<"th iteration **"<<endl<<flush;
		if(rank==root) cout<<"  * sampling *"<<endl<<flush;
		MPI::COMM_WORLD.Barrier();
		sampling_S();
		//flavorsymmetrize();

		if(rank==root) cout<<"  * stepsize (alpha) adjust *"<<endl<<flush;
		alpha_numerator = pow(static_cast<double>(iiteration),t);
		alpha_denomination += alpha_numerator;
		alpha = alpha_numerator/alpha_denomination;
		if(rank==root) cout<<"    alpha : "<<alpha<<endl;

		if(rank==root) cout<<"  * iterate Ww *"<<endl<<flush;
		iterate_W(alpha);

		if(rank==root) cout<<"  * fourier transform into Wt *"<<endl<<flush;
		interaction_t W2t;
		W2t.initWt2HubbardAtom();
		for(int iflavor=0; iflavor<Nflavor; iflavor++) for(int jflavor=0; jflavor<Nflavor; jflavor++){
			WwSum_prev[iflavor][jflavor] = WwSum[iflavor][jflavor];

			if(iflavor==jflavor){ WwSum[iflavor][iflavor] -= W2w; }
			else{ WwSum[iflavor][jflavor] -= Uw; }

			WwSum[iflavor][jflavor].fourier_ktox();
			Wt[0][iflavor][jflavor].initZero();
			Wt[0][iflavor][jflavor].setOtx(InverseTimeFourier(WwSum[iflavor][jflavor].getOwx()));

			if(iflavor==jflavor) Wt[0][iflavor][jflavor] += W2t;
		}

		//printiterativeoutput(global.getdir(),iiteration);
		//printoutput(global.getdir());
		if(rank==root) cout<<endl<<endl;
	}
	if(rank==root) printoutput_S(global.getdir());
};
void controller::iterativesampling(){
        //cout<<"At controller, sampling starts"<<endl<<flush;
        const int root = 0;
        const int rank = MPI::COMM_WORLD.Get_rank();

        //cout<<"accumulation quantities initialization"<<endl;
        //double norm_tmp;
        //double mistmp;
        complex<double> norm_acc;
        complex<double> miscellaneous[4]; // 0: sign, 1: Rvalid, 2: Rproposal, 3: Racceptance

        function_w Owtmp;
        Eigen::MatrixXcd M_Owx, M_Owk;
        Eigen::MatrixXcd **M_Owxacc, **M_Owxsqacc;
        Eigen::MatrixXcd **M_Owkacc, **M_Owksqacc;

        M_Owxacc = new Eigen::MatrixXcd*[Nflavor];
        M_Owkacc = new Eigen::MatrixXcd*[Nflavor];
        M_Owxsqacc = new Eigen::MatrixXcd*[Nflavor];
        M_Owksqacc = new Eigen::MatrixXcd*[Nflavor];
        for(int iflavor=0; iflavor<Nflavor; iflavor++){
                M_Owxacc[iflavor] = new Eigen::MatrixXcd[Nflavor];
                M_Owkacc[iflavor] = new Eigen::MatrixXcd[Nflavor];
                M_Owxsqacc[iflavor] = new Eigen::MatrixXcd[Nflavor];
                M_Owksqacc[iflavor] = new Eigen::MatrixXcd[Nflavor];
        }

        // 5 : normalization, sign, Rvalid, Rproposal, Racceptance
        // 4 : M_Owxacc, M_Owkacc, M_Owxsqacc, M_Owksqacc
        int Ndata = 5 + 4*Nflavor*Nflavor*Nw*Nx;
        complex<double> *source, *result, *output0, *output1, *output2;
        source  = new complex<double>[Ndata];
        result  = new complex<double>[Ndata];
        output0  = new complex<double>[5];
        output1  = new complex<double>[Nw*Nx];
        output2  = new complex<double>[Nw*Nx];

        //int position;

        diagram Diag;
        for(int iordermax=2; iordermax<Nordermax+1; iordermax++){
                norm_acc = 0.;
                for(int i=0; i<4; i++) miscellaneous[i] = 0.;
                for(int iflavor=0; iflavor<Nflavor; iflavor++) for(int jflavor=0; jflavor<Nflavor; jflavor++){
                        M_Owxacc[iflavor][jflavor] = Eigen::MatrixXcd::Zero(Nw,Nx);
                        M_Owkacc[iflavor][jflavor] = Eigen::MatrixXcd::Zero(Nw,Nx);
                        M_Owxsqacc[iflavor][jflavor] = Eigen::MatrixXcd::Zero(Nw,Nx);
                        M_Owksqacc[iflavor][jflavor] = Eigen::MatrixXcd::Zero(Nw,Nx);
                }

                Diag.init_P_diagram();
                Diag.setNordermax(iordermax);

                MPI::COMM_WORLD.Barrier();
                if(rank==root) cout<<"controller::"<<iordermax<<"th-order sampling : "<<flush;
                MPI::COMM_WORLD.Barrier();
                for(int i=0; i<Nerr; i++){
                        // sampling
                        cout<<"*"<<flush;
                        Diag.diagMC_P(iordermax);

			}
                MPI::COMM_WORLD.Barrier();
                if(rank==root) cout<<"|"<<endl<<flush;
	}
	}

        delete[] source;
        delete[] result;
        delete[] output1;
        delete[] output2;
        delete[] M_Owxacc;
        delete[] M_Owxsqacc;
        delete[] M_Owkacc;
        delete[] M_Owksqacc;
};
void controller::sampling(){
	const int root = 0;
	const int rank = MPI::COMM_WORLD.Get_rank();
	if(rank==root) cout<<"At controller, sampling starts"<<endl<<flush;

	if(rank==root) cout<<"accumulation quantities initialization"<<endl;
	double norm_tmp;
	complex<double> norm_acc[2][Nordermax+1];
	double mistmp;
	complex<double> miscellaneous[6]; 
	// 0: Rproposal, 1: Racceptance, 2: Rvalid, 3: Rirreducible, 4: RGcompact, 5: RWcompact 

	double Ploc_tmp;
	complex<double> Ploc_acc[Nordermax];

	selfenergy_w Swtmp;
	polarization_w Pwtmp;

	Eigen::MatrixXcd M_Swx, M_Swk;
	Eigen::MatrixXcd ***M_swxacc, ***M_swxsqacc;
	Eigen::MatrixXcd ***M_swkacc, ***M_swksqacc;
	Eigen::MatrixXcd ***M_Swxacc, ***M_Swxsqacc;
	Eigen::MatrixXcd ***M_Swkacc, ***M_Swksqacc;

	Eigen::MatrixXcd M_Pwx, M_Pwk;
	Eigen::MatrixXcd ***M_pwxacc, ***M_pwxsqacc;
	Eigen::MatrixXcd ***M_pwkacc, ***M_pwksqacc;
	Eigen::MatrixXcd ***M_Pwxacc, ***M_Pwxsqacc;
	Eigen::MatrixXcd ***M_Pwkacc, ***M_Pwksqacc;

	M_swxacc = new Eigen::MatrixXcd**[Nordermax];
	M_swkacc = new Eigen::MatrixXcd**[Nordermax];
	M_swxsqacc = new Eigen::MatrixXcd**[Nordermax];
	M_swksqacc = new Eigen::MatrixXcd**[Nordermax];

	M_Swxacc = new Eigen::MatrixXcd**[Nordermax];
	M_Swkacc = new Eigen::MatrixXcd**[Nordermax];
	M_Swxsqacc = new Eigen::MatrixXcd**[Nordermax];
	M_Swksqacc = new Eigen::MatrixXcd**[Nordermax];

	M_pwxacc = new Eigen::MatrixXcd**[Nordermax];
	M_pwkacc = new Eigen::MatrixXcd**[Nordermax];
	M_pwxsqacc = new Eigen::MatrixXcd**[Nordermax];
	M_pwksqacc = new Eigen::MatrixXcd**[Nordermax];

	M_Pwxacc = new Eigen::MatrixXcd**[Nordermax];
	M_Pwkacc = new Eigen::MatrixXcd**[Nordermax];
	M_Pwxsqacc = new Eigen::MatrixXcd**[Nordermax];
	M_Pwksqacc = new Eigen::MatrixXcd**[Nordermax];

	for(int iorder=0; iorder<Nordermax; iorder++){
		M_swxacc[iorder] = new Eigen::MatrixXcd*[Nflavor];
		M_swkacc[iorder] = new Eigen::MatrixXcd*[Nflavor];
		M_swxsqacc[iorder] = new Eigen::MatrixXcd*[Nflavor];
		M_swksqacc[iorder] = new Eigen::MatrixXcd*[Nflavor];

		M_Swxacc[iorder] = new Eigen::MatrixXcd*[Nflavor];
		M_Swkacc[iorder] = new Eigen::MatrixXcd*[Nflavor];
		M_Swxsqacc[iorder] = new Eigen::MatrixXcd*[Nflavor];
		M_Swksqacc[iorder] = new Eigen::MatrixXcd*[Nflavor];

		M_pwxacc[iorder] = new Eigen::MatrixXcd*[Nflavor];
		M_pwkacc[iorder] = new Eigen::MatrixXcd*[Nflavor];
		M_pwxsqacc[iorder] = new Eigen::MatrixXcd*[Nflavor];
		M_pwksqacc[iorder] = new Eigen::MatrixXcd*[Nflavor];

		M_Pwxacc[iorder] = new Eigen::MatrixXcd*[Nflavor];
		M_Pwkacc[iorder] = new Eigen::MatrixXcd*[Nflavor];
		M_Pwxsqacc[iorder] = new Eigen::MatrixXcd*[Nflavor];
		M_Pwksqacc[iorder] = new Eigen::MatrixXcd*[Nflavor];

		for(int iflavor=0; iflavor<Nflavor; iflavor++){
			M_swxacc[iorder][iflavor] = new Eigen::MatrixXcd[Nflavor];
			M_swkacc[iorder][iflavor] = new Eigen::MatrixXcd[Nflavor];
			M_swxsqacc[iorder][iflavor] = new Eigen::MatrixXcd[Nflavor];
			M_swksqacc[iorder][iflavor] = new Eigen::MatrixXcd[Nflavor];

			M_Swxacc[iorder][iflavor] = new Eigen::MatrixXcd[Nflavor];
			M_Swkacc[iorder][iflavor] = new Eigen::MatrixXcd[Nflavor];
			M_Swxsqacc[iorder][iflavor] = new Eigen::MatrixXcd[Nflavor];
			M_Swksqacc[iorder][iflavor] = new Eigen::MatrixXcd[Nflavor];

			M_pwxacc[iorder][iflavor] = new Eigen::MatrixXcd[Nflavor];
			M_pwkacc[iorder][iflavor] = new Eigen::MatrixXcd[Nflavor];
			M_pwxsqacc[iorder][iflavor] = new Eigen::MatrixXcd[Nflavor];
			M_pwksqacc[iorder][iflavor] = new Eigen::MatrixXcd[Nflavor];

			M_Pwxacc[iorder][iflavor] = new Eigen::MatrixXcd[Nflavor];
			M_Pwkacc[iorder][iflavor] = new Eigen::MatrixXcd[Nflavor];
			M_Pwxsqacc[iorder][iflavor] = new Eigen::MatrixXcd[Nflavor];
			M_Pwksqacc[iorder][iflavor] = new Eigen::MatrixXcd[Nflavor];
		}
	}


	//Diag.printoutput(global.getdir());
	//printinput(global.getdir());

	if(rank==root) cout<<"initialize accumulation quatities"<<endl;
	diagram Diag;

	for(int i=0; i<Nordermax+1; i++) for(int sp=0; sp<2; sp++) norm_acc[sp][i] = 0.;

	for(int i=0; i<6; i++) miscellaneous[i] = 0.;
	for(int iorder=0; iorder<Nordermax; iorder++){
		Ploc_acc[iorder] = 0.;
		for(int iflavor=0; iflavor<Nflavor; iflavor++) for(int jflavor=0; jflavor<Nflavor; jflavor++){
			M_Swxacc[iorder][iflavor][jflavor] = Eigen::MatrixXcd::Zero(Nw,Nx);
			M_Swkacc[iorder][iflavor][jflavor] = Eigen::MatrixXcd::Zero(Nw,Nx);
			M_Swxsqacc[iorder][iflavor][jflavor] = Eigen::MatrixXcd::Zero(Nw,Nx);
			M_Swksqacc[iorder][iflavor][jflavor] = Eigen::MatrixXcd::Zero(Nw,Nx);

			M_swxacc[iorder][iflavor][jflavor] = Eigen::MatrixXcd::Zero(Nw,Nx);
			M_swkacc[iorder][iflavor][jflavor] = Eigen::MatrixXcd::Zero(Nw,Nx);
			M_swxsqacc[iorder][iflavor][jflavor] = Eigen::MatrixXcd::Zero(Nw,Nx);
			M_swksqacc[iorder][iflavor][jflavor] = Eigen::MatrixXcd::Zero(Nw,Nx);


			M_Pwxacc[iorder][iflavor][jflavor] = Eigen::MatrixXcd::Zero(Nw,Nx);
			M_Pwkacc[iorder][iflavor][jflavor] = Eigen::MatrixXcd::Zero(Nw,Nx);
			M_Pwxsqacc[iorder][iflavor][jflavor] = Eigen::MatrixXcd::Zero(Nw,Nx);
			M_Pwksqacc[iorder][iflavor][jflavor] = Eigen::MatrixXcd::Zero(Nw,Nx);

			M_pwxacc[iorder][iflavor][jflavor] = Eigen::MatrixXcd::Zero(Nw,Nx);
			M_pwkacc[iorder][iflavor][jflavor] = Eigen::MatrixXcd::Zero(Nw,Nx);
			M_pwxsqacc[iorder][iflavor][jflavor] = Eigen::MatrixXcd::Zero(Nw,Nx);
			M_pwksqacc[iorder][iflavor][jflavor] = Eigen::MatrixXcd::Zero(Nw,Nx);
		}
	}

	if(rank==root) cout<<"actual sampling process"<<endl;
	//Diag.init_S_diagram();
	Diag.init_P_diagram();

	//Diag.thermalization();
	if(rank==root) cout<<"controller::sampling : "<<flush;
	MPI::COMM_WORLD.Barrier();
	for(int i=0; i<Nerr; i++){
		cout<<"*"<<flush;
		Diag.bold_diagMC_GW();
		Diag.normalize_S(0,0,U*n0);
		Diag.normalize_P(0,0,U*n0);

		//cout<<imag(Diag.getSw(1,1,1).getOwx(0,0))<<",";
		//cout<<real(Diag.getPw(2,1,1).getOwx(0,0))<<",";

		//mistmp = Diag.getAvesign();		miscellaneous[0] += complex<double>(mistmp, mistmp*mistmp);
		mistmp = Diag.getRproposal();		miscellaneous[0] += complex<double>(mistmp, mistmp*mistmp);
		mistmp = Diag.getRacceptance();		miscellaneous[1] += complex<double>(mistmp, mistmp*mistmp);
		mistmp = Diag.getRvalid();		miscellaneous[2] += complex<double>(mistmp, mistmp*mistmp);
		mistmp = Diag.getRirreducible();	miscellaneous[3] += complex<double>(mistmp, mistmp*mistmp);
		mistmp = Diag.getRGcompact();		miscellaneous[4] += complex<double>(mistmp, mistmp*mistmp);
		mistmp = Diag.getRWcompact();		miscellaneous[5] += complex<double>(mistmp, mistmp*mistmp);

		for(int n=0; n<Nordermax+1; n++) for(int sp=0; sp<2; sp++){
			norm_tmp = Diag.getnorm(sp,n)*(U*n0)/Diag.getnorm(0,0);
			norm_acc[sp][n] += complex<double>(norm_tmp,norm_tmp*norm_tmp);
		}

		for(int n=0; n<Nordermax; n++){
			Ploc_tmp = Diag.getPloc(n);
			Ploc_acc[n] += complex<double>(Ploc_tmp,Ploc_tmp*Ploc_tmp);
			for(int iflavor=0; iflavor<Nflavor; iflavor++) for(int jflavor=0; jflavor<Nflavor; jflavor++){
				Swtmp = Diag.getSwSum(n+1,iflavor,jflavor);
				M_Swx = Swtmp.getOwx();
				M_Swk = Swtmp.getOwk();
				M_Swxacc[n][iflavor][jflavor] += M_Swx;
				M_Swkacc[n][iflavor][jflavor] += M_Swk;
				M_Swxsqacc[n][iflavor][jflavor] += M_Swx.real().cwiseAbs2() + complex<double>(0.0,1.0)*M_Swx.imag().cwiseAbs2();
				M_Swksqacc[n][iflavor][jflavor] += M_Swk.real().cwiseAbs2() + complex<double>(0.0,1.0)*M_Swk.imag().cwiseAbs2();

				Swtmp = Diag.getSw(n,iflavor,jflavor);
				M_Swx = Swtmp.getOwx();
				M_Swk = Swtmp.getOwk();
				M_swxacc[n][iflavor][jflavor] += M_Swx;
				M_swkacc[n][iflavor][jflavor] += M_Swk;
				M_swxsqacc[n][iflavor][jflavor] += M_Swx.real().cwiseAbs2() + complex<double>(0.0,1.0)*M_Swx.imag().cwiseAbs2();
				M_swksqacc[n][iflavor][jflavor] += M_Swk.real().cwiseAbs2() + complex<double>(0.0,1.0)*M_Swk.imag().cwiseAbs2();


				Pwtmp = Diag.getPwSum(n,iflavor,jflavor);
				M_Pwx = Pwtmp.getOwx();
				M_Pwk = Pwtmp.getOwk();
				M_Pwxacc[n][iflavor][jflavor] += M_Pwx;
				M_Pwkacc[n][iflavor][jflavor] += M_Pwk;
				M_Pwxsqacc[n][iflavor][jflavor] += M_Pwx.real().cwiseAbs2() + complex<double>(0.0,1.0)*M_Pwx.imag().cwiseAbs2();
				M_Pwksqacc[n][iflavor][jflavor] += M_Pwk.real().cwiseAbs2() + complex<double>(0.0,1.0)*M_Pwk.imag().cwiseAbs2();

				Pwtmp = Diag.getPw(n,iflavor,jflavor);
				M_Pwx = Pwtmp.getOwx();
				M_Pwk = Pwtmp.getOwk();
				M_pwxacc[n][iflavor][jflavor] += M_Pwx;
				M_pwkacc[n][iflavor][jflavor] += M_Pwk;
				M_pwxsqacc[n][iflavor][jflavor] += M_Pwx.real().cwiseAbs2() + complex<double>(0.0,1.0)*M_Pwx.imag().cwiseAbs2();
				M_pwksqacc[n][iflavor][jflavor] += M_Pwk.real().cwiseAbs2() + complex<double>(0.0,1.0)*M_Pwk.imag().cwiseAbs2();
			}
		}
	}

	if(rank==root) cout<<endl<<flush;

	if(rank==root) cout<<"MPI communication starts"<<endl<<flush;
	if(rank==root) cout<<"MPI Pack"<<endl<<flush;
	
	// 2*(Nordermax+1) : 0) normalization[s], 1) normalization[p], 
	// 1 : Ploc
	// 6 : Rvalid, Rproposal, Racceptance
	// 8 : M_Swxacc, M_Swkacc, M_Swxsqacc, M_Swksqacc (S<-->s)
	// 8 : M_Pwxacc, M_Pwkacc, M_Pwxsqacc, M_Pwksqacc (P<-->p)
	int Ndata = 2*(Nordermax+1) + 1*(Nordermax) + 6 + 16*Nordermax*Nflavor*Nflavor*Nw*Nx;
	complex<double> *source, *result, *outputn, *output0, *output1, *output2;
	source  = new complex<double>[Ndata];
	result  = new complex<double>[Ndata];
	outputn  = new complex<double>[2*(Nordermax+1)];
	output0  = new complex<double>[6];
	output1  = new complex<double>[Nw*Nx];
	output2  = new complex<double>[Nw*Nx];

	int position = 0;
	MPI::DOUBLE_COMPLEX.Pack(norm_acc,2*(Nordermax+1),source,16*Ndata,position,MPI::COMM_WORLD);
	MPI::DOUBLE_COMPLEX.Pack(Ploc_acc,Nordermax,source,16*Ndata,position,MPI::COMM_WORLD);
	MPI::DOUBLE_COMPLEX.Pack(miscellaneous,6,source,16*Ndata,position,MPI::COMM_WORLD);
	for(int n=0; n<Nordermax; n++)
	for(int iflavor=0; iflavor<Nflavor; iflavor++) for(int jflavor=0; jflavor<Nflavor; jflavor++){
		MPI::DOUBLE_COMPLEX.Pack(M_Swxacc[n][iflavor][jflavor].data(),Nx*Nw,source,16*Ndata,position,MPI::COMM_WORLD);
		MPI::DOUBLE_COMPLEX.Pack(M_Swxsqacc[n][iflavor][jflavor].data(),Nx*Nw,source,16*Ndata,position,MPI::COMM_WORLD);
		MPI::DOUBLE_COMPLEX.Pack(M_Swkacc[n][iflavor][jflavor].data(),Nx*Nw,source,16*Ndata,position,MPI::COMM_WORLD);
		MPI::DOUBLE_COMPLEX.Pack(M_Swksqacc[n][iflavor][jflavor].data(),Nx*Nw,source,16*Ndata,position,MPI::COMM_WORLD);

		MPI::DOUBLE_COMPLEX.Pack(M_swxacc[n][iflavor][jflavor].data(),Nx*Nw,source,16*Ndata,position,MPI::COMM_WORLD);
		MPI::DOUBLE_COMPLEX.Pack(M_swxsqacc[n][iflavor][jflavor].data(),Nx*Nw,source,16*Ndata,position,MPI::COMM_WORLD);
		MPI::DOUBLE_COMPLEX.Pack(M_swkacc[n][iflavor][jflavor].data(),Nx*Nw,source,16*Ndata,position,MPI::COMM_WORLD);
		MPI::DOUBLE_COMPLEX.Pack(M_swksqacc[n][iflavor][jflavor].data(),Nx*Nw,source,16*Ndata,position,MPI::COMM_WORLD);
	}
	for(int n=0; n<Nordermax; n++)
	for(int iflavor=0; iflavor<Nflavor; iflavor++) for(int jflavor=0; jflavor<Nflavor; jflavor++){
		MPI::DOUBLE_COMPLEX.Pack(M_Pwxacc[n][iflavor][jflavor].data(),Nx*Nw,source,16*Ndata,position,MPI::COMM_WORLD);
		MPI::DOUBLE_COMPLEX.Pack(M_Pwxsqacc[n][iflavor][jflavor].data(),Nx*Nw,source,16*Ndata,position,MPI::COMM_WORLD);
		MPI::DOUBLE_COMPLEX.Pack(M_Pwkacc[n][iflavor][jflavor].data(),Nx*Nw,source,16*Ndata,position,MPI::COMM_WORLD);
		MPI::DOUBLE_COMPLEX.Pack(M_Pwksqacc[n][iflavor][jflavor].data(),Nx*Nw,source,16*Ndata,position,MPI::COMM_WORLD);

		MPI::DOUBLE_COMPLEX.Pack(M_pwxacc[n][iflavor][jflavor].data(),Nx*Nw,source,16*Ndata,position,MPI::COMM_WORLD);
		MPI::DOUBLE_COMPLEX.Pack(M_pwxsqacc[n][iflavor][jflavor].data(),Nx*Nw,source,16*Ndata,position,MPI::COMM_WORLD);
		MPI::DOUBLE_COMPLEX.Pack(M_pwkacc[n][iflavor][jflavor].data(),Nx*Nw,source,16*Ndata,position,MPI::COMM_WORLD);
		MPI::DOUBLE_COMPLEX.Pack(M_pwksqacc[n][iflavor][jflavor].data(),Nx*Nw,source,16*Ndata,position,MPI::COMM_WORLD);
	}

	if(rank==root) cout<<"MPI Reduce and Bcast"<<endl<<flush;
	MPI::COMM_WORLD.Reduce(source,result,Ndata,MPI::DOUBLE_COMPLEX,MPI::SUM,root);
	MPI::COMM_WORLD.Bcast(result,Ndata,MPI::DOUBLE_COMPLEX,root);

	if(rank==root) cout<<"MPI Unpack"<<endl<<flush;
	position = 0;
	if(rank==root) cout<<"normalization"<<endl<<flush;
	MPI::DOUBLE_COMPLEX.Unpack(result,16*Ndata,outputn,2*(Nordermax+1),position,MPI::COMM_WORLD);
	for(int sp=0; sp<2; sp++){
		for(int iorder=0; iorder<Nordermax+1; iorder++){
			normalization[sp][iorder] = outputn[sp*(Nordermax+1)+iorder].real()/NerrTot;
			dnormalization[sp][iorder] = 
				sqrt(abs(
					outputn[sp*(Nordermax+1)+iorder].imag()/NerrTot
					- normalization[sp][iorder]*normalization[sp][iorder]
				)/NerrTot);
		};
	}
	if(rank==root) cout<<"Ptxloc"<<endl<<flush;
	MPI::DOUBLE_COMPLEX.Unpack(result,16*Ndata,outputn,Nordermax,position,MPI::COMM_WORLD);
	for(int iorder=0; iorder<Nordermax; iorder++){
		Ploc[iorder] = outputn[iorder].real()/NerrTot;
		dPloc[iorder] = 
			sqrt(abs(
				outputn[iorder].imag()/NerrTot
				- Ploc[iorder]*Ploc[iorder]
			)/NerrTot);
	};
	if(rank==root) cout<<"miscellaneous"<<endl<<flush;
	MPI::DOUBLE_COMPLEX.Unpack(result,16*Ndata,output0,6,position,MPI::COMM_WORLD);
	Rproposal = output0[0].real()/NerrTot;
	dRproposal = sqrt(abs(
			output0[0].imag()/NerrTot
			- Rproposal*Rproposal
			)/NerrTot);
	Racceptance = output0[1].real()/NerrTot;
	dRacceptance = sqrt(abs(
			output0[1].imag()/NerrTot
			- Racceptance*Racceptance
			)/NerrTot);
	Rvalid = output0[2].real()/NerrTot;
	dRvalid = sqrt(abs(
			output0[2].imag()/NerrTot
			- Rvalid*Rvalid
			)/NerrTot);
	Rirreducible = output0[3].real()/NerrTot;
	dRirreducible = sqrt(abs(
			output0[3].imag()/NerrTot
			- Rirreducible*Rirreducible
			)/NerrTot);
	RGcompact = output0[4].real()/NerrTot;
	dRGcompact = sqrt(abs(
			output0[4].imag()/NerrTot
			- RGcompact*RGcompact
			)/NerrTot);
	RWcompact = output0[5].real()/NerrTot;
	dRWcompact = sqrt(abs(
			output0[5].imag()/NerrTot
			- RWcompact*RWcompact
			)/NerrTot);

	if(rank==root) cout<<"selfenergy"<<endl<<flush;
	for(int n=0; n<Nordermax; n++)
	for(int iflavor=0; iflavor<Nflavor; iflavor++) for(int jflavor=0; jflavor<Nflavor; jflavor++){
		//cout<<"("<<n<<","<<iflavor<<","<<jflavor<<")"<<flush;
		MPI::DOUBLE_COMPLEX.Unpack(result,16*Ndata,output1,Nx*Nw,position,MPI::COMM_WORLD);
		MPI::DOUBLE_COMPLEX.Unpack(result,16*Ndata,output2,Nx*Nw,position,MPI::COMM_WORLD);
		for(int i=0; i<Nw*Nx; i++){
			output1[i] /= static_cast<double>(NerrTot);
			output2[i] = complex<double>(
					sqrt(abs(
						output2[i].real()/static_cast<double>(NerrTot) 
						- output1[i].real()*output1[i].real()
					)/static_cast<double>(NerrTot)),

					sqrt(abs(
						output2[i].imag()/static_cast<double>(NerrTot) 
						- output1[i].imag()*output1[i].imag()
					)/static_cast<double>(NerrTot))
					);
		}
		Sw[n][iflavor][jflavor].setOwx(output1);
		dSw[n][iflavor][jflavor].setOwx(output2);


		MPI::DOUBLE_COMPLEX.Unpack(result,16*Ndata,output1,Nx*Nw,position,MPI::COMM_WORLD);
		MPI::DOUBLE_COMPLEX.Unpack(result,16*Ndata,output2,Nx*Nw,position,MPI::COMM_WORLD);
		for(int i=0; i<Nw*Nx; i++){
			output1[i] /= static_cast<double>(NerrTot);
			output2[i] = complex<double>(
					sqrt(abs(
						output2[i].real()/static_cast<double>(NerrTot) 
						- output1[i].real()*output1[i].real()
					)/static_cast<double>(NerrTot)),

					sqrt(abs(
						output2[i].imag()/static_cast<double>(NerrTot) 
						- output1[i].imag()*output1[i].imag()
					)/static_cast<double>(NerrTot))
					);
		}
		Sw[n][iflavor][jflavor].setOwk(output1);
		dSw[n][iflavor][jflavor].setOwk(output2);

		MPI::DOUBLE_COMPLEX.Unpack(result,16*Ndata,output1,Nx*Nw,position,MPI::COMM_WORLD);
		MPI::DOUBLE_COMPLEX.Unpack(result,16*Ndata,output2,Nx*Nw,position,MPI::COMM_WORLD);
		for(int i=0; i<Nw*Nx; i++){
			output1[i] /= static_cast<double>(NerrTot);
			output2[i] = complex<double>(
					sqrt(abs(
						output2[i].real()/static_cast<double>(NerrTot) 
						- output1[i].real()*output1[i].real()
					)/static_cast<double>(NerrTot)),

					sqrt(abs(
						output2[i].imag()/static_cast<double>(NerrTot) 
						- output1[i].imag()*output1[i].imag()
					)/static_cast<double>(NerrTot))
					);
		}
		sw[n][iflavor][jflavor].setOwx(output1);
		dsw[n][iflavor][jflavor].setOwx(output2);


		MPI::DOUBLE_COMPLEX.Unpack(result,16*Ndata,output1,Nx*Nw,position,MPI::COMM_WORLD);
		MPI::DOUBLE_COMPLEX.Unpack(result,16*Ndata,output2,Nx*Nw,position,MPI::COMM_WORLD);
		for(int i=0; i<Nw*Nx; i++){
			output1[i] /= static_cast<double>(NerrTot);
			output2[i] = complex<double>(
					sqrt(abs(
						output2[i].real()/static_cast<double>(NerrTot) 
						- output1[i].real()*output1[i].real()
					)/static_cast<double>(NerrTot)),

					sqrt(abs(
						output2[i].imag()/static_cast<double>(NerrTot) 
						- output1[i].imag()*output1[i].imag()
					)/static_cast<double>(NerrTot))
					);
		}
		sw[n][iflavor][jflavor].setOwk(output1);
		dsw[n][iflavor][jflavor].setOwk(output2);
	}

	if(rank==root) cout<<"polarization"<<endl<<flush;
	for(int n=0; n<Nordermax; n++)
	for(int iflavor=0; iflavor<Nflavor; iflavor++) for(int jflavor=0; jflavor<Nflavor; jflavor++){
		//cout<<"("<<n<<","<<iflavor<<","<<jflavor<<")"<<flush;
		MPI::DOUBLE_COMPLEX.Unpack(result,16*Ndata,output1,Nx*Nw,position,MPI::COMM_WORLD);
		MPI::DOUBLE_COMPLEX.Unpack(result,16*Ndata,output2,Nx*Nw,position,MPI::COMM_WORLD);
		for(int i=0; i<Nw*Nx; i++){
			output1[i] /= static_cast<double>(NerrTot);
			output2[i] = complex<double>(
					sqrt(abs(
						output2[i].real()/static_cast<double>(NerrTot) 
						- output1[i].real()*output1[i].real()
					)/static_cast<double>(NerrTot)),

					sqrt(abs(
						output2[i].imag()/static_cast<double>(NerrTot) 
						- output1[i].imag()*output1[i].imag()
					)/static_cast<double>(NerrTot))
					);
		}
		Pw[n][iflavor][jflavor].setOwx(output1);
		dPw[n][iflavor][jflavor].setOwx(output2);

		MPI::DOUBLE_COMPLEX.Unpack(result,16*Ndata,output1,Nx*Nw,position,MPI::COMM_WORLD);
		MPI::DOUBLE_COMPLEX.Unpack(result,16*Ndata,output2,Nx*Nw,position,MPI::COMM_WORLD);
		for(int i=0; i<Nw*Nx; i++){
			output1[i] /= static_cast<double>(NerrTot);
			output2[i] = complex<double>(
					sqrt(abs(
						output2[i].real()/static_cast<double>(NerrTot) 
						- output1[i].real()*output1[i].real()
					)/static_cast<double>(NerrTot)),

					sqrt(abs(
						output2[i].imag()/static_cast<double>(NerrTot) 
						- output1[i].imag()*output1[i].imag()
					)/static_cast<double>(NerrTot))
					);
		}
		Pw[n][iflavor][jflavor].setOwk(output1);
		dPw[n][iflavor][jflavor].setOwk(output2);


		MPI::DOUBLE_COMPLEX.Unpack(result,16*Ndata,output1,Nx*Nw,position,MPI::COMM_WORLD);
		MPI::DOUBLE_COMPLEX.Unpack(result,16*Ndata,output2,Nx*Nw,position,MPI::COMM_WORLD);
		for(int i=0; i<Nw*Nx; i++){
			output1[i] /= static_cast<double>(NerrTot);
			output2[i] = complex<double>(
					sqrt(abs(
						output2[i].real()/static_cast<double>(NerrTot) 
						- output1[i].real()*output1[i].real()
					)/static_cast<double>(NerrTot)),

					sqrt(abs(
						output2[i].imag()/static_cast<double>(NerrTot) 
						- output1[i].imag()*output1[i].imag()
					)/static_cast<double>(NerrTot))
					);
		}
		pw[n][iflavor][jflavor].setOwx(output1);
		dpw[n][iflavor][jflavor].setOwx(output2);

		MPI::DOUBLE_COMPLEX.Unpack(result,16*Ndata,output1,Nx*Nw,position,MPI::COMM_WORLD);
		MPI::DOUBLE_COMPLEX.Unpack(result,16*Ndata,output2,Nx*Nw,position,MPI::COMM_WORLD);
		for(int i=0; i<Nw*Nx; i++){
			output1[i] /= static_cast<double>(NerrTot);
			output2[i] = complex<double>(
					sqrt(abs(
						output2[i].real()/static_cast<double>(NerrTot) 
						- output1[i].real()*output1[i].real()
					)/static_cast<double>(NerrTot)),

					sqrt(abs(
						output2[i].imag()/static_cast<double>(NerrTot) 
						- output1[i].imag()*output1[i].imag()
					)/static_cast<double>(NerrTot))
					);
		}
		pw[n][iflavor][jflavor].setOwk(output1);
		dpw[n][iflavor][jflavor].setOwk(output2);
	}
	//for(int iflavor=0; iflavor<Nflavor; iflavor++) for(int jflavor=0; jflavor<Nflavor; jflavor++)
	//	dPwSum[iflavor][jflavor] = dPw[Nordermax][iflavor][jflavor];

	if(rank==root) cout<<"memory deallocation"<<endl<<flush;

	delete[] source;
	delete[] result;

	delete[] outputn;
	delete[] output0;
	delete[] output1;
	delete[] output2;

	for(int iorder=0; iorder<Nordermax; iorder++){
		for(int iflavor=0; iflavor<Nflavor; iflavor++){
			delete[] M_Swxacc[iorder][iflavor];
			delete[] M_Swxsqacc[iorder][iflavor];
			delete[] M_Swkacc[iorder][iflavor];
			delete[] M_Swksqacc[iorder][iflavor];

			delete[] M_swxacc[iorder][iflavor];
			delete[] M_swxsqacc[iorder][iflavor];
			delete[] M_swkacc[iorder][iflavor];
			delete[] M_swksqacc[iorder][iflavor];
		}
		delete[] M_Swxacc[iorder];
		delete[] M_Swxsqacc[iorder];
		delete[] M_Swkacc[iorder];
		delete[] M_Swksqacc[iorder];

		delete[] M_swxacc[iorder];
		delete[] M_swxsqacc[iorder];
		delete[] M_swkacc[iorder];
		delete[] M_swksqacc[iorder];
	}
	delete[] M_Swxacc;
	delete[] M_Swxsqacc;
	delete[] M_Swkacc;
	delete[] M_Swksqacc;

	delete[] M_swxacc;
	delete[] M_swxsqacc;
	delete[] M_swkacc;
	delete[] M_swksqacc;

	for(int iorder=0; iorder<Nordermax; iorder++){
		for(int iflavor=0; iflavor<Nflavor; iflavor++){
			delete[] M_Pwxacc[iorder][iflavor];
			delete[] M_Pwxsqacc[iorder][iflavor];
			delete[] M_Pwkacc[iorder][iflavor];
			delete[] M_Pwksqacc[iorder][iflavor];

			delete[] M_pwxacc[iorder][iflavor];
			delete[] M_pwxsqacc[iorder][iflavor];
			delete[] M_pwkacc[iorder][iflavor];
			delete[] M_pwksqacc[iorder][iflavor];
		}
		delete[] M_Pwxacc[iorder];
		delete[] M_Pwxsqacc[iorder];
		delete[] M_Pwkacc[iorder];
		delete[] M_Pwksqacc[iorder];

		delete[] M_pwxacc[iorder];
		delete[] M_pwxsqacc[iorder];
		delete[] M_pwkacc[iorder];
		delete[] M_pwksqacc[iorder];
	}
	delete[] M_Pwxacc;
	delete[] M_Pwxsqacc;
	delete[] M_Pwkacc;
	delete[] M_Pwksqacc;

	delete[] M_pwxacc;
	delete[] M_pwxsqacc;
	delete[] M_pwkacc;
	delete[] M_pwksqacc;

	if(rank==root) cout<<"At controller, sampling ends"<<endl<<flush;
};
/*void controller::sampling(){
        //cout<<"At controller, sampling starts"<<endl<<flush;
        const int root = 0;
        const int rank = MPI::COMM_WORLD.Get_rank();

        //cout<<"accumulation quantities initialization"<<endl;
        double norm_tmp;
        complex<double> norm_acc[Nordermax];
        double mistmp;
        complex<double> miscellaneous[4]; // 0: sign, 1: Rvalid, 2: Rproposal, 3: Racceptance
        for(int i=0; i<4; i++) miscellaneous[i] = 0.;

        selfenergy_w Swtmp;
        Eigen::MatrixXcd M_Swx, M_Swk;
        Eigen::MatrixXcd ***M_Swxacc, ***M_Swxsqacc;
        Eigen::MatrixXcd ***M_Swkacc, ***M_Swksqacc;
        M_Swxacc = new Eigen::MatrixXcd**[Nordermax];
        M_Swkacc = new Eigen::MatrixXcd**[Nordermax];
        M_Swxsqacc = new Eigen::MatrixXcd**[Nordermax];
        M_Swksqacc = new Eigen::MatrixXcd**[Nordermax];
        for(int iorder=0; iorder<Nordermax; iorder++){
                M_Swxacc[iorder] = new Eigen::MatrixXcd*[Nflavor];
                M_Swkacc[iorder] = new Eigen::MatrixXcd*[Nflavor];
                M_Swxsqacc[iorder] = new Eigen::MatrixXcd*[Nflavor];
                M_Swksqacc[iorder] = new Eigen::MatrixXcd*[Nflavor];
                for(int iflavor=0; iflavor<Nflavor; iflavor++){
                        M_Swxacc[iorder][iflavor] = new Eigen::MatrixXcd[Nflavor];
                        M_Swkacc[iorder][iflavor] = new Eigen::MatrixXcd[Nflavor];
                        M_Swxsqacc[iorder][iflavor] = new Eigen::MatrixXcd[Nflavor];
                        M_Swksqacc[iorder][iflavor] = new Eigen::MatrixXcd[Nflavor];
                        for(int jflavor=0; jflavor<Nflavor; jflavor++){
                                M_Swxacc[iorder][iflavor][jflavor] = Eigen::MatrixXcd::Zero(Nw,Nx);
                                M_Swkacc[iorder][iflavor][jflavor] = Eigen::MatrixXcd::Zero(Nw,Nx);
                                M_Swxsqacc[iorder][iflavor][jflavor] = Eigen::MatrixXcd::Zero(Nw,Nx);
                                M_Swksqacc[iorder][iflavor][jflavor] = Eigen::MatrixXcd::Zero(Nw,Nx);
                        }
                }
        }

        //cout<<"actual sampling process"<<endl;
        diagram Diag;
        Diag.init_S_diagram();
        MPI::COMM_WORLD.Barrier();
        if(rank==root) cout<<"   controller::sampling : "<<flush;
        MPI::COMM_WORLD.Barrier();
        for(int i=0; i<Nerr; i++){
                //cout<<"***** "<<i<<"th controller::sampling ****"<<endl;
                cout<<"*"<<flush;
                //Diag.bold_diagMC_G();
                //Diag.normalize_S(1,U*n0);
                Diag.iterative_bold_diagMC_G(1,U*n0);

                mistmp = Diag.getAvesign();     miscellaneous[0] += complex<double>(mistmp, mistmp*mistmp);
                mistmp = Diag.getRvalid();      miscellaneous[1] += complex<double>(mistmp, mistmp*mistmp);
                mistmp = Diag.getRproposal();   miscellaneous[2] += complex<double>(mistmp, mistmp*mistmp);
                mistmp = Diag.getRacceptance(); miscellaneous[3] += complex<double>(mistmp, mistmp*mistmp);
                for(int n=0; n<Nordermax; n++){
                        norm_tmp = Diag.getnorm(n+1)*(U*n0)/Diag.getnorm(1);
                        norm_acc[n] += complex<double>(norm_tmp,norm_tmp*norm_tmp);
                        for(int iflavor=0; iflavor<Nflavor; iflavor++) for(int jflavor=0; jflavor<Nflavor; jflavor++){
                                Swtmp = Diag.getSw(n,iflavor,jflavor);
                                M_Swx = Swtmp.getOwx();
                                M_Swk = Swtmp.getOwk();
                                M_Swxacc[n][iflavor][jflavor] += M_Swx;
                                M_Swkacc[n][iflavor][jflavor] += M_Swk;
                                M_Swxsqacc[n][iflavor][jflavor] += M_Swx.real().cwiseAbs2() + complex<double>(0.0,1.0)*M_Swx.imag().cwiseAbs2();
                                M_Swksqacc[n][iflavor][jflavor] += M_Swk.real().cwiseAbs2() + complex<double>(0.0,1.0)*M_Swk.imag().cwiseAbs2();
                                //if(n==1 && iflavor==0 && jflavor==0){ cout<<"[("<<M_Swx(0,0)<<"),"; }
                                //if(n==1 && iflavor==0 && jflavor==0){ cout<<M_Swx(0,0).real()<<","; }
                        }
                }
                //cout<<(M_Swx.real().cwiseAbs2() + complex<double>(0.0,1.0)*M_Swx.imag().cwiseAbs2())<<endl;
        }
        if(rank==root) cout<<endl<<flush;
	// Nordermax : normalization
        // 4 : sign, Rvalid, Rproposal, Racceptance
        // 4 : M_Swxacc, M_Swkacc, M_Swxsqacc, M_Swksqacc
        int Ndata = Nordermax + 4 + 4*Nordermax*Nflavor*Nflavor*Nw*Nx;
        complex<double> *source, *result, *outputn, *output0, *output1, *output2;
        source  = new complex<double>[Ndata];
        result  = new complex<double>[Ndata];
        outputn  = new complex<double>[Nordermax];
        output0  = new complex<double>[4];
        output1  = new complex<double>[Nw*Nx];
        output2  = new complex<double>[Nw*Nx];

        int position = 0;
        MPI::DOUBLE_COMPLEX.Pack(norm_acc,Nordermax,source,16*Ndata,position,MPI::COMM_WORLD);
        MPI::DOUBLE_COMPLEX.Pack(miscellaneous,4,source,16*Ndata,position,MPI::COMM_WORLD);
        for(int n=0; n<Nordermax; n++)
        for(int iflavor=0; iflavor<Nflavor; iflavor++) for(int jflavor=0; jflavor<Nflavor; jflavor++){
                MPI::DOUBLE_COMPLEX.Pack(M_Swxacc[n][iflavor][jflavor].data(),Nx*Nw,source,16*Ndata,position,MPI::COMM_WORLD);
                MPI::DOUBLE_COMPLEX.Pack(M_Swxsqacc[n][iflavor][jflavor].data(),Nx*Nw,source,16*Ndata,position,MPI::COMM_WORLD);
                MPI::DOUBLE_COMPLEX.Pack(M_Swkacc[n][iflavor][jflavor].data(),Nx*Nw,source,16*Ndata,position,MPI::COMM_WORLD);
                MPI::DOUBLE_COMPLEX.Pack(M_Swksqacc[n][iflavor][jflavor].data(),Nx*Nw,source,16*Ndata,position,MPI::COMM_WORLD);
        }

        //cout<<"MPI Reduce and Bcast"<<endl<<flush;
        MPI::COMM_WORLD.Reduce(source,result,Ndata,MPI::DOUBLE_COMPLEX,MPI::SUM,root);
        MPI::COMM_WORLD.Bcast(result,Ndata,MPI::DOUBLE_COMPLEX,root);

        // cout<<"MPI Unpack"<<endl<<flush;
        position = 0;
        // normalization
        MPI::DOUBLE_COMPLEX.Unpack(result,16*Ndata,outputn,Nordermax,position,MPI::COMM_WORLD);
        for(int iorder=0; iorder<Nordermax; iorder++){
                normalization[iorder] = outputn[iorder].real()/NerrTot;
                dnormalization[iorder] =
                        sqrt(abs(
                                outputn[iorder].imag()/NerrTot
                                - normalization[iorder]*normalization[iorder]
                        )/NerrTot);
        };
        // miscellaneous
        MPI::DOUBLE_COMPLEX.Unpack(result,16*Ndata,output0,4,position,MPI::COMM_WORLD);
        sign = output0[0].real()/NerrTot;
        dsign = sqrt(abs(
                        output0[0].imag()/NerrTot
                        - sign*sign
                        )/NerrTot);
        Rvalid = output0[1].real()/NerrTot;
        dRvalid = sqrt(abs(
                        output0[1].imag()/NerrTot
                        - Rvalid*Rvalid
                        )/NerrTot);
        Rproposal = output0[2].real()/NerrTot;
        dRproposal = sqrt(abs(
                        output0[2].imag()/NerrTot
                        - Rproposal*Rproposal
                        )/NerrTot);
        Racceptance = output0[3].real()/NerrTot;
        dRacceptance = sqrt(abs(
                        output0[3].imag()/NerrTot
                        - Racceptance*Racceptance
                        )/NerrTot);
        // selfenergy
        for(int n=0; n<Nordermax; n++)
        for(int iflavor=0; iflavor<Nflavor; iflavor++) for(int jflavor=0; jflavor<Nflavor; jflavor++){
                MPI::DOUBLE_COMPLEX.Unpack(result,16*Ndata,output1,Nx*Nw,position,MPI::COMM_WORLD);
                MPI::DOUBLE_COMPLEX.Unpack(result,16*Ndata,output2,Nx*Nw,position,MPI::COMM_WORLD);
                for(int i=0; i<Nw*Nx; i++){
                        output1[i] /= static_cast<double>(NerrTot);
                        output2[i] = complex<double>(
                                        sqrt(abs(
                                                output2[i].real()/static_cast<double>(NerrTot)
                                                - output1[i].real()*output1[i].real()
                                        )/static_cast<double>(NerrTot)),

                                        sqrt(abs(
                                                output2[i].imag()/static_cast<double>(NerrTot)
                                                - output1[i].imag()*output1[i].imag()
                                        )/static_cast<double>(NerrTot))
                                        );
                }

		                Sw[n][iflavor][jflavor].setOwx(output1);
                dSw[n][iflavor][jflavor].setOwx(output2);

                MPI::DOUBLE_COMPLEX.Unpack(result,16*Ndata,output1,Nx*Nw,position,MPI::COMM_WORLD);
                MPI::DOUBLE_COMPLEX.Unpack(result,16*Ndata,output2,Nx*Nw,position,MPI::COMM_WORLD);
                for(int i=0; i<Nw*Nx; i++){
                        output1[i] /= static_cast<double>(NerrTot);
                        output2[i] = complex<double>(
                                        sqrt(abs(
                                                output2[i].real()/static_cast<double>(NerrTot)
                                                - output1[i].real()*output1[i].real()
                                        )/static_cast<double>(NerrTot)),

                                        sqrt(abs(
                                                output2[i].imag()/static_cast<double>(NerrTot)
                                                - output1[i].imag()*output1[i].imag()
                                        )/static_cast<double>(NerrTot))
                                        );
                }
                Sw[n][iflavor][jflavor].setOwk(output1);
                dSw[n][iflavor][jflavor].setOwk(output2);
        }

        delete[] source;
        delete[] result;
        delete[] outputn;
        delete[] output0;
        delete[] output1;
        delete[] output2;
        delete[] M_Swxacc;
        delete[] M_Swxsqacc;
        delete[] M_Swkacc;
        delete[] M_Swksqacc;
};*/

void controller::sampling_S(){
	const int root = 0;
	const int rank = MPI::COMM_WORLD.Get_rank();
	if(rank==root) cout<<"At controller, Sigma sampling starts"<<endl<<flush;

	if(rank==root) cout<<"accumulation quantities initialization"<<endl;
	double norm_tmp;
	complex<double> norm_acc[Nordermax+1];
	double mistmp;
	complex<double> miscellaneous[6]; 
	// 0: Rproposal, 1: Racceptance, 2: Rvalid, 3: Rirreducible, 4: RGcompact, 5: RWcompact 

	selfenergy_w Swtmp;

	Eigen::MatrixXcd M_Swx, M_Swk;
	Eigen::MatrixXcd ***M_swxacc, ***M_swxsqacc;
	Eigen::MatrixXcd ***M_swkacc, ***M_swksqacc;
	Eigen::MatrixXcd ***M_Swxacc, ***M_Swxsqacc;
	Eigen::MatrixXcd ***M_Swkacc, ***M_Swksqacc;

	M_swxacc = new Eigen::MatrixXcd**[Nordermax];
	M_swkacc = new Eigen::MatrixXcd**[Nordermax];
	M_swxsqacc = new Eigen::MatrixXcd**[Nordermax];
	M_swksqacc = new Eigen::MatrixXcd**[Nordermax];

	M_Swxacc = new Eigen::MatrixXcd**[Nordermax];
	M_Swkacc = new Eigen::MatrixXcd**[Nordermax];
	M_Swxsqacc = new Eigen::MatrixXcd**[Nordermax];
	M_Swksqacc = new Eigen::MatrixXcd**[Nordermax];
	for(int iorder=0; iorder<Nordermax; iorder++){
		M_swxacc[iorder] = new Eigen::MatrixXcd*[Nflavor];
		M_swkacc[iorder] = new Eigen::MatrixXcd*[Nflavor];
		M_swxsqacc[iorder] = new Eigen::MatrixXcd*[Nflavor];
		M_swksqacc[iorder] = new Eigen::MatrixXcd*[Nflavor];

		M_Swxacc[iorder] = new Eigen::MatrixXcd*[Nflavor];
		M_Swkacc[iorder] = new Eigen::MatrixXcd*[Nflavor];
		M_Swxsqacc[iorder] = new Eigen::MatrixXcd*[Nflavor];
		M_Swksqacc[iorder] = new Eigen::MatrixXcd*[Nflavor];
		for(int iflavor=0; iflavor<Nflavor; iflavor++){
			M_swxacc[iorder][iflavor] = new Eigen::MatrixXcd[Nflavor];
			M_swkacc[iorder][iflavor] = new Eigen::MatrixXcd[Nflavor];
			M_swxsqacc[iorder][iflavor] = new Eigen::MatrixXcd[Nflavor];
			M_swksqacc[iorder][iflavor] = new Eigen::MatrixXcd[Nflavor];

			M_Swxacc[iorder][iflavor] = new Eigen::MatrixXcd[Nflavor];
			M_Swkacc[iorder][iflavor] = new Eigen::MatrixXcd[Nflavor];
			M_Swxsqacc[iorder][iflavor] = new Eigen::MatrixXcd[Nflavor];
			M_Swksqacc[iorder][iflavor] = new Eigen::MatrixXcd[Nflavor];
		}
	}

	//Diag.printoutput(global.getdir());
	//printinput(global.getdir());

	if(rank==root) cout<<"initialize accumulation quatities"<<endl;

	diagram Diag;

	for(int i=0; i<6; i++) miscellaneous[i] = 0.;
	for(int iorder=0; iorder<Nordermax; iorder++)
	for(int iflavor=0; iflavor<Nflavor; iflavor++) for(int jflavor=0; jflavor<Nflavor; jflavor++){
		M_Swxacc[iorder][iflavor][jflavor] = Eigen::MatrixXcd::Zero(Nw,Nx);
		M_Swkacc[iorder][iflavor][jflavor] = Eigen::MatrixXcd::Zero(Nw,Nx);
		M_Swxsqacc[iorder][iflavor][jflavor] = Eigen::MatrixXcd::Zero(Nw,Nx);
		M_Swksqacc[iorder][iflavor][jflavor] = Eigen::MatrixXcd::Zero(Nw,Nx);

		M_swxacc[iorder][iflavor][jflavor] = Eigen::MatrixXcd::Zero(Nw,Nx);
		M_swkacc[iorder][iflavor][jflavor] = Eigen::MatrixXcd::Zero(Nw,Nx);
		M_swxsqacc[iorder][iflavor][jflavor] = Eigen::MatrixXcd::Zero(Nw,Nx);
		M_swksqacc[iorder][iflavor][jflavor] = Eigen::MatrixXcd::Zero(Nw,Nx);
	}

	if(rank==root) cout<<"actual sampling process"<<endl;
	Diag.init_S_diagram();
	//Diag.init_P_diagram();
	//Diag.thermalization();
	if(rank==root) cout<<"controller::sampling : "<<flush;
	MPI::COMM_WORLD.Barrier();
	for(int i=0; i<Nerr; i++){
		cout<<"*"<<flush;
		Diag.bold_diagMC_GW_S();
		Diag.normalize_S(0,0,U*n0);

		//cout<<imag(Diag.getSw(1,1,1).getOwx(0,0))<<",";
		//cout<<real(Diag.getPw(2,1,1).getOwx(0,0))<<",";

		//mistmp = Diag.getAvesign();		miscellaneous[0] += complex<double>(mistmp, mistmp*mistmp);
		mistmp = Diag.getRproposal();		miscellaneous[0] += complex<double>(mistmp, mistmp*mistmp);
		mistmp = Diag.getRacceptance();		miscellaneous[1] += complex<double>(mistmp, mistmp*mistmp);
		mistmp = Diag.getRvalid();		miscellaneous[2] += complex<double>(mistmp, mistmp*mistmp);
		mistmp = Diag.getRirreducible();	miscellaneous[3] += complex<double>(mistmp, mistmp*mistmp);
		mistmp = Diag.getRGcompact();		miscellaneous[4] += complex<double>(mistmp, mistmp*mistmp);
		mistmp = Diag.getRWcompact();		miscellaneous[5] += complex<double>(mistmp, mistmp*mistmp);
		for(int n=0; n<Nordermax; n++)
		for(int iflavor=0; iflavor<Nflavor; iflavor++) for(int jflavor=0; jflavor<Nflavor; jflavor++){
			Swtmp = Diag.getSwSum(n+1,iflavor,jflavor);
			M_Swx = Swtmp.getOwx();
			M_Swk = Swtmp.getOwk();
			M_Swxacc[n][iflavor][jflavor] += M_Swx;
			M_Swkacc[n][iflavor][jflavor] += M_Swk;
			M_Swxsqacc[n][iflavor][jflavor] += M_Swx.real().cwiseAbs2() + complex<double>(0.0,1.0)*M_Swx.imag().cwiseAbs2();
			M_Swksqacc[n][iflavor][jflavor] += M_Swk.real().cwiseAbs2() + complex<double>(0.0,1.0)*M_Swk.imag().cwiseAbs2();

			Swtmp = Diag.getSw(n,iflavor,jflavor);
			M_Swx = Swtmp.getOwx();
			M_Swk = Swtmp.getOwk();
			M_swxacc[n][iflavor][jflavor] += M_Swx;
			M_swkacc[n][iflavor][jflavor] += M_Swk;
			M_swxsqacc[n][iflavor][jflavor] += M_Swx.real().cwiseAbs2() + complex<double>(0.0,1.0)*M_Swx.imag().cwiseAbs2();
			M_swksqacc[n][iflavor][jflavor] += M_Swk.real().cwiseAbs2() + complex<double>(0.0,1.0)*M_Swk.imag().cwiseAbs2();
		}
		for(int n=0; n<Nordermax+1; n++){
			norm_tmp = Diag.getnorm(0,n)*(U*n0)/Diag.getnorm(0,0);
			norm_acc[n] += complex<double>(norm_tmp,norm_tmp*norm_tmp);
		}

	}
	if(rank==root) cout<<endl<<flush;

	if(rank==root) cout<<"MPI communication starts"<<endl<<flush;
	if(rank==root) cout<<"MPI Pack"<<endl<<flush;
	
	// (Nordermax+1) : 0) normalization[s]
	// 6 : Rvalid, Rproposal, Racceptance
	// 8 : M_Swxacc, M_Swkacc, M_Swxsqacc, M_Swksqacc (S<-->s)
	int Ndata = (Nordermax+1) + 6 + 8*Nordermax*Nflavor*Nflavor*Nw*Nx;
	complex<double> *source, *result, *outputn, *output0, *output1, *output2;
	source  = new complex<double>[Ndata];
	result  = new complex<double>[Ndata];
	outputn  = new complex<double>[Nordermax+1];
	output0  = new complex<double>[6];
	output1  = new complex<double>[Nw*Nx];
	output2  = new complex<double>[Nw*Nx];

	int position = 0;
	MPI::DOUBLE_COMPLEX.Pack(norm_acc,Nordermax+1,source,16*Ndata,position,MPI::COMM_WORLD);
	MPI::DOUBLE_COMPLEX.Pack(miscellaneous,6,source,16*Ndata,position,MPI::COMM_WORLD);
	for(int n=0; n<Nordermax; n++)
	for(int iflavor=0; iflavor<Nflavor; iflavor++) for(int jflavor=0; jflavor<Nflavor; jflavor++){
		MPI::DOUBLE_COMPLEX.Pack(M_Swxacc[n][iflavor][jflavor].data(),Nx*Nw,source,16*Ndata,position,MPI::COMM_WORLD);
		MPI::DOUBLE_COMPLEX.Pack(M_Swxsqacc[n][iflavor][jflavor].data(),Nx*Nw,source,16*Ndata,position,MPI::COMM_WORLD);
		MPI::DOUBLE_COMPLEX.Pack(M_Swkacc[n][iflavor][jflavor].data(),Nx*Nw,source,16*Ndata,position,MPI::COMM_WORLD);
		MPI::DOUBLE_COMPLEX.Pack(M_Swksqacc[n][iflavor][jflavor].data(),Nx*Nw,source,16*Ndata,position,MPI::COMM_WORLD);

		MPI::DOUBLE_COMPLEX.Pack(M_swxacc[n][iflavor][jflavor].data(),Nx*Nw,source,16*Ndata,position,MPI::COMM_WORLD);
		MPI::DOUBLE_COMPLEX.Pack(M_swxsqacc[n][iflavor][jflavor].data(),Nx*Nw,source,16*Ndata,position,MPI::COMM_WORLD);
		MPI::DOUBLE_COMPLEX.Pack(M_swkacc[n][iflavor][jflavor].data(),Nx*Nw,source,16*Ndata,position,MPI::COMM_WORLD);
		MPI::DOUBLE_COMPLEX.Pack(M_swksqacc[n][iflavor][jflavor].data(),Nx*Nw,source,16*Ndata,position,MPI::COMM_WORLD);
	}

	if(rank==root) cout<<"MPI Reduce and Bcast"<<endl<<flush;
	MPI::COMM_WORLD.Reduce(source,result,Ndata,MPI::DOUBLE_COMPLEX,MPI::SUM,root);
	MPI::COMM_WORLD.Bcast(result,Ndata,MPI::DOUBLE_COMPLEX,root);

	if(rank==root) cout<<"MPI Unpack"<<endl<<flush;
	position = 0;
	if(rank==root) cout<<"normalization"<<endl<<flush;
	MPI::DOUBLE_COMPLEX.Unpack(result,16*Ndata,outputn,Nordermax+1,position,MPI::COMM_WORLD);
	for(int iorder=0; iorder<Nordermax+1; iorder++){
		normalization[0][iorder] = outputn[iorder].real()/NerrTot;
		dnormalization[0][iorder] = 
			sqrt(abs(
				outputn[iorder].imag()/NerrTot
				- normalization[0][iorder]*normalization[0][iorder]
			)/NerrTot);
	}
	if(rank==root) cout<<"miscellaneous"<<endl<<flush;
	MPI::DOUBLE_COMPLEX.Unpack(result,16*Ndata,output0,6,position,MPI::COMM_WORLD);
	Rproposal = output0[0].real()/NerrTot;
	dRproposal = sqrt(abs(
			output0[0].imag()/NerrTot
			- Rproposal*Rproposal
			)/NerrTot);
	Racceptance = output0[1].real()/NerrTot;
	dRacceptance = sqrt(abs(
			output0[1].imag()/NerrTot
			- Racceptance*Racceptance
			)/NerrTot);
	Rvalid = output0[2].real()/NerrTot;
	dRvalid = sqrt(abs(
			output0[2].imag()/NerrTot
			- Rvalid*Rvalid
			)/NerrTot);
	Rirreducible = output0[3].real()/NerrTot;
	dRirreducible = sqrt(abs(
			output0[3].imag()/NerrTot
			- Rirreducible*Rirreducible
			)/NerrTot);
	RGcompact = output0[4].real()/NerrTot;
	dRGcompact = sqrt(abs(
			output0[4].imag()/NerrTot
			- RGcompact*RGcompact
			)/NerrTot);
	RWcompact = output0[5].real()/NerrTot;
	dRWcompact = sqrt(abs(
			output0[5].imag()/NerrTot
			- RWcompact*RWcompact
			)/NerrTot);

	if(rank==root) cout<<"selfenergy"<<endl<<flush;
	for(int n=0; n<Nordermax; n++)
	for(int iflavor=0; iflavor<Nflavor; iflavor++) for(int jflavor=0; jflavor<Nflavor; jflavor++){
		MPI::DOUBLE_COMPLEX.Unpack(result,16*Ndata,output1,Nx*Nw,position,MPI::COMM_WORLD);
		MPI::DOUBLE_COMPLEX.Unpack(result,16*Ndata,output2,Nx*Nw,position,MPI::COMM_WORLD);
		for(int i=0; i<Nw*Nx; i++){
			output1[i] /= static_cast<double>(NerrTot);
			output2[i] = complex<double>(
					sqrt(abs(
						output2[i].real()/static_cast<double>(NerrTot) 
						- output1[i].real()*output1[i].real()
					)/static_cast<double>(NerrTot)),

					sqrt(abs(
						output2[i].imag()/static_cast<double>(NerrTot) 
						- output1[i].imag()*output1[i].imag()
					)/static_cast<double>(NerrTot))
					);
		}
		Sw[n][iflavor][jflavor].setOwx(output1);
		dSw[n][iflavor][jflavor].setOwx(output2);


		MPI::DOUBLE_COMPLEX.Unpack(result,16*Ndata,output1,Nx*Nw,position,MPI::COMM_WORLD);
		MPI::DOUBLE_COMPLEX.Unpack(result,16*Ndata,output2,Nx*Nw,position,MPI::COMM_WORLD);
		for(int i=0; i<Nw*Nx; i++){
			output1[i] /= static_cast<double>(NerrTot);
			output2[i] = complex<double>(
					sqrt(abs(
						output2[i].real()/static_cast<double>(NerrTot) 
						- output1[i].real()*output1[i].real()
					)/static_cast<double>(NerrTot)),

					sqrt(abs(
						output2[i].imag()/static_cast<double>(NerrTot) 
						- output1[i].imag()*output1[i].imag()
					)/static_cast<double>(NerrTot))
					);
		}
		Sw[n][iflavor][jflavor].setOwk(output1);
		dSw[n][iflavor][jflavor].setOwk(output2);

		MPI::DOUBLE_COMPLEX.Unpack(result,16*Ndata,output1,Nx*Nw,position,MPI::COMM_WORLD);
		MPI::DOUBLE_COMPLEX.Unpack(result,16*Ndata,output2,Nx*Nw,position,MPI::COMM_WORLD);
		for(int i=0; i<Nw*Nx; i++){
			output1[i] /= static_cast<double>(NerrTot);
			output2[i] = complex<double>(
					sqrt(abs(
						output2[i].real()/static_cast<double>(NerrTot) 
						- output1[i].real()*output1[i].real()
					)/static_cast<double>(NerrTot)),

					sqrt(abs(
						output2[i].imag()/static_cast<double>(NerrTot) 
						- output1[i].imag()*output1[i].imag()
					)/static_cast<double>(NerrTot))
					);
		}
		sw[n][iflavor][jflavor].setOwx(output1);
		dsw[n][iflavor][jflavor].setOwx(output2);


		MPI::DOUBLE_COMPLEX.Unpack(result,16*Ndata,output1,Nx*Nw,position,MPI::COMM_WORLD);
		MPI::DOUBLE_COMPLEX.Unpack(result,16*Ndata,output2,Nx*Nw,position,MPI::COMM_WORLD);
		for(int i=0; i<Nw*Nx; i++){
			output1[i] /= static_cast<double>(NerrTot);
			output2[i] = complex<double>(
					sqrt(abs(
						output2[i].real()/static_cast<double>(NerrTot) 
						- output1[i].real()*output1[i].real()
					)/static_cast<double>(NerrTot)),

					sqrt(abs(
						output2[i].imag()/static_cast<double>(NerrTot) 
						- output1[i].imag()*output1[i].imag()
					)/static_cast<double>(NerrTot))
					);
		}
		sw[n][iflavor][jflavor].setOwk(output1);
		dsw[n][iflavor][jflavor].setOwk(output2);
	}

	if(rank==root) cout<<"memory deallocation"<<endl<<flush;

	delete[] source;
	delete[] result;

	delete[] outputn;
	delete[] output0;
	delete[] output1;
	delete[] output2;

	for(int iorder=0; iorder<Nordermax; iorder++){
		for(int iflavor=0; iflavor<Nflavor; iflavor++){
			delete[] M_Swxacc[iorder][iflavor];
			delete[] M_Swxsqacc[iorder][iflavor];
			delete[] M_Swkacc[iorder][iflavor];
			delete[] M_Swksqacc[iorder][iflavor];

			delete[] M_swxacc[iorder][iflavor];
			delete[] M_swxsqacc[iorder][iflavor];
			delete[] M_swkacc[iorder][iflavor];
			delete[] M_swksqacc[iorder][iflavor];
		}
		delete[] M_Swxacc[iorder];
		delete[] M_Swxsqacc[iorder];
		delete[] M_Swkacc[iorder];
		delete[] M_Swksqacc[iorder];

		delete[] M_swxacc[iorder];
		delete[] M_swxsqacc[iorder];
		delete[] M_swkacc[iorder];
		delete[] M_swksqacc[iorder];
	}
	delete[] M_Swxacc;
	delete[] M_Swxsqacc;
	delete[] M_Swkacc;
	delete[] M_Swksqacc;

	delete[] M_swxacc;
	delete[] M_swxsqacc;
	delete[] M_swkacc;
	delete[] M_swksqacc;

	if(rank==root) cout<<"At controller, sampling ends"<<endl<<flush;
};

/*void controller::sampling(){
	//cout<<"At controller, sampling starts"<<endl<<flush;
	const int root = 0;
	const int rank = MPI::COMM_WORLD.Get_rank();

	//cout<<"accumulation quantities initialization"<<endl;
	double norm_tmp;
	complex<double> norm_acc[Nordermax];
	double mistmp;
	complex<double> miscellaneous[4]; // 0: sign, 1: Rvalid, 2: Rproposal, 3: Racceptance
	for(int i=0; i<4; i++) miscellaneous[i] = 0.;

	selfenergy_w Swtmp;
	Eigen::MatrixXcd M_Swx, M_Swk;
	Eigen::MatrixXcd ***M_Swxacc, ***M_Swxsqacc;
	Eigen::MatrixXcd ***M_Swkacc, ***M_Swksqacc;
	M_Swxacc = new Eigen::MatrixXcd**[Nordermax];
	M_Swkacc = new Eigen::MatrixXcd**[Nordermax];
	M_Swxsqacc = new Eigen::MatrixXcd**[Nordermax];
	M_Swksqacc = new Eigen::MatrixXcd**[Nordermax];
	for(int iorder=0; iorder<Nordermax; iorder++){
		M_Swxacc[iorder] = new Eigen::MatrixXcd*[Nflavor];
		M_Swkacc[iorder] = new Eigen::MatrixXcd*[Nflavor];
		M_Swxsqacc[iorder] = new Eigen::MatrixXcd*[Nflavor];
		M_Swksqacc[iorder] = new Eigen::MatrixXcd*[Nflavor];
		for(int iflavor=0; iflavor<Nflavor; iflavor++){
			M_Swxacc[iorder][iflavor] = new Eigen::MatrixXcd[Nflavor];
			M_Swkacc[iorder][iflavor] = new Eigen::MatrixXcd[Nflavor];
			M_Swxsqacc[iorder][iflavor] = new Eigen::MatrixXcd[Nflavor];
			M_Swksqacc[iorder][iflavor] = new Eigen::MatrixXcd[Nflavor];
			for(int jflavor=0; jflavor<Nflavor; jflavor++){
				M_Swxacc[iorder][iflavor][jflavor] = Eigen::MatrixXcd::Zero(Nw,Nx);
				M_Swkacc[iorder][iflavor][jflavor] = Eigen::MatrixXcd::Zero(Nw,Nx);
				M_Swxsqacc[iorder][iflavor][jflavor] = Eigen::MatrixXcd::Zero(Nw,Nx);
				M_Swksqacc[iorder][iflavor][jflavor] = Eigen::MatrixXcd::Zero(Nw,Nx);
			}
		}
	}

	//cout<<"actual sampling process"<<endl;
	diagram Diag;
	Diag.init_S_diagram();
	if(rank==root) cout<<"controller::sampling : "<<flush;
	MPI::COMM_WORLD.Barrier();
	for(int i=0; i<Nerr; i++){
		cout<<"*"<<flush;
		Diag.bold_diagMC_G();
		Diag.normalize_S(1,U*n0);

		mistmp = Diag.getAvesign();	miscellaneous[0] += complex<double>(mistmp, mistmp*mistmp);
		mistmp = Diag.getRvalid();	miscellaneous[1] += complex<double>(mistmp, mistmp*mistmp);
		mistmp = Diag.getRproposal();	miscellaneous[2] += complex<double>(mistmp, mistmp*mistmp);
		mistmp = Diag.getRacceptance();	miscellaneous[3] += complex<double>(mistmp, mistmp*mistmp);
		for(int n=0; n<Nordermax; n++){
			norm_tmp = Diag.getnorm(n)*(U*n0)/Diag.getnorm(0);
			norm_acc[n] += complex<double>(norm_tmp,norm_tmp*norm_tmp);
			for(int iflavor=0; iflavor<Nflavor; iflavor++) for(int jflavor=0; jflavor<Nflavor; jflavor++){
				Swtmp = Diag.getSw(n,iflavor,jflavor);
				M_Swx = Swtmp.getOwx();
				M_Swk = Swtmp.getOwk();
				M_Swxacc[n][iflavor][jflavor] += M_Swx;
				M_Swkacc[n][iflavor][jflavor] += M_Swk;
				M_Swxsqacc[n][iflavor][jflavor] += M_Swx.real().cwiseAbs2() + complex<double>(0.0,1.0)*M_Swx.imag().cwiseAbs2();
				M_Swksqacc[n][iflavor][jflavor] += M_Swk.real().cwiseAbs2() + complex<double>(0.0,1.0)*M_Swk.imag().cwiseAbs2();
				//if(n==1 && iflavor==0 && jflavor==0){ cout<<"[("<<M_Swx(0,0)<<"),"; }
				//if(n==1 && iflavor==0 && jflavor==0){ cout<<M_Swx(0,0).real()<<","; }
			}
		}
		//cout<<(M_Swx.real().cwiseAbs2() + complex<double>(0.0,1.0)*M_Swx.imag().cwiseAbs2())<<endl;
	}
	if(rank==root) cout<<endl<<flush;

	//cout<<"MPI communication starts"<<endl<<flush;
	//cout<<"MPI Pack"<<endl<<flush;
	
	// Nordermax : normalization
	// 4 : sign, Rvalid, Rproposal, Racceptance
	// 4 : M_Swxacc, M_Swkacc, M_Swxsqacc, M_Swksqacc
	int Ndata = Nordermax + 4 + 4*Nordermax*Nflavor*Nflavor*Nw*Nx;
	complex<double> *source, *result, *outputn, *output0, *output1, *output2;
	source  = new complex<double>[Ndata];
	result  = new complex<double>[Ndata];
	outputn  = new complex<double>[Nordermax];
	output0  = new complex<double>[4];
	output1  = new complex<double>[Nw*Nx];
	output2  = new complex<double>[Nw*Nx];

	int position = 0;
	MPI::DOUBLE_COMPLEX.Pack(norm_acc,Nordermax,source,16*Ndata,position,MPI::COMM_WORLD);
	MPI::DOUBLE_COMPLEX.Pack(miscellaneous,4,source,16*Ndata,position,MPI::COMM_WORLD);
	for(int n=0; n<Nordermax; n++)
	for(int iflavor=0; iflavor<Nflavor; iflavor++) for(int jflavor=0; jflavor<Nflavor; jflavor++){
		MPI::DOUBLE_COMPLEX.Pack(M_Swxacc[n][iflavor][jflavor].data(),Nx*Nw,source,16*Ndata,position,MPI::COMM_WORLD);
		MPI::DOUBLE_COMPLEX.Pack(M_Swxsqacc[n][iflavor][jflavor].data(),Nx*Nw,source,16*Ndata,position,MPI::COMM_WORLD);
		MPI::DOUBLE_COMPLEX.Pack(M_Swkacc[n][iflavor][jflavor].data(),Nx*Nw,source,16*Ndata,position,MPI::COMM_WORLD);
		MPI::DOUBLE_COMPLEX.Pack(M_Swksqacc[n][iflavor][jflavor].data(),Nx*Nw,source,16*Ndata,position,MPI::COMM_WORLD);
	}

	//cout<<"MPI Reduce and Bcast"<<endl<<flush;
	MPI::COMM_WORLD.Reduce(source,result,Ndata,MPI::DOUBLE_COMPLEX,MPI::SUM,root);
	MPI::COMM_WORLD.Bcast(result,Ndata,MPI::DOUBLE_COMPLEX,root);

	// cout<<"MPI Unpack"<<endl<<flush;
	position = 0;
	// normalization
	MPI::DOUBLE_COMPLEX.Unpack(result,16*Ndata,outputn,Nordermax,position,MPI::COMM_WORLD);
	for(int iorder=0; iorder<Nordermax; iorder++){
		normalization[iorder] = outputn[iorder].real()/NerrTot;
		dnormalization[iorder] = 
			sqrt(abs(
				outputn[iorder].imag()/NerrTot
				- normalization[iorder]*normalization[iorder]
			)/NerrTot);
	};
	// miscellaneous
	MPI::DOUBLE_COMPLEX.Unpack(result,16*Ndata,output0,4,position,MPI::COMM_WORLD);
	sign = output0[0].real()/NerrTot;
	dsign = sqrt(abs(
			output0[0].imag()/NerrTot
			- sign*sign
			)/NerrTot);
	Rvalid = output0[1].real()/NerrTot;
	dRvalid = sqrt(abs(
			output0[1].imag()/NerrTot
			- Rvalid*Rvalid
			)/NerrTot);
	Rproposal = output0[2].real()/NerrTot;
	dRproposal = sqrt(abs(
			output0[2].imag()/NerrTot
			- Rproposal*Rproposal
			)/NerrTot);
	Racceptance = output0[3].real()/NerrTot;
	dRacceptance = sqrt(abs(
			output0[3].imag()/NerrTot
			- Racceptance*Racceptance
			)/NerrTot);
	// selfenergy
	for(int n=0; n<Nordermax; n++)
	for(int iflavor=0; iflavor<Nflavor; iflavor++) for(int jflavor=0; jflavor<Nflavor; jflavor++){
		MPI::DOUBLE_COMPLEX.Unpack(result,16*Ndata,output1,Nx*Nw,position,MPI::COMM_WORLD);
		MPI::DOUBLE_COMPLEX.Unpack(result,16*Ndata,output2,Nx*Nw,position,MPI::COMM_WORLD);
		for(int i=0; i<Nw*Nx; i++){
			output1[i] /= static_cast<double>(NerrTot);
			output2[i] = complex<double>(
					sqrt(abs(
						output2[i].real()/static_cast<double>(NerrTot) 
						- output1[i].real()*output1[i].real()
					)/static_cast<double>(NerrTot)),

					sqrt(abs(
						output2[i].imag()/static_cast<double>(NerrTot) 
						- output1[i].imag()*output1[i].imag()
					)/static_cast<double>(NerrTot))
					);
		}
		Sw[n][iflavor][jflavor].setOwx(output1);
		dSw[n][iflavor][jflavor].setOwx(output2);

		MPI::DOUBLE_COMPLEX.Unpack(result,16*Ndata,output1,Nx*Nw,position,MPI::COMM_WORLD);
		MPI::DOUBLE_COMPLEX.Unpack(result,16*Ndata,output2,Nx*Nw,position,MPI::COMM_WORLD);
		for(int i=0; i<Nw*Nx; i++){
			output1[i] /= static_cast<double>(NerrTot);
			output2[i] = complex<double>(
					sqrt(abs(
						output2[i].real()/static_cast<double>(NerrTot) 
						- output1[i].real()*output1[i].real()
					)/static_cast<double>(NerrTot)),

					sqrt(abs(
						output2[i].imag()/static_cast<double>(NerrTot) 
						- output1[i].imag()*output1[i].imag()
					)/static_cast<double>(NerrTot))
					);
		}
		Sw[n][iflavor][jflavor].setOwk(output1);
		dSw[n][iflavor][jflavor].setOwk(output2);
	}

	delete[] source;
	delete[] result;
	delete[] outputn;
	delete[] output0;
	delete[] output1;
	delete[] output2;
	delete[] M_Swxacc;
	delete[] M_Swxsqacc;
	delete[] M_Swkacc;
	delete[] M_Swksqacc;
};*/
/*void controller::iterativesampling(){
	//cout<<"At controller, sampling starts"<<endl<<flush;
	const int root = 0;
	const int rank = MPI::COMM_WORLD.Get_rank();

	//cout<<"accumulation quantities initialization"<<endl;
	double norm_tmp;
	complex<double> norm_acc;
	double mistmp;
	complex<double> miscellaneous[4]; // 0: sign, 1: Rvalid, 2: Rproposal, 3: Racceptance

	function_w Owtmp;
	Eigen::MatrixXcd M_Owx, M_Owk;
	Eigen::MatrixXcd **M_Owxacc, **M_Owxsqacc;
	Eigen::MatrixXcd **M_Owkacc, **M_Owksqacc;

	M_Owxacc = new Eigen::MatrixXcd*[Nflavor];
	M_Owkacc = new Eigen::MatrixXcd*[Nflavor];
	M_Owxsqacc = new Eigen::MatrixXcd*[Nflavor];
	M_Owksqacc = new Eigen::MatrixXcd*[Nflavor];
	for(int iflavor=0; iflavor<Nflavor; iflavor++){
		M_Owxacc[iflavor] = new Eigen::MatrixXcd[Nflavor];
		M_Owkacc[iflavor] = new Eigen::MatrixXcd[Nflavor];
		M_Owxsqacc[iflavor] = new Eigen::MatrixXcd[Nflavor];
		M_Owksqacc[iflavor] = new Eigen::MatrixXcd[Nflavor];
	}

	// 5 : normalization, sign, Rvalid, Rproposal, Racceptance
	// 4 : M_Owxacc, M_Owkacc, M_Owxsqacc, M_Owksqacc
	int Ndata = 5 + 4*Nflavor*Nflavor*Nw*Nx;
	complex<double> *source, *result, *output0, *output1, *output2;
	source  = new complex<double>[Ndata];
	result  = new complex<double>[Ndata];
	output0  = new complex<double>[5];
	output1  = new complex<double>[Nw*Nx];
	output2  = new complex<double>[Nw*Nx];

	int position;

	diagram Diag;
	for(int iordermax=2; iordermax<Nordermax+1; iordermax++){
		norm_acc = 0.;
		for(int i=0; i<4; i++) miscellaneous[i] = 0.;
		for(int iflavor=0; iflavor<Nflavor; iflavor++) for(int jflavor=0; jflavor<Nflavor; jflavor++){
			M_Owxacc[iflavor][jflavor] = Eigen::MatrixXcd::Zero(Nw,Nx);
			M_Owkacc[iflavor][jflavor] = Eigen::MatrixXcd::Zero(Nw,Nx);
			M_Owxsqacc[iflavor][jflavor] = Eigen::MatrixXcd::Zero(Nw,Nx);
			M_Owksqacc[iflavor][jflavor] = Eigen::MatrixXcd::Zero(Nw,Nx);
		}

		Diag.init_P_diagram();
		Diag.setNordermax(iordermax);

		MPI::COMM_WORLD.Barrier();
		if(rank==root) cout<<"controller::"<<iordermax<<"th-order sampling : "<<flush;
		MPI::COMM_WORLD.Barrier();
		for(int i=0; i<Nerr; i++){
			// sampling
			cout<<"*"<<flush;
			Diag.diagMC_P(iordermax);
			//Diag.normalize_P(iordermax-1,normalization[iordermax-2],iordermax);

			// accumulation
			//norm_tmp = Diag.getnorm(iordermax-1)*normalization[iordermax-2]/Diag.getnorm(iordermax-2);
			//norm_acc += complex<double>(norm_tmp,norm_tmp*norm_tmp);
			//mistmp = Diag.getAvesign();	miscellaneous[0] += complex<double>(mistmp, mistmp*mistmp);
			//mistmp = Diag.getRvalid();	miscellaneous[1] += complex<double>(mistmp, mistmp*mistmp);
			//mistmp = Diag.getRproposal();	miscellaneous[2] += complex<double>(mistmp, mistmp*mistmp);
			//mistmp = Diag.getRacceptance();	miscellaneous[3] += complex<double>(mistmp, mistmp*mistmp);

			//for(int iflavor=0; iflavor<Nflavor; iflavor++) for(int jflavor=0; jflavor<Nflavor; jflavor++){
			//	Owtmp = Diag.getPw(iordermax,iflavor,jflavor);
			//	M_Owx = Owtmp.getOwx();
			//	M_Owk = Owtmp.getOwk();
			//	M_Owxacc[iflavor][jflavor] += M_Owx;
			//	M_Owkacc[iflavor][jflavor] += M_Owk;
			//	M_Owxsqacc[iflavor][jflavor] += M_Owx.real().cwiseAbs2() + complex<double>(0.0,1.0)*M_Owx.imag().cwiseAbs2();
			//	M_Owksqacc[iflavor][jflavor] += M_Owk.real().cwiseAbs2() + complex<double>(0.0,1.0)*M_Owk.imag().cwiseAbs2();

			//}
		}
		MPI::COMM_WORLD.Barrier();
		if(rank==root) cout<<"|"<<endl<<flush;

		//for(int i=0; i<Nerr; i++){
		//	Diag.init_S_diagram();
		//	// sampling
		//	cout<<"*"<<flush;
		//	Diag.diagMC_S(iordermax);
		//	Diag.normalize_S(iordermax-1,normalization[iordermax-2],iordermax);

		//	// accumulation
		//	norm_tmp = Diag.getnorm(iordermax-1)*normalization[iordermax-2]/Diag.getnorm(iordermax-2);
		//	norm_acc += complex<double>(norm_tmp,norm_tmp*norm_tmp);
		//	mistmp = Diag.getAvesign();	miscellaneous[0] += complex<double>(mistmp, mistmp*mistmp);
		//	mistmp = Diag.getRvalid();	miscellaneous[1] += complex<double>(mistmp, mistmp*mistmp);
		//	mistmp = Diag.getRproposal();	miscellaneous[2] += complex<double>(mistmp, mistmp*mistmp);
		//	mistmp = Diag.getRacceptance();	miscellaneous[3] += complex<double>(mistmp, mistmp*mistmp);

		//	for(int iflavor=0; iflavor<Nflavor; iflavor++) for(int jflavor=0; jflavor<Nflavor; jflavor++){
		//		Owtmp = Diag.getSw(iordermax-1,iflavor,jflavor);
		//		M_Owx = Owtmp.getOwx();
		//		M_Owk = Owtmp.getOwk();
		//		M_Owxacc[iflavor][jflavor] += M_Owx;
		//		M_Owkacc[iflavor][jflavor] += M_Owk;
		//		M_Owxsqacc[iflavor][jflavor] += M_Owx.real().cwiseAbs2() + complex<double>(0.0,1.0)*M_Owx.imag().cwiseAbs2();
		//		M_Owksqacc[iflavor][jflavor] += M_Owk.real().cwiseAbs2() + complex<double>(0.0,1.0)*M_Owk.imag().cwiseAbs2();

		//	}

		//	//cout<<"("<<norm_tmp<<")";
		//}
		//MPI::COMM_WORLD.Barrier();
		//if(rank==root) cout<<"|"<<endl<<flush;

		//cout<<"MPI communication starts"<<endl<<flush;
		//cout<<"MPI Pack"<<endl<<flush;
		position = 0;
		MPI::DOUBLE_COMPLEX.Pack(&norm_acc,1,source,16*Ndata,position,MPI::COMM_WORLD);
		MPI::DOUBLE_COMPLEX.Pack(miscellaneous,4,source,16*Ndata,position,MPI::COMM_WORLD);
		for(int iflavor=0; iflavor<Nflavor; iflavor++) for(int jflavor=0; jflavor<Nflavor; jflavor++){
			MPI::DOUBLE_COMPLEX.Pack(M_Owxacc[iflavor][jflavor].data(),Nx*Nw,source,16*Ndata,position,MPI::COMM_WORLD);
			MPI::DOUBLE_COMPLEX.Pack(M_Owxsqacc[iflavor][jflavor].data(),Nx*Nw,source,16*Ndata,position,MPI::COMM_WORLD);
			MPI::DOUBLE_COMPLEX.Pack(M_Owkacc[iflavor][jflavor].data(),Nx*Nw,source,16*Ndata,position,MPI::COMM_WORLD);
			MPI::DOUBLE_COMPLEX.Pack(M_Owksqacc[iflavor][jflavor].data(),Nx*Nw,source,16*Ndata,position,MPI::COMM_WORLD);
		}

		//cout<<"MPI Reduce and Bcast"<<endl<<flush;
		MPI::COMM_WORLD.Reduce(source,result,Ndata,MPI::DOUBLE_COMPLEX,MPI::SUM,root);
		MPI::COMM_WORLD.Bcast(result,Ndata,MPI::DOUBLE_COMPLEX,root);

		//cout<<"MPI Unpack"<<endl<<flush;
		position = 0;
		MPI::DOUBLE_COMPLEX.Unpack(result,16*Ndata,output0,5,position,MPI::COMM_WORLD);

		output0[0].real(output0[0].real()/static_cast<double>(NerrTot));
		normalization[iordermax-1] = output0[0].real();
		dnormalization[iordermax-1] = sqrt(abs(
						output0[0].imag()/static_cast<double>(NerrTot) 
						- output0[0].real()*output0[0].real()
						)/static_cast<double>(NerrTot));
		sign = output0[1].real()/NerrTot;
		dsign = sqrt(abs(
				output0[1].imag()/NerrTot
				- sign*sign
				)/NerrTot);
		Rvalid = output0[2].real()/NerrTot;
		dRvalid = sqrt(abs(
				output0[2].imag()/NerrTot
				- Rvalid*Rvalid
				)/NerrTot);
		Rproposal = output0[3].real()/NerrTot;
		dRproposal = sqrt(abs(
				output0[3].imag()/NerrTot
				- Rproposal*Rproposal
				)/NerrTot);
		Racceptance = output0[4].real()/NerrTot;
		dRacceptance = sqrt(abs(
				output0[4].imag()/NerrTot
				- Racceptance*Racceptance
				)/NerrTot);

		for(int iflavor=0; iflavor<Nflavor; iflavor++) for(int jflavor=0; jflavor<Nflavor; jflavor++){
			MPI::DOUBLE_COMPLEX.Unpack(result,16*Ndata,output1,Nx*Nw,position,MPI::COMM_WORLD);
			MPI::DOUBLE_COMPLEX.Unpack(result,16*Ndata,output2,Nx*Nw,position,MPI::COMM_WORLD);
			for(int i=0; i<Nw*Nx; i++){
				output1[i] /= static_cast<double>(NerrTot);
				output2[i] = complex<double>(
						sqrt(abs(
							output2[i].real()/static_cast<double>(NerrTot) 
							- output1[i].real()*output1[i].real()
						)/static_cast<double>(NerrTot)),

						sqrt(abs(
							output2[i].imag()/static_cast<double>(NerrTot) 
							- output1[i].imag()*output1[i].imag()
						)/static_cast<double>(NerrTot))
						);
			}
			Sw[iordermax-1][iflavor][jflavor].setOwx(output1);
			dSw[iordermax-1][iflavor][jflavor].setOwx(output2);
			//if(n==1 && iflavor==0 && jflavor==0) cout<<output1[0]<<" "<<output2[0]<<endl;

			MPI::DOUBLE_COMPLEX.Unpack(result,16*Ndata,output1,Nx*Nw,position,MPI::COMM_WORLD);
			MPI::DOUBLE_COMPLEX.Unpack(result,16*Ndata,output2,Nx*Nw,position,MPI::COMM_WORLD);
			for(int i=0; i<Nw*Nx; i++){
				output1[i] /= static_cast<double>(NerrTot);
				output2[i] = complex<double>(
						sqrt(abs(
							output2[i].real()/static_cast<double>(NerrTot) 
							- output1[i].real()*output1[i].real()
						)/static_cast<double>(NerrTot)),

						sqrt(abs(
							output2[i].imag()/static_cast<double>(NerrTot) 
							- output1[i].imag()*output1[i].imag()
						)/static_cast<double>(NerrTot))
						);
			}
			Sw[iordermax-1][iflavor][jflavor].setOwk(output1);
			dSw[iordermax-1][iflavor][jflavor].setOwk(output2);
		}
		//construct_Gw(iordermax);

		//printiterativeoutput(parameters::global.getdir(), iordermax);

	}

	delete[] source;
	delete[] result;
	delete[] output1;
	delete[] output2;
	delete[] M_Owxacc;
	delete[] M_Owxsqacc;
	delete[] M_Owkacc;
	delete[] M_Owksqacc;
};*/
void controller::fourier_Gw2t(){
        Eigen::MatrixXcd Gwxtmp;
        Eigen::MatrixXd Gtxtmp;
        for(int iflavor=0; iflavor<Nflavor; iflavor++){
                Gwxtmp = Gw[Nordermax-1][iflavor][iflavor].getOwx();
                for(int n=0; n<Nw; n++) Gwxtmp(n,0) -= complex<double>(0.0,-1.0/w(n));
                Gtxtmp = FermionicInverseTimeFourier(Gwxtmp);
                for(int it=0; it<Ntime; it++) Gtxtmp(it,0) += -0.5;

                Gtxtmp.conservativeResize(Ntime+1,Nx);
                Gtxtmp(Ntime,0) =  -1.0 - Gtxtmp(0,0);
                for(int xi=1; xi<Nx; xi++) Gtxtmp(Ntime,xi) = -Gtxtmp(0,xi);

                Gt[0][iflavor][iflavor].setOtx( Gtxtmp );
        }
};
void controller::fourier_Gt2w(){
        Eigen::MatrixXcd Gwxtmp;
        Eigen::MatrixXd Gtxtmp;
        for(int iflavor=0; iflavor<Nflavor; iflavor++){
                Gtxtmp = Gt[0][iflavor][iflavor].getOtx();
                for(int ti=0; ti<Ntime; ti++) Gtxtmp(ti,0) -= -0.5;
                Gwxtmp = FermionicTimeFourier(Gtxtmp);
                for(int n=0; n<Nw; n++) Gwxtmp(n,0) += complex<double>(0.0,-1.0/w(n));

                Gw[Nordermax-1][iflavor][iflavor].setOwx( Gwxtmp );
                Gw[Nordermax-1][iflavor][iflavor].setOwk( Gwxtmp );
        }
};

void controller::selfconsistentloop(const double& alpha){
        const int root = 0;
        const int rank = MPI::COMM_WORLD.Get_rank();
        for(int iscl=0; iscl<Nscl; iscl++){
                if(rank==root) cout<<"** "<<iscl<<"th iteration **"<<endl<<flush;
                sampling();
                if(rank==root) cout<<" * applying Dyson equation *"<<endl<<flush;
                dyson_equation();
                if(rank==root) cout<<" * update freq. Green function *: "<<flush;
                update_Gw(alpha);
                if(rank==root) cout<<" * fourier transform time. Green function *"<<endl<<flush;
                fourier_Gw2t();
                if(rank==root) cout<<" * print iteration *"<<endl<<flush;
                if(rank==root) printoutput(global.getdir(),iscl);
        }
};
void controller::selfconsistentloop(const double& lambda, const double& alpha){
        const int root = 0;
        const int rank = MPI::COMM_WORLD.Get_rank();
        for(int iscl=0; iscl<Nscl; iscl++){
                if(rank==root) cout<<"** "<<iscl<<"th iteration **"<<endl<<flush;
                sampling();
                if(rank==root) cout<<" * applying modified Dyson equation *"<<endl<<flush;
                dyson_modified(lambda);
                if(rank==root) cout<<" * update freq. Green function *"<<endl<<flush;
                update_Gw(alpha);
                if(rank==root) cout<<" * fourier transform time. Green function *"<<endl<<flush;
                fourier_Gw2t();
                if(rank==root) cout<<" * print iteration *"<<endl<<flush;
                if(rank==root) printoutput(global.getdir(),iscl);
        }
};
void controller::construct_Gw(const int& iorder){
	//cout<<"* construct_Gw module *"<<endl<<flush;
	green_w Gw_tmp;
	Gw_tmp.init();
	for(int iflavor=0; iflavor<Nflavor; iflavor++){
		Gw[iorder-1][iflavor][iflavor].init();
		for(auto& it0: order_combination[iorder-1]){
			Gw_tmp = Gw[0][iflavor][iflavor];
			for(auto& it1: it0) Gw_tmp *= Sw[(it1)-1][iflavor][iflavor]*Gw[0][iflavor][iflavor];
			Gw_tmp.fourier_ktox();
			Gw[iorder-1][iflavor][iflavor] += Gw_tmp*static_cast<double>(count_combination(it0));
		}
	}
	//cout<<"***********************"<<endl<<flush;
};
void controller::printstatistics(std::ostream& ostr){
	ostr<<setw(25)<<left<<"# Rproposal"
		<<setw(25)<<left<<"dRproposal"
		<<setw(25)<<left<<"Racceptance"
		<<setw(25)<<left<<"dRacceptance"
		<<setw(25)<<left<<"Rvalid"
		<<setw(25)<<left<<"dRvalid"
		<<setw(25)<<left<<"Rirreducible"
		<<setw(25)<<left<<"dRirreducible"
		<<setw(25)<<left<<"RGcompact"
		<<setw(25)<<left<<"dRGcompact"
		<<setw(25)<<left<<"RWcompact"
		<<setw(25)<<left<<"dRWcompact"<<endl;
		
	ostr<<setw(25)<<left<<Rproposal
		<<setw(25)<<left<<dRproposal
		<<setw(25)<<left<<Racceptance
		<<setw(25)<<left<<dRacceptance
		<<setw(25)<<left<<Rvalid
		<<setw(25)<<left<<dRvalid
		<<setw(25)<<left<<Rirreducible
		<<setw(25)<<left<<dRirreducible
		<<setw(25)<<left<<RGcompact
		<<setw(25)<<left<<dRGcompact
		<<setw(25)<<left<<RWcompact
		<<setw(25)<<left<<dRWcompact<<endl;
};
void controller::printnormalization(const int& sp, std::ostream& ostr){
	ostr<<setw(25)<<left<<"# order"
		<<setw(25)<<right<<"normalizaed_value"
		<<setw(25)<<right<<"dnormalizaed_value"<<endl;
	for(int i=0; i<normalization[sp].size(); i++)
		ostr<<setw(25)<<left<<i
			<<setw(25)<<right<<scientific<<normalization[sp][i]
			<<setw(25)<<right<<scientific<<dnormalization[sp][i]<<endl;
};
void controller::printPloc(std::ostream& ostr){
	ostr<<setw(25)<<left<<"# order"
		<<setw(25)<<right<<"Ploc_value"
		<<setw(25)<<right<<"dPloc_value"<<endl;
	for(int i=0; i<Ploc.size(); i++)
		ostr<<setw(25)<<left<<i
			<<setw(25)<<right<<scientific<<Ploc[i]
			<<setw(25)<<right<<scientific<<dPloc[i]<<endl;
};
void controller::printSwx(ostream& ostr, const int& iorder, const int& iflavor, const int& jflavor){
	stringstream ss;
	ostr<<setw(25)<<left<<"# w_n";
	for(int ix=0; ix<Nx; ix++){
		ss.str(""); ss<<"("<<Lattice.getxgrid(ix)(0)<<","<<Lattice.getxgrid(ix)(1)<<")";
		ostr<<setw(50)<<right<<ss.str();
	}
	ostr<<endl;

	ostr<<setprecision(15)<<scientific;
	for(int n=0; n<Nw; n++){
		ostr<<setw(25)<<left<<w(n);
		for(int ix=0; ix<Nx; ix++){
			ostr<<setw(25)<<right<<Sw[iorder-1][iflavor][jflavor].getOwx(n,ix).real()
				<<setw(25)<<right<<Sw[iorder-1][iflavor][jflavor].getOwx(n,ix).imag();
		}
		ostr<<endl;
	}
};
void controller::printSwk(ostream& ostr, const int& iorder, const int& iflavor, const int& jflavor){
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
			ostr<<setw(25)<<right<<Sw[iorder-1][iflavor][jflavor].getOwk(n,ik).real()
				<<setw(25)<<right<<Sw[iorder-1][iflavor][jflavor].getOwk(n,ik).imag();
		}
		ostr<<endl;
	}
};
void controller::printdSwx(ostream& ostr, const int& iorder, const int& iflavor, const int& jflavor){
	stringstream ss;
	ostr<<setw(25)<<left<<"# w_n";
	for(int ix=0; ix<Nx; ix++){
		ss.str(""); ss<<"("<<Lattice.getxgrid(ix)(0)<<","<<Lattice.getxgrid(ix)(1)<<")";
		ostr<<setw(50)<<right<<ss.str();
	}
	ostr<<endl;

	ostr<<setprecision(15)<<scientific;
	for(int n=0; n<Nw; n++){
		ostr<<setw(25)<<left<<w(n);
		for(int ix=0; ix<Nx; ix++){
			ostr<<setw(25)<<right<<dSw[iorder-1][iflavor][jflavor].getOwx(n,ix).real()
				<<setw(25)<<right<<dSw[iorder-1][iflavor][jflavor].getOwx(n,ix).imag();
		}
		ostr<<endl;
	}
};
void controller::printdSwk(ostream& ostr, const int& iorder, const int& iflavor, const int& jflavor){
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
			ostr<<setw(25)<<right<<dSw[iorder-1][iflavor][jflavor].getOwk(n,ik).real()
				<<setw(25)<<right<<dSw[iorder-1][iflavor][jflavor].getOwk(n,ik).imag();
		}
		ostr<<endl;
	}
};
void controller::printPwx(ostream& ostr, const int& iorder, const int& iflavor, const int& jflavor){
	stringstream ss;
	ostr<<setw(25)<<left<<"# q_n";
	for(int ik=0; ik<Nx; ik++){
		ss.str(""); ss<<"("<<Lattice.getkgrid(ik)(0)<<","<<Lattice.getkgrid(ik)(1)<<")";
		ostr<<setw(50)<<right<<ss.str();
	}
	ostr<<endl;

	ostr.precision(15);
	ostr<<scientific;
	for(int n=0; n<Nw; n++){
		ostr<<setw(25)<<left<<wB(n);
		for(int ix=0; ix<Nx; ix++){
			ostr<<setw(25)<<right<<Pw[iorder][iflavor][jflavor].getOwx(n,ix).real()
				<<setw(25)<<right<<Pw[iorder][iflavor][jflavor].getOwx(n,ix).imag();
		}
		ostr<<endl;
	}
};
void controller::printPwk(ostream& ostr, const int& iorder, const int& iflavor, const int& jflavor){
	stringstream ss;
	ostr<<setw(25)<<left<<"# q_n";
	for(int ik=0; ik<Nx; ik++){
		ss.str(""); ss<<"("<<Lattice.getkgrid(ik)(0)<<","<<Lattice.getkgrid(ik)(1)<<")";
		ostr<<setw(50)<<right<<ss.str();
	}
	ostr<<endl;

	ostr.precision(15);
	ostr<<scientific;
	for(int n=0; n<Nw; n++){
		ostr<<setw(25)<<left<<wB(n);
		for(int ik=0; ik<Nx; ik++){
			ostr<<setw(25)<<right<<Pw[iorder][iflavor][jflavor].getOwk(n,ik).real()
				<<setw(25)<<right<<Pw[iorder][iflavor][jflavor].getOwk(n,ik).imag();
		}
		ostr<<endl;
	}
};
void controller::printdPwx(ostream& ostr, const int& iorder, const int& iflavor, const int& jflavor){
	stringstream ss;
	ostr<<setw(25)<<left<<"# q_n";
	for(int ik=0; ik<Nx; ik++){
		ss.str(""); ss<<"("<<Lattice.getkgrid(ik)(0)<<","<<Lattice.getkgrid(ik)(1)<<")";
		ostr<<setw(50)<<right<<ss.str();
	}
	ostr<<endl;

	ostr.precision(15);
	ostr<<scientific;
	for(int n=0; n<Nw; n++){
		ostr<<setw(25)<<left<<wB(n);
		for(int ix=0; ix<Nx; ix++){
			ostr<<setw(25)<<right<<dPw[iorder][iflavor][jflavor].getOwx(n,ix).real()
				<<setw(25)<<right<<dPw[iorder][iflavor][jflavor].getOwx(n,ix).imag();
		}
		ostr<<endl;
	}
};
void controller::printdPwk(ostream& ostr, const int& iorder, const int& iflavor, const int& jflavor){
	stringstream ss;
	ostr<<setw(25)<<left<<"# q_n";
	for(int ik=0; ik<Nx; ik++){
		ss.str(""); ss<<"("<<Lattice.getkgrid(ik)(0)<<","<<Lattice.getkgrid(ik)(1)<<")";
		ostr<<setw(50)<<right<<ss.str();
	}
	ostr<<endl;

	ostr.precision(15);
	ostr<<scientific;
	for(int n=0; n<Nw; n++){
		ostr<<setw(25)<<left<<wB(n);
		for(int ik=0; ik<Nx; ik++){
			ostr<<setw(25)<<right<<dPw[iorder][iflavor][jflavor].getOwk(n,ik).real()
				<<setw(25)<<right<<dPw[iorder][iflavor][jflavor].getOwk(n,ik).imag();
		}
		ostr<<endl;
	}
};

void controller::printswx(ostream& ostr, const int& iorder, const int& iflavor, const int& jflavor){
	stringstream ss;
	ostr<<setw(25)<<left<<"# w_n";
	for(int ix=0; ix<Nx; ix++){
		ss.str(""); ss<<"("<<Lattice.getxgrid(ix)(0)<<","<<Lattice.getxgrid(ix)(1)<<")";
		ostr<<setw(50)<<right<<ss.str();
	}
	ostr<<endl;

	ostr<<setprecision(15)<<scientific;
	for(int n=0; n<Nw; n++){
		ostr<<setw(25)<<left<<w(n);
		for(int ix=0; ix<Nx; ix++){
			ostr<<setw(25)<<right<<sw[iorder-1][iflavor][jflavor].getOwx(n,ix).real()
				<<setw(25)<<right<<sw[iorder-1][iflavor][jflavor].getOwx(n,ix).imag();
		}
		ostr<<endl;
	}
};
void controller::printswk(ostream& ostr, const int& iorder, const int& iflavor, const int& jflavor){
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
			ostr<<setw(25)<<right<<sw[iorder-1][iflavor][jflavor].getOwk(n,ik).real()
				<<setw(25)<<right<<sw[iorder-1][iflavor][jflavor].getOwk(n,ik).imag();
		}
		ostr<<endl;
	}
};
void controller::printdswx(ostream& ostr, const int& iorder, const int& iflavor, const int& jflavor){
	stringstream ss;
	ostr<<setw(25)<<left<<"# w_n";
	for(int ix=0; ix<Nx; ix++){
		ss.str(""); ss<<"("<<Lattice.getxgrid(ix)(0)<<","<<Lattice.getxgrid(ix)(1)<<")";
		ostr<<setw(50)<<right<<ss.str();
	}
	ostr<<endl;

	ostr<<setprecision(15)<<scientific;
	for(int n=0; n<Nw; n++){
		ostr<<setw(25)<<left<<w(n);
		for(int ix=0; ix<Nx; ix++){
			ostr<<setw(25)<<right<<dsw[iorder-1][iflavor][jflavor].getOwx(n,ix).real()
				<<setw(25)<<right<<dsw[iorder-1][iflavor][jflavor].getOwx(n,ix).imag();
		}
		ostr<<endl;
	}
};
void controller::printdswk(ostream& ostr, const int& iorder, const int& iflavor, const int& jflavor){
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
			ostr<<setw(25)<<right<<dsw[iorder-1][iflavor][jflavor].getOwk(n,ik).real()
				<<setw(25)<<right<<dsw[iorder-1][iflavor][jflavor].getOwk(n,ik).imag();
		}
		ostr<<endl;
	}
};
void controller::printpwx(ostream& ostr, const int& iorder, const int& iflavor, const int& jflavor){
	stringstream ss;
	ostr<<setw(25)<<left<<"# q_n";
	for(int ik=0; ik<Nx; ik++){
		ss.str(""); ss<<"("<<Lattice.getkgrid(ik)(0)<<","<<Lattice.getkgrid(ik)(1)<<")";
		ostr<<setw(50)<<right<<ss.str();
	}
	ostr<<endl;

	ostr.precision(15);
	ostr<<scientific;
	for(int n=0; n<Nw; n++){
		ostr<<setw(25)<<left<<wB(n);
		for(int ix=0; ix<Nx; ix++){
			ostr<<setw(25)<<right<<pw[iorder][iflavor][jflavor].getOwx(n,ix).real()
				<<setw(25)<<right<<pw[iorder][iflavor][jflavor].getOwx(n,ix).imag();
		}
		ostr<<endl;
	}
};
void controller::printpwk(ostream& ostr, const int& iorder, const int& iflavor, const int& jflavor){
	stringstream ss;
	ostr<<setw(25)<<left<<"# q_n";
	for(int ik=0; ik<Nx; ik++){
		ss.str(""); ss<<"("<<Lattice.getkgrid(ik)(0)<<","<<Lattice.getkgrid(ik)(1)<<")";
		ostr<<setw(50)<<right<<ss.str();
	}
	ostr<<endl;

	ostr.precision(15);
	ostr<<scientific;
	for(int n=0; n<Nw; n++){
		ostr<<setw(25)<<left<<wB(n);
		for(int ik=0; ik<Nx; ik++){
			ostr<<setw(25)<<right<<pw[iorder][iflavor][jflavor].getOwk(n,ik).real()
				<<setw(25)<<right<<pw[iorder][iflavor][jflavor].getOwk(n,ik).imag();
		}
		ostr<<endl;
	}
};
void controller::printdpwx(ostream& ostr, const int& iorder, const int& iflavor, const int& jflavor){
	stringstream ss;
	ostr<<setw(25)<<left<<"# q_n";
	for(int ik=0; ik<Nx; ik++){
		ss.str(""); ss<<"("<<Lattice.getkgrid(ik)(0)<<","<<Lattice.getkgrid(ik)(1)<<")";
		ostr<<setw(50)<<right<<ss.str();
	}
	ostr<<endl;

	ostr.precision(15);
	ostr<<scientific;
	for(int n=0; n<Nw; n++){
		ostr<<setw(25)<<left<<wB(n);
		for(int ix=0; ix<Nx; ix++){
			ostr<<setw(25)<<right<<dpw[iorder][iflavor][jflavor].getOwx(n,ix).real()
				<<setw(25)<<right<<dpw[iorder][iflavor][jflavor].getOwx(n,ix).imag();
		}
		ostr<<endl;
	}
};
void controller::printdpwk(ostream& ostr, const int& iorder, const int& iflavor, const int& jflavor){
	stringstream ss;
	ostr<<setw(25)<<left<<"# q_n";
	for(int ik=0; ik<Nx; ik++){
		ss.str(""); ss<<"("<<Lattice.getkgrid(ik)(0)<<","<<Lattice.getkgrid(ik)(1)<<")";
		ostr<<setw(50)<<right<<ss.str();
	}
	ostr<<endl;

	ostr.precision(15);
	ostr<<scientific;
	for(int n=0; n<Nw; n++){
		ostr<<setw(25)<<left<<wB(n);
		for(int ik=0; ik<Nx; ik++){
			ostr<<setw(25)<<right<<dpw[iorder][iflavor][jflavor].getOwk(n,ik).real()
				<<setw(25)<<right<<dpw[iorder][iflavor][jflavor].getOwk(n,ik).imag();
		}
		ostr<<endl;
	}
};


void controller::printSwkSum(ostream& ostr, const int& iflavor, const int& jflavor){
	stringstream ss;
	ostr<<setw(25)<<left<<"# w_n";
	for(int ik=0; ik<Nx; ik++){
		ss.str(""); ss<<"("<<Lattice.getkgrid(ik)(0)<<","<<Lattice.getkgrid(ik)(1)<<")";
		ostr<<setw(50)<<right<<ss.str();
	}
	ostr<<endl;

	ostr.precision(15);
	ostr<<scientific;
	for(int n=0; n<Nw; n++){
		ostr<<setw(25)<<left<<w(n);
		for(int ik=0; ik<Nx; ik++){
			ostr<<setw(25)<<right<<SwSum[iflavor][jflavor].getOwk(n,ik).real()
				<<setw(25)<<right<<SwSum[iflavor][jflavor].getOwk(n,ik).imag();
		}
		ostr<<endl;
	}
};
void controller::printPwkSum(ostream& ostr, const int& iflavor, const int& jflavor){
	stringstream ss;
	ostr<<setw(25)<<left<<"# q_n";
	for(int ik=0; ik<Nx; ik++){
		ss.str(""); ss<<"("<<Lattice.getkgrid(ik)(0)<<","<<Lattice.getkgrid(ik)(1)<<")";
		ostr<<setw(50)<<right<<ss.str();
	}
	ostr<<endl;

	ostr.precision(15);
	ostr<<scientific;
	for(int n=0; n<Nw; n++){
		ostr<<setw(25)<<left<<wB(n);
		for(int ik=0; ik<Nx; ik++){
			ostr<<setw(25)<<right<<PwSum[iflavor][jflavor].getOwk(n,ik).real()
				<<setw(25)<<right<<PwSum[iflavor][jflavor].getOwk(n,ik).imag();
		}
		ostr<<endl;
	}
};
void controller::printWwkSum_prev(ostream& ostr, const int& iflavor, const int& jflavor){
	stringstream ss;
	ostr<<setw(25)<<left<<"# w_n";
	for(int ik=0; ik<Nx; ik++){
		ss.str(""); ss<<"("<<Lattice.getkgrid(ik)(0)<<","<<Lattice.getkgrid(ik)(1)<<")";
		ostr<<setw(50)<<right<<ss.str();
	}
	ostr<<endl;

	ostr.precision(15);
	ostr<<scientific;
	for(int n=0; n<Nw; n++){
		ostr<<setw(25)<<left<<wB(n);
		for(int ik=0; ik<Nx; ik++){
			ostr<<setw(25)<<right<<WwSum_prev[iflavor][jflavor].getOwk(n,ik).real()
				<<setw(25)<<right<<WwSum_prev[iflavor][jflavor].getOwk(n,ik).imag();
		}
		ostr<<endl;
	}
};
void controller::printWwkSum(ostream& ostr, const int& iflavor, const int& jflavor){
	stringstream ss;
	ostr<<setw(25)<<left<<"# w_n";
	for(int ik=0; ik<Nx; ik++){
		ss.str(""); ss<<"("<<Lattice.getkgrid(ik)(0)<<","<<Lattice.getkgrid(ik)(1)<<")";
		ostr<<setw(50)<<right<<ss.str();
	}
	ostr<<endl;

	ostr.precision(15);
	ostr<<scientific;
	for(int n=0; n<Nw; n++){
		ostr<<setw(25)<<left<<wB(n);
		for(int ik=0; ik<Nx; ik++){
			ostr<<setw(25)<<right<<WwSum[iflavor][jflavor].getOwk(n,ik).real()
				<<setw(25)<<right<<WwSum[iflavor][jflavor].getOwk(n,ik).imag();
		}
		ostr<<endl;
	}
};
void controller::printWtx(ostream& ostr, const int& dress, const int& iflavor, const int& jflavor){
	stringstream ss;
	ostr<<setw(25)<<left<<"# t";
	for(int x=0; x<Nx; x++){
		ss.str(""); ss<<"("<<Lattice.getkgrid(x)(0)<<","<<Lattice.getkgrid(x)(1)<<")";
		ostr<<setw(50)<<right<<ss.str();
	}
	ostr<<endl;

	Eigen::MatrixXd Wt_tmp = Wt[dress][iflavor][jflavor].getOtx();

	ostr.precision(15);
	ostr<<scientific;
	for(int ti=0; ti<Ntime; ti++){
		ostr<<setw(25)<<left<<ti*dt;
		for(int x=0; x<Nx; x++){
			ostr<<setw(25)<<right<<Wt_tmp(ti,x);
		}
		ostr<<endl;
	}
};
void controller::printinput(std::string foldername){
	ofstream ofstr;
	stringstream ss;

	ss.str(""); ss<<foldername<<"/Gtx_up.in";
	ofstr.open(ss.str().c_str());
	Gt[0][1][1].printX(ofstr); 
	ofstr.close();

	ss.str(""); ss<<foldername<<"/Gtx_dn.in";
	ofstr.open(ss.str().c_str());
	Gt[0][0][0].printX(ofstr); 
	ofstr.close();

	ss.str(""); ss<<foldername<<"/Pwk_Sum_up_up.in";
	ofstr.open(ss.str().c_str());
	printPwkSum(ofstr,1,1);
	ofstr.close();

	ss.str(""); ss<<foldername<<"/Pwk_Sum_up_dn.in";
	ofstr.open(ss.str().c_str());
	printPwkSum(ofstr,1,0);
	ofstr.close();

	ss.str(""); ss<<foldername<<"/Wwk_up_dn.in";
	ofstr.open(ss.str().c_str());
	printWwkSum_prev(ofstr,1,0);
	ofstr.close();

	ss.str(""); ss<<foldername<<"/Wwk_up_up.in";
	ofstr.open(ss.str().c_str());
	printWwkSum_prev(ofstr,1,1);
	ofstr.close();

	ss.str(""); ss<<foldername<<"/Wtx_up_up.in";
	ofstr.open(ss.str().c_str());
	Wt[0][1][1].printX(ofstr); 
	ofstr.close();

	ss.str(""); ss<<foldername<<"/Wtx_up_dn.in";
	ofstr.open(ss.str().c_str());
	Wt[0][1][0].printX(ofstr); 
	ofstr.close();
};
void controller::printoutput(std::string foldername){
	ofstream ofstr;
	stringstream ss;

	ss.str(""); ss<<foldername<<"/statistics.log";
	ofstr.open(ss.str().c_str());
	printstatistics(ofstr);
	ofstr.close();

	ss.str(""); ss<<foldername<<"/Snormalization.dat";
	ofstr.open(ss.str().c_str());
	printnormalization(0,ofstr);
	ofstr.close();

	ss.str(""); ss<<foldername<<"/Pnormalization.dat";
	ofstr.open(ss.str().c_str());
	printnormalization(1,ofstr);
	ofstr.close();

	for(int i=1; i<Nordermax+1; i++){
		//self-energy
		ss.str(""); ss<<foldername<<"/Swx_o"<<i<<"_up.dat";
		ofstr.open(ss.str().c_str());
		printSwx(ofstr,i,1,1);
		ofstr.close();

		ss.str(""); ss<<foldername<<"/Swx_o"<<i<<"_dn.dat";
		ofstr.open(ss.str().c_str());
		printSwx(ofstr,i,0,0);
		ofstr.close();

		ss.str(""); ss<<foldername<<"/Swk_o"<<i<<"_up.dat";
		ofstr.open(ss.str().c_str());
		printSwk(ofstr,i,1,1);
		ofstr.close();

		ss.str(""); ss<<foldername<<"/Swk_o"<<i<<"_dn.dat";
		ofstr.open(ss.str().c_str());
		printSwk(ofstr,i,0,0);
		ofstr.close();

		//self-energy errorbar
		ss.str(""); ss<<foldername<<"/dSwx_o"<<i<<"_up.dat";
		ofstr.open(ss.str().c_str());
		printdSwx(ofstr,i,1,1);
		ofstr.close();

		ss.str(""); ss<<foldername<<"/dSwx_o"<<i<<"_dn.dat";
		ofstr.open(ss.str().c_str());
		printdSwx(ofstr,i,0,0);
		ofstr.close();

		ss.str(""); ss<<foldername<<"/dSwk_o"<<i<<"_up.dat";
		ofstr.open(ss.str().c_str());
		printdSwk(ofstr,i,1,1);
		ofstr.close();

		ss.str(""); ss<<foldername<<"/dSwk_o"<<i<<"_dn.dat";
		ofstr.open(ss.str().c_str());
		printdSwk(ofstr,i,0,0);
		ofstr.close();


		//self-energy
		ss.str(""); ss<<foldername<<"/swx_n"<<i<<"_up.dat";
		ofstr.open(ss.str().c_str());
		printswx(ofstr,i,1,1);
		ofstr.close();

		ss.str(""); ss<<foldername<<"/swx_n"<<i<<"_dn.dat";
		ofstr.open(ss.str().c_str());
		printswx(ofstr,i,0,0);
		ofstr.close();

		ss.str(""); ss<<foldername<<"/swk_n"<<i<<"_up.dat";
		ofstr.open(ss.str().c_str());
		printswk(ofstr,i,1,1);
		ofstr.close();

		ss.str(""); ss<<foldername<<"/swk_n"<<i<<"_dn.dat";
		ofstr.open(ss.str().c_str());
		printswk(ofstr,i,0,0);
		ofstr.close();

		//self-energy errorbar
		ss.str(""); ss<<foldername<<"/dswx_n"<<i<<"_up.dat";
		ofstr.open(ss.str().c_str());
		printdswx(ofstr,i,1,1);
		ofstr.close();

		ss.str(""); ss<<foldername<<"/dswx_n"<<i<<"_dn.dat";
		ofstr.open(ss.str().c_str());
		printdswx(ofstr,i,0,0);
		ofstr.close();

		ss.str(""); ss<<foldername<<"/dswk_n"<<i<<"_up.dat";
		ofstr.open(ss.str().c_str());
		printdswk(ofstr,i,1,1);
		ofstr.close();

		ss.str(""); ss<<foldername<<"/dswk_n"<<i<<"_dn.dat";
		ofstr.open(ss.str().c_str());
		printdswk(ofstr,i,0,0);
		ofstr.close();

	}

	ss.str(""); ss<<foldername<<"/Ptxloc.dat";
	ofstr.open(ss.str().c_str());
	printPloc(ofstr);
	ofstr.close();

	for(int i=0; i<Nordermax; i++){
		// polarization
		ss.str(""); ss<<foldername<<"/Pwx_o"<<i<<"_up_up.dat";
		ofstr.open(ss.str().c_str());
		printPwx(ofstr,i,1,1);
		ofstr.close();

		ss.str(""); ss<<foldername<<"/Pwx_o"<<i<<"_up_dn.dat";
		ofstr.open(ss.str().c_str());
		printPwx(ofstr,i,1,0);
		ofstr.close();

		ss.str(""); ss<<foldername<<"/Pwx_o"<<i<<"_dn_dn.dat";
		ofstr.open(ss.str().c_str());
		printPwx(ofstr,i,0,0);
		ofstr.close();

		ss.str(""); ss<<foldername<<"/Pwx_o"<<i<<"_dn_up.dat";
		ofstr.open(ss.str().c_str());
		printPwx(ofstr,i,0,1);
		ofstr.close();

		ss.str(""); ss<<foldername<<"/Pwk_o"<<i<<"_up_up.dat";
		ofstr.open(ss.str().c_str());
		printPwk(ofstr,i,1,1);
		ofstr.close();

		ss.str(""); ss<<foldername<<"/Pwk_o"<<i<<"_up_dn.dat";
		ofstr.open(ss.str().c_str());
		printPwk(ofstr,i,1,0);
		ofstr.close();

		ss.str(""); ss<<foldername<<"/Pwk_o"<<i<<"_dn_dn.dat";
		ofstr.open(ss.str().c_str());
		printPwk(ofstr,i,0,0);
		ofstr.close();

		ss.str(""); ss<<foldername<<"/Pwk_o"<<i<<"_dn_up.dat";
		ofstr.open(ss.str().c_str());
		printPwk(ofstr,i,0,1);
		ofstr.close();

		// polarization errorbar
		ss.str(""); ss<<foldername<<"/dPwx_o"<<i<<"_up_up.dat";
		ofstr.open(ss.str().c_str());
		printdPwx(ofstr,i,1,1);
		ofstr.close();

		ss.str(""); ss<<foldername<<"/dPwx_o"<<i<<"_up_dn.dat";
		ofstr.open(ss.str().c_str());
		printdPwx(ofstr,i,1,0);
		ofstr.close();

		ss.str(""); ss<<foldername<<"/dPwx_o"<<i<<"_dn_dn.dat";
		ofstr.open(ss.str().c_str());
		printdPwx(ofstr,i,0,0);
		ofstr.close();

		ss.str(""); ss<<foldername<<"/dPwx_o"<<i<<"_dn_up.dat";
		ofstr.open(ss.str().c_str());
		printdPwx(ofstr,i,0,1);
		ofstr.close();

		ss.str(""); ss<<foldername<<"/dPwk_o"<<i<<"_up_up.dat";
		ofstr.open(ss.str().c_str());
		printdPwk(ofstr,i,1,1);
		ofstr.close();

		ss.str(""); ss<<foldername<<"/dPwk_o"<<i<<"_up_dn.dat";
		ofstr.open(ss.str().c_str());
		printdPwk(ofstr,i,1,0);
		ofstr.close();

		ss.str(""); ss<<foldername<<"/dPwk_o"<<i<<"_dn_dn.dat";
		ofstr.open(ss.str().c_str());
		printdPwk(ofstr,i,0,0);
		ofstr.close();

		ss.str(""); ss<<foldername<<"/dPwk_o"<<i<<"_dn_up.dat";
		ofstr.open(ss.str().c_str());
		printdPwk(ofstr,i,0,1);
		ofstr.close();



		// polarization
		ss.str(""); ss<<foldername<<"/pwx_n"<<i<<"_up_up.dat";
		ofstr.open(ss.str().c_str());
		printpwx(ofstr,i,1,1);
		ofstr.close();

		ss.str(""); ss<<foldername<<"/pwx_n"<<i<<"_up_dn.dat";
		ofstr.open(ss.str().c_str());
		printpwx(ofstr,i,1,0);
		ofstr.close();

		ss.str(""); ss<<foldername<<"/pwx_n"<<i<<"_dn_dn.dat";
		ofstr.open(ss.str().c_str());
		printpwx(ofstr,i,0,0);
		ofstr.close();

		ss.str(""); ss<<foldername<<"/pwx_n"<<i<<"_dn_up.dat";
		ofstr.open(ss.str().c_str());
		printpwx(ofstr,i,0,1);
		ofstr.close();

		ss.str(""); ss<<foldername<<"/pwk_n"<<i<<"_up_up.dat";
		ofstr.open(ss.str().c_str());
		printpwk(ofstr,i,1,1);
		ofstr.close();

		ss.str(""); ss<<foldername<<"/pwk_n"<<i<<"_up_dn.dat";
		ofstr.open(ss.str().c_str());
		printpwk(ofstr,i,1,0);
		ofstr.close();

		ss.str(""); ss<<foldername<<"/pwk_n"<<i<<"_dn_dn.dat";
		ofstr.open(ss.str().c_str());
		printpwk(ofstr,i,0,0);
		ofstr.close();

		ss.str(""); ss<<foldername<<"/pwk_n"<<i<<"_dn_up.dat";
		ofstr.open(ss.str().c_str());
		printpwk(ofstr,i,0,1);
		ofstr.close();

		// polarization errorbar
		ss.str(""); ss<<foldername<<"/dpwx_n"<<i<<"_up_up.dat";
		ofstr.open(ss.str().c_str());
		printdpwx(ofstr,i,1,1);
		ofstr.close();

		ss.str(""); ss<<foldername<<"/dpwx_n"<<i<<"_up_dn.dat";
		ofstr.open(ss.str().c_str());
		printdpwx(ofstr,i,1,0);
		ofstr.close();

		ss.str(""); ss<<foldername<<"/dpwx_n"<<i<<"_dn_dn.dat";
		ofstr.open(ss.str().c_str());
		printdpwx(ofstr,i,0,0);
		ofstr.close();

		ss.str(""); ss<<foldername<<"/dpwx_n"<<i<<"_dn_up.dat";
		ofstr.open(ss.str().c_str());
		printdpwx(ofstr,i,0,1);
		ofstr.close();

		ss.str(""); ss<<foldername<<"/dpwk_n"<<i<<"_up_up.dat";
		ofstr.open(ss.str().c_str());
		printdpwk(ofstr,i,1,1);
		ofstr.close();

		ss.str(""); ss<<foldername<<"/dpwk_n"<<i<<"_up_dn.dat";
		ofstr.open(ss.str().c_str());
		printdpwk(ofstr,i,1,0);
		ofstr.close();

		ss.str(""); ss<<foldername<<"/dpwk_n"<<i<<"_dn_dn.dat";
		ofstr.open(ss.str().c_str());
		printdpwk(ofstr,i,0,0);
		ofstr.close();

		ss.str(""); ss<<foldername<<"/dpwk_n"<<i<<"_dn_up.dat";
		ofstr.open(ss.str().c_str());
		printdpwk(ofstr,i,0,1);
		ofstr.close();
	}
};
void controller::printoutput_S(std::string foldername){
	ofstream ofstr;
	stringstream ss;

	ss.str(""); ss<<foldername<<"/statistics.log";
	ofstr.open(ss.str().c_str());
	printstatistics(ofstr);
	ofstr.close();

	ss.str(""); ss<<foldername<<"/Snormalization.dat";
	ofstr.open(ss.str().c_str());
	printnormalization(0,ofstr);
	ofstr.close();

	for(int i=1; i<Nordermax+1; i++){
		//self-energy
		ss.str(""); ss<<foldername<<"/Swx_o"<<i<<"_up.dat";
		ofstr.open(ss.str().c_str());
		printSwx(ofstr,i,1,1);
		ofstr.close();

		ss.str(""); ss<<foldername<<"/Swx_o"<<i<<"_dn.dat";
		ofstr.open(ss.str().c_str());
		printSwx(ofstr,i,0,0);
		ofstr.close();

		ss.str(""); ss<<foldername<<"/Swk_o"<<i<<"_up.dat";
		ofstr.open(ss.str().c_str());
		printSwk(ofstr,i,1,1);
		ofstr.close();

		ss.str(""); ss<<foldername<<"/Swk_o"<<i<<"_dn.dat";
		ofstr.open(ss.str().c_str());
		printSwk(ofstr,i,0,0);
		ofstr.close();

		//self-energy errorbar
		ss.str(""); ss<<foldername<<"/dSwx_o"<<i<<"_up.dat";
		ofstr.open(ss.str().c_str());
		printdSwx(ofstr,i,1,1);
		ofstr.close();

		ss.str(""); ss<<foldername<<"/dSwx_o"<<i<<"_dn.dat";
		ofstr.open(ss.str().c_str());
		printdSwx(ofstr,i,0,0);
		ofstr.close();

		ss.str(""); ss<<foldername<<"/dSwk_o"<<i<<"_up.dat";
		ofstr.open(ss.str().c_str());
		printdSwk(ofstr,i,1,1);
		ofstr.close();

		ss.str(""); ss<<foldername<<"/dSwk_o"<<i<<"_dn.dat";
		ofstr.open(ss.str().c_str());
		printdSwk(ofstr,i,0,0);
		ofstr.close();


		//self-energy
		ss.str(""); ss<<foldername<<"/swx_n"<<i<<"_up.dat";
		ofstr.open(ss.str().c_str());
		printswx(ofstr,i,1,1);
		ofstr.close();

		ss.str(""); ss<<foldername<<"/swx_n"<<i<<"_dn.dat";
		ofstr.open(ss.str().c_str());
		printswx(ofstr,i,0,0);
		ofstr.close();

		ss.str(""); ss<<foldername<<"/swk_n"<<i<<"_up.dat";
		ofstr.open(ss.str().c_str());
		printswk(ofstr,i,1,1);
		ofstr.close();

		ss.str(""); ss<<foldername<<"/swk_n"<<i<<"_dn.dat";
		ofstr.open(ss.str().c_str());
		printswk(ofstr,i,0,0);
		ofstr.close();

		//self-energy errorbar
		ss.str(""); ss<<foldername<<"/dswx_n"<<i<<"_up.dat";
		ofstr.open(ss.str().c_str());
		printdswx(ofstr,i,1,1);
		ofstr.close();

		ss.str(""); ss<<foldername<<"/dswx_n"<<i<<"_dn.dat";
		ofstr.open(ss.str().c_str());
		printdswx(ofstr,i,0,0);
		ofstr.close();

		ss.str(""); ss<<foldername<<"/dswk_n"<<i<<"_up.dat";
		ofstr.open(ss.str().c_str());
		printdswk(ofstr,i,1,1);
		ofstr.close();

		ss.str(""); ss<<foldername<<"/dswk_n"<<i<<"_dn.dat";
		ofstr.open(ss.str().c_str());
		printdswk(ofstr,i,0,0);
		ofstr.close();

	}
};

void controller::printiterativeoutput(std::string foldername, const int& index){
	ofstream ofstr;
	stringstream ss;

	ss.str(""); ss<<foldername<<"/Pwk_Sum_up_up."<<index;
	ofstr.open(ss.str().c_str());
	printPwkSum(ofstr,1,1);
	ofstr.close();

	ss.str(""); ss<<foldername<<"/Pwk_Sum_up_dn."<<index;
	ofstr.open(ss.str().c_str());
	printPwkSum(ofstr,1,0);
	ofstr.close();

	ss.str(""); ss<<foldername<<"/Wwk_up_dn."<<index;
	ofstr.open(ss.str().c_str());
	printWwkSum_prev(ofstr,1,0);
	ofstr.close();

	ss.str(""); ss<<foldername<<"/Wwk_up_up."<<index;
	ofstr.open(ss.str().c_str());
	printWwkSum_prev(ofstr,1,1);
	ofstr.close();

	/*ss.str(""); ss<<foldername<<"/statistics.log";
	ofstr.open(ss.str().c_str(), ofstream::app);
	ofstr<<"# Nordermax : "<<iordermax<<endl;
	printstatistics(ofstr);
	ofstr.close();

	ss.str(""); ss<<foldername<<"/normalization.dat";
	ofstr.open(ss.str().c_str(), ofstream::app);
	ofstr<<setw(25)<<left<<iordermax
		<<setw(25)<<right<<scientific<<normalization[iordermax-1]
		<<setw(25)<<right<<scientific<<dnormalization[iordermax-1]<<endl;
	ofstr.close();

	ss.str(""); ss<<foldername<<"/Swx_"<<iordermax<<"th_up.dat";
	ofstr.open(ss.str().c_str());
	printSwx(ofstr,iordermax,1,1);
	ofstr.close();

	ss.str(""); ss<<foldername<<"/Swx_"<<iordermax<<"th_dn.dat";
	ofstr.open(ss.str().c_str());
	printSwx(ofstr,iordermax,0,0);
	ofstr.close();

	ss.str(""); ss<<foldername<<"/Swk_"<<iordermax<<"th_up.dat";
	ofstr.open(ss.str().c_str());
	printSwk(ofstr,iordermax,1,1);
	ofstr.close();

	ss.str(""); ss<<foldername<<"/Swk_"<<iordermax<<"th_dn.dat";
	ofstr.open(ss.str().c_str());
	printSwk(ofstr,iordermax,0,0);
	ofstr.close();


	ss.str(""); ss<<foldername<<"/dSwx_"<<iordermax<<"th_up.dat";
	ofstr.open(ss.str().c_str());
	printdSwx(ofstr,iordermax,1,1);
	ofstr.close();

	ss.str(""); ss<<foldername<<"/dSwx_"<<iordermax<<"th_dn.dat";
	ofstr.open(ss.str().c_str());
	printdSwx(ofstr,iordermax,0,0);
	ofstr.close();

	ss.str(""); ss<<foldername<<"/dSwk_"<<iordermax<<"th_up.dat";
	ofstr.open(ss.str().c_str());
	printdSwk(ofstr,iordermax,1,1);
	ofstr.close();

	ss.str(""); ss<<foldername<<"/dSwk_"<<iordermax<<"th_dn.dat";
	ofstr.open(ss.str().c_str());
	printdSwk(ofstr,iordermax,0,0);
	ofstr.close();*/

	/*ss.str(""); ss<<foldername<<"/Gwx_"<<iordermax<<"th_up.dat";
	ofstr.open(ss.str().c_str());
	Gw[iordermax-1][1][1].printOwk(ofstr);
	ofstr.close();

	ss.str(""); ss<<foldername<<"/Gwx_"<<iordermax<<"th_dn.dat";
	ofstr.open(ss.str().c_str());
	Gw[iordermax-1][0][0].printOwk(ofstr);
	ofstr.close();*/
};

/*void controller::test(){
	//cout<<"At controller, sampling starts"<<endl<<flush;
	//const int root = 0;
	//const int rank = MPI::COMM_WORLD.Get_rank();

	//cout<<"accumulation quantities initialization"<<endl;
	//cout<<"actual sampling process"<<endl;
	diagram Diag;
	Diag.debug();
};
void controller::sampling(){
	//cout<<"At controller, sampling starts"<<endl<<flush;
	const int root = 0;
	const int rank = MPI::COMM_WORLD.Get_rank();

	//cout<<"accumulation quantities initialization"<<endl;
	selfenergy_w Swtmp;
	Eigen::MatrixXcd M_Swx, M_Swk;
	Eigen::MatrixXcd ***M_Swxacc, ***M_Swxsqacc;
	Eigen::MatrixXcd ***M_Swkacc, ***M_Swksqacc;
	M_Swxacc = new Eigen::MatrixXcd**[Nordermax];
	M_Swkacc = new Eigen::MatrixXcd**[Nordermax];
	M_Swxsqacc = new Eigen::MatrixXcd**[Nordermax];
	M_Swksqacc = new Eigen::MatrixXcd**[Nordermax];
	for(int iorder=0; iorder<Nordermax; iorder++){
		M_Swxacc[iorder] = new Eigen::MatrixXcd*[Nflavor];
		M_Swkacc[iorder] = new Eigen::MatrixXcd*[Nflavor];
		M_Swxsqacc[iorder] = new Eigen::MatrixXcd*[Nflavor];
		M_Swksqacc[iorder] = new Eigen::MatrixXcd*[Nflavor];
		for(int iflavor=0; iflavor<Nflavor; iflavor++){
			M_Swxacc[iorder][iflavor] = new Eigen::MatrixXcd[Nflavor];
			M_Swkacc[iorder][iflavor] = new Eigen::MatrixXcd[Nflavor];
			M_Swxsqacc[iorder][iflavor] = new Eigen::MatrixXcd[Nflavor];
			M_Swksqacc[iorder][iflavor] = new Eigen::MatrixXcd[Nflavor];
			for(int jflavor=0; jflavor<Nflavor; jflavor++){
				M_Swxacc[iorder][iflavor][jflavor] = Eigen::MatrixXcd::Zero(Nw,Nx);
				M_Swkacc[iorder][iflavor][jflavor] = Eigen::MatrixXcd::Zero(Nw,Nx);
				M_Swxsqacc[iorder][iflavor][jflavor] = Eigen::MatrixXcd::Zero(Nw,Nx);
				M_Swksqacc[iorder][iflavor][jflavor] = Eigen::MatrixXcd::Zero(Nw,Nx);
			}
		}
	}

	//cout<<"actual sampling process"<<endl;
	diagram Diag;
	Diag.thermalization();
	if(rank==root) cout<<"controller::sampling : "<<flush;
	MPI::COMM_WORLD.Barrier();
	for(int i=0; i<Nerr; i++){
		cout<<"*"<<flush;
		Diag.diagMC();
		//Diag.statistics();

		for(int n=0; n<Nordermax; n++)
		for(int iflavor=0; iflavor<Nflavor; iflavor++) for(int jflavor=0; jflavor<Nflavor; jflavor++){
			Swtmp = Diag.getOw(n,iflavor,jflavor);
			M_Swx = Swtmp.getOwx();
			M_Swk = Swtmp.getOwk();
			M_Swxacc[n][iflavor][jflavor] += M_Swx;
			M_Swkacc[n][iflavor][jflavor] += M_Swk;
			M_Swxsqacc[n][iflavor][jflavor] += M_Swx.real().cwiseAbs2() + complex<double>(0.0,1.0)*M_Swx.imag().cwiseAbs2();
			M_Swksqacc[n][iflavor][jflavor] += M_Swk.real().cwiseAbs2() + complex<double>(0.0,1.0)*M_Swk.imag().cwiseAbs2();

			//if(n==1 && iflavor==0 && jflavor==0){ cout<<"[("<<M_Swx(0,0)<<"),"; }
		}
		//cout<<(M_Swx.real().cwiseAbs2() + complex<double>(0.0,1.0)*M_Swx.imag().cwiseAbs2())<<endl;
	}
	if(rank==root) cout<<endl<<flush;

	//cout<<"MPI communication starts"<<endl<<flush;
	//cout<<"MPI Pack"<<endl<<flush;
	
	// 4 : M_Swxacc, M_Swkacc, M_Swxsqacc, M_Swksqacc
	int Ndata = 4*Nordermax*Nflavor*Nflavor*Nw*Nx;
	complex<double> *source, *result, *output1, *output2;
	source  = new complex<double>[Ndata];
	result  = new complex<double>[Ndata];
	output1  = new complex<double>[Nw*Nx];
	output2  = new complex<double>[Nw*Nx];

	int position = 0;
	for(int n=0; n<Nordermax; n++)
	for(int iflavor=0; iflavor<Nflavor; iflavor++) for(int jflavor=0; jflavor<Nflavor; jflavor++){
		MPI::DOUBLE_COMPLEX.Pack(M_Swxacc[n][iflavor][jflavor].data(),Nx*Nw,source,16*Ndata,position,MPI::COMM_WORLD);
		MPI::DOUBLE_COMPLEX.Pack(M_Swxsqacc[n][iflavor][jflavor].data(),Nx*Nw,source,16*Ndata,position,MPI::COMM_WORLD);
		MPI::DOUBLE_COMPLEX.Pack(M_Swkacc[n][iflavor][jflavor].data(),Nx*Nw,source,16*Ndata,position,MPI::COMM_WORLD);
		MPI::DOUBLE_COMPLEX.Pack(M_Swksqacc[n][iflavor][jflavor].data(),Nx*Nw,source,16*Ndata,position,MPI::COMM_WORLD);
	}

	//cout<<"MPI Reduce and Bcast"<<endl<<flush;
	MPI::COMM_WORLD.Reduce(source,result,Ndata,MPI::DOUBLE_COMPLEX,MPI::SUM,root);
	MPI::COMM_WORLD.Bcast(result,Ndata,MPI::DOUBLE_COMPLEX,root);

	//cout<<"MPI Unpack"<<endl<<flush;
	position = 0;
	for(int n=0; n<Nordermax; n++)
	for(int iflavor=0; iflavor<Nflavor; iflavor++) for(int jflavor=0; jflavor<Nflavor; jflavor++){
		MPI::DOUBLE_COMPLEX.Unpack(result,16*Ndata,output1,Nx*Nw,position,MPI::COMM_WORLD);
		MPI::DOUBLE_COMPLEX.Unpack(result,16*Ndata,output2,Nx*Nw,position,MPI::COMM_WORLD);
		for(int i=0; i<Nw*Nx; i++){
			output1[i] /= static_cast<double>(NerrTot);
			output2[i] = complex<double>(
					sqrt(abs(
						output2[i].real()/static_cast<double>(NerrTot) 
						- output1[i].real()*output1[i].real()
					)/static_cast<double>(NerrTot)),

					sqrt(abs(
						output2[i].imag()/static_cast<double>(NerrTot) 
						- output1[i].imag()*output1[i].imag()
					)/static_cast<double>(NerrTot))
					);
		}
		Sw[n][iflavor][jflavor].setOwx(output1);
		dSw[n][iflavor][jflavor].setOwx(output2);
		//if(n==1 && iflavor==0 && jflavor==0) cout<<output1[0]<<" "<<output2[0]<<endl;

		MPI::DOUBLE_COMPLEX.Unpack(result,16*Ndata,output1,Nx*Nw,position,MPI::COMM_WORLD);
		MPI::DOUBLE_COMPLEX.Unpack(result,16*Ndata,output2,Nx*Nw,position,MPI::COMM_WORLD);
		for(int i=0; i<Nw*Nx; i++){
			output1[i] /= static_cast<double>(NerrTot);
			output2[i] = complex<double>(
					sqrt(abs(
						output2[i].real()/static_cast<double>(NerrTot) 
						- output1[i].real()*output1[i].real()
					)/static_cast<double>(NerrTot)),

					sqrt(abs(
						output2[i].imag()/static_cast<double>(NerrTot) 
						- output1[i].imag()*output1[i].imag()
					)/static_cast<double>(NerrTot))
					);
		}
		Sw[n][iflavor][jflavor].setOwk(output1);
		dSw[n][iflavor][jflavor].setOwk(output2);
	}

	delete[] source;
	delete[] result;
	delete[] output1;
	delete[] output2;
	delete[] M_Swxacc;
	delete[] M_Swxsqacc;
	delete[] M_Swkacc;
	delete[] M_Swksqacc;
};*/
/*void constructG(const int& iorder){
	selfenergy_w SG;
	double coefficient = 1.;
	for(int k=0; k<iorder; k++){
		for(int iflavor=0; iflavor<Nflavor; iflavor++)
			SG += coefficient*Sw[iorder-k-1][iflavor][iflavor]*G[k][iflavor][iflavor]
		coefficient *= static_cast<double>(iorder-k-1)/(k+1);
	}
	for(int iflavor=0; iflavor<Nflavor; iflavor++){
		G[iorder][iflavor][iflavor] = G[0][iflavor][iflavor]*SG;
	}
};
void controller::printSwx(ostream& ostr, const int& iorder, const int& iflavor, const int& jflavor){
	stringstream ss;
	ostr<<setw(25)<<left<<"# w_n";
	for(int ix=0; ix<Nx; ix++){
		ss.str(""); ss<<"("<<Lattice.getxgrid(ix)(0)<<","<<Lattice.getxgrid(ix)(1)<<")";
		ostr<<setw(50)<<right<<ss.str();
	}
	ostr<<endl;

	ostr<<setprecision(15)<<scientific;
	for(int n=0; n<Nw; n++){
		ostr<<setw(25)<<left<<w(n);
		for(int ix=0; ix<Nx; ix++){
			ostr<<setw(25)<<right<<Sw[iorder-1][iflavor][jflavor].getOwx(n,ix).real()
				<<setw(25)<<right<<Sw[iorder-1][iflavor][jflavor].getOwx(n,ix).imag();
		}
		ostr<<endl;
	}
};
void controller::printSwk(ostream& ostr, const int& iorder, const int& iflavor, const int& jflavor){
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
			ostr<<setw(25)<<right<<Sw[iorder-1][iflavor][jflavor].getOwk(n,ik).real()
				<<setw(25)<<right<<Sw[iorder-1][iflavor][jflavor].getOwk(n,ik).imag();
		}
		ostr<<endl;
	}
};
void controller::printdSwx(ostream& ostr, const int& iorder, const int& iflavor, const int& jflavor){
	stringstream ss;
	ostr<<setw(25)<<left<<"# w_n";
	for(int ix=0; ix<Nx; ix++){
		ss.str(""); ss<<"("<<Lattice.getxgrid(ix)(0)<<","<<Lattice.getxgrid(ix)(1)<<")";
		ostr<<setw(50)<<right<<ss.str();
	}
	ostr<<endl;

	ostr<<setprecision(15)<<scientific;
	for(int n=0; n<Nw; n++){
		ostr<<setw(25)<<left<<w(n);
		for(int ix=0; ix<Nx; ix++){
			ostr<<setw(25)<<right<<dSw[iorder-1][iflavor][jflavor].getOwx(n,ix).real()
				<<setw(25)<<right<<dSw[iorder-1][iflavor][jflavor].getOwx(n,ix).imag();
		}
		ostr<<endl;
	}
};
void controller::printdSwk(ostream& ostr, const int& iorder, const int& iflavor, const int& jflavor){
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
			ostr<<setw(25)<<right<<dSw[iorder-1][iflavor][jflavor].getOwk(n,ik).real()
				<<setw(25)<<right<<dSw[iorder-1][iflavor][jflavor].getOwk(n,ik).imag();
		}
		ostr<<endl;
	}
};

void controller::printoutput(std::string foldername){
	ofstream ofstr;
	stringstream ss;
	for(int i=2; i<Nordermax+1; i++){
		ss.str(""); ss<<foldername<<"/Swx_"<<i<<"th_up.dat";
		ofstr.open(ss.str().c_str());
		printSwx(ofstr,i,1,1);
		ofstr.close();

		ss.str(""); ss<<foldername<<"/Swx_"<<i<<"th_dn.dat";
		ofstr.open(ss.str().c_str());
		printSwx(ofstr,i,0,0);
		ofstr.close();

		ss.str(""); ss<<foldername<<"/Swk_"<<i<<"th_up.dat";
		ofstr.open(ss.str().c_str());
		printSwk(ofstr,i,1,1);
		ofstr.close();

		ss.str(""); ss<<foldername<<"/Swk_"<<i<<"th_dn.dat";
		ofstr.open(ss.str().c_str());
		printSwk(ofstr,i,0,0);
		ofstr.close();


		ss.str(""); ss<<foldername<<"/dSwx_"<<i<<"th_up.dat";
		ofstr.open(ss.str().c_str());
		printdSwx(ofstr,i,1,1);
		ofstr.close();

		ss.str(""); ss<<foldername<<"/dSwx_"<<i<<"th_dn.dat";
		ofstr.open(ss.str().c_str());
		printdSwx(ofstr,i,0,0);
		ofstr.close();

		ss.str(""); ss<<foldername<<"/dSwk_"<<i<<"th_up.dat";
		ofstr.open(ss.str().c_str());
		printdSwk(ofstr,i,1,1);
		ofstr.close();

		ss.str(""); ss<<foldername<<"/dSwk_"<<i<<"th_dn.dat";
		ofstr.open(ss.str().c_str());
		printdSwk(ofstr,i,0,0);
		ofstr.close();

	}

	//ss.str(""); ss<<foldername<<"/Swx_"<<3<<"th_up.dat";
	//ofstr.open(ss.str().c_str());
	//printSwx(ofstr,3,1,1);
	//ofstr.close();

};*/
