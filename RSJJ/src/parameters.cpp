#include "parameters.hpp"
#include <mpi.h>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <algorithm>
#include <chrono>

namespace parameters{
	nlohmann::json	input;

	const int Dspace = 1;
	const int Nptvertex = 3;
	const int Nflavor = 1;

	int Nordermax;
	//int Nordermax, NselfloopMax;
	//int Nauxdress = 1;

	double beta;
	double dt, dtW;
	int Ntau;
	int NWtau;
	//int Nout = 10;
	//int Nw = 100;
	//int Nx, Lx, dim;

	double epsilon;
	double B_field;
	double g_coupling;
	double Omega11, omega_c;
	int Nkx;
	//double w_boson;
	//double W;
	//double n0;

	bool is_mixing;
	bool does_printTt, does_printQt;

	//MtRng64 mt;
	//int site_generate_radius = 8;
	std::mt19937_64 mt;
	std::uniform_real_distribution<double> unidist;
	//std::binomial_distribution<int> bidist;
	//std::vector<double> biprobability;

	std::vector<double> Ri_N;

	std::vector<double> lattice_para;
	//lattice Lattice;
	//hubbardatom Lattice;
	//squarelattice Lattice(2);

	//pseudogreen_t Gt;
	//interaction_t Wt;

	int Nupdate;
	int Nth;
	const int Nmeasureinterval = 10;
	const int NthMax = 1e4;
	long long Nmc, NmcTot;
	int Nerr, NerrTot;
	int Nscl = 10;

	Folder global;

	void MtInit(){
		std::seed_seq seed = {0x12345ULL, 0x23456ULL, 0x34567ULL, 0x45678ULL};
		mt.seed(seed);
	};
	void MtInit(const unsigned long long& seed_i){
		std::seed_seq seed = {seed_i+0x12345ULL, seed_i+0x23456ULL, seed_i+0x34567ULL, seed_i+0x45678ULL};
		mt.seed(seed);
	};
	void MtInit_systime(const unsigned long long& seed_i){
		unsigned seed_sys = std::chrono::system_clock::now().time_since_epoch().count();
		std::seed_seq seed = {seed_i+seed_sys+0x12345ULL, seed_i+seed_sys+0x23456ULL, seed_i+seed_sys+0x34567ULL, seed_i+seed_sys+0x45678ULL};
		mt.seed(seed);
	};
	int getSimpsonFactor(const int& ix, const int& Nx){
		if(ix==0 || ix==Nx) return 1;
		else if(ix%2==1) return 4;
		else return 2;
	};

	/*void select_one_excluding_one(const int& N_i, const int& e_i, int& i_o){
		i_o = (N_i-1)*unidist(mt);
		if(i_o>e_i-1) i_o++;
	};*/
	/*void select_one_excluding_two(const int& N_i, const int& e_i, const int& e_j, int& i_o){
		i_o = N_i*unidist(mt);
		if(i_o>std::min(e_i,e_j)-1) i_o++;
		if(i_o>std::max(e_i,e_j)-1) i_o++;
	};
	int gen_index_merge(const int& m1, const int& m2, const int& o){
		if(o<std::max(m1,m2)) return o;
		else if(o==std::max(m1,m2)) return std::min(m1,m2);
		else return o-1;
	};
	int gen_index_split(const int& m, const int& o){
		if(o>m) return o+1;
		else return o;
	}
	double w(const int& n_i){ return M_PI*(2.*n_i+1.)/beta; };
	double wB(const int& n_i){ return M_PI*2.*n_i/beta; };
	double nF(const double& ep){ return 1./(exp(beta*ep)+1); };
	d+ouble nB(const double& ep){ return 1./(exp(beta*ep)-1); };*/
	void ReadData(char* argv[]){
		//std::cout<<"* ReadData *"<<std::endl<<std::flush;
		std::ifstream para(argv[1]);
	       	para >> input;
	       	para.close(); 

		const int root = 0;
		const int rank = MPI::COMM_WORLD.Get_rank();
		if( input.find("is_mt_random")!=input.end() ){
			if( input["is_mt_random"] ) MtInit_systime(rank);
			else MtInit(rank);
		}
		else MtInit(rank);
		//MtInit(rank);

		//std::cout<<"check point -1"<<std::endl<<std::flush;

		Nordermax = input["Nordermax"];
	       	//NselfloopMax = input["NselfloopMax"];;

		beta = input["beta"];
		Ntau = input["Ntau"];
		dt = beta/Ntau;
		NWtau = input["NWtau"];
		dtW = beta/Ntau;

		//std::cout<<"check point 0"<<std::endl<<std::flush;

		//W = input["W"];
		epsilon = input["epsilon"];
		B_field = input["B_field"];
		//w_boson = input["w_boson"];
		g_coupling = input["g_coupling"];
		Omega11 = input["Omega11"];
		omega_c = input["omega_c"];
		Nkx = input["Nkx"];

		//std::cout<<"check point 1"<<std::endl<<std::flush;

		if(input.find("is_mixing")!=input.end()) is_mixing = input["is_mixing"];
		else is_mixing = true;
		if(input.find("does_printTt")!=input.end()) does_printTt = input["does_printTt"];
		else does_printTt = false;
		if(input.find("does_printQt")!=input.end()) does_printQt = input["does_printQt"];
		else does_printQt = false;

		//if(input.find("Nw")!=input.end()) Nw = input["Nw"];
		//if(input.find("Lx")!=input.end()) Lattice->setLgrid(input["Lx"]);

		Ri_N.resize(Nordermax+1);
		if(input.find("InverseReweight")!=input.end()){
			if(input["InverseReweight"].size()!=Nordermax+1){
				std::cout<<"** THE NUMBER OF REWEIGHT ELEMENTS DOESN'T MATCH WITH NORDERMAX **"<<std::endl;
				exit(EXIT_FAILURE);
			}
			double Ri_N0 = input["InverseReweight"][0];
			for(int i=0; i<Nordermax+1; i++) Ri_N[i] = static_cast<double>(input["InverseReweight"][i])/Ri_N0;
		}
		else for(int i=0; i<Nordermax+1; i++) Ri_N[i] = 1.;

		//Nx = Lattice->getNgrid();
		//Lx = Lattice->getLgrid();
		//dim = Lattice->getdim();

		//if(site_generate_radius>Lx-1) site_generate_radius = Lx-1;
		//bidist = std::binomial_distribution<int>(site_generate_radius,0.5);
		//biprobability = construct_binormal_probability(site_generate_radius);
		
		//U = input["U"];
		//n0 = input["n0"];
		//mu = input["mu"];

		//Gt.initGt_atomic();
		//Wt.initWt_spinBoson();

		/*
		lattice_para.push_back(mu);
		Lattice->set_para(lattice_para);

		Gt = new green_t[Nflavor];
		Wt = new interaction_t[Nflavor];
		for(int iflavor=0; iflavor<Nflavor; iflavor++){
			//Gt[iflavor].initGt0(*Lattice);
			//Wt[iflavor].set_const(-1.0);
			Gt[iflavor].initGt_spinBoson();
			Wt[iflavor].initWt_spinBoson();
		}
		*/

		Nupdate = 5;

		Nth = input["Nth"];
		NmcTot = input["Nmc"];
		if(input.find("Nscl")!=input.end()) Nscl = static_cast<int>(input["Nscl"]);

		int nprocs = MPI::COMM_WORLD.Get_size();
                Nerr    = static_cast<int>( ceil( 10./nprocs) );
                NerrTot = Nerr*nprocs;
                Nmc     = NmcTot/NerrTot;
                //Nth     = std::min(Nmc,static_cast<long long>(NthMax));

		int global_size;
		std::string global_folder_name;
		if(rank==root){
			global.open(FolderName());
			global.mkdir();
			global_folder_name = global.getdir();
			global_size = global_folder_name.size();
		}
		MPI::COMM_WORLD.Bcast(&global_size,1,MPI::INT,root);
		char* c_str_folder_name = new char[global_size+1];
		if(rank==root)
			for(int i=0; i<global_size+1; i++) c_str_folder_name[i] = global_folder_name.c_str()[i];
		MPI::COMM_WORLD.Bcast(c_str_folder_name,global_size+1,MPI::CHAR,root);
		if(rank!=root){
			global.open(std::string(c_str_folder_name));
		}
		delete[] c_str_folder_name;

		//std::stringstream ss;
		//ss<<global.getdir()<<"/log.dat";
		//std::ofstream ofstr(ss.str().c_str());
		//print(ofstr);
		//ofstr.close();
	}
	std::string FolderName(){
		std::stringstream ss(std::stringstream::in|std::stringstream::out);
		ss.setf(std::ios::fixed, std::ios::floatfield); 
		ss<<"dir_Nordermax"<<Nordermax;
		//ss<<"_Lx"<<Lx;
		//ss.precision(3); ss<<"_n0_"<<n0;
		ss.precision(4); ss<<"_beta"<<beta;
		ss.precision(4); ss<<"_ep"<<epsilon;
		ss.precision(4); ss<<"_B"<<B_field;
		ss.precision(4); ss<<"_g"<<std::scientific<<g_coupling;
		ss.precision(4); ss<<"_Omega11_"<<Omega11;
		ss.precision(4); ss<<"_omega_c"<<omega_c;
		ss<<"_Ntau"<<Ntau;
		ss<<"_NWtau"<<NWtau;
		ss<<"_Nkx"<<Nkx;
		ss<<"_Nmc"<<NmcTot;
		ss<<"_Nscl"<<Nscl;

		return	ss.str();
	};
	void print(std::ostream& ostr){
		ostr<<"** information of parameters **"<<std::endl
			//<<"on "<<Lx<<" x "<<Lx<<" lattice"<<std::endl
			<<"Nordermax = "<<Nordermax<<std::endl;
			//<<"inverse reweight factor by order: "<<"{";
		//for(int i=0; i<Nordermax+1; i++) ostr<<"order "<<i<<": "<<Ri_N[i]<<"; "; ostr<<"}"<<std::endl;
		//ostr<<"NselfloopMax = "<<NselfloopMax<<std::endl
			
		ostr<<"beta = "<<beta<<std::endl
			<<"epsilon = "<<epsilon<<std::endl
			<<"B_field = "<<B_field<<std::endl
			<<"g_coupling = "<<g_coupling<<std::endl
			<<"Omega11 = "<<Omega11<<std::endl
			<<"omega_c = "<<omega_c<<std::endl
			<<"Nkx = "<<Nkx<<std::endl
			<<"Ntau = "<<Ntau<<std::endl
			<<"NWtau = "<<NWtau<<std::endl
			<<"Nth = "<<Nth<<std::endl
			<<"Nmc = "<<NmcTot<<" ("<<Nmc<<")"<<std::endl
			<<"Nscl = "<<Nscl<<std::endl
			<<"*******************************"<<std::endl;
			//<<"W = "<<W<<std::endl
			//<<"n0 = "<<n0<<std::endl
	};
	/*std::vector<double> construct_binormal_probability(const int& t_i){
		std::vector<double> probvec(t_i+1);

		double sum = 0.;
		probvec[0] = 1.;
		sum += probvec[0];
		for(int n=1; n<t_i+1; n++){
			probvec[n] = probvec[n-1]*(t_i+1-n)/n;
			sum += probvec[n];
		}
		for(int n=0; n<t_i+1; n++) probvec[n] /= sum;

		//std::cout<<"* probability distribution for site selection *"<<std::endl;
		//for(int n=0; n<t_i+1; n++) std::cout<<probvec[n]<<" ";
		//std::cout<<std::endl;
		//sum = 0.;
		//for(int n=0; n<t_i+1; n++){ sum += probvec[n]; }
		//std::cout<<"probability normalization : "<<std::setprecision(10)<<sum<<std::endl;

		return probvec;
	};*/
};
