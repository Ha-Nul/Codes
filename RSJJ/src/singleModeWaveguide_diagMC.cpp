#include <mpi.h>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include "parameters.hpp"
#include "diagram.hpp"
#include "dyson.hpp"
//#include "configuration.hpp"
//#include "green.hpp"
//#include "controller.hpp"
using namespace std;
using namespace parameters;

int main(int argc, char* argv[]){
	MPI::Init(argc,argv);

       	const int root = 0;
	int rank = MPI::COMM_WORLD.Get_rank();
	//int nprocs = MPI::COMM_WORLD.Get_size();

	switch( argc ){
		case 1:
			if(rank==root) cout<<"Usage: "<<argv[0]<<" input.json"<<endl;
			return 0;
			break;
		case 2:
			break;
		default:
			if(rank==root) cout<<"wrong input!"<<endl;
			return 0;
			break;
	}

	if(rank==root){
		cout<<"************************"<<endl;
		cout<<"** DiagMC calculation **"<<endl;
		cout<<"************************"<<endl;
	}

	parameters::ReadData(argv);
	
	if(rank==root){
		parameters::print(cout);
	}

	dyson DysonEquation(parameters::beta, parameters::Ntau, parameters::epsilon, parameters::B_field);
	//DysonEquation.initWt(parameters::g_coupling,parameters::w_boson);
	DysonEquation.initWtSingleModeWaveguide(parameters::g_coupling,parameters::Omega11,parameters::omega_c,parameters::Nkx);
	DysonEquation.initGt();
	DysonEquation.normalizeGt();

	if( parameters::input.find("Gt_input")!=input.end() ){
		string Gt_in_str = input["Gt_input"];
		ifstream ifstr( Gt_in_str.c_str() );
		DysonEquation.readGt( ifstr );
		ifstr.close();
	}

	ofstream ofstr, ofstr_iter;
	stringstream ss;

	const int Nprint = 5;
	if(rank==root){
		ss.str(""); ss<<parameters::global.getdir()<<"/Wt.in";
		ofstr.open( ss.str().c_str() );
		DysonEquation.printWt(ofstr);
		ofstr.close();

		ss.str(""); ss<<parameters::global.getdir()<<"/Gt.in";
		ofstr.open( ss.str().c_str() );
		DysonEquation.printGt(ofstr);
		ofstr.close();

		ss.str(""); ss<<parameters::global.getdir()<<"/G05beta.iter";
		ofstr_iter.open(ss.str().c_str());
		cout<<"iteration ("<<Nscl<<"): "<<flush;
	}

	for(int i=0; i<Nscl; i++){
		if(rank==root) cout<<i<<flush;
		//DysonEquation.initQt2ndOrder();
		DysonEquation.diagMC_Qt(Nordermax);

		DysonEquation.Tupdate();
		DysonEquation.updateSt();
		DysonEquation.updateGt_volterra();
		//DysonEquation.updateGt_dyson(epsilon);

		DysonEquation.normalizeGt();
		if( is_mixing ) DysonEquation.mixingGt();

		if(rank==root){
			DysonEquation.printG05beta(ofstr_iter,i);

		if(rank==root){
			DysonEquation.printG05beta(ofstr_iter,i);

			ss.str(""); ss<<parameters::global.getdir()<<"/Gt.dat"<<i;
			ofstr.open(ss.str().c_str());
			DysonEquation.printGt(ofstr);
			ofstr.close();

			DysonEquation.calculateStS0();
			ss.str(""); ss<<parameters::global.getdir()<<"/StS0.dat"<<i;
			ofstr.open(ss.str().c_str());
			DysonEquation.printStS0(ofstr);
			ofstr.close();

			if(i>Nprint-1){
				ss.str(""); ss<<parameters::global.getdir()<<"/Gt.dat"<<i-Nprint;
				remove( ss.str().c_str() );

				ss.str(""); ss<<parameters::global.getdir()<<"/StS0.dat"<<i-Nprint;
				remove( ss.str().c_str() );
			}

		}
		}
	}
	if(rank==root){
		cout<<endl<<flush;
		ofstr_iter.close();
	}


	if(rank==root){
		cout<<endl<<flush;
		ofstr_iter.close();

		if( does_printQt ){
			int Nflavor = 2;
			for(int i=0; i<Nflavor*Nflavor; i++) for(int j=0; j<Nflavor*Nflavor; j++){
				for(int itn=0; itn<Ntau+1; itn++){
					ss.str(""); ss<<parameters::global.getdir()<<"/Qt"<<itn<<"_"<<i<<j<<".dat";
					ofstr.open(ss.str().c_str());
					DysonEquation.printQt(ofstr,itn,i,j);
					ofstr.close();
				}
			}
		}
		if( does_printTt ){
			for(int i=0; i<Nflavor; i++) for(int j=0; j<Nflavor; j++){
				ss.str(""); ss<<parameters::global.getdir()<<"/Tt"<<i<<j<<".dat";
				ofstr.open(ss.str().c_str());
				DysonEquation.printTt(ofstr,i,j);
				ofstr.close();
			}
		}

		ss.str(""); ss<<parameters::global.getdir()<<"/Qseries.dat";
		ofstr.open(ss.str().c_str());
		DysonEquation.printQseries(ofstr,Nordermax);
		ofstr.close();

		ss.str(""); ss<<parameters::global.getdir()<<"/Qseries_err.dat";
		ofstr.open(ss.str().c_str());
		DysonEquation.printQseries_err(ofstr,Nordermax);
		ofstr.close();

		ss.str(""); ss<<parameters::global.getdir()<<"/St.dat";
		ofstr.open(ss.str().c_str());
		DysonEquation.printSt(ofstr);
		ofstr.close();

		ss.str(""); ss<<parameters::global.getdir()<<"/Gt.dat";
		ofstr.open(ss.str().c_str());
		DysonEquation.printGt(ofstr);
		ofstr.close();

		DysonEquation.calculateStS0();
		ss.str(""); ss<<parameters::global.getdir()<<"/StS0.dat";
		ofstr.open(ss.str().c_str());
		DysonEquation.printStS0(ofstr);
		ofstr.close();
	}

	MPI::Finalize();

	return 0;
};

