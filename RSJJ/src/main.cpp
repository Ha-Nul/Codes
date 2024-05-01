#include <mpi.h>
#include <iostream>
#include <fstream>
#include "parameters.hpp"
#include "diagram.hpp"
#include "green.hpp"
//#include "controller.hpp"
using namespace std;
using namespace parameters;

namespace parameters{
	//squarelattice Lattice(2);
	lattice* Lattice = new spinlessatom;
};

int main(int argc, char* argv[]){
	MPI::Init(argc,argv);
       	const int root = 0;
	int rank = MPI::COMM_WORLD.Get_rank();
	int nprocs = MPI::COMM_WORLD.Get_size();

	switch( argc ){
		case 1:
			cout<<"Usage: "<<argv[0]<<" input.json"<<endl;
			return 0;
			break;
		case 2:
			break;
		default:
			cout<<"wrong input!"<<endl;
			return 0;
			break;
	}

	parameters::ReadData(argv);
	if(rank==root) parameters::print(cout);

	stringstream ss; ss<<parameters::global.getdir()<<"/Gtx.in";
	ofstream ofstr( ss.str().c_str() );
	Gt[0].printX(ofstr);
	ofstr.close();

	ss.str(""); ss<<parameters::global.getdir()<<"/Wtx.in";
	ofstr.open( ss.str().c_str() );
	Wt[0].printX(ofstr);
	ofstr.close();

	diagram Diag;
	//Diag.debug();

	Diag.initDiagram();
	Diag.MCsampling();

	ss.str(""); ss<<parameters::global.getdir()<<"/normalization.dat";
	ofstr.open( ss.str().c_str() );
	Diag.print_normalization(ofstr);
	ofstr.close();

	for(int iorder=1; iorder<Nordermax+1; iorder++){
		for(int it=0; it<Nout+1; it++){
			ss.str(""); ss<<parameters::global.getdir()<<"/Q_o"<<iorder<<"_it"<<it<<".dat";
			ofstr.open( ss.str().c_str() );
			Diag.printQ(ofstr,iorder,it);
			ofstr.close();
		}
	}

	MPI::Finalize();

	return 0;
}
