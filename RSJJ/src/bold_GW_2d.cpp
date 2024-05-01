#include <mpi.h>
#include <iostream>
#include <fstream>
#include "parameters.hpp"
#include "controller.hpp"
#include "green.hpp"
using namespace std;
using namespace parameters;

namespace parameters{
	lattice* Lattice = new squarelattice;
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

	controller Manager;
	//Manager.debug();

	Manager.bolditeration();
	//Manager.bolditeration_S();

	MPI::Finalize();

	return 0;
}
