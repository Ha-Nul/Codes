#pragma once

#include <iostream>
#include <string>
#include <random>
#include "json.hpp"
#include "dirmanip.hpp"
//#include "lattice.hpp"
//#include "green.hpp"

namespace parameters{
	extern nlohmann::json input;

	extern const int Dspace;
	extern const int Nptvertex;
	extern const int Nflavor;

	extern int Nordermax;
	//extern int Nordermax, NselfloopMax;
	//extern int Nauxdress;

	extern double beta;
	extern double dt, dtW;
	extern int Ntau;
	extern int NWtau;
	//extern int Nw;
	//extern int Nx, Lx;

	extern double epsilon;
	extern double B_field;
	extern double g_coupling;
	extern double Omega11, omega_c;
	extern int Nkx;
	//extern double w_boson;
	//extern double W;
	//extern double n0;

	extern bool is_mixing;
	extern bool does_printTt, does_printQt;

	//extern int site_generate_radius;
	extern std::mt19937_64 mt;
	extern std::uniform_real_distribution<double> unidist;
	//extern std::binomial_distribution<int> bidist;
	//extern std::vector<double> biprobability;

	extern std::vector<double> Ri_N;

	//extern lattice* Lattice;

	//extern lattice Lattice;
	//extern hubbardatom Lattice;
	//extern squarelattice Lattice;

	//extern pseudogreen_t Gt;
	//extern interaction_t Wt;

	extern int Nupdate;
	extern int Nth;
	extern const int Nmeasureinterval;
	extern const int NthMax;
	extern long long Nmc, NmcTot;
	extern int Nerr, NerrTot;
	extern int Nscl;

	extern Folder global;

	void MtInit();
	void MtInit(const unsigned long long& seed_i);
	void MtInit_systime(const unsigned long long& seed_i);
	int getSimpsonFactor(const int& ix, const int& Nx);
	//void select_two_index(const int& N_i, int& i_o, int& j_o);
	//void select_one_excluding_one(const int& N_i, const int& e_i, int& i_o);
	//void select_one_excluding_two(const int& N_i, const int& e_i, const int& e_j, int& i_o);
	//int gen_index_merge(const int& m1, const int& m2, const int& o);
	//int gen_index_split(const int& m, const int& o);
	//double w(const int& n_i);
	//double wB(const int& n_i);
	//double nF(const double& ep);
	//double nB(const double& ep);
	void ReadData(char* argv[]);
	std::string FolderName();
	void print(std::ostream& ostr);
	//std::vector<double> construct_binormal_probability(const int& t_i);
};
