#pragma once

#include <iostream>
#include <string>
#include <list>
#include <vector>
#include <Eigen/Dense>
#include <iterator>
#include <chrono>
#include "configuration.hpp"
#include "element.hpp"
#include "green.hpp"

class diagram{
	int Norder, Nordermax;

	int NtG, NtW;
	double dtau;

	//interaction_t mWt;
	//pseudogreen_t Gt;

	int sign;
	//int Nselfloop;

	configuration Conf;
	//configuration Conf_new;

	//coordinate x_b;
	//std::vector<int> ext_flavor;

	//std::vector<vertex> vertices;
	//std::vector<Eigen::MatrixXd> PGtop, PGbottom;
	//std::vector<double> Ws;

	//Eigen::MatrixXd fullPGtop, fullPGbottom;
	//std::vector<Eigen::MatrixXd>>::iterator PGiter;

	//std::list<vertex> vertices;
	//std::list<line> Gs;
	//std::list<line> Ws;

	////int miflavor, mjflavor;
	//std::list<line>::iterator ref_line, measuring_line;

	////std::list<fermionicloop> fermionicloops;
	////setUWinter SetUWInter;

	//std::vector<std::vector<int> > GadjL, GadjLtmp;
	//std::vector<std::vector<int> > WadjL, WadjLtmp;
	//std::vector<int> WmadjL, WmadjLtmp;

	//// 0: self-energy, 1: polarization
	std::vector<double> normalization;
	////std::vector<long> avesign[2];
	////std::vector<double> Plocacc;
	////double Rproposal, Racceptance;
	////double Rvalid, Rirreducible, RGcompact, RWcompact;

	//Eigen::VectorXcd Owx;

	vertex_4pt Qsum;
	std::vector<vertex_4pt> Qacc;
	//selfenergy_w*** Swxacc;
	//polarization_w*** Pwxacc;

	// monitoring quantities
	double NorderAv, NorderAcc;

	std::chrono::high_resolution_clock::time_point t1;
	std::chrono::high_resolution_clock::time_point t2;
	std::chrono::duration<double> time_span;

	public:
	diagram();
	//~diagram();

	void initDiagram(const std::vector<double>& mWt_i, const std::vector<Eigen::Matrix2d>& Gt_i);
	////void initWt(const std::vector<double>& mWt_i);
	////void initGt(const std::vector<Eigen::Matrix2d>& Gt_i);
	//void debug();

	////double calculate_fullPGtop();
	////double calculate_fullPGbottom();
	////double getInsertionRatio();
	////double getRemovalRatio();

	////int gen_index_removal(const int& r1, const int& r2, const int& o);
	////void initDiagram();
	//////void init_S_diagram();
	//////void init_P_diagram();
	////void setNordermax(const int& Nordermax_i);
	void select_two_index(const int& N_i, int& i_o, int& j_o);
	void select_one_excluding_one(const int& N_i, const int& e_i, int& i_o);
	void select_one_excluding_two(const int& N_i, const int& e_i, const int& e_j, int& i_o);
	int insertVertices();
	int removeVertices();
	//int swapWline();
	int shiftIm2();
	int shiftVertex();
	int flipExtFlavor();
	////int swapMeasuringLine();
	//////int reconnect();
	//////int intraswap();
	//////int interswap();
	//////int transformUW();
	//////int loopspinflip();
	//////int absorbbubble();
	//////int ejectbubble();
	//////int spinflip();
	//////void setSwx();
	//////void setPwx();
	void measure();
	//////void measure_S();
	//////void measure_S(const int& iorder);
	//////void measure_P();
	//////void measure_P(const int& iorder);
	////void debug();
	void thermalization(const int& Nth_i);
	void MCsampling(const int& Nmc_i);
	//////void diagMC_S();
	//////void diagMC_S(int iorder);
	//////void diagMC_P();
	//////void diagMC_P(int iorder);
	//////void diagMC();
	//////void diagMC(int iorder);
	//////void bold_diagMC_G();
	//////void bold_diagMC_GW();
	//////void bold_diagMC_GW_S();
	void normalize();
	void partialSum();
	//////void normalize_S(const int& reftype, const int& reforder, const double& norm_i);
	//////void normalize_S(const int& reftype, const int& reforder, const double& norm_i, const int& iorder);
	//////void normalize_P(const int& reftype, const int& reforder, const double& norm_i);
	//////void normalize_P(const int& reftype, const int& reforder, const double& norm_i, const int& iorder);
	//////double getAvesign();
	//////double getRproposal();
	//////double getRacceptance();
	//////double getRvalid();
	//////double getRirreducible();
	//////double getRGcompact();
	//////double getRWcompact();
	//////long getnorm(const int& type_index, const int& order_index);
	//////double getPloc(const int& order_index);
	//////double getPlocSum(const int& order_index);
	//////selfenergy_w getSw(const int& order_index, const int& iflavor_i, const int& jflavor_i);
	//////selfenergy_w getSwSum(const int& order_index, const int& iflavor_i, const int& jflavor_i);
	//////polarization_w getPw(const int& order_index, const int& iflavor_i, const int& jflavor_i);
	//////polarization_w getPwSum(const int& order_index, const int& iflavor_i, const int& jflavor_i);
	double trapezoidalWeight(const int& i, const int& );
	std::vector<double> getQseries(const int& iout_i, const int& iin_i);
	vertex_4pt getQt();
	////Eigen::MatrixXd getQt(const int& n);
	////void printlog(std::ostream& ostr);
	void printconfiguration(std::ostream& ostr);
	////void print_normalization(std::ostream& ostr);
	////void printQ(std::ostream& ostr, const int& iorder, const int& it);
	//////void printSwx(std::ostream& ostr, const int& iorder, const int& iflavor, const int& jflavor);
	//////void printSwk(std::ostream& ostr, const int& iorder, const int& iflavor, const int& jflavor);
	//////void printPwx(std::ostream& ostr, const int& iorder, const int& iflavor, const int& jflavor);
	//////void printPwk(std::ostream& ostr, const int& iorder, const int& iflavor, const int& jflavor);
	//////void printoutput(std::string foldername);
	//////void printoutput_S(std::string foldername);
	//////void printoutput_P(std::string foldername);
	//////void updatebubbles();
	//////void updatefermionicloops();
	//////void updateSetUWInter();
	//////bool checkconnection();
	//////bool checkconnectivity();
	//////bool checkconnectivity(const std::vector<std::vector<int> >& GadjL_i, const std::vector<std::vector<int> >& WadjL_i);
	//////bool checkconnectivityforremove(std::list<line>::iterator Wremoval);
	//////bool checkconnectivityforreconnect(std::list<line>::iterator G_i, std::list<line>::iterator G_j);
	//////bool checkconnectivityforinterswap(std::list<line>::iterator target);
	////bool checkirreducibility();
	//////bool checkirreducibility(std::vector<std::vector<int> > GadjL_i, std::vector<std::vector<int> > WadjL_i);
	//////bool brutecheckirreducibility();
	////bool check_irreducibility_for_remove(std::list<line>::iterator Wremoval);
	//////bool checkirreducibilityforreconnect(std::list<line>::iterator G_i, std::list<line>::iterator G_j);
	////bool check_irreducibility_for_swap(std::list<line>::iterator target);
	//////bool checkirreducibility(const bool& isG_i);
	////bool checkcompactness(const bool& isG_i);
	//////bool checkGcompactness(std::vector<std::vector<int> > GadjL_i, std::vector<int> WmadjL_i, std::vector<std::vector<int> > WadjL_i);
	//////bool checkWcompactness(std::vector<std::vector<int> > GadjL_i, std::vector<int> WmadjL_i, std::vector<std::vector<int> > WadjL_i);
	//////bool brutecheckcompactness(const bool& isG);
	////bool check_compactness_for_remove(std::list<line>::iterator Wremoval);
	//////bool checkcompactnessforreconnect(std::list<line>::iterator G_i, std::list<line>::iterator G_j);
	////bool check_compactness_for_swap(std::list<line>::iterator target);
	//////bool checkbubble();
	//////bool checkbubble(const std::vector<std::vector<int> >& GadjL_i);
	//////bool checkdiagramsign();
	//////bool checkNorder();
	//////bool checkUWlines();
	//////bool checksetUWinters();
	//////bool checkdress();
	//////bool checkNselfloop();
	//////bool checkfermionicloop();
	////bool check_ref_line();
	////bool check_measuring_line();
	//////bool checkdetailedbalance();
	bool checkSign();
	bool checkDetailedBalanceInsertRemove();
	bool checkDetailedBalanceRemoveInsert();
	bool checkDetailedBalanceShiftIm2();
	bool checkDetailedBalanceShiftVertex();
	bool checkDetailedBalanceFlipExtFlavor();
	////bool check_detailed_balance_remove_insert();
	////bool check_detailed_balance_shift();
	//////bool checkdetailedbalance_reconnect();
	//////bool checkdetailedbalance_interswap();
	//////bool checkdetailedbalance_intraswap();
	//////bool checkdetailedbalance_transformUW();
	//////bool checkdetailedbalance_eject_absorb_bubble();
	//////bool checkdetailedbalance_absorb_eject_bubble();
	//////bool checkdetailedbalance_loopspinflip();
	//////bool checkUlines();
	//////bool checkWlines();
};

