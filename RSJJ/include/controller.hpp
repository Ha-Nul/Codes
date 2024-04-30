#ifndef __controller_hpp__included__
#define __controller_hpp__included__

#include <iostream>
#include <string>
#include <vector>
#include <Eigen/Dense>
#include "green.hpp"

class controller{
	double sign, dsign;
	double Rproposal, dRproposal;
	double Racceptance, dRacceptance;
	double Rvalid, dRvalid;
	double Rirreducible, dRirreducible;
	double RGcompact, dRGcompact;
	double RWcompact, dRWcompact;

	// 0: self-energy, 1: polarization
	std::vector<double> normalization[2];
	std::vector<double> dnormalization[2];

	std::vector<double> Ploc;
	std::vector<double> dPloc;

	std::vector<std::vector<std::vector<int> > > order_combination;

	green_w ***Gw, ***dGw;
	selfenergy_w ***Sw, ***dSw;
	selfenergy_w ***sw, ***dsw;
	interaction_w ***Ww, ***dWw;
	polarization_w ***Pw, ***dPw;
	polarization_w ***pw, ***dpw;

	green_w **GwSum, **dGwSum;
	selfenergy_w **SwSum, **dSwSum;

	interaction_w Uw, W2w;
	interaction_w **WwSum_prev, **dWwSum_prev;
	interaction_w **WwSum, **dWwSum;
	polarization_w **PwSum, **dPwSum;

	std::vector<green_w> Gw_memory;
	std::vector<interaction_w> Ww_memory;


	public:
	controller();
	void init_order_combination();
	void debug();
	void construct_W2();
	Eigen::MatrixXcd construct_Sw1(const Eigen::MatrixXd& Gt, const Eigen::MatrixXd& Wt);
	void dysonequation_W();
	void iterate_W(const double& alpha_i);
	Eigen::MatrixXd InverseTimeFourier(const Eigen::MatrixXcd& Owx);
	void bolditeration();
	void bolditeration_S();
	void sampling();
	void sampling_S();
	//void iterativesampling();
	void construct_Gw(const int& iorder);
	void selfconsistentloop(const double& alpha);
	void selfconsistentloop(const double& lambda, const double& alpha);

	void printstatistics(std::ostream& ostr);
	void printnormalization(const int& sp, std::ostream& ostr);
	void printPloc(std::ostream& ostr);

	void printSwx(std::ostream& ostr, const int& iorder, const int& iflavor, const int& jflavor);
	void printSwk(std::ostream& ostr, const int& iorder, const int& iflavor, const int& jflavor);
	void printdSwx(std::ostream& ostr, const int& iorder, const int& iflavor, const int& jflavor);
	void printdSwk(std::ostream& ostr, const int& iorder, const int& iflavor, const int& jflavor);
	void printPwx(std::ostream& ostr, const int& iorder, const int& iflavor, const int& jflavor);
	void printPwk(std::ostream& ostr, const int& iorder, const int& iflavor, const int& jflavor);
	void printdPwx(std::ostream& ostr, const int& iorder, const int& iflavor, const int& jflavor);
	void printdPwk(std::ostream& ostr, const int& iorder, const int& iflavor, const int& jflavor);

	void printswx(std::ostream& ostr, const int& iorder, const int& iflavor, const int& jflavor);
	void printswk(std::ostream& ostr, const int& iorder, const int& iflavor, const int& jflavor);
	void printdswx(std::ostream& ostr, const int& iorder, const int& iflavor, const int& jflavor);
	void printdswk(std::ostream& ostr, const int& iorder, const int& iflavor, const int& jflavor);
	void printpwx(std::ostream& ostr, const int& iorder, const int& iflavor, const int& jflavor);
	void printpwk(std::ostream& ostr, const int& iorder, const int& iflavor, const int& jflavor);
	void printdpwx(std::ostream& ostr, const int& iorder, const int& iflavor, const int& jflavor);
	void printdpwk(std::ostream& ostr, const int& iorder, const int& iflavor, const int& jflavor);

	void printSwkSum(std::ostream& ostr, const int& iflavor, const int& jflavor);
	void printPwkSum(std::ostream& ostr, const int& iflavor, const int& jflavor);
	void printWwkSum_prev(std::ostream& ostr, const int& iflavor, const int& jflavor);
	void printWwkSum(std::ostream& ostr, const int& iflavor, const int& jflavor);
	void printWtx(std::ostream& ostr, const int& dress, const int& iflavor, const int& jflavor);
	void printinput(std::string foldername);
	void printoutput(std::string foldername);
	void printoutput_S(std::string foldername);
	void printiterativeoutput(std::string foldername, const int& iordermax);
};

#endif
