#include <iostream>
#include <fstream>
#include <vector>
#include <Eigen/Dense>
#include "green.hpp"

class dyson{
	const int Nflavor;

	int Ntau;
	double beta;
	double dtau;

	//int n,m;
	Eigen::MatrixXd VOp, Hloc;
	Eigen::MatrixXd ZeroM, ZeroV, IdentityM;
	//Eigen::Matrix2d Sz, Hatom, Zero2x2;
	//Eigen::MatrixXd Zero4x1;

	std::vector<double> Wt;

	double lambda, lambda_prev;
	//std::vector<Eigen::Matrix2d> G0t;
	//std::vector<Eigen::Matrix2d> Gt_prev;
	std::vector<Eigen::MatrixXd> Gt, Gt_prev;

	std::vector<Eigen::MatrixXd> St;

	std::vector<double> StS0;

	//std::vector<std::vector<Eigen::Matrix2d> > Tt_amp;
	//std::vector<std::vector<Eigen::Matrix2d> > Tt_amp_up, Tt_amp_dn, Tt_amp_W;

	std::vector<std::vector<Eigen::MatrixXd> > Tt;

	//std::vector<std::vector<std::vector<Eigen::Matrix4d> > > Qt_reducible;
	//std::vector<std::vector<std::vector<Eigen::Matrix4d> > > Qt;

	////vertex_4pt Qt, Qt_err;

	public:
	dyson(const int& Nflavor_i, const double& beta_i, const int& Ntau_i);
	void readHloc(std::ifstream& ifstr);
	void readVOp(std::ifstream& ifstr);
	//int gen4index(const int& a, const int& b);
	//std::vector<int> gen22index(const int& ab);
	//Eigen::MatrixXd transformTo2x2(const Eigen::MatrixXd& V4x1);
	//Eigen::MatrixXd transformTo4x1(const Eigen::MatrixXd& M2x2);
	//Eigen::Matrix4d tensorProd(const Eigen::MatrixXd& M2x2_1, const Eigen::MatrixXd& M2x2_2);
	//void initWt(const double& gCouple, const double& wBoson);
	void readWt(std::ifstream& ifstr);
	//int getSimpsonFactor(const int& ix, const int& Nx);
	//void initWtSingleModeWaveguide(const double& g0, const double& Omega11, const double& omega_c, const int& Nkx);
	void initGt();
	void readGt(std::ifstream& ifstr);
	////void setGtGlobal();
	//////void initGtexact(const double& epsilon, const double& gCouple, const double wBoson, const int& nBcut);
	////void setSt(const double& s_value);
	//void initQtReducible();
	//void initQt2ndOrder();
	//void initQt3rdOrder();
	////void diagMC_Qt(const int& Nordermax);
	//void initTt();
	//void initTt_amp();
	double weight(const int& n, const int& m);
	//double calculateBulkM_TCA(const int& n);
	//double calculateBulkM_Q2(const int& n);
	//void calculateBulkTtAmputatedGUp(const int& m);
	//void calculateBulkTtAmputatedGDn(const int& nm);
	//void calculateBulkTtAmputatedW(const int& m);
	//void calculateBulkQtReducible(const int& n);
	//void calculateBulkTtAmputated(const int& n, const int& m);
	//Eigen::Matrix2d calculateBulkB(const int& n, const int& m);
	Eigen::MatrixXd calculateBulkB_OCA(const int& n, const int& m);
	////double calculateBulkB_Ladder(const int& n, const int& m);
	//Eigen::Matrix2d calculateBulkB_TCA(const int& n, const int& m);
	std::vector<Eigen::MatrixXd> calculateBoundaryB_OCA(const int& n);
	//std::vector<Eigen::Matrix2d> calculateBoundaryB(const int& n);
	//Eigen::Matrix2d calculateBoundaryM(const int& n);
	//void calculateBoundaryTtAmputated(const int& n);
	//void Tupdate();
	////void TupdateBruteForce_OCA();
	////void TupdateBruteForce_Ladder();
	////void TupdateBruteForce_TCA();
	void Tupdate_OCA();
	void Tupdate_OCA(const Eigen::MatrixXd& X);
	//void Tupdate_TCA0();
	//void Tupdate_Q2_0();
	////void Tupdate_Ladder();
	void Tupdate_TOA();
	void Tupdate_TOA(const Eigen::MatrixXd& X);
	//void Tupdate_TCA();
	////void TupdateDebug();
	void updateSt_NCA();
	void updateSt();
	void updateGt_volterra();
	////void updateGt_dyson(const double& epsilon, const int& Niter);
	void normalizeGt();
	//void mixingGt();
	void calculateStS0(const Eigen::MatrixXd& X);
	void calculateStS0_NCA(const Eigen::MatrixXd& X);
	////void calculateStS0_OCA();
	void printWt(std::ostream& ostr);
	void printHloc(std::ostream& ostr);
	void printVOp(std::ostream& ostr);
	void printGt(std::ostream& ostr);
	//void printQt(std::ostream& ostr, const int& n, const int& iflavor, const int& jflavor);
	////void printQt_err(std::ostream& ostr, const int& n);
	//void printTt(std::ostream& ostr, const int& iflavor, const int& jflavor);
	void printSt(std::ostream& ostr);
	void printStS0(std::ostream& ostr);
	void printG05beta(std::ostream& ofstr, const int& iter);
};
