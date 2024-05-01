#pragma once
#include <vector>
#include <Eigen/Dense>
#include "element.hpp"

class interaction_t{
	static const double eps;
	
	int NtW;
	double dt;

	std::vector<double> mWt;

	public:
	interaction_t();
	void init(const std::vector<double>& mWt_i);
	//void initWt_spinBoson();
	void initWt_spinBoson(const int& NtW_i, const double& gCouple, const double& wBoson);
	void initWtSingleModeWaveguide(const int& Ntw_i, const double& g0, const double& Omega11, const double& omega_c, const int& Nkx);
	//double getWt(const double& t_i);
	double getWt(const coordinate& xi, const coordinate& xj) const;
	void print(std::ostream& ostr);
};

class pseudogreen_t{
	static const double eps;

	int NtG;
	double dt;

	Eigen::Matrix2d Sz;
	std::vector<Eigen::Matrix2d> Gt;

	public:
	pseudogreen_t();
	void initOp( const Eigen::MatrixXd& Op_i);
	void init(const std::vector<Eigen::Matrix2d>& Gt_i);
	//void initGt_atomic();
	//Eigen::Matrix2d getG(const double& t_i);
	//Eigen::Matrix2d getPG(const double& t_i);
	Eigen::Matrix2d getGtSz(const coordinate& xi, const coordinate& xj) const;
	Eigen::Matrix2d getSzGtSz(const coordinate& xi, const coordinate& xj) const;
	void print(std::ostream& ostr);
};

class vertex_4pt{
	std::vector<std::vector<std::vector<Eigen::Matrix4d> > > Qt;

	public:
	vertex_4pt();
	//~vertex_4pt();
	void initZero(const int& Nt_i);
	//void accumulate(const int& it, const int& it2, const int& it1, const double& value_i);
	void accumulate(const std::vector<int>& iflavor_i, const std::vector<int>& it_i, const double& value_i);
	//Eigen::MatrixXd getQtx();
	//Eigen::MatrixXd getQtk();
	//void setQtx(const Eigen::MatrixXd& Qtx_i);
	//void setQtk(const Eigen::MatrixXd& Qtk_i);
	//void operator+=(const function_t& ft_i);
	//void setQt(const std::vector<Eigen::MatrixXd>& Qt_i);
	//void setQt(const int& n, const Eigen::MatrixXd& Qt_i);
	//void setQt(const int& n, const double* Qt_i);
	void setQt(const int& it, const int& it2, const int& it1, const Eigen::MatrixXd& Qt_tt2t1);
	void setQt(const int& it, const int& it2, const int& it1, const double* Qt_tt2t1);
	Eigen::Matrix4d getQt(const int& it, const int& it2, const int& it1) const;
	double getQt(const int& it, const int& it2, const int& it1, const int& iout, const int& iin) const;
	void print(std::ostream& ostr,const int& it, const int& iout, const int& iin);
	void operator/=(const double& denominator_i);
	void operator+=(const vertex_4pt& Qt_i);
	vertex_4pt operator*(const vertex_4pt& Qt_i);
};


/*
#ifndef __Green_hpp__included__
#define __Green_hpp__included__

#include <iostream>
#include <vector>
#include <list>
#include <iterator>
#include <complex>
#include <Eigen/Dense>
#include "element.hpp"
#include "lattice.hpp"

class function_t{
	protected:
	static const double eps;
	Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor> Otk;
	Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor> Otx;

	public:
	function_t();
	void initZero();
	Eigen::MatrixXd getOtx();
	Eigen::MatrixXd getOtk();
	void setOtx(const Eigen::MatrixXd& Otx_i);
	void setOtk(const Eigen::MatrixXd& Otk_i);
	void operator+=(const function_t& ft_i);
	void printX(std::ostream& ostr);
};
class green_t: public function_t {
	public:
	green_t();
	void initGt0(const lattice& Lattice_i);
	void initGt_spinBoson();
	void initGtSpinlessAtom(const lattice& Lattice_i);
	//void initGtHubbardAtom(const lattice& Lattice_i);
	double getXvalue(const coordinate& xi, const coordinate& xj);
};

class interaction_t: public function_t {
	public:
	interaction_t();
	void set_const(const double& W_i);
	void initWt_spinBoson();
	//void initWt0();
	//void initWt1(const lattice& Lattice_i);
	//void initWt2HubbardAtom();
	//void initWt_upup_exactHubbardAtom();
	//void initWt_updn_exactHubbardAtom();
	//void initWt_upup_rpaHubbardAtom();
	//void initWt_updn_rpaHubbardAtom();
	double getXvalue(const coordinate& xi, const coordinate& xj);
};

class vertex_4pt{
	std::vector<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> > Qt;

	public:
	vertex_4pt();
	//~vertex_4pt();
	void initZero();
	void accumulate(const int& it, const int& it2, const int& it1, const double& value_i);
	//Eigen::MatrixXd getQtx();
	//Eigen::MatrixXd getQtk();
	//void setQtx(const Eigen::MatrixXd& Qtx_i);
	//void setQtk(const Eigen::MatrixXd& Qtk_i);
	//void operator+=(const function_t& ft_i);
	void setQt(const std::vector<Eigen::MatrixXd>& Qt_i);
	void setQt(const int& n, const Eigen::MatrixXd& Qt_i);
	void setQt(const int& n, const double* Qt_i);
	Eigen::MatrixXd getQt(const int& n) const;
	double getQt(const int& n, const int& i, const int& j) const;
	void print(std::ostream& ostr,const int& it);
	void operator/=(const double& denominator_i);
	void operator+=(const vertex_4pt& vertex_i);
	vertex_4pt operator*(const vertex_4pt& Q1);
};
*/

/*class function_w{
	protected:
	Eigen::Matrix<std::complex<double>,Eigen::Dynamic,Eigen::Dynamic> Owx;
	Eigen::Matrix<std::complex<double>,Eigen::Dynamic,Eigen::Dynamic> Owk;

	public:
	function_w();
	void init();
	void setZero();
	Eigen::Matrix<std::complex<double>,Eigen::Dynamic,Eigen::Dynamic> get();
	Eigen::MatrixXcd getOwx();
	Eigen::MatrixXcd getOwk();
	std::complex<double> getOwx(const int& n_i, const int& x_i);
	std::complex<double> getOwk(const int& n_i, const int& k_i);
	void* getbufferOwx();
	void* getbufferOwk();
	void setOwx(const std::complex<double>* Owx_i);
	void setOwk(const std::complex<double>* Owk_i);
	void setOwx(const Eigen::MatrixXcd& Owx_i);
	void setOwk(const Eigen::MatrixXcd& Owk_i);
	void setOwx(const int n, const int x, const std::complex<double> Owx_i);
	void setOwk(const int n, const int k, const std::complex<double> Owk_i);
	void accumulateOwx(const int& ix, const Eigen::Matrix<std::complex<double>,Eigen::Dynamic,1>& Owx_i);
	void operator+= (const function_w& Ow_i);
	void operator-= (const function_w& Ow_i);
	void operator/= (const double& denominator_i);
	function_w operator* (const double& scalar_i);
	function_w operator* (const function_w& Ow_i);
	void operator*= (const function_w& Ow_i);
	void fourier_xtok();
	void fourier_ktox();
	void print_x(std::ostream& ostr);
	void print_k(std::ostream& ostr);
	void printOwk(std::ostream& ostr);
	void printOwx(std::ostream& ostr);
};
class green_w: public function_w {
	public:
	green_w();
	void initGw0(const lattice& Lattice_i);
	void initGwHubbardAtomexact();
};
class interaction_w: public function_w {
	public:
	interaction_w();
	void initW1();
	void initW2();
};
class selfenergy_w: public function_w{
	public:
	selfenergy_w();
	void init2order(const lattice& Lattice_i);
};
class polarization_w: public function_w{
	public:
	polarization_w();
	void init0orderHubbardAtom();
};*/

