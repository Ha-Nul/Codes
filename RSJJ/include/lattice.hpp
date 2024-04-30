#ifndef __lattice_hpp__
#define __lattice_hpp__

#include <iostream>
#include <vector>
#include <Eigen/Dense>

class lattice {
	protected:
	int dim;
	int Lgrid, Ngrid;

	std::vector< Eigen::VectorXi > xunit;
	std::vector< Eigen::VectorXi > xgrid;
	std::vector< Eigen::VectorXd > kunit;
	std::vector< Eigen::VectorXd > kgrid;

	Eigen::MatrixXd kernel;
	std::vector<int> XindexTable;

	std::vector< double > energy;
	std::vector< double > weight;

	public:
	lattice();
	lattice( const int& dim_i, const int& Lgrid_i );
	void set( const int& dim_i, const int& Lgrid_i );
	void setdim( const int& d_i );
	//void setLgrid( const int& N_i );
	virtual void setLgrid( const int& Lgrid_i );
	void setNgrid( const int& N_i );
	int getdim() const;
	int getLgrid() const;
	int getNgrid() const;
	int getXindex(const Eigen::VectorXi& dx_i) const;
	std::vector< Eigen::VectorXd > getkunit() const;
	std::vector< Eigen::VectorXi > getxgrid() const;
	std::vector< Eigen::VectorXd > getkgrid() const;
	Eigen::VectorXi getxgrid(const int& ix) const;
	Eigen::VectorXd getkgrid(const int& ik) const;
	std::vector< double > getenergy() const;
	std::vector< double > getweight() const;
	void pushxunit( const Eigen::VectorXi& xunit_i );
	void pushxgrid( const Eigen::VectorXi& xgrid_i );
	void pushkunit( const Eigen::VectorXd& kunit_i );
	void pushkgrid( const Eigen::VectorXd& kgrid_i );
	void pushenergy( const double& energy_i );
	void pushweight( const double& weight_i );
	void setweight( std::vector<double> weight_vi );
	void setkernel();
	void setindex_table();
	Eigen::MatrixXd getkernel() const;
	double genP0norm(double beta_i);
	void printkgrid( std::ostream& ostr );
	void printenergy( std::ostream& ostr );
	virtual double dispersion(Eigen::VectorXd k) const;
	virtual void set_para(std::vector<double>& para_i);
};

class spinlessatom: public lattice{
	double eploc;

	public:
	spinlessatom();
	void set_para(const std::vector<double>& para_i);
	void setLgrid(const int& lgrid_i);
};

/*class hubbardatom: public lattice{
	double eploc;

	public:
	hubbardatom();
	void setLgrid(const int& lgrid_i);
};


class squarelattice: public lattice{
	double	t_nn;
	double	t_nnn;

	public:
	squarelattice();
	squarelattice(const int& lgrid_i);
	void setLgrid(const int& Lgrid_i);
	void sett_nnn(const int& t_nnn_i);
	double dispersion(Eigen::VectorXd k) const;
};*/

#endif
