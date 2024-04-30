#include "lattice.hpp"
#include <iomanip>
#include <cmath>
using namespace std;
	
////////////////////////////////////
// general lattice implementation //
///////////////////////////////////

lattice::lattice(){};
lattice::lattice( const int& dim_i, const int& Lgrid_i ) : dim(dim_i), Lgrid(Lgrid_i), Ngrid(pow(Lgrid_i,dim_i)) {};
void lattice::set( const int& dim_i, const int& Lgrid_i ){
	dim = dim_i;
	Lgrid = Lgrid_i;
	Ngrid = pow(Lgrid,dim);
};
void lattice::setdim( const int& d_i ){ dim = d_i; };
void lattice::setLgrid( const int& Lgrid_i ){ 
	Lgrid = Lgrid_i; 
	Ngrid = pow(Lgrid,dim);
};
void lattice::setNgrid( const int& N_i ){ Ngrid = N_i; };
int lattice::getdim() const { return dim; };
int lattice::getLgrid() const { return Lgrid; };
int lattice::getNgrid() const { return Ngrid; };
int lattice::getXindex(const Eigen::VectorXi& dx_i) const {
	int ix = 0, factor = 1, dx_abs;
	for(int i=dx_i.size()-1; i>-1; i--){
		dx_abs = (dx_i(i) + Lgrid)%Lgrid;
		ix += dx_abs*factor;
		factor *= Lgrid;
	}
	return ix;
};
std::vector< Eigen::VectorXd > lattice::getkunit() const { return kunit; };
std::vector< Eigen::VectorXi > lattice::getxgrid() const { return xgrid; };
std::vector< Eigen::VectorXd > lattice::getkgrid() const { return kgrid; };
Eigen::VectorXi lattice::getxgrid(const int& ix) const { return xgrid[ix]; };
Eigen::VectorXd lattice::getkgrid(const int& ik) const { return kgrid[ik]; };
std::vector<double> lattice::getenergy() const { return energy; };
std::vector<double> lattice::getweight() const { return weight; };
void lattice::pushxunit( const Eigen::VectorXi& xunit_i ){ xunit.push_back( xunit_i ); };
void lattice::pushxgrid( const Eigen::VectorXi& xgrid_i ){ xgrid.push_back( xgrid_i ); };
void lattice::pushkunit( const Eigen::VectorXd& kunit_i ){ kunit.push_back( kunit_i ); };
void lattice::pushkgrid( const Eigen::VectorXd& kgrid_i ){ kgrid.push_back( kgrid_i ); };
void lattice::pushenergy( const double& energy_i ){ energy.push_back( energy_i ); };
void lattice::setweight( std::vector<double> weight_vi ){ weight = weight_vi; };
void lattice::setkernel(){
	int Nx = xgrid.size(); int Nk = kgrid.size();

	std::vector<Eigen::VectorXd> xgridtmp(xgrid.size());
	for(int ix=0; ix<xgrid.size(); ix++) xgridtmp[ix] = xgrid[ix].cast<double>();

	kernel.resize(Nx,Nk);
	for(int i=0; i<Nx; i++) for(int j=0; j<Nk; j++) 
		kernel(i,j) = cos(kgrid[j].dot(xgridtmp[i]));

	kernel /= static_cast<double>(Nk);
};
Eigen::MatrixXd lattice::getkernel() const { return kernel; };
double lattice::genP0norm(double beta_i){
	double tmp, exptmp;

	tmp = 0.;
	for(int ik=0; ik<energy.size(); ik++){
		exptmp = exp(beta_i*energy[ik]);
		tmp += exptmp/(exptmp+1.)/(exptmp+1.);
	}
	return beta_i*tmp/energy.size();
};
void lattice::printkgrid( std::ostream& ostr ){
	for(vector< Eigen::VectorXd >::iterator iter=kgrid.begin(); iter!=kgrid.end(); ++iter){
		ostr << "(";
		for(int i=0; i<dim; i++)
			ostr << (*iter)(i) << ",";
		ostr << ")" << endl;
	}
};
void lattice::printenergy( std::ostream& ostr ){
	for(int i=0; i<energy.size(); i++){
		ostr << setw(15) << energy[i] << endl;
		//ostr << setw(15) << energy[i] 
		//	<< setw(15) << weight[i] <<endl;
	}
};
void lattice::set_para(std::vector<double>& para_i){};
double lattice::dispersion(Eigen::VectorXd k) const { return 0.; };


//////////////////////////////////
// spinless atom implementation //
//////////////////////////////////

spinlessatom::spinlessatom() : eploc(0.){
	set(1,1);
	pushxgrid( Eigen::Vector2i(0,0) );
	pushkgrid( Eigen::Vector2d(0.,0.) );
	pushenergy(eploc);
	setkernel();
};
void spinlessatom::set_para(const std::vector<double>& para_i){
	// para_i[0] = mu
	energy[0] = -para_i[0];
};
void spinlessatom::setLgrid( const int& Lgrid_i ){};


/*
/////////////////////////////////
// Hubbard atom implementation //
/////////////////////////////////

hubbardatom::hubbardatom()
: eploc(0.){
	set(1,1);
	pushxgrid( Eigen::Vector2i(0,0) );
	pushkgrid( Eigen::Vector2d(0.,0.) );
	pushenergy(eploc);
	setkernel();
};
void hubbardatom::setLgrid( const int& Lgrid_i ){};




///////////////////////////////////
// square lattice implementation //
///////////////////////////////////

squarelattice::squarelattice()
: t_nn(1.0), t_nnn(0.0) {
	int lgrid = 10;

	setdim(2);
	Lgrid = lgrid;
	Ngrid = pow(Lgrid,dim);

	Eigen::Vector2d k1(2.*M_PI/lgrid,0.0);
	Eigen::Vector2d k2(0.0,2.*M_PI/lgrid);
	
	pushkunit( Eigen::Vector2d (2.*M_PI,0.0) );
	pushkunit( Eigen::Vector2d (0.0,2.*M_PI) );

	Eigen::Vector2d k;
	for(int i=0; i<lgrid; i++) for(int j=0; j<lgrid; j++){
		pushxgrid( Eigen::Vector2i(i,j) );

		k = (-lgrid/2+i)*k1 + (-lgrid/2+j)*k2;
		pushkgrid(k);
		pushenergy(dispersion(k));
	}
	setkernel();
};
squarelattice::squarelattice(const int& lgrid_i)
: t_nn(1.0), t_nnn(0.0) {
	int lgrid = lgrid_i + (lgrid_i&1);

	setdim(2);
	Lgrid = lgrid;
	Ngrid = pow(Lgrid,dim);

	Eigen::Vector2d k1(2.*M_PI/lgrid,0.0);
	Eigen::Vector2d k2(0.0,2.*M_PI/lgrid);
	
	pushkunit( Eigen::Vector2d (2.*M_PI,0.0) );
	pushkunit( Eigen::Vector2d (0.0,2.*M_PI) );

	Eigen::Vector2d k;
	for(int i=0; i<lgrid; i++) for(int j=0; j<lgrid; j++){
		pushxgrid( Eigen::Vector2i(i,j) );

		k = (-lgrid/2+i)*k1 + (-lgrid/2+j)*k2;
		pushkgrid(k);
		pushenergy(dispersion(k));
	}
	setkernel();
};
void squarelattice::setLgrid(const int& Lgrid_i){
	int lgrid = Lgrid_i + (Lgrid_i&1);

	setdim(2);
	Lgrid = lgrid;
	Ngrid = pow(Lgrid,dim);

	Eigen::Vector2d k1(2.*M_PI/lgrid,0.0);
	Eigen::Vector2d k2(0.0,2.*M_PI/lgrid);
	
	xunit.clear();
	xgrid.clear();
	kunit.clear();
	kgrid.clear();
	energy.clear();
	weight.clear();

	pushkunit( Eigen::Vector2d (2.*M_PI,0.0) );
	pushkunit( Eigen::Vector2d (0.0,2.*M_PI) );

	Eigen::Vector2d k;
	for(int i=0; i<lgrid; i++) for(int j=0; j<lgrid; j++){
		pushxgrid( Eigen::Vector2i(i,j) );

		k = (-lgrid/2+i)*k1 + (-lgrid/2+j)*k2;
		pushkgrid(k);
		pushenergy(dispersion(k));
	}
	setkernel();
};
void squarelattice::set_t_nnn(const int& t_nnn_i){
	t_nnn = t_nnn_i;
};
double squarelattice::dispersion(Eigen::VectorXd k) const {
	//return -2.*t_nn*(cos(k(0)) + cos(k(1)));
	cosk0 = cos(k(0)); cosk1 = cos(k(1));
	return -2.*t_nn*(cosk0 + cosk1) - 4.t_nnn*cosk0*cosk1;
};
*/
