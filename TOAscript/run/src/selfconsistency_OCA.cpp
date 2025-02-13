#include <fstream>
#include <sstream>
#include "dyson.hpp"
#include "json.hpp"
using namespace std;

int main(int argc, char* argv[]){
	cout<<"*********************"<<endl;
	cout<<"** OCA calculation **"<<endl;
	cout<<"*********************"<<endl;
	nlohmann::json input;

	ifstream para(argv[1]);
	para >> input;
	para.close();

	int Ntau = input["Ntau"];
	double beta = input["beta"];
	string Winput = input["Winput"];
	string HlocInput = input["HlocInput"];
	string NmatrixInput = input["NmatrixInput"];
	int Nscl = input["Nscl"];
	bool is_mixing = input["is_mixing"];

	cout<<"** parameters **"<<endl;
	cout<<"   Ntau = "<<Ntau<<endl;
	cout<<"   beta = "<<beta<<endl;
	cout<<"   Nscl = "<<Nscl<<endl;
	cout<<"****************"<<endl;

	dyson DysonEquation(3, beta, Ntau);

	ifstream ifstr( Winput.c_str() );
	DysonEquation.readWt( ifstr );
	ifstr.close();

	ofstream ofstr;
	ofstr.open("Wt.in");
	DysonEquation.printWt(ofstr);
	ofstr.close();

	ifstr.open( HlocInput.c_str() );
	DysonEquation.readHloc( ifstr );
	ifstr.close();

	ofstr.open("Hloc.in");
	DysonEquation.printHloc(ofstr);
	ofstr.close();

	ifstr.open( NmatrixInput.c_str() );
	DysonEquation.readVOp( ifstr );
	ifstr.close();

	ofstr.open("VOp.in");
	DysonEquation.printVOp(ofstr);
	ofstr.close();


	//DysonEquation.initWt(g_coupling,w_boson);
	DysonEquation.initGt();
	DysonEquation.normalizeGt();

	ofstr.open("Gt.in");
	DysonEquation.printGt(ofstr);
	ofstr.close();

	stringstream ss;
	ss.str(""); ss<<"G05beta.iter";
	ofstream ofstr_iter(ss.str().c_str());
	cout<<"iteration ("<<Nscl<<"): "<<flush;
	for(int i=0; i<Nscl; i++){
		cout<<i<<flush;
		DysonEquation.Tupdate_OCA();

		DysonEquation.updateSt();
		DysonEquation.updateGt_volterra();

		DysonEquation.normalizeGt();
		//if( is_mixing ) DysonEquation.mixingGt();

		DysonEquation.printG05beta(ofstr_iter,i);
	}
	cout<<endl<<flush;
	ofstr_iter.close();

	ofstr.open("St.dat");
	DysonEquation.printSt(ofstr);
	ofstr.close();

	ofstr.open("Gt.dat");
	DysonEquation.printGt(ofstr);
	ofstr.close();

    Eigen::MatrixXd lambda1 = Eigen::MatrixXd::Zero(3,3);
    lambda1(0,1) = 1.0; lambda1(1,0) = 1.0;

	DysonEquation.Tupdate_OCA(lambda1);
	DysonEquation.calculateStS0(lambda1);
	ofstr.open("StS0.dat");
	DysonEquation.printStS0(ofstr);
	ofstr.close();


	/*DysonEquation.calculateStS0_NCA();
	ofstr.open("StS0.dat");
	DysonEquation.printStS0(ofstr);
	ofstr.close();*/

	/*DysonEquation.initGtexact(epsilon, g_coupling, w_boson, 100);
	ofstr.open("Gt.exact");
	DysonEquation.printGt(ofstr);
	ofstr.close();*/

	return 0;
};
