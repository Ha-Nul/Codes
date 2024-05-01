#include <fstream>
#include <sstream>
#include "dyson.hpp"
#include "json.hpp"
using namespace std;

int main(int argc, char* argv[]){
	cout<<"**********************"<<endl;
	cout<<"** Q2v calculation **"<<endl;
	cout<<"**********************"<<endl;
	nlohmann::json input;

	ifstream para(argv[1]);
	para >> input;
	para.close();

	int Ntau = input["Ntau"];
	double beta = input["beta"];
	double epsilon = input["epsilon"];
	double B_field = input["B_field"];
	double g_coupling = input["g_coupling"];
	double w_boson = input["w_boson"];
	int Nscl = input["Nscl"];
	bool is_mixing = input["is_mixing"];

	cout<<"** parameters **"<<endl;
	cout<<"   Ntau = "<<Ntau<<endl;
	cout<<"   beta = "<<beta<<endl;
	cout<<"   epsilon = "<<epsilon<<endl;
	cout<<"   B_field = "<<B_field<<endl;
	cout<<"   g_coupling = "<<g_coupling<<endl;
	cout<<"   w_boson = "<<w_boson<<endl;
	cout<<"   Nscl = "<<Nscl<<endl;
	cout<<"****************"<<endl;

	dyson DysonEquation(beta, Ntau, epsilon, B_field);

	DysonEquation.initWt(g_coupling,w_boson);
	DysonEquation.initGt();
	DysonEquation.normalizeGt();

	if( input.find("Gt_input")!=input.end() ){
		string Gt_in_str = input["Gt_input"];
		ifstream ifstr( Gt_in_str.c_str() );
		DysonEquation.readGt( ifstr );
		ifstr.close();
	}

	ofstream ofstr;
	ofstr.open("Wt.in");
	DysonEquation.printWt(ofstr);
	ofstr.close();

	ofstr.open("Gt.in");
	DysonEquation.printGt(ofstr);
	ofstr.close();

	stringstream ss;
	ss.str(""); ss<<"G05beta.iter";
	ofstream ofstr_iter(ss.str().c_str());
	cout<<"iteration ("<<Nscl<<"): "<<flush;
	for(int i=0; i<Nscl; i++){
		cout<<i<<flush;
		DysonEquation.initQt2ndOrder();

		DysonEquation.Tupdate();
		DysonEquation.updateSt();
		DysonEquation.updateGt_volterra();
		//DysonEquation.updateGt_dyson(epsilon);

		DysonEquation.normalizeGt();
		if( is_mixing ) DysonEquation.mixingGt();

		DysonEquation.printG05beta(ofstr_iter,i);
	}
	cout<<endl<<flush;
	ofstr_iter.close();

	int Nflavor = 2;
	for(int i=0; i<Nflavor; i++) for(int j=0; j<Nflavor; j++){
		/*for(int itn=0; itn<Ntau+1; itn++){
			ss.str(""); ss<<"Qt"<<itn<<"_"<<i<<j<<".dat";
			ofstr.open(ss.str().c_str());
			DysonEquation.printQt(ofstr,itn,i,j);
			ofstr.close();
		}*/

		ss.str(""); ss<<"Tt"<<i<<j<<".dat";
		ofstr.open(ss.str().c_str());
		DysonEquation.printTt(ofstr,i,j);
		ofstr.close();
	}

	ofstr.open("St.dat");
	DysonEquation.printSt(ofstr);
	ofstr.close();

	ofstr.open("Gt.dat");
	DysonEquation.printGt(ofstr);
	ofstr.close();

	DysonEquation.calculateStS0();
	ofstr.open("StS0.dat");
	DysonEquation.printStS0(ofstr);
	ofstr.close();

	return 0;
};
