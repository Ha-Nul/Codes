#include <fstream>
#include <sstream>
#include "dyson.hpp"
#include "json.hpp"
using namespace std;

int main(int argc, char* argv[]){
	nlohmann::json input;

	ifstream para(argv[1]);
	para >> input;
	para.close();

	int Ntau = input["Ntau"];
	double beta = input["beta"];
	double g_coupling = input["g_coupling"];
	double omega = input["omega"];
	double epsilon = input["epsilon"];

	cout<<"** parameters **"<<endl;
	cout<<"   Ntau = "<<Ntau<<endl;
	cout<<"   beta = "<<beta<<endl;
	cout<<"   g_coupling = "<<g_coupling<<endl;
	cout<<"   omega = "<<omega<<endl;
	cout<<"   epsilon = "<<epsilon<<endl;
	cout<<"****************"<<endl;

	dyson DysonEquation(beta, Ntau);

	DysonEquation.initWt(g_coupling,omega);
	DysonEquation.initGt(epsilon);
	DysonEquation.initQt2ndOrder();

	DysonEquation.initTt();


	ofstream ofstr;

	ofstr.open("Wt.in");
	DysonEquation.printWt(ofstr);
	ofstr.close();

	ofstr.open("Gt.in");
	DysonEquation.printGt(ofstr);
	ofstr.close();

	stringstream ss;
	for(int i=0; i<12; i++){
		ss.str(""); ss<<"Qt_"<<i<<".in";
		ofstr.open(ss.str().c_str());
		DysonEquation.printQt(ofstr,i);
		ofstr.close();
	}


	int Nscl = 1;
	for(int i=0; i<Nscl; i++){
		//DysonEquation.TupdateDebug();
		DysonEquation.Tupdate();
		DysonEquation.updateSt();
		DysonEquation.updateGt(0.0);
		DysonEquation.normalizeGt();
	}

	ofstr.open("Tt.dat");
	DysonEquation.printTt(ofstr);
	ofstr.close();

	ofstr.open("St.dat");
	DysonEquation.printSt(ofstr);
	ofstr.close();

	ofstr.open("Gt.dat");
	DysonEquation.printGt(ofstr);
	ofstr.close();

	DysonEquation.calculateStS0();
	ofstr.open("StS0.exact");
	DysonEquation.printStS0(ofstr);
	ofstr.close();

	/*DysonEquation.initGtexact(epsilon, g_coupling, omega, 100);
	ofstr.open("Gt.exact");
	DysonEquation.printGt(ofstr);
	ofstr.close();*/

	return 0;
};

