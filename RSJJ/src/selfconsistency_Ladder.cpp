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
	int Nscl = input["Nscl"];

	cout<<"** parameters **"<<endl;
	cout<<"   Ntau = "<<Ntau<<endl;
	cout<<"   beta = "<<beta<<endl;
	cout<<"   g_coupling = "<<g_coupling<<endl;
	cout<<"   omega = "<<omega<<endl;
	cout<<"   epsilon = "<<epsilon<<endl;
	cout<<"   Nscl = "<<Nscl<<endl;
	cout<<"****************"<<endl;

	dyson DysonEquation(beta, Ntau);

	DysonEquation.initWt(g_coupling,omega);
	DysonEquation.initGt(epsilon);

	DysonEquation.initTt();


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
	for(int i=0; i<Nscl; i++){
		DysonEquation.Tupdate_Ladder();
		DysonEquation.updateSt();
		DysonEquation.updateGt_volterra(epsilon);
		//DysonEquation.updateGt_dyson(epsilon);

		DysonEquation.normalizeGt();
		DysonEquation.mixingGt();

		DysonEquation.printG05beta(ofstr_iter,i);
	}
	ofstr_iter.close();

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
	ofstr.open("StS0.dat");
	DysonEquation.printStS0(ofstr);
	ofstr.close();

	/*
	for(int i=0; i<Nscl; i++){
		DysonEquation.TupdateBruteForce_Ladder();
		//DysonEquation.updateSt();
		////DysonEquation.updateGt_dyson(epsilon);
		//DysonEquation.updateGt_volterra(epsilon);
		//DysonEquation.normalizeGt();

		//DysonEquation.printG05beta(ofstr_iter,i);
	}

	ofstr.open("Tt.bruteforce");
	DysonEquation.printTt(ofstr);
	ofstr.close();

	//ofstr.open("St.bruteforce");
	//DysonEquation.printSt(ofstr);
	//ofstr.close();

	//ofstr.open("Gt.bruteforce");
	//DysonEquation.printGt(ofstr);
	//ofstr.close();

	//DysonEquation.calculateStS0();
	//ofstr.open("StS0.bruteforce");
	//DysonEquation.printStS0(ofstr);
	//ofstr.close();
	*/

	return 0;
};

