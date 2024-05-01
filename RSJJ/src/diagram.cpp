#include <mpi.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <iterator>
#include "parameters.hpp"
#include "element.hpp"
#include "diagram.hpp"
//#include "green.hpp"
//#include "graph.hpp"
using namespace std;
using namespace parameters;

////////////////////////////
// diagram implementation //
////////////////////////////
diagram::diagram(){};
void diagram::initDiagram(const vector<double>& mWt_i,const vector<Eigen::Matrix2d>& Gt_i){
	this->Nordermax = parameters::Nordermax;
	Norder = 2;
	//cout<<"Nordermax in diagram: "<<(this->Nordermax)<<endl;

	sign = 1;
	dtau = parameters::beta/(Gt_i.size()-1);
	//cout<<endl
	//	<<"dtau = "<<dtau<<endl;
	//cout<<"beta/dtau = "<<beta/dtau<<endl;

	NtW = mWt_i.size();
	NtG = Gt_i.size();

	//mWt.init(mWt_i);
	//Gt.init(Gt_i);

	//cout<<"mWt_i.size(): "<<mWt_i.size()<<endl;
	//cout<<"Gt_i.size(): "<<Gt_i.size()<<endl;
	//mWt.init(mWt_i);
	//Gt.init(Gt_i);

	Conf.init(Gt_i);
	//Conf.init(mWt_i,Gt_i);
};
//void diagram::initWt(const vector<double>& mWt_i){
//	NtW = mWt_i.size();
//	cout<<"NtW: "<<NtW<<endl;
//	mWt.init(mWt_i);
//};
//void diagram::initGt(const vector<Eigen::Matrix2d>& Gt_i){
//	NtG = Gt_i.size();
//	cout<<"NtG: "<<NtG<<endl;
//	Gt.init(Gt_i);
//};

//void diagram::debug(){
//	normalization.resize(this->Nordermax+1);
//	for(int i=0; i<normalization.size(); i++) normalization[i] = 0;
//
//	Qacc.resize(this->Nordermax+1);
//	for(int iorder=0; iorder<Nordermax; iorder++){
//		cout<<"order = "<<iorder<<endl;
//		Qacc[iorder].initZero(NtG);
//		//for(int it=0; it<Nout+1; it++){ Qacc[iorder].print(cout,it); }
//	}
//
//	for(int imc=0; imc<1; imc++){
//		measure();
//	}
//	//Conf.calculate_WProduct();
//	//Conf.calculate_GtTopProduct();
//	//Conf.calculate_GtBottomProduct();
//};
//
////void diagram::select_two_index(const int& N_i, int& i_o, int& j_o){
////	i_o = N_i*unidist(mt);
////	j_o = (N_i-1)*unidist(mt);
////	if(j_o>i_o-1) j_o++;
////	else{
////		int i_o_tmp = i_o;
////		i_o = j_o;
////		j_o = i_o_tmp;
////	}
////};
////
////int diagram::gen_index_removal(const int& r1, const int& r2, const int& o){
////	if(o<std::min(r1,r2)) return o;
////	else if(o<std::max(r1,r2)) return o-1;
////	else return o-2;
////};
///*void diagram::initDiagram(){
//	//int Ndress = 0;
//
//	coordinate xi_n, xj_n;
//	xi_n.setref();
//	xj_n.genrand();
//
//	int iflavor_n, jflavor_n;
//	iflavor_n = 0; jflavor_n = 0;
//
//	vertices.emplace_back(vertex(iflavor_n,xi_n));
//	list<vertex>::iterator vi_iter = vertices.begin();
//	vertices.emplace_back(vertex(jflavor_n,xj_n));
//	list<vertex>::iterator vj_iter = next(vi_iter);
//
//	Ws.emplace_back(line(true,1,iflavor_n,jflavor_n,xi_n,xj_n));
//	list<line>::iterator W_iter = Ws.begin();
//	W_iter->setconnected_vertices(vi_iter,vj_iter);
//
//	list<line>::iterator Gij_iter, Gji_iter;
//	Gs.emplace_back(line(false,0,iflavor_n,jflavor_n,xi_n,xj_n));
//	Gij_iter = Gs.begin();
//	Gs.emplace_back(line(false,0,jflavor_n,iflavor_n,xj_n,xi_n));
//	Gji_iter = next(Gij_iter);
//
//	//miflavor = jflavor_n;
//	//mjflavor = jflavor_n;
//	ref_line = Gji_iter;
//	measuring_line = Gij_iter;
//
//	Gij_iter->setconnected_vertices(vi_iter,vj_iter);
//	Gji_iter->setconnected_vertices(vj_iter,vi_iter);
//
//	vi_iter->setconnected_lines(W_iter,Gji_iter,Gij_iter);
//	vj_iter->setconnected_lines(W_iter,Gij_iter,Gji_iter);
//
//	//fermionicloops.emplace_front(vi_iter);
//	//list<fermionicloop>::iterator loopi_iter = fermionicloops.begin();
//	//fermionicloops.emplace_front(vj_iter);
//	//list<fermionicloop>::iterator loopj_iter = fermionicloops.begin();
//
//	//loopi_iter->calvalue();
//	//loopj_iter->calvalue();
//
//	//vi_iter->setconnected_loop(loopi_iter);
//	//vj_iter->setconnected_loop(loopj_iter);
//
//	//updateSetUWInter();
//	//SetUWInter.calvalue();
//
//	Norder = 1;
//	if( W_iter->getvalueij()>0.0 ) sign = 1;
//	else sign = -1;
//	//Nselfloop = 2;
//
//	this->Nordermax = parameters::Nordermax;
//};*/
///*void diagram::init_P_diagram(){
//	int Ndress = 0;
//
//	coordinate xi_n, xj_n;
//	xi_n.genrand(); xj_n.genrand();
//
//	int iflavor_n, jflavor_n;
//	iflavor_n = static_cast<int>(Nflavor*unidist(mt));
//	jflavor_n = iflavor_n;
//
//	vertices.emplace_front(vertex(iflavor_n,xi_n));
//	list<vertex>::iterator vi_iter = vertices.begin();
//	vertices.emplace_front(vertex(jflavor_n,xj_n));
//	list<vertex>::iterator vj_iter = vertices.begin();
//
//	Ws.emplace_front(line(false,2,Ndress,iflavor_n,jflavor_n,xi_n,xj_n));
//	list<line>::iterator W_iter = Ws.begin();
//	W_iter->setconnected_vertices(vi_iter,vj_iter);
//	list<line>::iterator Gi_iter, Gj_iter;
//
//	Gs.emplace_front(line(true,0,0,iflavor_n,iflavor_n,xi_n,xj_n));
//	Gi_iter = Gs.begin();
//	Gs.emplace_front(line(true,0,0,jflavor_n,jflavor_n,xj_n,xi_n));
//	Gj_iter = Gs.begin();
//
//	miflavor = iflavor_n;
//	mjflavor = jflavor_n;
//	measuringline = W_iter;
//
//	Gi_iter->setconnected_vertices(vi_iter,vj_iter);
//	Gj_iter->setconnected_vertices(vj_iter,vi_iter);
//
//	vi_iter->setconnected_lines(W_iter,Gj_iter,Gi_iter);
//	vj_iter->setconnected_lines(W_iter,Gi_iter,Gj_iter);
//
//	fermionicloops.emplace_front(vi_iter);
//	fermionicloops.begin()->push_back(vj_iter);
//	list<fermionicloop>::iterator loop_iter = fermionicloops.begin();
//	loop_iter->calvalue();
//
//	vi_iter->setconnected_loop(loop_iter);
//	vj_iter->setconnected_loop(loop_iter);
//
//	updateSetUWInter();
//	SetUWInter.calvalue();
//
//	Norder = 0;
//	sign = 1;
//	Nselfloop = 0;
//
//	this->Nordermax = parameters::Nordermax;
//};*/
//
////void diagram::setNordermax(const int& Nordermax_i){ this->Nordermax = Nordermax_i; };
//
//// updates, meaning of return values; 
//// 0: invalid proposal, 1: valid proposal,not accepted, 
//// 2: valid proposal, accepted, measuringline NOT involved
//// 3: valid proposal, accepted, measuringline involved
void diagram::select_two_index(const int& N_i, int& i_o, int& j_o){
	i_o = N_i*unidist(mt);
	j_o = (N_i-1)*unidist(mt);
	if(j_o>i_o-1) j_o++;
	else{
		int i_o_tmp = i_o;
		i_o = j_o;
		j_o = i_o_tmp;
	}
};
void diagram::select_one_excluding_one(const int& N_i, const int& e_i, int& i_o){
	i_o = (N_i-1)*unidist(mt);
	if(i_o>e_i-1) i_o++;
};
void diagram::select_one_excluding_two(const int& N_i, const int& e_i, const int& e_j, int& i_o){
	i_o = N_i*unidist(mt);
	if(i_o>std::min(e_i,e_j)-1) i_o++;
	if(i_o>std::max(e_i,e_j)-1) i_o++;
};

int diagram::insertVertices(){
	if( Conf.getWsSize()==(this->Nordermax)){ return 0; }
	else{
		int Nvs = Conf.getVerticesSize();
		int NWs = Conf.getWsSize();
		int iv, jv;
		select_two_index(Nvs,iv,jv);

		coordinate xi_n, xj_n, xb;

		xb.settime(parameters::beta);

		//printconfiguration(cout);
		//cout<<"Nvs = "<<Nvs<<endl<<flush;
		//cout<<"NWs = "<<NWs<<endl<<flush;
		//cout<<"iv = "<<iv<<endl<<flush;
		//cout<<"jv = "<<jv<<endl<<flush;

		//int iflavor_n, jflavor_n;

		/*list<line>::iterator Gi_iter = next(Gs.begin(),iG);
		list<line>::iterator Gj_iter = next(Gs.begin(),jG);

		list<vertex>::iterator vi_Gi = Gi_iter->getconnected_vi();
		list<vertex>::iterator vj_Gi = Gi_iter->getconnected_vj();

		list<vertex>::iterator vi_Gj = Gj_iter->getconnected_vi();
		list<vertex>::iterator vj_Gj = Gj_iter->getconnected_vj();

		iflavor_n = Gi_iter->getjflavor();
		jflavor_n = Gj_iter->getjflavor();*/


		double ti_interval, tj_interval;
		//xi_n.genrand( Conf.getVertex(iv),Conf.getVertex(iv+1) );
		ti_interval = Conf.getVertex(iv+1).time - Conf.getVertex(iv).time;
		if(jv==Nvs-1){
			//xj_n.genrand(Conf.getVertex(jv),xb);
			tj_interval = xb.time - Conf.getVertex(jv).time;
		}
		else{
			//xj_n.genrand(Conf.getVertex(jv),Conf.getVertex(jv+1));
			tj_interval = Conf.getVertex(jv+1).time - Conf.getVertex(jv).time;
		}

		//double proposal = static_cast<double>(Nvs*(Nvs-1))/2./NWs
		//		* (vj_Gi->getx().gettime() - vi_Gi->getx().gettime()) 
		//	 	* fmod(vj_Gj->getx().gettime() - vi_Gj->getx().gettime() + beta,beta);

		double proposal = static_cast<double>(Nvs*(Nvs-1))/2./NWs
				* ti_interval * tj_interval;

		// for probability check
		//Conf_new.copy( Conf );

		if( iv==Conf.getim2()-1 || jv==Conf.getim2()-1 ){
			proposal *= 2.;
		}

		//Conf.print(cout);
		//Conf.printProposed(cout);

		//cout<<"Maybe here"<<endl<<flush;
		Conf.proposeInsertVertices(iv, jv);
		//cout<<"Isn't it"<<endl<<flush;


		/*
		line Wij(true,1,iflavor_n,jflavor_n,xi_n,xj_n);
		line Gin_m = (*Gi_iter); Gin_m.setj(iflavor_n,xi_n);
		line Gjn_m = (*Gj_iter); Gjn_m.setj(jflavor_n,xj_n);
		line Gin_c(true,0,iflavor_n,Gi_iter->getjflavor(),xi_n,Gi_iter->getxj());
		line Gjn_c(true,0,jflavor_n,Gj_iter->getjflavor(),xj_n,Gj_iter->getxj());

		// measuring line selection
		bool measuring_line_involved = false;
		bool relocate_measuring_line = false;

		if( Gi_iter==measuring_line ){
			measuring_line_involved = true;
			if( Gi_iter!=Gs.begin() ){
				proposal *= 2.;
				if(unidist(mt)<0.5){
					relocate_measuring_line = true;
					Gin_m.setphysical(true);
					Gin_c.setphysical(false);
				}
			}
			else{
				relocate_measuring_line = true;
				Gin_m.setphysical(true);
				Gin_c.setphysical(false);
			}
		}

		bool ref_line_involved = false;
		bool relocate_ref_line = false;
		if( Gj_iter!=ref_line ){
			if( Gj_iter == measuring_line ){
				measuring_line_involved = true;
				proposal *= 2.;
				if(unidist(mt)<0.5){
					relocate_measuring_line = true;
					Gjn_m.setphysical(true);
					Gjn_c.setphysical(false);
				}
			}
		}
		else{
			ref_line_involved = true;
			relocate_ref_line = true;
			Gjn_m.setphysical(true);
			Gjn_c.setphysical(false);
		}

		double probability 
			= proposal
			* Ri_N[Norder] / Ri_N[Norder+1]
			* Wij.getvalueij()
			* Gin_m.getvalueij()*Gjn_m.getvalueij()
			* Gin_c.getvalueij()*Gjn_c.getvalueij()
			/ Gi_iter->getvalueij() / Gj_iter->getvalueij();
		*/

		//Conf.print(cout);
		//Conf.printProposed(cout);

		double probability 
			= proposal
			* Ri_N[Norder] / Ri_N[Norder+1]
			* Conf.getProposedWeight()
			/ Conf.getWeight();

		int signfactor;
		if(probability<0.) signfactor = -1;
		else signfactor = 1;

		if( unidist(mt)<abs(probability) ){
			//cout<<"* insertVertices info *"<<endl;
			//cout<<"  (iG,jG) = ("<<iG<<","<<jG<<")"<<endl;
			//cout<<"  xi_n = "<<xi_n<<endl;
			//cout<<"  xj_n = "<<xj_n<<endl;
			//cout<<"  proposal = "<<proposal<<endl;
			//cout<<"  probability = "<<probability<<endl;
			//cout<<"  (vj_Gi->getx().gettime() - vi_Gi->getx().gettime()) = "
			//	<<(vj_Gi->getx().gettime() - vi_Gi->getx().gettime())<<endl;
			//cout<<"  fmod(vj_Gj->getx().gettime() - vi_Gj->getx().gettime() + beta,beta) = "
			//	<<fmod(vj_Gj->getx().gettime() - vi_Gj->getx().gettime() + beta,beta)<<endl;
			//cout<<"***********************"<<endl;

			Conf.updateInsertVertices();

			//Conf.print(cout);
			//Conf.printProposed(cout);


			/*
			vertices.emplace(next(vi_Gi),vertex(iflavor_n,xi_n));
			vertices.emplace(next(vi_Gj),vertex(jflavor_n,xj_n));
			list<vertex>::iterator vi_iter = next(vi_Gi);
			list<vertex>::iterator vj_iter = next(vi_Gj);

			Ws.push_back(Wij); list<line>::iterator W_iter = prev(Ws.end());
			W_iter->setconnected_vertices(vi_iter,vj_iter);
			//cout<<*W_iter<<endl;

			Gs.insert(next(Gi_iter),Gin_c); list<line>::iterator Gi_next = next(Gi_iter);
			Gs.insert(next(Gj_iter),Gjn_c); list<line>::iterator Gj_next = next(Gj_iter);

			Gi_iter->setj(iflavor_n,xi_n); Gi_iter->setphysical(Gin_m.getphysical());
			Gj_iter->setj(jflavor_n,xj_n); Gj_iter->setphysical(Gjn_m.getphysical());

			list<vertex>::iterator vi_next = Gi_iter->getconnected_vj();
			list<vertex>::iterator vj_next = Gj_iter->getconnected_vj();

			Gi_next->setconnected_vertices(vi_iter,vi_next);
			Gj_next->setconnected_vertices(vj_iter,vj_next);
			vi_next->setconnected_Gin(Gi_next);
			vj_next->setconnected_Gin(Gj_next);

			Gi_iter->setconnected_vj(vi_iter);
			Gj_iter->setconnected_vj(vj_iter);

			vi_iter->setconnected_lines(W_iter,Gi_iter,Gi_next);
			vj_iter->setconnected_lines(W_iter,Gj_iter,Gj_next);

			if( relocate_measuring_line ){
				if( !Gi_next->getphysical() )
					measuring_line = Gi_next;
				else
					measuring_line = Gj_next;
			}
			if( relocate_ref_line ){
				ref_line = prev(Gs.end());
			}
			*/

			Norder += 1;
			sign *= signfactor;

			//cout<<"-->update accepted"<<endl;
			//if(measuring_line_involved || ref_line_involved) return 3;
			//else return 2;
			return 3;
		}
		//cout<<"-->update rejected"<<endl;
		return 1;
	}
};
int diagram::removeVertices(){
	int NWs= Conf.getWsSize();
	if(NWs<3) return 0;

	int Nvs = Conf.getVerticesSize();
	int iW = 1+(NWs-1)*unidist(mt);
	int iv_new = Conf.getiv(iW);
	int jv_new = Conf.getjv(iW);
	int im2 = Conf.getim2();


	//printconfiguration(cout);
	//cout<<"Nvs = "<<Nvs<<endl;
	//cout<<"iW = "<<iW<<endl;
	//cout<<"iv_new = "<<iv_new<<endl;
	//cout<<"jv_new = "<<jv_new<<endl;
	//cout<<"im2 = "<<im2<<endl;
	//cout<<flush;

	/*
	list<line>::iterator W_iter = next(Ws.begin(),1+iW);

	list<vertex>::iterator vi_iter = W_iter->getconnected_vi();
	list<vertex>::iterator vj_iter = W_iter->getconnected_vj();

	list<line>::iterator Giin_iter = vi_iter->getconnected_Gin();
	list<line>::iterator Giout_iter = vi_iter->getconnected_Gout();

	list<line>::iterator Gjin_iter = vj_iter->getconnected_Gin();
	list<line>::iterator Gjout_iter = vj_iter->getconnected_Gout();

	list<vertex>::iterator vi_next = Giout_iter->getconnected_vj();
	list<vertex>::iterator vj_next = Gjout_iter->getconnected_vj();
	*/

	if( (iv_new==0 || iv_new==1) && im2==2) return 0;
	if( (jv_new==Nvs-1 || jv_new==Nvs-2) && im2==Nvs-2) return 0;


	/*
	bool PhyGiin = Giin_iter->getphysical(), PhyGiout = Giout_iter->getphysical();
	bool PhyGjin = Gjin_iter->getphysical(), PhyGjout = Gjout_iter->getphysical();
	bool unphysical_line_involved = !(PhyGiin && PhyGiout && PhyGjin && PhyGjout);

	if( NWs!=2 
		&& 
		(
		 (
		  ( vj_iter==prev(vertices.end(),1) || vj_iter==prev(vertices.end(),2) ) 
		  && !(prev(Gs.end(),3)->getphysical())
		  )
		 || 
		 ( 
		  ( vi_iter==vertices.begin() || vi_iter==next(vertices.begin()) ) 
		  && !(next(Gs.begin())->getphysical())
		  )
		 )
		) {
		return 0;
	}
	*/

	/*if( NWs!=2 
		&& 
		( Giin_iter==Gs.begin() && !PhyGiout ) )
		|| ( Gjout_iter==prev(ref_line) && !PhyGjin ) {
		return 0;
	}*/

	/*
	bool measuring_line_choice = false;
	bool relocate_measuring_line = false;
	bool relocate_ref_line = false;
	if( Giin_iter==measuring_line || Giout_iter==measuring_line ){
		if( Giin_iter!=Gs.begin() ){
			measuring_line_choice = true;
		}
		else{
			measuring_line_choice = false;
		}
		if( Giout_iter==measuring_line ){ relocate_measuring_line = true; }
	}

	if( Gjout_iter!=ref_line ){
		if( Gjin_iter==measuring_line || Gjout_iter==measuring_line ){
			measuring_line_choice = true;
			if( Gjout_iter==measuring_line ){ relocate_measuring_line = true; }
		}
	}
	else{
		relocate_ref_line = true;
	}
	*/

	//cout<<"Maybe here"<<endl<<flush;
	Conf.proposeRemoveVertices(iW);

	//cout<<"* after proposal *"<<endl<<flush;
	//Conf.print(cout);
	//Conf.printProposed(cout);

	//cout<<"Isn't it?"<<endl<<flush;

	//Conf_new.copy( Conf );
	//Conf_new.remove( iW );

	/*
	line Gi_n((Giin_iter->getphysical()&&Giout_iter->getphysical()), 0,
			Giin_iter->getiflavor(),Giout_iter->getjflavor(),
			Giin_iter->getxi(),Giout_iter->getxj());
	line Gj_n((Gjin_iter->getphysical()&&Gjout_iter->getphysical()), 0,
			Gjin_iter->getiflavor(),Gjout_iter->getjflavor(),
			Gjin_iter->getxi(),Gjout_iter->getxj());
	*/

	double proposal = 2.*static_cast<double>(NWs-1)/(Nvs-2)/(Nvs-3)
		/ Conf.get_tinterval(iv_new) / Conf.get_tinterval(jv_new);

	if( iv_new==im2-1 || iv_new==im2 ) proposal /= 2.;
	else if( jv_new==im2-1 || jv_new==im2 ) proposal /= 2.;

	double probability 
		= proposal
		* Ri_N[Norder] / Ri_N[Norder-1]
		* Conf.getProposedWeight()
		/ Conf.getWeight();

	/*double probability 
		= proposal
		* Ri_N[Norder] / Ri_N[Norder-1]
		* Gi_n.getvalueij()*Gj_n.getvalueij()
		/ Giin_iter->getvalueij() / Giout_iter->getvalueij()
		/ Gjin_iter->getvalueij() / Gjout_iter->getvalueij()
		/ W_iter->getvalueij();*/

	//cout<<"* removeVertices info *"<<endl;
	//cout<<"  proposal = "<<proposal<<endl;
	//cout<<"  probability = "<<probability<<endl;
	//cout<<"***********************"<<endl;

	int signfactor;
	if(probability<0.) signfactor = -1;
	else signfactor = 1;

	if( unidist(mt)<abs(probability) ){
		//cout<<endl<<"2PI check for configuration:"<<endl<<flush;
		//Conf.print(cout);
		//Conf.printProposed(cout);
		//cout<<"Wremove: "<<*W_iter<<endl<<flush;
		//cout<<"relocate_measuring_line: "<<relocate_measuring_line<<endl<<flush;
		if( Conf.checkProposedIrreducibility() ){
			//cout<<"which one?"<<endl<<flush;
			if( !(Conf.checkProposedCompactness()) ){ return 0; }
		}
		else return 0;
		/*if(check_irreducibility_for_remove(W_iter)){ 
			if(!(check_compactness_for_remove(W_iter))){ return 0; }
		}
		else return 0;*/
		//cout<<"passed"<<endl<<flush;

		Conf.updateRemoveVertices();

		//Conf.print(cout);
		//Conf.printProposed(cout);
		//Conf = Conf_new;

		/*vi_next->setconnected_Gin(Giin_iter);
		Giin_iter->setphysical(Giin_iter->getphysical()&&Giout_iter->getphysical());
		Giin_iter->setj(Giout_iter->getjflavor(),Giout_iter->getxj());
		Giin_iter->setconnected_vj(Giout_iter->getconnected_vj());

		vj_next = Gjout_iter->getconnected_vj();
		vj_next->setconnected_Gin(Gjin_iter);
		Gjin_iter->setphysical(Gjin_iter->getphysical()&&Gjout_iter->getphysical());
		Gjin_iter->setj(Gjout_iter->getjflavor(),Gjout_iter->getxj());
		Gjin_iter->setconnected_vj(Gjout_iter->getconnected_vj());

		vertices.erase(vi_iter); vertices.erase(vj_iter);
		Ws.erase(W_iter);
		Gs.erase(Giout_iter); Gs.erase(Gjout_iter);

		if( relocate_measuring_line ){
			if( !Giin_iter->getphysical() )
				measuring_line = Giin_iter;
			else
				measuring_line = Gjin_iter;
		}

		if(relocate_ref_line){
			ref_line = Gjin_iter;
		}*/

		Norder -= 1;
		sign *= signfactor;

		//cout<<"-->update accepted"<<endl;
		//if( unphysical_line_involved ) return 3;
		//else return 2;
		return 3;
	}
	//cout<<"-->update rejected"<<endl;
	return 1;
};
/*int diagram::swapWline(){
	int Nvs= Conf.getVerticesSize();
	if(Nvs<5) return 0;

	int iv = Nvs*unidist(mt);
	int iv_partner = Conf.getiv_partner(iv);
	int jv;
	select_one_excluding_two(Nvs, iv, iv_partner, jv);

	Conf.proposeSwapWline(iv,jv);	

	double probability = Conf.getProposedWeight() / Conf.getWeight();

	int signfactor;
	if(probability<0.) signfactor = -1;
	else signfactor = 1;

	if( unidist(mt)<abs(probability) ){
		//cout<<"(accepted)"<<endl;
		if( Conf.checkProposedIrreducibility() ){
			if( !(Conf.checkProposedCompactness()) ){ return 0; }
		}
		else return 0;

		Conf.updateSwapWline();
		sign *= signfactor;

		//G_iter->setxj(x_n);
		//G_next->setxi(x_n);
		//if( W_iter->getconnected_vi()==G_iter->getconnected_vj() )
		//	W_iter->setxi(x_n);
		//else
		//	W_iter->setxj(x_n);

		//v_iter->setx(x_n);
		

		//if( !(G_iter->getphysical() && G_next->getphysical()) ){
		//	return 3;
		//}
		return 3;
	}
	//cout<<"(rejected)"<<endl;
	return 1;
};*/
int diagram::shiftIm2(){
	//cout<<"* shiftIm2 update *"<<endl;
	//Conf_new.copy( Conf );
	//Conf_new.swapMeasuringLine();
	int Nvs = Conf.getVerticesSize();

	if(Nvs<5) return 0;

	int im2 = Conf.getim2();

	int dim2_new;
	select_one_excluding_one(Nvs-3,im2-2,dim2_new);
	int im2_new = 2 + dim2_new;

	//cout<<"im2 : "<<im2<<endl;
	//cout<<"im2_new input: "<< im2_new<<endl;
	
	Conf.proposeShiftIm2(im2_new);

	//Conf.print(cout);
	//Conf.printProposed(cout);

	/*int NWs = Ws.size();
	if(NWs<3) return 0;

	int NGs = Gs.size();
	int iG = distance(Gs.begin(),measuring_line);
	int jG;
	select_one_excluding_one(NGs-3,iG-1,jG);
	jG++;

	line Gi_new = *measuring_line;
	Gi_new.setphysical(true);
	list<line>::iterator Gj_iter = next(Gs.begin(),jG);

	double probability = Gi_new.getvalueij() / Gj_iter->getvalueij();*/

	double probability = Conf.getProposedWeight() / Conf.getWeight();

	int signfactor;
	if(probability<0.) signfactor = -1;
	else signfactor = 1;

	//cout<<"** shiftIm2 info **"<<endl;
	//cout<<"   probability = "<<probability<<endl;
	//cout<<"**********************"<<endl;

	if( unidist(mt)<abs(probability) ){
		//cout<<"checking connectivity"<<endl<<flush;
		if( Conf.checkProposedIrreducibility() ){
			//cout<<"which one"<<endl<<flush;
			if( !(Conf.checkProposedCompactness()) ){ return 0; }
		}
		else return 0;
		//cout<<"valid update"<<endl<<flush;

		Conf.updateShiftIm2();

		sign *= signfactor;

		return 3;
	}
	return 1;
};

int diagram::shiftVertex(){
	//cout<<endl<<"* shiftVertex *"<<endl;
	//Conf_new.copy( Conf );
	//Conf_new.randomShift();

	int Nvs = Conf.getVerticesSize();
	int iv = 1 + (Nvs-1)*unidist(mt);

	//Conf.initProposal();
	Conf.proposeShift(iv);
	//Conf.updateShift();
	//Conf.print(cout);
	//Conf.printProposed(cout);

	double proposedWeight = Conf.getProposedWeight();
	double currentWeight = Conf.getWeight();

	double probability = proposedWeight / currentWeight;

	int signfactor;
	if(probability<0.) signfactor = -1;
	else signfactor = 1;

	if( unidist(mt)<abs(probability) ){
		//cout<<"(accepted)"<<endl;
		Conf.updateShift();
		sign *= signfactor;

		//G_iter->setxj(x_n);
		//G_next->setxi(x_n);
		//if( W_iter->getconnected_vi()==G_iter->getconnected_vj() )
		//	W_iter->setxi(x_n);
		//else
		//	W_iter->setxj(x_n);

		//v_iter->setx(x_n);
		

		//if( !(G_iter->getphysical() && G_next->getphysical()) ){
		//	return 3;
		//}
		return 3;
	}
	//cout<<"(rejected)"<<endl;
	return 1;

	/*int NGs = Gs.size();
	int iG = (NGs-1)*unidist(mt);

	list<line>::iterator G_iter = next(Gs.begin(),iG);
	list<line>::iterator G_next = next(G_iter);
	list<line>::iterator W_iter = G_iter->getconnected_vj()->getconnected_W();

	list<vertex>::iterator v_iter = G_iter->getconnected_vj();

	coordinate x_i = G_iter->getconnected_vi()->getx();
	coordinate x_j = G_next->getconnected_vj()->getx();
	//cout<<"x_i: "<<x_i<<endl
	//	<<"x_j: "<<x_j<<endl;
	//double dtau = x_j.time - x_i.time;
	//if( dtau<1.0e-15 )
	//	dtau += beta;

	coordinate x_n;
	if( G_iter->getconnected_vi()==G_next->getconnected_vj() ){
		x_n.genrand();
	}
	else{
		x_n.genrand(x_i,x_j);
	}

	line Gn = *G_iter; 
	line Gn_next  = *G_next;
	line Wn  = *W_iter;

	Gn.setxj(x_n);
	Gn_next.setxi(x_n);
	if( W_iter->getconnected_vi()==G_iter->getconnected_vj() )
		Wn.setxi(x_n);
	else
		Wn.setxj(x_n);

	double probability = Gn.getvalueij() * Gn_next.getvalueij() * Wn.getvalueij()
			   / G_iter->getvalueij() / G_next->getvalueij() / W_iter->getvalueij();
	*/

	//cout<<"** shiftVertex info **"<<endl;
	//cout<<"   probability = "<<probability<<endl;
	//cout<<"   Gn.getvalueij() = "<<Gn.getvalueij()<<endl;
	//cout<<"   Gn_next.getvalueij() = "<<Gn_next.getvalueij()<<endl;
	//cout<<"   Wn.getvalueij() = "<<Wn.getvalueij()<<endl;
	//cout<<"   G_iter->getvalueij() = "<<G_iter->getvalueij()<<endl;
	//cout<<"   G_next->getvalueij() = "<<G_next->getvalueij()<<endl;
	//cout<<"   W_iter->getvalueij() = "<<W_iter->getvalueij()<<endl;
	//cout<<"**********************"<<endl;

};

int diagram::flipExtFlavor(){
	//cout<<endl<<"* flipExtFlavor *"<<endl;
	//Conf_new.copy( Conf );
	//Conf_new.randomShift();

	int NExtV = 4;
	int iv = NExtV*unidist(mt);

	//cout<<"Before proposal"<<endl;
	//Conf.print(cout);
	//Conf.initProposal();
	Conf.proposeFlip(iv);
	//Conf.print(cout);
	//cout<<"After proposal"<<endl;
	//Conf.printProposed(cout);
	//Conf.updateFlip();

	//double proposedWeight = Conf.getProposedWeight();
	//double currentWeight = Conf.getWeight();

	double probability = Conf.getProposedWeight() / Conf.getWeight();

	int signfactor;
	if(probability<0.) signfactor = -1;
	else signfactor = 1;

	if( unidist(mt)<abs(probability) ){
		//cout<<"(accepted)"<<endl;
		Conf.updateFlip();
		sign *= signfactor;

		/*G_iter->setxj(x_n);
		G_next->setxi(x_n);
		if( W_iter->getconnected_vi()==G_iter->getconnected_vj() )
			W_iter->setxi(x_n);
		else
			W_iter->setxj(x_n);

		v_iter->setx(x_n);*/
		

		//if( !(G_iter->getphysical() && G_next->getphysical()) ){
		//	return 3;
		//}
		return 2;
	}
	//cout<<"(rejected)"<<endl;
	return 1;

	/*int NGs = Gs.size();
	int iG = (NGs-1)*unidist(mt);

	list<line>::iterator G_iter = next(Gs.begin(),iG);
	list<line>::iterator G_next = next(G_iter);
	list<line>::iterator W_iter = G_iter->getconnected_vj()->getconnected_W();

	list<vertex>::iterator v_iter = G_iter->getconnected_vj();

	coordinate x_i = G_iter->getconnected_vi()->getx();
	coordinate x_j = G_next->getconnected_vj()->getx();
	//cout<<"x_i: "<<x_i<<endl
	//	<<"x_j: "<<x_j<<endl;
	//double dtau = x_j.time - x_i.time;
	//if( dtau<1.0e-15 )
	//	dtau += beta;

	coordinate x_n;
	if( G_iter->getconnected_vi()==G_next->getconnected_vj() ){
		x_n.genrand();
	}
	else{
		x_n.genrand(x_i,x_j);
	}

	line Gn = *G_iter; 
	line Gn_next  = *G_next;
	line Wn  = *W_iter;

	Gn.setxj(x_n);
	Gn_next.setxi(x_n);
	if( W_iter->getconnected_vi()==G_iter->getconnected_vj() )
		Wn.setxi(x_n);
	else
		Wn.setxj(x_n);

	double probability = Gn.getvalueij() * Gn_next.getvalueij() * Wn.getvalueij()
			   / G_iter->getvalueij() / G_next->getvalueij() / W_iter->getvalueij();
	*/

	//cout<<"** shiftVertex info **"<<endl;
	//cout<<"   probability = "<<probability<<endl;
	//cout<<"   Gn.getvalueij() = "<<Gn.getvalueij()<<endl;
	//cout<<"   Gn_next.getvalueij() = "<<Gn_next.getvalueij()<<endl;
	//cout<<"   Wn.getvalueij() = "<<Wn.getvalueij()<<endl;
	//cout<<"   G_iter->getvalueij() = "<<G_iter->getvalueij()<<endl;
	//cout<<"   G_next->getvalueij() = "<<G_next->getvalueij()<<endl;
	//cout<<"   W_iter->getvalueij() = "<<W_iter->getvalueij()<<endl;
	//cout<<"**********************"<<endl;

};

////int diagram::shuffleExtFlavor(){
////	Conf_new.copy( Conf );
////	Conf_new.swapMeasuringLine();
////
////	double probability = Conf_new.getWeight() / Conf.getWeight()
////
////	int signfactor;
////	if(probability<0.) signfactor = -1;
////	else signfactor = 1;
////
////	if( unidist(mt)<abs(probability) ){
////		Conf = Conf_new;
////
////		return 3;
////	}
////	return 1;
////}
//
///*int diagram::reconnect(){
//	int Nvs = vertices.size();
//	//if(Nvs<3) return 0;
//
//	int iv, jv;
//	select_two_index(Nvs, iv, jv);
//
//	list<vertex>::iterator vi_iter = next(vertices.begin(),iv);
//	list<vertex>::iterator vj_iter = next(vertices.begin(),jv);
//
//	//cout<<"selected vertices"<<endl;
//	//cout<<*vi_iter<<*vj_iter<<endl;
//
//	list<line>::iterator Giin_iter = vi_iter->getconnected_Gin();
//	list<line>::iterator Gjin_iter = vj_iter->getconnected_Gin();
//
//	line Giin_new(*Gjin_iter);
//	line Gjin_new(*Giin_iter);
//	Giin_new.setj(Giin_iter->getjflavor(),Giin_iter->getxj());
//	Gjin_new.setj(Gjin_iter->getjflavor(),Gjin_iter->getxj());
//	Giin_new.setconnected_vertices(Gjin_iter->getconnected_vi(),vi_iter);
//	Gjin_new.setconnected_vertices(Giin_iter->getconnected_vi(),vj_iter);
//
//	//line Giin_new(*Giin_iter);
//	//line Gjin_new(*Gjin_iter);
//	//Giin_new.setj(Gjin_iter->getjflavor(),Gjin_iter->getxj());
//	//Gjin_new.setj(Giin_iter->getjflavor(),Giin_iter->getxj());
//	//Giin_new.setconnected_vertices(Giin_iter->getconnected_vi(),vj_iter);
//	//Gjin_new.setconnected_vertices(Gjin_iter->getconnected_vi(),vi_iter);
//
//
//	int dNselfloop = 0;
//	if(Giin_iter==vi_iter->getconnected_Gout()) dNselfloop--;
//	if(Gjin_iter==vj_iter->getconnected_Gout()) dNselfloop--;
//	if(Giin_new.getconnected_vi()==Giin_new.getconnected_vj()) dNselfloop++;
//	if(Gjin_new.getconnected_vi()==Gjin_new.getconnected_vj()) dNselfloop++;
//	if(Ws.size()!=1 && Nselfloop+dNselfloop>NselfloopMax) return 0;
//	//if((Us.size()!=1 || Ws.size()!=0) && Nselfloop+dNselfloop>NselfloopMax) return 0;
//
//	// probability check
//	double probability 
//		= -1.; // Fermionic sign factor
//
//		// * Giin_new.getvalue() * Gjin_new.getvalue() 
//		/// Giin_iter->getvalue() / Gjin_iter->getvalue();
//
//	//cout<<"loop update proposal"<<endl<<flush;
//
//	list<fermionicloop>::iterator loopi_iter = vi_iter->getconnected_loop();
//	list<fermionicloop>::iterator loopj_iter = vj_iter->getconnected_loop();
//
//	bool sameloop = (loopi_iter==loopj_iter);
//
//	setUWinter trial;
//	list<line>::iterator UW_iter;
//	list<fermionicloop>::iterator loop1_iter;
//	list<fermionicloop>::iterator loop2_iter;
//
//	list<list<vertex>::iterator> vset_split1, vset_split2, vset_combined;
//
//	double loop1_valuerh_n, loop1_valuelh_n, loop2_valuerh_n, loop2_valuelh_n;
//	double loop_valuerh_n, loop_valuelh_n;
//
//	if( sameloop ){
//		//cout<<"loop split case"<<endl<<flush;
//		vset_split1 = loopi_iter->getpath(vi_iter,vj_iter);
//		vset_split2 = loopi_iter->getpath(vj_iter,vi_iter);
//
//		int iloop = distance(fermionicloops.begin(),loopi_iter);
//
//		fermionicloop loop1_pseudo(vset_split1);
//		fermionicloop loop2_pseudo(vset_split2);
//
//		loop1_pseudo.calvalue();
//		loop2_pseudo.calvalue();
//
//		loop1_valuerh_n = loop1_pseudo.getvaluerh() * Giin_new.getvalueij() / Gjin_iter->getvalueij();
//		loop1_valuelh_n = loop1_pseudo.getvaluelh() * Giin_new.getvalueji() / Gjin_iter->getvalueji();
//		loop2_valuerh_n = loop2_pseudo.getvaluerh() * Gjin_new.getvalueij() / Giin_iter->getvalueij();
//		loop2_valuelh_n = loop2_pseudo.getvaluelh() * Gjin_new.getvalueji() / Giin_iter->getvalueji();
//
//		trial.setNfermionicloops(fermionicloops.size()+1);
//		bool v1_in, v2_in;
//		list<list<vertex>::iterator>::iterator v_iiter;
//		list<vertex>::iterator v1_iter;
//		list<vertex>::iterator v2_iter;
//		for(UW_iter=Ws.begin(); UW_iter!=Ws.end(); ++UW_iter){
//			v1_iter = UW_iter->getconnected_vi();
//			v2_iter = UW_iter->getconnected_vj();
//			loop1_iter = v1_iter->getconnected_loop();
//			loop2_iter = v2_iter->getconnected_loop();
//
//			if( (loop1_iter==loopi_iter) && (loop2_iter==loopi_iter) ){
//				v1_in = false; v2_in = false;
//				for(v_iiter=vset_split1.begin(); v_iiter!=vset_split1.end(); ++v_iiter){
//					if(v1_iter==*v_iiter) v1_in = true;
//					else if(v2_iter==*v_iiter) v2_in = true;
//				}
//				if(v1_in^v2_in){
//					probability /= UW_iter->getvalueij();
//					if(v1_in) trial.push_back(*UW_iter,iloop,iloop+1);
//					else trial.push_back(*UW_iter,iloop+1,iloop);
//				}
//			}
//			if(loop1_iter!=loop2_iter){
//				if(loop1_iter==loopi_iter){
//					v1_in = false;
//					for(v_iiter=vset_split1.begin(); v_iiter!=vset_split1.end(); ++v_iiter){
//						if(v1_iter==*v_iiter){
//							v1_in = true;
//							break;
//						}
//					}
//					if(v1_in){
//						trial.push_back(
//							*UW_iter,
//							iloop,
//							gen_index_split(iloop,distance(fermionicloops.begin(),loop2_iter))
//						);
//					}
//					else{
//						trial.push_back(
//							*UW_iter,
//							iloop+1,
//							gen_index_split(iloop,distance(fermionicloops.begin(),loop2_iter))
//						);
//					}
//				}
//				else if(loop2_iter==loopi_iter){
//					v1_in = false;
//					for(v_iiter=vset_split1.begin(); v_iiter!=vset_split1.end(); ++v_iiter){
//						if(v2_iter==*v_iiter){
//							v1_in = true;
//							break;
//						}
//					}
//					if(v1_in){
//						trial.push_back(
//							*UW_iter,
//							gen_index_split(iloop,distance(fermionicloops.begin(),loop1_iter)),
//							iloop
//						);
//					}
//					else{
//						trial.push_back(
//							*UW_iter,
//							gen_index_split(iloop,distance(fermionicloops.begin(),loop1_iter)),
//							iloop+1
//						);
//					}
//				}
//				else{
//					trial.push_back(
//						*UW_iter,
//						gen_index_split(iloop,distance(fermionicloops.begin(),loop1_iter)),
//						gen_index_split(iloop,distance(fermionicloops.begin(),loop2_iter))
//					);
//				}
//			}
//		}
//
//		v1_iter = measuringline->getconnected_vi();
//		v2_iter = measuringline->getconnected_vj();
//		loop1_iter = v1_iter->getconnected_loop();
//		loop2_iter = v2_iter->getconnected_loop();
//		bool vi_1in = false; 
//		bool vj_1in = false;
//		for(v_iiter=vset_split1.begin(); v_iiter!=vset_split1.end(); ++v_iiter){
//			if(v1_iter==*v_iiter){ vi_1in = true; }
//			if(v2_iter==*v_iiter){ vj_1in = true; }
//		}
//		bool vi_2in = false; 
//		bool vj_2in = false;
//		for(v_iiter=vset_split2.begin(); v_iiter!=vset_split2.end(); ++v_iiter){
//			if(v1_iter==*v_iiter){ vi_2in = true; }
//			if(v2_iter==*v_iiter){ vj_2in = true; }
//		}
//
//		if( measuringline->gettype()==0 ){
//			if(vi_1in){
//				trial.fixloopflavor(iloop, miflavor);
//				probability *= loop1_valuerh_n * 0.5*(loop2_valuerh_n+loop2_valuelh_n)
//						/ loopi_iter->getvaluerh();
//				//probability *= loop1_valuerh_n * (loop2_valuerh_n+loop2_valuelh_n)
//				//		/ loopi_iter->getvaluerh();
//			}
//			else if(vi_2in){
//				trial.fixloopflavor(iloop+1, miflavor);
//				probability *= 0.5*(loop1_valuerh_n+loop1_valuelh_n) * loop2_valuerh_n
//						/ loopi_iter->getvaluerh();
//				//probability *= (loop1_valuerh_n+loop1_valuelh_n) * loop2_valuerh_n
//				//		/ loopi_iter->getvaluerh();
//			}
//			else{
//				trial.fixloopflavor(
//					gen_index_split(iloop,distance(fermionicloops.begin(),loop1_iter)),
//					miflavor
//				);
//				probability *= 0.5
//						* (loop1_valuerh_n+loop1_valuelh_n) * (loop2_valuerh_n+loop2_valuelh_n)
//						/ (loopi_iter->getvaluerh() + loopi_iter->getvaluelh());
//				//probability *= (loop1_valuerh_n+loop1_valuelh_n) * (loop2_valuerh_n+loop2_valuelh_n)
//				//		/ (loopi_iter->getvaluerh() + loopi_iter->getvaluelh());
//			}
//		}
//		else{
//			if(vi_1in) trial.fixloopflavor(iloop, miflavor);
//			else if(vi_2in) trial.fixloopflavor(iloop+1, miflavor);
//			else{
//				trial.fixloopflavor(
//					gen_index_split(iloop,distance(fermionicloops.begin(),loop1_iter)),
//					miflavor
//				);
//			}
//			if(vj_1in) trial.fixloopflavor(iloop, mjflavor);
//			else if(vj_2in) trial.fixloopflavor(iloop+1, mjflavor);
//			else{
//				trial.fixloopflavor(
//					gen_index_split(iloop,distance(fermionicloops.begin(),loop2_iter)),
//					mjflavor
//				);
//			}
//			probability *= 0.5
//					* (loop1_valuerh_n+loop1_valuelh_n) * (loop2_valuerh_n+loop2_valuelh_n)
//					/ (loopi_iter->getvaluerh() + loopi_iter->getvaluelh());
//			//probability *= (loop1_valuerh_n+loop1_valuelh_n) * (loop2_valuerh_n+loop2_valuelh_n)
//			//		/ (loopi_iter->getvaluerh() + loopi_iter->getvaluelh());
//		}
//
//		trial.calvalue();
//		probability *= trial.getvalue() / SetUWInter.getvalue();
//
//		//printconfiguration(cout);
//		//cout<<"trial.calvalue()"<<endl<<flush;
//		//cout<<trial<<endl;
//		//cout<<"symmetrized UWinter value: "<<trial.getvalue()<<endl;
//	}
//	else{
//		//cout<<"loop merge case"<<endl<<flush;
//		int iloop = distance(fermionicloops.begin(),loopi_iter);
//		int jloop = distance(fermionicloops.begin(),loopj_iter);
//
//		if( 
//			(SetUWInter.getisloopflavorfixed(iloop) && SetUWInter.getisloopflavorfixed(jloop))
//			&& 
//			(SetUWInter.getloopflavor(iloop)!=SetUWInter.getloopflavor(jloop))
//		){ return 0; }
//
//
//		loop_valuerh_n = loopi_iter->getvaluerh()*loopj_iter->getvaluerh()
//					* Giin_new.getvalueij() * Gjin_new.getvalueij() 
//					/ Giin_iter->getvalueij() / Gjin_iter->getvalueij();
//		loop_valuelh_n = loopi_iter->getvaluelh()*loopj_iter->getvaluelh()
//					* Giin_new.getvalueji() * Gjin_new.getvalueji() 
//					/ Giin_iter->getvalueji() / Gjin_iter->getvalueji();
//
//		//cout<<"loopi_iter->getvaluerh() = "<<loopi_iter->getvaluerh()<<endl;
//		//cout<<"loopi_iter->getvaluelh() = "<<loopi_iter->getvaluelh()<<endl;
//		//cout<<"Giin_new.getvalueij()    = "<<Giin_new.getvalueij()  <<endl;
//		//cout<<"Gjin_new.getvalueij()    = "<<Gjin_new.getvalueij()  <<endl;
//		//cout<<"Giin_iter->getvalueij()  = "<<Giin_iter->getvalueij()<<endl; 
//		//cout<<"Gjin_iter->getvalueij()  = "<<Gjin_iter->getvalueij()<<endl;
//		//cout<<"construct trial"<<endl;
//		trial.setNfermionicloops(fermionicloops.size()-1);
//		int i1loop, i2loop;
//		list<line>::iterator UW_iter;
//		for(UW_iter=Ws.begin(); UW_iter!=Ws.end(); ++UW_iter){
//			loop1_iter = UW_iter->getconnected_vi()->getconnected_loop();
//			loop2_iter = UW_iter->getconnected_vj()->getconnected_loop();
//
//			i1loop = gen_index_merge(iloop,jloop,distance(fermionicloops.begin(),loop1_iter));
//			i2loop = gen_index_merge(iloop,jloop,distance(fermionicloops.begin(),loop2_iter));
//
//			//cout<<"("<<i1loop<<";"<<distance(fermionicloops.begin(),loop1_iter)<<")"
//			//	<<"("<<i2loop<<";"<<distance(fermionicloops.begin(),loop2_iter)<<")"<<endl;
//
//			if(loop1_iter!=loop2_iter){
//				if(i1loop==i2loop){ probability *= UW_iter->getvalueij(); }
//				else{ trial.push_back(*UW_iter,i1loop,i2loop); }
//			}
//		}
//		//cout<<"mark measuringline"<<endl;
//		loop1_iter = measuringline->getconnected_vi()->getconnected_loop();
//		loop2_iter = measuringline->getconnected_vj()->getconnected_loop();
//		i1loop = gen_index_merge(iloop,jloop,distance(fermionicloops.begin(),loop1_iter));
//		i2loop = gen_index_merge(iloop,jloop,distance(fermionicloops.begin(),loop2_iter));
//		trial.fixloopflavor( i1loop, miflavor );
//		trial.fixloopflavor( i2loop, mjflavor );
//
//		if(measuringline->gettype()==0){
//			if(loop1_iter==loopi_iter){
//				probability *= loop_valuerh_n 
//						/ loopi_iter->getvaluerh() 
//						/ 0.5/(loopj_iter->getvaluerh()+loopj_iter->getvaluelh());
//				//probability *= loop_valuerh_n 
//				//		/ loopi_iter->getvaluerh() 
//				//		/ (loopj_iter->getvaluerh()+loopj_iter->getvaluelh());
//			}
//			else if(loop1_iter==loopj_iter){
//				probability *= loop_valuerh_n 
//						/ 0.5/(loopi_iter->getvaluerh()+loopi_iter->getvaluelh())
//						/ loopj_iter->getvaluerh();
//				//probability *= loop_valuerh_n 
//				//		/ (loopi_iter->getvaluerh()+loopi_iter->getvaluelh())
//				//		/ loopj_iter->getvaluerh();
//			}
//			else{
//				probability *= (loop_valuerh_n+loop_valuelh_n) 
//						/ 0.5
//						/ (loopi_iter->getvaluerh()+loopi_iter->getvaluelh())
//						/ (loopj_iter->getvaluerh()+loopj_iter->getvaluelh());
//				//probability *= (loop_valuerh_n+loop_valuelh_n) 
//				//		/ (loopi_iter->getvaluerh()+loopi_iter->getvaluelh())
//				//		/ (loopj_iter->getvaluerh()+loopj_iter->getvaluelh());
//			}
//		}
//		else{
//			probability *= (loop_valuerh_n+loop_valuelh_n) 
//					/ 0.5
//					/ (loopi_iter->getvaluerh()+loopi_iter->getvaluelh())
//					/ (loopj_iter->getvaluerh()+loopj_iter->getvaluelh());
//			//probability *= (loop_valuerh_n+loop_valuelh_n) 
//			//		/ (loopi_iter->getvaluerh()+loopi_iter->getvaluelh())
//			//		/ (loopj_iter->getvaluerh()+loopj_iter->getvaluelh());
//		}
//		//cout<<"calculate trial value"<<endl;
//		//cout<<"symmetrized UWinter value: "<<trial.getvalue()<<endl;
//		trial.calvalue();
//
//		probability *= trial.getvalue() / SetUWInter.getvalue();
//
//	}
//
//	//printconfiguration(cout);
//	//cout<<"* reconnect info *"<<endl;
//	//if(sameloop){
//	//	cout<<"  loop split case"<<endl;
//	//	cout<<"    loop1_valuerh_n = "<<loop1_valuerh_n<<endl;
//	//	cout<<"    loop1_valuelh_n = "<<loop1_valuelh_n<<endl;
//	//	cout<<"    loop2_valuerh_n = "<<loop2_valuerh_n<<endl;
//	//	cout<<"    loop2_valuelh_n = "<<loop2_valuelh_n<<endl;
//	//}
//	//else{
//	//	cout<<"  loop merge case"<<endl;
//	//	cout<<"    loop_valuerh_n = "<<loop_valuerh_n<<endl;
//	//	cout<<"    loop_valuelh_n = "<<loop_valuelh_n<<endl;
//	//}
//	//cout<<"  probability = "<<probability<<endl;
//	//cout<<"  trial info. "<<endl<<trial<<endl;
//	//cout<<"  SetUWInter.getvalue() = "<<SetUWInter.getvalue()<<endl;
//	//cout<<"*******************"<<endl;
//
//	int signfactor;
//	if(probability<0.) signfactor = -1;
//	else signfactor = 1;
//
//	if(unidist(mt)<abs(probability)){
//		int vi, vj;
//		GadjLtmp.clear();
//		GadjLtmp.resize(vertices.size());
//		for(list<line>::iterator giter=Gs.begin(); giter!=Gs.end(); ++giter){
//			if(giter->getphysical()){
//				if(giter==Giin_iter){
//					vi = distance(vertices.begin(),giter->getconnected_vi());
//					vj = distance(vertices.begin(),Gjin_iter->getconnected_vj());;
//				}
//				else if(giter==Gjin_iter){
//					vi = distance(vertices.begin(),giter->getconnected_vi());
//					vj = distance(vertices.begin(),Giin_iter->getconnected_vj());;
//				}
//				else{
//					vi = distance(vertices.begin(),giter->getconnected_vi());
//					vj = distance(vertices.begin(),giter->getconnected_vj());
//				}
//				if(vi!=vj){
//					GadjLtmp[vi].push_back(vj);
//					GadjLtmp[vj].push_back(vi);
//				}
//			}
//		}
//
//		WmadjLtmp.clear();
//		WadjLtmp.clear();
//		WadjLtmp.resize(vertices.size());
//		for(list<line>::iterator witer=Ws.begin(); witer!=Ws.end(); ++witer){
//			vi = distance(vertices.begin(),witer->getconnected_vi());
//			vj = distance(vertices.begin(),witer->getconnected_vj());
//			if(vi!=vj){
//				if(witer->getphysical()){
//					WadjLtmp[vi].push_back(vj);
//					WadjLtmp[vj].push_back(vi);
//				}
//				else{
//					WmadjLtmp.push_back(vi);
//					WmadjLtmp.push_back(vj);
//				}
//			}
//		}
//
//		//cout<<"* GadjLtmp *"<<endl;
//		//for(auto i: GadjLtmp){
//		//	for(auto j: i) cout<<j;
//		//	cout<<endl;
//		//}
//		//cout<<"* WadjLtmp *"<<endl;
//		//for(auto i: WadjLtmp){
//		//	for(auto j: i) cout<<j;
//		//	cout<<endl;
//		//}
//
//		//cout<<"compactness check"<<endl;
//		//if(checkconnectivity(GadjLtmp,WadjLtmp)!=checkconnectivityforreconnect(Giin_iter,Gjin_iter)){
//		//	cout<<"reconnect update; connectivity check problem"<<endl;
//		//	exit(EXIT_FAILURE);
//		//}
//		//if(checkirreducibility(GadjLtmp,WadjLtmp)!=checkirreducibilityforreconnect(Giin_iter,Gjin_iter)){
//		//	cout<<"reconnect update; irreducibility check problem"<<endl;
//		//	printconfiguration(cout);
//		//	cout<<"is connected after update? "<<checkconnectivity(GadjLtmp,WadjLtmp)<<endl;
//		//	exit(EXIT_FAILURE);
//		//}
//		//if( (checkGcompactness(GadjLtmp,WadjLtmp)&&checkWcompactness(GadjLtmp,WmadjLtmp,WadjLtmp))!=checkcompactnessforreconnect(Giin_iter,Gjin_iter) ){
//		//	cout<<"reconnect update; compactness check problem"<<endl;
//		//	exit(EXIT_FAILURE);
//		//}
//
//		if(!checkconnectivity(GadjLtmp,WadjLtmp)) return 0;
//		if(!checkirreducibility(GadjLtmp,WadjLtmp)) return 0;
//		//if(!checkGcompactness(GadjLtmp,WmadjLtmp,WadjLtmp)) return 0;
//		//if(!checkWcompactness(GadjLtmp,WmadjLtmp,WadjLtmp)) return 0;
//		if(checkbubble(GadjLtmp)) return 0;
//
//		//cout<<"is irreducible: "<<checkirreducibility(GadjLtmp,WadjLtmp);
//		//cout<<"is compact in G: "<<checkGcompactness(GadjLtmp,WadjLtmp);
//		//cout<<"is compact in W: "<<checkWcompactness(GadjLtmp,WadjLtmp);
//		//if(!(checkconnectivityforreconnect(Giin_iter,Gjin_iter))) return 0;
//		//if(!(checkirreducibilityforreconnect(Giin_iter,Gjin_iter))) return 0;
//		//if(!(checkcompactnessforreconnect(Giin_iter,Gjin_iter))) return 0;
//
//		//cout<<"start to update quantities"<<endl;
//		bool measuring_line_involved = false;
//		if( Giin_iter->getphysical()^Gjin_iter->getphysical() )
//			measuring_line_involved = true;
//
//		vi_iter->setconnected_Gin(Gjin_iter);
//		vj_iter->setconnected_Gin(Giin_iter);
//		Giin_iter->setxj(vj_iter->getx());
//		Gjin_iter->setxj(vi_iter->getx());
//		Giin_iter->setconnected_vj(vj_iter);
//		Gjin_iter->setconnected_vj(vi_iter);
//
//		//cout<<"loop quantities update"<<endl;
//		list<list<vertex>::iterator>::iterator v_iiter;
//		if( sameloop ){
//			//cout<<"loop split case"<<endl;
//			loopi_iter->set(vset_split1);
//			fermionicloops.insert(next(loopi_iter),fermionicloop(vset_split2));
//			loopj_iter = next(loopi_iter);
//
//			loopi_iter->setvalue(loop1_valuerh_n,loop1_valuelh_n);
//			loopj_iter->setvalue(loop2_valuerh_n,loop2_valuelh_n);
//
//			for(v_iiter=vset_split1.begin(); v_iiter!=vset_split1.end(); ++v_iiter)
//				(*v_iiter)->setconnected_loop(loopi_iter);
//			for(v_iiter=vset_split2.begin(); v_iiter!=vset_split2.end(); ++v_iiter)
//				(*v_iiter)->setconnected_loop(loopj_iter);
//		}
//		else{
//			//cout<<"loop merge case"<<endl;
//			//cout<<*loopi_iter<<*loopj_iter<<endl;
//			vset_combined = loopi_iter->getpath(vi_iter,vi_iter);
//			vset_combined.splice(vset_combined.end(),loopj_iter->getpath(vj_iter,vj_iter));
//			int iloop = distance(fermionicloops.begin(),loopi_iter);
//			int jloop = distance(fermionicloops.begin(),loopj_iter);
//			//cout<<"vset_combined"<<endl;
//			//for(v_iiter=vset_combined.begin(); v_iiter!=vset_combined.end(); ++v_iiter)
//			//	cout<<(**v_iiter);
//			//cout<<endl;
//			if(iloop<jloop){
//				loopi_iter->set(vset_combined);
//				fermionicloops.erase(loopj_iter);
//
//				loopi_iter->setvalue(loop_valuerh_n,loop_valuelh_n);
//
//				for(v_iiter=vset_combined.begin(); v_iiter!=vset_combined.end(); ++v_iiter)
//					(*v_iiter)->setconnected_loop(loopi_iter);
//			}
//			else{
//				loopj_iter->set(vset_combined);
//				fermionicloops.erase(loopi_iter);
//
//				loopj_iter->setvalue(loop_valuerh_n,loop_valuelh_n);
//
//				for(v_iiter=vset_combined.begin(); v_iiter!=vset_combined.end(); ++v_iiter)
//					(*v_iiter)->setconnected_loop(loopj_iter);
//			}
//		}
//		SetUWInter = trial;
//		//updateSetUWInter();
//		//SetUWInter.setvalue(trial.getvalue());
//
//		sign *= signfactor;
//		Nselfloop += dNselfloop;
//
//		//printconfiguration(cout);
//		//exit(EXIT_FAILURE);
//
//		//cout<<"-->update accepted"<<endl;
//		if(measuring_line_involved) return 3;
//		else return 2;
//	}
//	//cout<<"-->update rejected"<<endl;
//	return 1;
//};*/
///*int diagram::intraSwap(){
//	//cout<<"* intraswap module *"<<endl<<flush;
//	int index_i, index_j;
//	list<line>::iterator line_j_iter;
//	bool ismeasuringlineGline = (measuringline->gettype()==0);
//	if(ismeasuringlineGline){
//		index_i = distance(Gs.begin(),measuringline);
//		select_one_excluding_one(Gs.size(),index_i,index_j);
//		line_j_iter = next(Gs.begin(),index_j);
//	}
//	else{
//		index_i = distance(Ws.begin(),measuringline);
//		select_one_excluding_one(Ws.size(),index_i,index_j);
//		line_j_iter = next(Ws.begin(),index_j);
//	}
//
//	line line_i_new(*measuringline);
//	line_i_new.setphysical(true);
//
//	double proposal = 1.0;
//	double probability = 1.0;
//	
//	int miflavor_n, mjflavor_n;
//	int iloop, jloop;
//	list<fermionicloop>::iterator loopi_iter, loopj_iter;
//
//	double loopi_valuerh_n, loopi_valuelh_n;
//	double loopj_valuerh_n, loopj_valuelh_n;
//
//	vector<bool> isloopflavorfixed_n(fermionicloops.size(),false);
//	vector<int> loopflavor_n(fermionicloops.size(),0);
//	double valtrial;
//
//	bool islinei_inter, islinej_inter;
//	if(ismeasuringlineGline){
//		loopi_iter = measuringline->getconnected_vi()->getconnected_loop();
//		loopj_iter = line_j_iter->getconnected_vi()->getconnected_loop();
//		if(loopi_iter==loopj_iter){
//			loopi_valuerh_n = loopi_iter->getvaluerh() 
//					* line_i_new.getvalueij() 
//					/ line_j_iter->getvalueij();
//			loopi_valuelh_n = loopi_iter->getvaluelh() 
//					* line_i_new.getvalueji() 
//					/ line_j_iter->getvalueji();
//
//			miflavor_n = Nflavor*unidist(mt);
//			mjflavor_n = miflavor_n;
//
//			jloop = distance(fermionicloops.begin(),loopj_iter);
//			isloopflavorfixed_n[jloop] = true;
//			loopflavor_n[jloop] = miflavor_n;
//			valtrial = SetUWInter.getvalue();
//
//			probability *= loopi_valuerh_n / loopi_iter->getvaluerh();
//		}
//		else{
//			loopi_valuerh_n = loopi_iter->getvaluerh() 
//					* line_i_new.getvalueij();
//			loopi_valuelh_n = loopi_iter->getvaluelh() 
//					* line_i_new.getvalueji();
//			loopj_valuerh_n = loopj_iter->getvaluerh() 
//					/ line_j_iter->getvalueij();
//			loopj_valuelh_n = loopj_iter->getvaluelh() 
//					/ line_j_iter->getvalueji();
//
//			miflavor_n = Nflavor*unidist(mt);
//			mjflavor_n = miflavor_n;
//
//			jloop = distance(fermionicloops.begin(),loopj_iter);
//			isloopflavorfixed_n[jloop] = true;
//			loopflavor_n[jloop] = miflavor_n;
//			valtrial = SetUWInter.calvalueafterinterswap(isloopflavorfixed_n,loopflavor_n);
//
//			probability *= 0.5*(loopi_valuerh_n + loopi_valuelh_n) / loopi_iter->getvaluerh() 
//					* loopj_valuerh_n / 0.5/(loopj_iter->getvaluerh()+loopj_iter->getvaluelh())
//					* valtrial / SetUWInter.getvalue();
//		}
//	}
//	else{
//		islinei_inter = (measuringline->getconnected_vi()->getconnected_loop()
//					!= measuringline->getconnected_vj()->getconnected_loop());
//		islinej_inter = (line_j_iter->getconnected_vi()->getconnected_loop()
//					!= line_j_iter->getconnected_vj()->getconnected_loop());
//
//		if(!islinej_inter){
//			miflavor_n = Nflavor*unidist(mt);
//			mjflavor_n = miflavor_n;
//
//			if(!islinei_inter){
//				probability *= line_i_new.getvalueij() / line_j_iter->getvalueij();
//				valtrial = SetUWInter.getvalue();
//			}
//			else{
//				proposal /= 2.;
//
//				loopi_iter = line_j_iter->getconnected_vi()->getconnected_loop();
//
//				iloop = distance(fermionicloops.begin(),loopi_iter);
//				isloopflavorfixed_n[iloop] = true;
//				loopflavor_n[iloop] = miflavor_n;
//				valtrial = SetUWInter.calvalueafterinterswap(*measuringline,isloopflavorfixed_n,loopflavor_n);
//
//				probability *= 1. / line_j_iter->getvalueij()
//						* valtrial / SetUWInter.getvalue();
//			}
//		}
//		else{
//			miflavor_n = Nflavor*unidist(mt);
//			mjflavor_n = Nflavor*unidist(mt);
//
//			loopi_iter = line_j_iter->getconnected_vi()->getconnected_loop();
//			loopj_iter = line_j_iter->getconnected_vj()->getconnected_loop();
//
//			iloop = distance(fermionicloops.begin(),loopi_iter);
//			jloop = distance(fermionicloops.begin(),loopj_iter);
//			isloopflavorfixed_n[iloop] = true;
//			isloopflavorfixed_n[jloop] = true;
//			loopflavor_n[iloop] = miflavor_n;
//			loopflavor_n[jloop] = mjflavor_n;
//
//			if(!islinei_inter){
//				proposal *= 2.;
//				valtrial = SetUWInter.calvalueafterinterswap(*line_j_iter,isloopflavorfixed_n,loopflavor_n);
//			}
//			else{
//				valtrial = SetUWInter.calvalueafterintraswap(*measuringline,*line_j_iter,isloopflavorfixed_n,loopflavor_n);
//			}
//		}
//	}
//
//	probability *= proposal;
//
//	int signfactor;
//	if(probability<0.) signfactor = -1;
//	else signfactor = 1;
//
//	//if(NdressChoice==3 && measuringline->getisG() && Ndress==2){
//	//printconfiguration(cout);
//	//cout<<"* interswap info *"<<endl;
//	//cout<<"measuringline = "<<*measuringline<<endl;
//	//cout<<"line_i_new = "<<line_i_new<<endl;
//	//cout<<"line_j_iter = "<<*line_j_iter<<endl;
//	//cout<<"Norder = "<<Norder<<endl;
//	//cout<<"dNorder = "<<dNorder<<endl<<flush;
//	//cout<<"probability = "<<probability<<endl;
//	//cout<<"signfactor = "<<signfactor<<endl;
//	//cout<<"line_i_new.getvalue() = "<<line_i_new.getvalue()<<endl;
//	//cout<<"line_j_iter->getvalue() = "<<line_j_iter->getvalue()<<endl;
//	//cout<<endl<<endl<<flush;
//	//cout<<"type of measuringline = "<<measuringline->gettype()<<endl;
//	//cout<<"(index_i,index_j) = ("<<index_i<<","<<index_j<<")"<<endl;
//	//cout<<"NdressChoice = "<<NdressChoice<<endl;
//	//cout<<"NdressChoiceInverse = "<<(1 + Nordermax - Norder + line_j_iter->getdress())<<endl;
//	//cout<<"Ndress = "<<Ndress<<endl<<flush;
//	//}
//
//	if(unidist(mt)<abs(probability)){
//		int vi, vj;
//		GadjLtmp.clear();
//		GadjLtmp.resize(vertices.size());
//		for(list<line>::iterator giter=Gs.begin(); giter!=Gs.end(); ++giter){
//			if(giter!=line_j_iter){
//				vi = distance(vertices.begin(),giter->getconnected_vi());
//				vj = distance(vertices.begin(),giter->getconnected_vj());
//				GadjLtmp[vi].push_back(vj);
//				GadjLtmp[vj].push_back(vi);
//			}
//		}
//		WmadjLtmp.clear();
//		WadjLtmp.clear();
//		WadjLtmp.resize(vertices.size());
//		for(list<line>::iterator witer=Ws.begin(); witer!=Ws.end(); ++witer){
//			vi = distance(vertices.begin(),witer->getconnected_vi());
//			vj = distance(vertices.begin(),witer->getconnected_vj());
//			if(witer!=line_j_iter){
//				WadjLtmp[vi].push_back(vj);
//				WadjLtmp[vj].push_back(vi);
//			}
//			else{
//				WmadjLtmp.push_back(vi);
//				WmadjLtmp.push_back(vj);
//			}
//		}
//
//		if(!checkconnectivity(GadjLtmp,WadjLtmp)) return 0;
//		if(!checkirreducibility(GadjLtmp,WadjLtmp)) return 0;
//		//if(!checkGcompactness(GadjLtmp,WmadjLtmp,WadjLtmp)) return 0;
//		//if(!checkWcompactness(GadjLtmp,WmadjLtmp,WadjLtmp)) return 0;
//		if(checkbubble(GadjLtmp)) return 0;
//
//		//if(!(checkconnectivityforinterswap(line_j_iter))){ return 0; }
//		//if(!(checkirreducibilityforinterswap(line_j_iter))){ return 0; }
//		//if(!(checkcompactnessforinterswap(line_j_iter))){ return 0; }
//
//
//		//SetUWInter = trial;
//		SetUWInter.setisloopflavorfixed(isloopflavorfixed_n);
//		SetUWInter.setloopflavor(loopflavor_n);
//		if( ismeasuringlineGline ){
//			if(loopi_iter==loopj_iter)
//				loopi_iter->setvalue(loopi_valuerh_n,loopi_valuelh_n);
//			else{
//				loopi_iter->setvalue(loopi_valuerh_n,loopi_valuelh_n);
//				loopj_iter->setvalue(loopj_valuerh_n,loopj_valuelh_n);
//			}
//		}
//		else{
//			if(islinei_inter)
//				SetUWInter.setphysical(*measuringline,true);
//			if(islinej_inter)
//				SetUWInter.setphysical(*line_j_iter,false);
//		}
//		SetUWInter.setvalue(valtrial);
//
//		measuringline->setphysical(true);
//		line_j_iter->setphysical(false);
//		measuringline = line_j_iter;
//
//		miflavor = miflavor_n;
//		mjflavor = mjflavor_n;
//
//		sign *= signfactor;
//
//		//cout<<"-->update accepted"<<endl;
//		return 3;
//	}
//
//	//cout<<"-->update rejected"<<endl;
//	return 1;
//};*/
///*int diagram::interswap(){
//	//cout<<"* interswap module *"<<endl<<flush;
//	int index_j;
//	int baredNorder;
//	list<line>::iterator line_j_iter;
//	if(measuringline->gettype()==0){
//		index_j = static_cast<int>( Ws.size()*unidist(mt) );
//		line_j_iter = next(Ws.begin(),index_j);
//		baredNorder = -1;
//	}
//	else{
//		index_j = static_cast<int>( Gs.size()*unidist(mt) );
//		line_j_iter = next(Gs.begin(),index_j);
//		baredNorder = +1;
//	}
//
//	int dNorder = baredNorder;
//	if(Norder+dNorder>(this->Nordermax)){ return 0; }
//
//	line line_i_new(*measuringline);
//	line_i_new.setphysical(true);
//
//	double proposal = 1.0;
//
//	//Wow if I turn off below two lines, it works! Super strange... Orz...
//	//It make sence.  It is artifact of closed loop MC
//	//if(measuringline->gettype()==0) proposal *= 0.5;
//	//else proposal *= 2.0;
//
//	double probability = 
//		-1 // Fermionic sign factor
//		* Ri_N[Norder] / Ri_N[Norder+dNorder];
//		// line_i_new.getvalue() / line_j_iter->getvalue();
//	
//	int miflavor_n, mjflavor_n;
//	int iloop, jloop;
//	list<vertex>::iterator vi_iter, vj_iter;
//	list<fermionicloop>::iterator loop_iter;
//	list<fermionicloop>::iterator loopi_iter, loopj_iter;
//
//	double loop_valuerh_n, loop_valuelh_n;
//
//	bool ismeasuringlineGline = (measuringline->gettype()==0);
//	vector<bool> isloopflavorfixed_n(fermionicloops.size(),false);
//	vector<int> loopflavor_n(fermionicloops.size(),0);
//	double valtrial;
//	//setUWinter trial = SetUWInter;
//	if(ismeasuringlineGline){
//		loop_iter = measuringline->getconnected_vi()->getconnected_loop();
//		loop_valuerh_n = loop_iter->getvaluerh()*line_i_new.getvalueij();
//		loop_valuelh_n = loop_iter->getvaluelh()*line_i_new.getvalueji();
//
//		vi_iter = line_j_iter->getconnected_vi();
//		vj_iter = line_j_iter->getconnected_vj();
//		loopi_iter = vi_iter->getconnected_loop();
//		loopj_iter = vj_iter->getconnected_loop();
//		if(loopi_iter==loopj_iter){
//			iloop  = distance(fermionicloops.begin(), loopi_iter);
//
//			miflavor_n = Nflavor*unidist(mt);
//			mjflavor_n = miflavor_n;
//
//			isloopflavorfixed_n[iloop] = true;
//			loopflavor_n[iloop] = miflavor_n;
//			valtrial = SetUWInter.calvalueafterinterswap(isloopflavorfixed_n,loopflavor_n);
//			//trial.releaseloopflavor();
//			//trial.fixloopflavor(iloop,miflavor_n);
//			//trial.calvalue();
//
//			probability *= 0.5*(loop_valuerh_n + loop_valuelh_n) / loop_iter->getvaluerh() 
//					/ line_j_iter->getvalueij()
//					* valtrial / SetUWInter.getvalue();
//					// * trial.getvalue() / SetUWInter.getvalue();
//			//probability *= (loop_valuerh_n + loop_valuelh_n) / loop_iter->getvaluerh() 
//			//		/ line_j_iter->getvalueij();
//		}
//		else{
//			iloop  = distance(fermionicloops.begin(), loopi_iter);
//			jloop  = distance(fermionicloops.begin(), loopj_iter);
//
//			if(line_j_iter->gettype()==1){
//				miflavor_n = Nflavor*unidist(mt);
//				mjflavor_n = (miflavor_n^1);
//			}
//			else{
//				proposal *= 2.;
//				miflavor_n = Nflavor*unidist(mt);
//				mjflavor_n = Nflavor*unidist(mt);
//			}
//
//			isloopflavorfixed_n[iloop] = true;
//			isloopflavorfixed_n[jloop] = true;
//			loopflavor_n[iloop] = miflavor_n;
//			loopflavor_n[jloop] = mjflavor_n;
//			valtrial = SetUWInter.calvalueafterinterswap(*line_j_iter,isloopflavorfixed_n,loopflavor_n);
//			//trial.setphysical(*line_j_iter,false);
//			//trial.releaseloopflavor();
//			//trial.fixloopflavor(iloop,miflavor_n);
//			//trial.fixloopflavor(jloop,mjflavor_n);
//			//trial.calvalue();
//
//			probability *= 0.5*(loop_valuerh_n + loop_valuelh_n) / loop_iter->getvaluerh() 
//				* valtrial / SetUWInter.getvalue();
//				// * trial.getvalue() / SetUWInter.getvalue();
//			//probability *= (loop_valuerh_n + loop_valuelh_n) / loop_iter->getvaluerh() 
//			//	* trial.getvalue() / SetUWInter.getvalue();
//		}
//	}
//	else{
//		loop_iter = line_j_iter->getconnected_vi()->getconnected_loop();
//		loop_valuerh_n = loop_iter->getvaluerh() / line_j_iter->getvalueij();
//		loop_valuelh_n = loop_iter->getvaluelh() / line_j_iter->getvalueji();
//
//		miflavor_n = Nflavor*unidist(mt);
//		mjflavor_n = miflavor_n;
//
//		vi_iter = line_j_iter->getconnected_vi();
//		loopi_iter = vi_iter->getconnected_loop();
//		iloop  = distance(fermionicloops.begin(), loopi_iter);
//
//		isloopflavorfixed_n[iloop] = true;
//		loopflavor_n[iloop] = miflavor_n;
//		if( measuringline->getconnected_vi()->getconnected_loop()
//			== measuringline->getconnected_vj()->getconnected_loop() ){
//			probability *= line_i_new.getvalueij();
//			valtrial = SetUWInter.calvalueafterinterswap(isloopflavorfixed_n,loopflavor_n);
//		}
//		else{
//			if(measuringline->gettype()==2){
//				proposal /= 2.;
//			}
//			//trial.setphysical(*measuringline,true);
//			valtrial = SetUWInter.calvalueafterinterswap(*measuringline,isloopflavorfixed_n,loopflavor_n);
//		}
//
//		//trial.releaseloopflavor();
//		//trial.fixloopflavor(iloop,miflavor_n);
//		//trial.calvalue();
//
//		probability *= loop_valuerh_n / 0.5/(loop_iter->getvaluerh()+loop_iter->getvaluelh())
//				* valtrial / SetUWInter.getvalue();
//				// * trial.getvalue() / SetUWInter.getvalue();
//		//probability *= loop_valuerh_n / (loop_iter->getvaluerh()+loop_iter->getvaluelh())
//		//		* trial.getvalue() / SetUWInter.getvalue();
//	}
//	probability *= proposal;
//
//	int signfactor;
//	if(probability<0.) signfactor = -1;
//	else signfactor = 1;
//
//	//if(NdressChoice==3 && measuringline->getisG() && Ndress==2){
//	//printconfiguration(cout);
//	//cout<<"* interswap info *"<<endl;
//	//cout<<"measuringline = "<<*measuringline<<endl;
//	//cout<<"line_i_new = "<<line_i_new<<endl;
//	//cout<<"line_j_iter = "<<*line_j_iter<<endl;
//	//cout<<"Norder = "<<Norder<<endl;
//	//cout<<"dNorder = "<<dNorder<<endl<<flush;
//	//cout<<"probability = "<<probability<<endl;
//	//cout<<"signfactor = "<<signfactor<<endl;
//	//cout<<"line_i_new.getvalue() = "<<line_i_new.getvalue()<<endl;
//	//cout<<"line_j_iter->getvalue() = "<<line_j_iter->getvalue()<<endl;
//	//cout<<endl<<endl<<flush;
//	//cout<<"type of measuringline = "<<measuringline->gettype()<<endl;
//	//cout<<"(index_i,index_j) = ("<<index_i<<","<<index_j<<")"<<endl;
//	//cout<<"NdressChoice = "<<NdressChoice<<endl;
//	//cout<<"NdressChoiceInverse = "<<(1 + Nordermax - Norder + line_j_iter->getdress())<<endl;
//	//cout<<"Ndress = "<<Ndress<<endl<<flush;
//	//}
//
//	if(unidist(mt)<abs(probability)){
//		int vi, vj;
//		GadjLtmp.clear();
//		GadjLtmp.resize(vertices.size());
//		for(list<line>::iterator giter=Gs.begin(); giter!=Gs.end(); ++giter){
//			if(giter!=line_j_iter){
//				vi = distance(vertices.begin(),giter->getconnected_vi());
//				vj = distance(vertices.begin(),giter->getconnected_vj());
//				if(vi!=vj){
//					GadjLtmp[vi].push_back(vj);
//					GadjLtmp[vj].push_back(vi);
//				}
//			}
//		}
//		WmadjLtmp.clear();
//		WadjLtmp.clear();
//		WadjLtmp.resize(vertices.size());
//		//cout<<"* WadjLtmp before *"<<endl;
//		//for(auto i: WadjLtmp){
//		//	for(auto j: i) cout<<j;
//		//	cout<<endl;
//		//}
//
//		for(list<line>::iterator witer=Ws.begin(); witer!=Ws.end(); ++witer){
//			vi = distance(vertices.begin(),witer->getconnected_vi());
//			vj = distance(vertices.begin(),witer->getconnected_vj());
//			if(witer!=line_j_iter){
//				WadjLtmp[vi].push_back(vj);
//				WadjLtmp[vj].push_back(vi);
//			}
//			else{
//				WmadjLtmp.push_back(vi);
//				WmadjLtmp.push_back(vj);
//			}
//		}
//
//		//printconfiguration(cout);
//		//cout<<"* WmadjLtmp *"<<endl;
//		//for(auto i: WmadjLtmp) cout<<i;
//		//cout<<endl;
//		//cout<<"* WadjLtmp *"<<endl;
//		//for(auto i: WadjLtmp){
//		//	for(auto j: i) cout<<j;
//		//	cout<<endl;
//		//}
//		//if(checkconnectivity(GadjLtmp,WadjLtmp)!=checkconnectivityforreconnect(Giin_iter,Gjin_iter)){
//		//	cout<<"reconnect update; connectivity check problem"<<endl;
//		//	exit(EXIT_FAILURE);
//		//}
//		//if(checkirreducibility(GadjLtmp,WadjLtmp)!=checkirreducibilityforreconnect(Giin_iter,Gjin_iter)){
//		//	cout<<"reconnect update; irreducibility check problem"<<endl;
//		//	printconfiguration(cout);
//		//	cout<<"is connected after update? "<<checkconnectivity(GadjLtmp,WadjLtmp)<<endl;
//		//	exit(EXIT_FAILURE);
//		//}
//		//if( (checkGcompactness(GadjLtmp,WmadjLtmp,WadjLtmp)&&checkWcompactness(GadjLtmp,WmadjLtmp,WadjLtmp))!=checkcompactnessforinterswap(line_j_iter) ){
//		//	cout<<"interswap update; compactness check problem"<<endl;
//		//	//printconfiguration(cout);
//		//	cout<<"* GadjLtmp *"<<endl;
//		//	for(auto i: GadjLtmp){
//		//		for(auto j: i) cout<<j;
//		//		cout<<endl;
//		//	}
//		//	cout<<"* WmadjLtmp *"<<endl;
//		//	for(auto i: WmadjLtmp) cout<<i;
//		//	cout<<endl;
//		//	cout<<"* WadjLtmp *"<<endl;
//		//	for(auto i: WadjLtmp){
//		//		for(auto j: i) cout<<j;
//		//		cout<<endl;
//		//	}
//		//	cout<<flush;
//		//	exit(EXIT_FAILURE);
//		//}
//
//		if(!checkconnectivity(GadjLtmp,WadjLtmp)) return 0;
//		if(!checkirreducibility(GadjLtmp,WadjLtmp)) return 0;
//		//if(!checkGcompactness(GadjLtmp,WmadjLtmp,WadjLtmp)) return 0;
//		//if(!checkWcompactness(GadjLtmp,WmadjLtmp,WadjLtmp)) return 0;
//		if(checkbubble(GadjLtmp)) return 0;
//		//if(!(checkconnectivityforinterswap(line_j_iter))){ return 0; }
//		//if(!(checkirreducibilityforinterswap(line_j_iter))){ return 0; }
//		//if(!(checkcompactnessforinterswap(line_j_iter))){ return 0; }
//
//		loop_iter->setvalue(loop_valuerh_n,loop_valuelh_n);
//
//		//SetUWInter = trial;
//		SetUWInter.setisloopflavorfixed(isloopflavorfixed_n);
//		SetUWInter.setloopflavor(loopflavor_n);
//		if( ismeasuringlineGline
//			&& ( line_j_iter->getconnected_vi()->getconnected_loop() 
//				!= line_j_iter->getconnected_vj()->getconnected_loop() )
//		){ SetUWInter.setphysical(*line_j_iter,false); }
//		else if( (!ismeasuringlineGline)
//			&& ( measuringline->getconnected_vi()->getconnected_loop() 
//				!= measuringline->getconnected_vj()->getconnected_loop() )
//	        ){ SetUWInter.setphysical(*measuringline,true); }
//		SetUWInter.setvalue(valtrial);
//
//		measuringline->setphysical(true);
//		line_j_iter->setphysical(false);
//		measuringline = line_j_iter;
//
//		miflavor = miflavor_n;
//		mjflavor = mjflavor_n;
//
//		Norder += dNorder;
//		sign *= signfactor;
//
//		//cout<<"-->update accepted"<<endl;
//		return 3;
//	}
//
//	//cout<<"-->update rejected"<<endl;
//	return 1;
//};*/
///*int diagram::transformUW(){
//	list<list<line>::iterator> UWinters;
//	list<line>::iterator UW_iter;
//	list<fermionicloop>::iterator loopi_iter, loopj_iter;
//	for(UW_iter=Ws.begin(); UW_iter!=Ws.end(); ++UW_iter){
//		loopi_iter = UW_iter->getconnected_vi()->getconnected_loop();
//		loopj_iter = UW_iter->getconnected_vj()->getconnected_loop();
//		if( loopi_iter!=loopj_iter && UW_iter->getphysical() ) UWinters.push_back(UW_iter);
//	}
//	if(UWinters.size()==0) return 0;
//	else{
//		int NUWinters = UWinters.size();
//		int iUWinters = static_cast<int>(NUWinters*unidist(mt));
//		int ijselection = static_cast<int>(2*unidist(mt));
//
//		UW_iter = (*next(UWinters.begin(),iUWinters));
//		int UWtype = UW_iter->gettype();
//		list<vertex>::iterator v_iter, v_oppo;
//		if(ijselection==0){
//			v_iter = UW_iter->getconnected_vi();
//			v_oppo = UW_iter->getconnected_vj();
//		}
//		else{
//			v_iter = UW_iter->getconnected_vj();
//			v_oppo = UW_iter->getconnected_vi();
//		}
//
//		double proposal;
//		double probability;
//
//		coordinate x_n;
//		line UW_n = *UW_iter;
//		if(UWtype==1){
//			x_n.genrand(v_iter->getx());
//
//			proposal = beta/(v_iter->getx().getsitegenprobability(x_n));
//			if(ijselection==0) 
//				UW_n = line(UW_iter->getphysical(),2,UW_iter->getdress(),UW_iter->getiflavor(),UW_iter->getjflavor(),x_n,UW_iter->getxj());
//			else 
//				UW_n = line(UW_iter->getphysical(),2,UW_iter->getdress(),UW_iter->getiflavor(),UW_iter->getjflavor(),UW_iter->getxi(),x_n);
//		}
//		else{
//			x_n = v_oppo->getx();
//
//			proposal = (v_iter->getx().getsitegenprobability(x_n))/beta;
//			if(ijselection==0)
//				UW_n = line(UW_iter->getphysical(),1,UW_iter->getdress(),UW_iter->getiflavor(),UW_iter->getjflavor(),x_n,UW_iter->getxj());
//			else 
//				UW_n = line(UW_iter->getphysical(),1,UW_iter->getdress(),UW_iter->getiflavor(),UW_iter->getjflavor(),UW_iter->getxi(),x_n);
//		}
//
//		list<fermionicloop>::iterator loop_iter = v_iter->getconnected_loop();
//		bool ismeasuringlineinloop = ( measuringline->getconnected_vi()->getconnected_loop()==loop_iter );
//
//		//setUWinter trial(SetUWInter);;
//		//trial.transform(*UW_iter,UW_n);
//		//trial.calvalue();
//
//		//printconfiguration(cout);
//		double valtrial = SetUWInter.calvalueaftertransformUW(*UW_iter,UW_n);
//		//cout<<"*UW_iter: "<<*UW_iter<<endl;
//		//cout<<"UW_n: "<<UW_n<<endl;
//		//cout<<"  valtrial = "<<valtrial<<endl;
//
//		list<line>::iterator Gin_iter = v_iter->getconnected_Gin();
//		list<line>::iterator Gout_iter = v_iter->getconnected_Gout();
//		bool is_selfloop = (Gin_iter==Gout_iter);
//		line Gin_n = *Gin_iter;
//		line Gout_n = *Gout_iter;
//
//		double loop_valuerh_n, loop_valuelh_n;
//
//		if(is_selfloop){
//			Gin_n.setx(x_n,x_n);
//
//			probability = proposal 
//				* valtrial / SetUWInter.getvalue();
//			//probability = proposal 
//			//	* trial.getvalue() / SetUWInter.getvalue();
//		}
//		else{
//			Gin_n.setxj(x_n);
//			Gout_n.setxi(x_n);
//			loop_valuerh_n = Gin_n.getvalueij() * Gout_n.getvalueij() 
//					/ Gin_iter->getvalueij() / Gout_iter->getvalueij() 
//					* loop_iter->getvaluerh();
//			loop_valuelh_n = Gin_n.getvalueji() * Gout_n.getvalueji() 
//					/ Gin_iter->getvalueji() / Gout_iter->getvalueji() 
//					* loop_iter->getvaluelh();
//
//			if(ismeasuringlineinloop){
//				probability = proposal
//						* loop_valuerh_n / loop_iter->getvaluerh()
//						* valtrial / SetUWInter.getvalue();
//				//probability = proposal
//				//		* loop_valuerh_n / loop_iter->getvaluerh()
//				//		* trial.getvalue() / SetUWInter.getvalue();
//			}
//			else{
//				probability = proposal
//						* (loop_valuerh_n+loop_valuelh_n) 
//						/ (loop_iter->getvaluerh()+loop_iter->getvaluelh())
//						* valtrial / SetUWInter.getvalue();
//				//probability = proposal
//				//		* (loop_valuerh_n+loop_valuelh_n) 
//				//		/ (loop_iter->getvaluerh()+loop_iter->getvaluelh())
//				//		* trial.getvalue() / SetUWInter.getvalue();
//			}
//		}
//
//		int signfactor;
//		if(probability<0.) signfactor = -1;
//		else signfactor = 1;
//
//		if(unidist(mt)<abs(probability)){
//			SetUWInter.transform(*UW_iter,UW_n);
//			SetUWInter.setvalue(valtrial);
//
//			v_iter->setx(x_n);
//			UW_iter->assign(UW_n);
//			if(is_selfloop) 
//				Gin_iter->assign(Gin_n);
//			else{
//				Gin_iter->assign(Gin_n);
//				Gout_iter->assign(Gout_n);
//				loop_iter->setvalue(loop_valuerh_n,loop_valuelh_n);
//			}
//
//			//SetUWInter = trial;
//
//			sign *= signfactor;
//
//			bool measuring_line_involved = true;
//			measuring_line_involved &= UW_iter->getphysical();
//			measuring_line_involved &= v_iter->getconnected_Gin()->getphysical();
//			measuring_line_involved &= v_iter->getconnected_Gout()->getphysical();
//
//			if(measuring_line_involved) return 3;
//			else return 2;
//		}
//		return 1;
//	}
//};*/
///*int diagram::intraswap(){
//	int index_i, index_j;
//	list<line>::iterator line_j_iter;
//	if(measuringline->gettype()==0){
//		index_i = distance(Gs.begin(),measuringline);
//		select_one_excluding_one(Gs.size(),index_i,index_j);
//		line_j_iter = next(Gs.begin(),index_j);
//	}
//	else{
//		int NWs = Ws.size();
//		//int NWs = Us.size()+Ws.size();
//		if(NWs==1) return 0;
//		if(measuringline->gettype()!=0)
//			index_i = distance(Ws.begin(),measuringline);
//		//if(measuringline->gettype()==1)
//		//	index_i = distance(Us.begin(),measuringline);
//		//else 
//		//	index_i = Us.size()+distance(Ws.begin(),measuringline);
//		select_one_excluding_one(NWs,index_i,index_j);
//
//		line_j_iter = next(Ws.begin(),index_j);
//		//if(index_j<Us.size()) line_j_iter = next(Us.begin(),index_j);
//		//else line_j_iter = next(Ws.begin(),index_j-Us.size());
//	}
//
//	int dNorder = measuringline->getdress() - line_j_iter->getdress();
//	if(Norder+dNorder>(this->Nordermax)){ return 0; }
//
//	line line_i_new(*measuringline);
//	line_i_new.setphysical(true);
//
//	double probability = 
//		Ri_N[Norder] / Ri_N[Norder+dNorder]
//		* line_i_new.getvalue() / line_j_iter->getvalue();
//
//	int signfactor;
//	if(probability<0.) signfactor = -1;
//	else signfactor = 1;
//
//	//cout<<"* intraswap info *"<<endl;
//	//cout<<"does measuringline indicate the G-line? "<<measuringline->getisG()<<endl;
//	//cout<<"(index_i,index_j) = ("<<index_i<<","<<index_j<<")"<<endl;
//	//cout<<"Norder = "<<Norder<<endl;
//	//cout<<"probability = "<<probability<<endl;
//	//cout<<"signfactor = "<<signfactor<<endl;
//	//printcionfiguration(cout);
//
//	if(unidist(mt)<abs(probability)){
//		measuringline->setphysical(true);
//		line_j_iter->setphysical(false);
//		measuringline = line_j_iter;
//
//		Norder += dNorder;
//		sign *= signfactor;
//
//		//cout<<"-->update accepted"<<endl;
//		return 3;
//	}
//
//	//cout<<"-->update rejected"<<endl;
//	return 1;
//};*/
///*int diagram::absorbbubble(){
//	//cout<<"* absorbbubble module *"<<endl;
//	updatebubbles();
//	if(bubbles.size()==0 || Us.size()==1) return 0;
//
//	int Nbubble = bubbles.size();
//	int ibubble = static_cast<int>(Nbubble*unidist(mt));
//	list<bubble>::iterator biter = next(bubbles.begin(),ibubble);
//
//	list<vertex>::iterator vi_outer, vj_outer, vi_inner, vj_inner;
//	list<line>::iterator Ui_iter, Uj_iter;
//
//	vi_inner = biter->getconnected_vi();
//	vj_inner = biter->getconnected_vj();
//	Ui_iter = vi_inner->getconnected_W();
//	Uj_iter = vj_inner->getconnected_W();
//
//	if(Ui_iter->gettype()==2 || Uj_iter->gettype()==2){ return 0; }
//
//	vi_outer = ((Ui_iter->getconnected_vi()==vi_inner) ? (Ui_iter->getconnected_vj()) : (Ui_iter->getconnected_vi()));
//	vj_outer = ((Uj_iter->getconnected_vi()==vj_inner) ? (Uj_iter->getconnected_vj()) : (Uj_iter->getconnected_vi()));
//
//	int Wdress = 0;
//	int W_type = 2;
//	bool bothphysical = (Ui_iter->getphysical() && Uj_iter->getphysical());
//	if(!bothphysical) return 0;
//	//if(!bothphysical && Wdress>Nauxdress) return 0;
//	line Wnew(bothphysical, W_type, Wdress, vi_outer->getflavor(), vj_outer->getflavor(), vi_outer->getx(), vj_outer->getx());
//
//
//	int dNorder = -1;
//	if(Norder+dNorder>(this->Nordermax)){ return 0; }
//
//	double proposal = 1.0;
//	if(!bothphysical) proposal /= 2.; 
//	proposal *= static_cast<double>(Nbubble) / (Ws.size()+1);
//	double probability = proposal 
//			* Ri_N[Norder] / Ri_N[Norder+dNorder]
//			* Wnew.getvalue() 
//			/ Ui_iter->getvalue() / biter->getvalue() / Uj_iter->getvalue();
//
//	int signfactor;
//	if(probability<0.) signfactor = -1;
//	else signfactor = 1;
//
//	//cout<<"* absorbbubble module *"<<endl;
//	//cout<<" selected bubble : "<<endl<<*biter<<endl;
//	//cout<<" probability = "<<probability<<endl;
//	//cout<<" proposal = "<<proposal<<endl;
//	//cout<<" Nbubble = "<<Nbubble<<endl;
//	//cout<<" Norder = "<<Norder<<endl;
//	//cout<<" dNorder = "<<dNorder<<endl;
//	//cout<<" Ri_N[Norder] = "<<Ri_N[Norder]<<endl;
//	//cout<<" Ri_N[Norder+dNorder] = "<<Ri_N[Norder+dNorder]<<endl;
//	//cout<<" Ws.size() = "<<Ws.size()<<endl;
//	//cout<<" Wnew.getvalue() = "<<Wnew.getvalue()<<endl;
//	//cout<<" Ui_iter->getvalue() = "<<Ui_iter->getvalue()<<endl;
//	//cout<<" biter->getvalue() = "<<biter->getvalue()<<endl;
//	//cout<<" Uj_iter->getvalue() = "<<Uj_iter->getvalue()<<endl;
//	//printconfiguration(cout);
//
//	if(unidist(mt)<abs(probability)){
//		Ws.push_front(Wnew); list<line>::iterator W_iter = Ws.begin();
//		//W_iter->setphysical( Ui_iter->getphysical() && Uj_iter->getphysical() );
//		W_iter->setconnected_vertices(vi_outer, vj_outer);
//		vi_outer->setconnected_W(W_iter);
//		vj_outer->setconnected_W(W_iter);
//		if(!bothphysical) measuringline = W_iter;
//
//		Us.erase(Ui_iter);
//		Us.erase(Uj_iter);
//		Gs.erase(biter->getGi());
//		Gs.erase(biter->getGj());
//		vertices.erase(vi_inner);
//		vertices.erase(vj_inner);
//		bubbles.erase(biter);
//
//		Norder += dNorder;
//		sign *= signfactor;
//
//		//cout<<"--> update accepted"<<endl;
//		return 2;
//	}
//
//	//cout<<"--> update rejected"<<endl;
//	return 1;
//};
//int diagram::ejectbubble(){
//	//cout<<"* ejectbubble module *"<<endl;
//	updatebubbles();
//	int Nbubble = bubbles.size();
//	int NW = Ws.size();
//	if(NW<1) return 0;
//	int iW = static_cast<int>(NW*unidist(mt));
//
//	list<line>::iterator W_iter = next(Ws.begin(),iW);
//	list<vertex>::iterator vi_outer = W_iter->getconnected_vi();
//	list<vertex>::iterator vj_outer = W_iter->getconnected_vj();
//
//	if(vi_outer->getflavor()!=vj_outer->getflavor()) return 0;
//
//	int Worder = 1 + W_iter->getdress();
//	int iflavor = (vi_outer->getflavor()^1);
//
//	double proposal = 1.;
//	coordinate xi_n, xj_n;
//
//	bool Ui_physical, Uj_physical;
//	bool measuring_line_involved = (!(W_iter->getphysical()));
//	int dNorder = 1;
//	if(Norder+dNorder>(this->Nordermax)) return 0;
//	if(measuring_line_involved){
//		return 0;
//		//if(unidist(mt)<0.5){ 
//		//	Ui_physical = true;
//		//       	Uj_physical = false; 
//		//}
//		//else{ 
//		//	Ui_physical = false;
//		//       	Uj_physical = true; 
//		//}
//		//proposal *= 2.;
//	}
//	else{
//		Ui_physical = true;
//		Uj_physical = true;
//	}
//
//	xi_n = vi_outer->getx();
//	xj_n = vj_outer->getx();
//
//	line Ui_new(Ui_physical, 1, 0, vi_outer->getflavor(), iflavor, vi_outer->getx(), xi_n);
//	line Gi_new(true, 0, 0, iflavor, iflavor, xi_n, xj_n);
//	line Gj_new(true, 0, 0, iflavor, iflavor, xj_n, xi_n);
//	line Uj_new(Uj_physical, 1, 0, iflavor, vj_outer->getflavor(), xj_n, vj_outer->getx());
//
//	proposal *= static_cast<double>(NW)/(Nbubble+1);
//	double probability = proposal 
//			* Ri_N[Norder] / Ri_N[Norder+dNorder]
//			* Ui_new.getvalue() * (-1.) * Gi_new.getvalue() * Gj_new.getvalue() * Uj_new.getvalue() 
//			/ W_iter->getvalue();
//
//	int signfactor;
//	if(probability<0.) signfactor = -1;
//	else signfactor = 1;
//
//	//cout<<"* ejectbubble module *"<<endl;
//	//cout<<"selected W-line: "<<(*W_iter)<<endl;
//	//cout<<"proposal = "<<proposal<<endl;
//	//cout<<"probability = "<<probability<<endl;
//	//cout<<"Nbubble = "<<Nbubble<<endl;
//	//cout<<"NW = "<<NW<<endl;
//	//cout<<"(Ui_physical,Uj_physical) = ("<<Ui_physical<<","<<Uj_physical<<")"<<endl;
//	//cout<<"Ui_new.getvalue() = "<<Ui_new.getvalue()<<endl;
//	//cout<<"Gi_new.getvalue() = "<<Gi_new.getvalue()<<endl;
//	//cout<<"Gj_new.getvalue() = "<<Gj_new.getvalue()<<endl;
//	//cout<<"Uj_new.getvalue() = "<<Uj_new.getvalue()<<endl;
//	//printconfiguration(cout);
//
//	if(unidist(mt)<abs(probability)){
//		Ws.erase(W_iter);
//		Us.push_front(Ui_new); list<line>::iterator Ui_iter = Us.begin();
//		Us.push_front(Uj_new); list<line>::iterator Uj_iter = Us.begin();
//		Gs.push_front(Gi_new); list<line>::iterator Gi_iter = Gs.begin();
//		Gs.push_front(Gj_new); list<line>::iterator Gj_iter = Gs.begin();
//
//		vertices.emplace_front(iflavor,xi_n); list<vertex>::iterator vi_iter = vertices.begin();
//		vertices.emplace_front(iflavor,xj_n); list<vertex>::iterator vj_iter = vertices.begin();
//
//		vi_iter->setconnected_lines(Ui_iter,Gj_iter,Gi_iter);
//		vj_iter->setconnected_lines(Uj_iter,Gi_iter,Gj_iter);
//		vi_outer->setconnected_W(Ui_iter);
//		vj_outer->setconnected_W(Uj_iter);
//		Gi_iter->setconnected_vertices(vi_iter,vj_iter);
//		Gj_iter->setconnected_vertices(vj_iter,vi_iter);
//		Ui_iter->setconnected_vertices(vi_outer,vi_iter);
//		Uj_iter->setconnected_vertices(vj_iter,vj_outer);
//
//		if(measuring_line_involved){
//			if(!Ui_physical) measuringline = Ui_iter;
//			else measuringline = Uj_iter;
//		}
//
//		bubbles.emplace_front(Gi_iter,Gj_iter);
//
//		Norder += dNorder;
//		sign *= signfactor;
//
//		//cout<<"--> update accepted"<<endl;
//		return 2;
//	}
//
//	//cout<<"--> update rejected"<<endl;
//	return 1;
//};
//int diagram::spinflip(){
//	list<vertex>::iterator viter;
//	list<line>::iterator liter;
//	int iflavor, jflavor;
//
//	for(viter=vertices.begin(); viter!=vertices.end(); ++viter){
//		iflavor = viter->getflavor();
//		viter->setflavor((iflavor^1));
//	}
//	for(liter=Gs.begin(); liter!=Gs.end(); ++liter){
//		iflavor = liter->getiflavor();
//		jflavor = liter->getjflavor();
//		liter->setflavor((iflavor^1),(jflavor^1));
//	}
//	for(liter=Us.begin(); liter!=Us.end(); ++liter){
//		iflavor = liter->getiflavor();
//		jflavor = liter->getjflavor();
//		liter->setflavor((iflavor^1),(jflavor^1));
//	}
//	for(liter=Ws.begin(); liter!=Ws.end(); ++liter){
//		iflavor = liter->getiflavor();
//		jflavor = liter->getjflavor();
//		liter->setflavor((iflavor^1),(jflavor^1));
//	}
//
//	return 2;
//};*/
///*int diagram::loopspinflip(){
//	//cout<<"** loop spin flip update **"<<endl<<flush;
//	//cout<<"   updatefermionicloops"<<endl<<flush;
//	updatefermionicloops();
//
//	//for(auto fl:fermionicloops) cout<<fl;
//	
//	int Nloop = fermionicloops.size();
//	int iloop = static_cast<int>(Nloop*unidist(mt));
//
//	list<fermionicloop>::iterator fl_loop = next(fermionicloops.begin(),iloop);
//	fl_loop->proposeloopspinflip();
//
//	double probability = fl_loop->get_loopspinflip_probability();
//	//cout<<"probability = "<<probability<<endl;
//
//	int signfactor;
//	if(probability<0.) signfactor = -1;
//	else signfactor = 1;
//
//	//cout<<"* loopspinflip module *"<<endl;
//	//cout<<"probability = "<<probability<<endl;
//	//printconfiguration(cout);
//
//	if(unidist(mt)<abs(probability)){
//		fl_loop->loopspinflip();
//		sign *= signfactor;
//
//		//cout<<"--> update accepted"<<endl;
//		//printconfiguration(cout);
//		return 3;
//	}
//
//	//for(auto fl:fermionicloops) cout<<fl;
//	//exit(EXIT_FAILURE);
//
//	return 0;
//};*/
//
///*void diagram::setSwx(){
//	double tji = measuringline->getxi().time - measuringline->getxj().time;
//
//	double w0t = w(0)*tji;
//	double dwt = 2.*M_PI*tji/beta;
//
//	complex<double> Swx_tmp = complex<double>(cos(w0t),sin(w0t));
//	complex<double> factor = complex<double>(cos(dwt),sin(dwt));
//
//	for(int n=0; n<Owx.size(); n++){
//		Owx(n) = Swx_tmp;
//		Swx_tmp *= factor;
//	}
//};
//void diagram::setPwx(){
//	double tji = measuringline->getxi().time - measuringline->getxj().time;
//
//	double w0t = wB(0)*tji;
//	double dwt = 2.*M_PI*tji/beta;
//
//	complex<double> Pwx_tmp = complex<double>(cos(w0t),sin(w0t));
//	complex<double> factor = complex<double>(cos(dwt),sin(dwt));
//
//	for(int n=0; n<Owx.size(); n++){
//		Owx(n) = Pwx_tmp;
//		Pwx_tmp *= factor;
//	}
//};*/
//
void diagram::measure(){
	//cout<<"* measure *"<<endl;
	//int iflavor = measuringline->getjflavor();
	//int jflavor = measuringline->getiflavor();

	//int measuringlineisselfloop = 0;
	//if(measuringline->getconnected_vi()==measuringline->getconnected_vj()) measuringlineisselfloop++;

	//int m_type = measuringline->gettype();
	//if(m_type==0){ 
	//	if(Norder==1 && (Ws.front().gettype()==1)) normalization[0][0]++; 
	//	else normalization[0][Norder]++; 
	//}
	//else if(m_type==2) normalization[1][Norder]++;

	normalization[Norder] += 1;
	if(Norder==2)
		//normalization[0] += 6.0 / Conf.getWeight() / pow(beta,3);
		//normalization[0] += 6.0 / Conf.getWeight() / pow(beta,3) / pow(2,4);
		normalization[0] += 6.0 / Conf.getWeight() / pow(parameters::beta,3) / pow(2,4);

	//NorderAcc += Norder;

	vector<coordinate> vertexExt = Conf.getExtVertices();
	vector<int> flavorExt = Conf.getFlavorExt();
	vector<int> itList(3);
	itList[2] = static_cast<int>( (vertexExt[0].gettime()+0.5*dtau)/dtau );
	itList[1] = static_cast<int>( (vertexExt[1].gettime()+0.5*dtau)/dtau );
	itList[0] = static_cast<int>( (vertexExt[2].gettime()+0.5*dtau)/dtau );

	//cout<<"dtau = "<<dtau<<endl;
	//cout<<"vertexExt: "<<endl;
	//for(auto& it: vertexExt) cout<<it<<" ";
	//cout<<endl;
	//cout<<"flavorExt: "<<endl;
	//for(auto& it: flavorExt) cout<<it<<" ";
	//cout<<endl;
	//cout<<"indexTime: "<<endl;
	//for(auto& it: itList) cout<<it<<" ";
	//cout<<endl;

	//double accValue = sign;
	double accValue = sign/pow(dtau,3);

	if(itList[0]==0 || itList[0]==NtG-1) accValue *= 2.0;
	if(itList[1]==0 || itList[1]==NtG-1) accValue *= 2.0;
	if(itList[2]==0 || itList[2]==NtG-1) accValue *= 2.0;

	if(itList[0]==itList[1] && itList[1]!=itList[2]) accValue *= 2.0;
	if(itList[0]!=itList[1] && itList[1]==itList[2]) accValue *= 2.0;
	if(itList[0]==itList[1] && itList[1]==itList[2]) accValue *= 6.0;

	//cout<<(*measuring_line)<<endl;

	/*if(Norder==2){
		printconfiguration(cout);
		cout<<"("<<itList[0]<<","<<itList[1]<<","<<itList[2]<<")"<<endl;
		for(auto& iter: vertexExt) cout<<iter;
		cout<<endl;
	}*/

	Qacc[Norder-1].accumulate(flavorExt,itList,accValue);
	
	//cout<<"after measurement Qcc value: "<<Qacc[Norder-1].getQt(itList[0],itList[1],itList[2],flavorExt[0],flavorExt[1])<<endl<<endl;

	//if(Nselfloop-measuringlineisselfloop==0){
	//	if(m_type==0){
	//		Swxacc[Norder-1][miflavor][mjflavor].accumulateOwx(
	//			Lattice.getXindex(measuringline->getxi().space - measuringline->getxj().space),
	//			Owx*static_cast<double>(sign*Ri_N[Norder])
	//			);
	//	}
	//	else if(m_type==1){
	//		Plocacc[Norder] += static_cast<double>(sign*Ri_N[Norder]);
	//	}
	//	else{
	//		Pwxacc[Norder][miflavor][mjflavor].accumulateOwx(
	//			Lattice.getXindex(measuringline->getxi().space - measuringline->getxj().space),
	//			Owx*static_cast<double>(sign*Ri_N[Norder])
	//			);
	//	}
	//}
};
///*void diagram::measure_S(){
//	int iflavor = measuringline->getjflavor();
//	int jflavor = measuringline->getiflavor();
//
//	int measuringlineisselfloop = 0;
//	if(measuringline->getconnected_vi()==measuringline->getconnected_vj()) measuringlineisselfloop++;
//
//	if(Norder==1 && Us.size()==1) normalization[0][0]++; 
//	else normalization[0][Norder]++; 
//	if(Nselfloop-measuringlineisselfloop==0){
//		Swxacc[Norder-1][iflavor][jflavor].accumulateOwx(
//			Lattice.getXindex(measuringline->getxi().space - measuringline->getxj().space),
//			Owx*static_cast<double>(sign*Ri_N[Norder])
//			);
//		//Swxacc[Norder-1][iflavor][jflavor] += (Owx*static_cast<double>(sign));
//	}
//};
//void diagram::measure_S(const int& iorder){
//	int iflavor = measuringline->getjflavor();
//	int jflavor = measuringline->getiflavor();
//
//	int measuringlineisselfloop = 0;
//	if(measuringline->getconnected_vi()==measuringline->getconnected_vj()) measuringlineisselfloop++;
//
//	if(Norder==1 && Us.size()==1) normalization[0][0]++; 
//	else normalization[0][Norder]++; 
//
//	if((Nselfloop-measuringlineisselfloop==0) && (Norder==iorder)){
//		Swxacc[Norder-1][iflavor][jflavor].accumulateOwx(
//			Lattice.getXindex(measuringline->getxi().space - measuringline->getxj().space),
//			Owx*static_cast<double>(sign*Ri_N[Norder])
//			);
//	}
//};
//void diagram::measure_P(){
//	int iflavor = measuringline->getjflavor();
//	int jflavor = measuringline->getiflavor();
//
//	int measuringlineisselfloop = 0;
//	if(measuringline->getconnected_vi()==measuringline->getconnected_vj()) measuringlineisselfloop++;
//
//	normalization[1][Norder]++;
//	if(Nselfloop-measuringlineisselfloop==0){
//		Pwxacc[Norder][iflavor][jflavor].accumulateOwx(
//			Lattice.getXindex(measuringline->getxi().space - measuringline->getxj().space),
//			Owx*static_cast<double>(sign*Ri_N[Norder])
//			);
//		//printconfiguration(cout);
//	}
//};
//void diagram::measure_P(const int& iorder){
//	int iflavor = measuringline->getjflavor();
//	int jflavor = measuringline->getiflavor();
//
//	int measuringlineisselfloop = 0;
//	if(measuringline->getconnected_vi()==measuringline->getconnected_vj()) measuringlineisselfloop++;
//
//	normalization[1][Norder]++;
//	if((Nselfloop-measuringlineisselfloop==0) && (Norder==iorder)){
//		Pwxacc[Norder][iflavor][jflavor].accumulateOwx(
//			Lattice.getXindex(measuringline->getxi().space - measuringline->getxj().space),
//			Owx*static_cast<double>(sign*Ri_N[Norder])
//			);
//		//printconfiguration(cout);
//	}
//};*/
////void diagram::debug(){
////	initDiagram();
////	printconfiguration(cout);
////
////	int ox;
////	/*for(int ith=0; ith<Nth; ith++){
////		//cout<<"*** "<<ith<<"th iterator ***"<<endl;
////		//cout<<"update: ";
////		switch( static_cast<int>(Nupdate*unidist(mt)) ){
////			case 0:
////				//cout<<"insertVertices"<<endl;
////				insertVertices();
////				break;
////			case 1:
////				//cout<<"removeVertices"<<endl;
////				removeVertices();
////				break;
////			case 2:
////				//cout<<"shiftVertex"<<endl;
////				shiftVertex();
////				break;
////			default:
////				//cout<<"swapMeasuringLine"<<endl;
////				swapMeasuringLine();
////				break;
////		}
////	}*/
////
////	//cout<<"shiftVertex"<<endl;
////	//ox = shiftVertex();
////	//printconfiguration(cout);
////
////	for(long long imc=0; imc<Nmc; imc++){
////		cout<<"*** "<<imc<<"th iterator ***"<<endl<<"update: ";
////		switch( static_cast<int>(Nupdate*unidist(mt)) ){
////			case 0:
////				cout<<"addvertice"<<endl;
////				ox = insertVertices();
////				break;
////			case 1:
////				cout<<"removevertice"<<endl;
////				ox = removeVertices();
////				break;
////			case 2:
////				cout<<"shiftVertex"<<endl;
////				ox = shiftVertex();
////				break;
////			default:
////				cout<<"swapMeasuringLine"<<endl;
////				ox = swapMeasuringLine();
////				break;
////		}
////		printconfiguration(cout);
////		if(ox>1){
////			// checking modules
////			//printconfiguration(cout);
////			if( sign<0.0 ){
////				cout<<"*** "<<imc<<"th iterator ***"<<endl;
////				cout<<"sign is negative"<<endl;
////				printconfiguration(cout);
////				exit(EXIT_FAILURE);;
////			}
////			//if(!checkNorder()){
////			//	cout<<"Noder check wrong"<<endl;
////			//	cout<<"*** "<<imc<<"th iterator ***"<<endl;
////			//	exit(EXIT_FAILURE);;
////			//}
////			//if(!checkfermionicloop()){
////			//	cout<<"fermionicloop check wrong"<<endl;
////			//	cout<<"*** "<<imc<<"th iterator ***"<<endl;
////			//	exit(EXIT_FAILURE);;
////			//}
////			//if(!checkUWlines()){
////			//	cout<<"Ulines check wrong"<<endl;
////			//	cout<<"*** "<<imc<<"th iterator ***"<<endl;
////			//	exit(EXIT_FAILURE);;
////			//}
////			//if(!checksetUWinters()){
////			//	cout<<"SetUWInter check wrong"<<endl;
////			//	cout<<"*** "<<imc<<"th iterator ***"<<endl;
////			//	exit(EXIT_FAILURE);;
////			//}
////			//if(!checkdress()){
////			//	cout<<"dress check wrong"<<endl;
////			//	cout<<"*** "<<imc<<"th iterator ***"<<endl;
////			//	exit(EXIT_FAILURE);;
////			//}
////			//if(!checkNselfloop()){
////			//	cout<<"*** "<<imc<<"th iterator ***"<<endl;
////			//	cout<<"Nselfloop check wrong"<<endl;
////			//	exit(EXIT_FAILURE);;
////			//}
////			//if(!checkdiagramsign()){
////			//	cout<<"*** "<<imc<<"th iterator ***"<<endl;
////			//	cout<<"sign check wrong"<<endl;
////			//	exit(EXIT_FAILURE);;
////			//}
////			if(!check_measuring_line()){
////				cout<<"*** "<<imc<<"th iterator ***"<<endl;
////				cout<<"measuring_line check wrong"<<endl;
////				exit(EXIT_FAILURE);;
////			}
////			if(!check_ref_line()){
////				cout<<"*** "<<imc<<"th iterator ***"<<endl;
////				cout<<"ref_line check wrong"<<endl;
////				exit(EXIT_FAILURE);;
////			}
////			//if(!checkconnection()){
////			//	cout<<"*** "<<imc<<"th iterator ***"<<endl;
////			//	cout<<"connection check wrong"<<endl;
////			//	exit(EXIT_FAILURE);;
////			//}
////			//if( !checkconnectivity() ){
////                        //        cout<<"connectivity check wrong!"<<endl;
////			//	cout<<"*** "<<imc<<"th iterator ***"<<endl;
////                        //        exit(EXIT_FAILURE);
////                        //}
////			if( !checkirreducibility() ){
////                                cout<<"irreducibility check wrong!"<<endl;
////				cout<<"*** "<<imc<<"th iterator ***"<<endl;
////				exit(EXIT_FAILURE);
////			}
////			if( !checkcompactness(true) ){
////                                cout<<"G-compactness check wrong!"<<endl;
////				cout<<"*** "<<imc<<"th iterator ***"<<endl;
////				exit(EXIT_FAILURE);
////			}
////			//if( checkbubble() ){
////                        //        cout<<"bubble check wrong!"<<endl;
////			//	cout<<"*** "<<imc<<"th iterator ***"<<endl;
////			//	exit(EXIT_FAILURE);
////			//}
////			//if( !checkcompactness(false) ){
////                        //        cout<<"W-compactness check wrong!"<<endl;
////			//	cout<<"*** "<<imc<<"th iterator ***"<<endl;
////			//	exit(EXIT_FAILURE);
////			//}
////			//if( !(checkcompactness(true)&&checkcompactness(false)) ){
////                        //        cout<<"irreducibility check wrong!"<<endl;
////			//	cout<<"*** "<<imc<<"th iterator ***"<<endl;
////			//	exit(EXIT_FAILURE);
////			//}
////			//for(int i=0; i<10; i++){
////			//	if(!checkdetailedbalance()){
////			//		cout<<"*** "<<imc<<"th iterator ***"<<endl;
////			//		cout<<" ** "<<i<<"th update test **"<<endl;
////			//		cout<<"detalied balance check wrong"<<endl;
////			//		printconfiguration(cout);
////			//		exit(EXIT_FAILURE);;
////			//	}
////			//}
////			//if( Norder==(this->Nordermax) && fermionicloops.size()>1 )
////			//	exit(EXIT_FAILURE);
////			for(int i=0; i<10; i++){
////				printconfiguration(cout);
////				if(!check_detailed_balance_insert_remove()){
////					cout<<"*** "<<imc<<"th iterator ***"<<endl;
////					cout<<" ** "<<i<<"th update test **"<<endl;
////					cout<<"detalied balance check for insert-remove vertices wrong"<<endl;
////					printconfiguration(cout);
////					exit(EXIT_FAILURE);;
////				}
////			}
////			for(int i=0; i<10; i++){
////				if(!check_detailed_balance_remove_insert()){
////					cout<<"*** "<<imc<<"th iterator ***"<<endl;
////					cout<<" ** "<<i<<"th update test **"<<endl;
////					cout<<"detalied balance check for remove-insert vertices wrong"<<endl;
////					printconfiguration(cout);
////					exit(EXIT_FAILURE);;
////				}
////			}
////			for(int i=0; i<10; i++){
////				if(!check_detailed_balance_shift()){
////					cout<<"*** "<<imc<<"th iterator ***"<<endl;
////					cout<<" ** "<<i<<"th update test **"<<endl;
////					cout<<"detalied balance check for shift vertex wrong"<<endl;
////					printconfiguration(cout);
////					exit(EXIT_FAILURE);;
////				}
////			}
////			//for(int i=0; i<10; i++){
////			//	if(!checkdetailedbalance_reconnect()){
////			//		cout<<"*** "<<imc<<"th iterator ***"<<endl;
////			//		cout<<" ** "<<i<<"th update test **"<<endl;
////			//		cout<<"detalied balance check for reconnect wrong"<<endl;
////			//		printconfiguration(cout);
////			//		exit(EXIT_FAILURE);;
////			//	}
////			//}
////			//for(int i=0; i<10; i++){
////			//	if(!checkdetailedbalance_interswap()){
////			//		cout<<"*** "<<imc<<"th iterator ***"<<endl;
////			//		cout<<" ** "<<i<<"th update test **"<<endl;
////			//		cout<<"detalied balance check for interswap wrong"<<endl;
////			//		printconfiguration(cout);
////			//		exit(EXIT_FAILURE);;
////			//	}
////			//}
////			//for(int i=0; i<10; i++){
////			//	if(!checkdetailedbalance_intraswap()){
////			//		cout<<"*** "<<imc<<"th iterator ***"<<endl;
////			//		cout<<" ** "<<i<<"th update test **"<<endl;
////			//		cout<<"detalied balance check for intraswap wrong"<<endl;
////			//		printconfiguration(cout);
////			//		exit(EXIT_FAILURE);;
////			//	}
////			//}
////			//for(int i=0; i<10; i++){
////			//	if(!checkdetailedbalance_loopspinflip()){
////			//		cout<<"*** "<<imc<<"th iterator ***"<<endl;
////			//		cout<<" ** "<<i<<"th update test **"<<endl;
////			//		cout<<"detalied balance check wrong"<<endl;
////			//		printconfiguration(cout);
////			//		exit(EXIT_FAILURE);;
////			//	}
////			//}
////                        //if( brutecheckirreducibility() ^ checkirreducibility() ){
////                        //        cout<<"irreducibility check wrong!"<<endl;
////			//	cout<<"*** "<<imc<<"th iterator ***"<<endl;
////                        //        exit(EXIT_FAILURE);
////                        //}
////                        //if( checkirreducibility() ){
////                        //        if( brutecheckcompactness(true) ^ checkcompactness(true) ){
////                        //                cout<<"compactness check wrong!"<<endl;
////                        //                cout<<"brute force check: "<<brutecheckcompactness(true)<<endl;
////                        //                cout<<"algorithmic check: "<<checkcompactness(true)<<endl;
////                        //                printconfiguration(cout);
////                        //                exit(EXIT_FAILURE);
////                        //        }
////			//	else if( brutecheckcompactness(false) ^ checkcompactness(false) ){
////                        //                cout<<"compactness check wrong!"<<endl;
////                        //                cout<<"brute force check: "<<brutecheckcompactness(false)<<endl;
////                        //                cout<<"algorithmic check: "<<checkcompactness(false)<<endl;
////                        //                printconfiguration(cout);
////                        //                exit(EXIT_FAILURE);
////                        //        }
////                        //}
////		}
////	}
////};
void diagram::thermalization(const int& Nth_i){
	//cout<<"* thermalization started *"<<endl<<flush;
	for(int ith=0; ith<Nth_i; ith++){
		//cout<<"*** "<<ith<<"th iterator ***"<<endl;
		//cout<<"update: ";
		switch( static_cast<int>(parameters::Nupdate*unidist(mt)) ){
			case 0:
				//cout<<"insertVertices"<<endl;
				insertVertices();
				break;
			case 1:
				//cout<<"removeVertices"<<endl;
				removeVertices();
				break;
			case 2:
				//cout<<"shiftVertex"<<endl;
				shiftVertex();
				break;
			case 3:
				//cout<<"shiftIm2"<<endl;
				shiftIm2();
				break;
			default:
				//cout<<"flipExtFlavor"<<endl;
				flipExtFlavor();
				break;
		}

		//if(ith%100==0) printconfiguration(cout);
		if( Conf.getim2()!=Conf.getGtTopSize()+1 ){
			cout<<"im2 treatment is wrong"<<endl;
			exit(EXIT_FAILURE);
		}
	}

	//cout<<"**************************"<<endl<<flush;
};
void diagram::MCsampling(const int& Nmc_i){
	int ox;
       	bool valid;

	//cout<<"thermalization"<<endl;
	//thermalization();

	// initialize quantities
	//long Sumsign = 0, Nmeasure = 0, Nvalid = 0, Nproposed = 0, Naccepted = 0;
	//long Nmeasure = 0, Nvalid = 0, Nproposed = 0, Naccepted = 0;

	normalization.resize(this->Nordermax+1);
	for(int i=0; i<normalization.size(); i++) normalization[i] = 0;

	Qacc.resize(this->Nordermax);
	for(int iorder=0; iorder<Nordermax; iorder++){
		//cout<<"order = "<<iorder<<endl;
		Qacc[iorder].initZero(NtG);
		//for(int it=0; it<Nout+1; it++){ Qacc[iorder].print(cout,it); }
	}
	//Owx.resize(Nw);
	//setSwx();

	//Qacc = new selfenergy_w**[parameters::Nordermax];
	//for(int iorder=0; iorder<parameters::Nordermax; iorder++){
	//	Swxacc[iorder] = new selfenergy_w*[Nflavor];
	//	for(int iflavor=0; iflavor<Nflavor; iflavor++)
	//		Swxacc[iorder][iflavor] = new selfenergy_w[Nflavor];
	//}

	//for(int iorder=0; iorder<Nordermax; iorder++)
	//for(int iflavor=0; iflavor<Nflavor; iflavor++) for(int jflavor=0; jflavor<Nflavor; jflavor++)
	//	Swxacc[iorder][iflavor][jflavor].init();

	//cout<<"* MCsampling started *"<<endl<<flush;
	//printconfiguration(cout);
	NorderAcc = 0.0;
	t1 = chrono::high_resolution_clock::now();
	for(int imc=0; imc<Nmc_i; imc++){
		//cout<<"*** "<<imc<<"th iterator ***"<<endl<<"update: ";
		switch( static_cast<int>(parameters::Nupdate*unidist(mt)) ){
			case 0:
				//cout<<"insertVertices"<<endl;
				ox = insertVertices();
				break;
			case 1:
				//cout<<"removeVertices"<<endl;
				ox = removeVertices();
				break;
			case 2:
				//cout<<"shiftVertex"<<endl;
				ox = shiftVertex();
				break;
			case 3:
				//cout<<"shiftIm2"<<endl;
				ox = shiftIm2();
				break;
			default:
				//cout<<"flipExtFlavor"<<endl;
				ox = flipExtFlavor();
				break;
		}
		//printconfiguration(cout);
		measure();

		// checking modules
		/*
		if( !Conf.checkGtSign() ){
			printconfiguration(cout);
			cout<<"Warning: Gt sign structure is wrong"<<endl;
			exit(EXIT_FAILURE);
		}
		if( Norder==2 && sign<0 ){
			printconfiguration(cout);
			cout<<"Warning: overall sign is wrong"<<endl;
			exit(EXIT_FAILURE);
		}
		if( Conf.getWeight()*sign<0.0 ){
			printconfiguration(cout);
			cout<<"Warning: overall sign is wrong"<<endl;
			exit(EXIT_FAILURE);
		}
		if( Conf.getiv(0)==0 && Conf.getjv(0)==Conf.getVerticesSize()-1){
			printconfiguration(cout);
			cout<<"Warning: Uncompact diagram"<<endl;
			exit(EXIT_FAILURE);
		}
		*/
		/*
		if( abs(Conf.getWeight()-1.0)>1.0e-10 ){
			printconfiguration(cout);
			cout<<"Warning: MCweight is wrong"<<endl;
			exit(EXIT_FAILURE);
		}
		if( Conf.getiW(0)!=0 ){
			printconfiguration(cout);
			cout<<"Warning: W index of vertex 0 is wrong"<<endl;
			exit(EXIT_FAILURE);
		}
		if( !checkSign() ){
			printconfiguration(cout);
			cout<<"Warning: Sign is wrong"<<endl;
			exit(EXIT_FAILURE);
		}
		if( !Conf.checkWProduct() ){
			printconfiguration(cout);
			cout<<"Warning: WProduct is not equal do product of Ws"<<endl;
			exit(EXIT_FAILURE);
		}
		if( !Conf.checkGProduct() ){
			printconfiguration(cout);
			cout<<"Warning: GProduct is not equal do product of Gs"<<endl;
			exit(EXIT_FAILURE);
		}
		if( !checkDetailedBalanceInsertRemove() ){
			printconfiguration(cout);
			cout<<"Warning: Detailed balance problem between insert and remove"<<endl;
			exit(EXIT_FAILURE);
		}
		if( !checkDetailedBalanceRemoveInsert() ){
			printconfiguration(cout);
			cout<<"Warning: Detailed balance problem between remove and insert"<<endl;
			exit(EXIT_FAILURE);
		}
		if( !checkDetailedBalanceShiftIm2() ){
			printconfiguration(cout);
			cout<<"Warning: Detailed balance problem between shiftIm2 updates"<<endl;
			exit(EXIT_FAILURE);
		}
		if( !checkDetailedBalanceShiftVertex() ){
			printconfiguration(cout);
			cout<<"Warning: Detailed balance problem between shiftVertex updates"<<endl;
			exit(EXIT_FAILURE);
		}
		if( !checkDetailedBalanceFlipExtFlavor() ){
			printconfiguration(cout);
			cout<<"Warning: Detailed balance problem between flipExtFlavor updates"<<endl;
			exit(EXIT_FAILURE);
		}
		*/

		/*if(sign<0){
			cout<<imc<<"th MC step"<<endl;
			printconfiguration(cout);
			exit(EXIT_FAILURE);
		}*/

		//if(ox>1){
		//	Nproposed++; Naccepted++;
		//	if(Norder>1 && ox==3) setSwx();
		//}
		//else if(ox==1) Nproposed++;

		//if(imc%Nmeasureinterval==0){
		//	Nmeasure++;
		//	valid = (checkirreducibility()&&checkcompactness(true))&&checkcompactness(false);
		//	if(valid){
		//		Nvalid++;
		//		//Sumsign += sign;
		//		measure_S();
		//	}
		//}
	}
	//cout<<endl<<"normalization"<<endl;
	//for(auto& iterRef: normalization) cout<<iterRef<<" ";
	//cout<<endl;

	t2 = chrono::high_resolution_clock::now();
	time_span = chrono::duration_cast<chrono::duration<double>>(t2 - t1);

	//Q = Qacc;
	//printconfiguration(cout);

	normalize();

	//Qacc[1].print(cout,10,1,3);
	//Qacc[1].print(cout,38,1,3);

	//NorderAv = NorderAcc / Nmc;
	//Avesign = static_cast<double>(Sumsign)/Nvalid;
	//Rvalid = static_cast<double>(Nvalid)/Nmeasure;
	//Rproposal = static_cast<double>(Nproposed)/Nmc;
	//Racceptance = static_cast<double>(Naccepted)/Nproposed;

	//int rank = MPI::COMM_WORLD.Get_rank();
	//stringstream ss;
	//ss<<parameters::global.getdir()<<"/MClog"<<rank<<".dat";
	//ofstream ofstr(ss.str().c_str());
	//printlog(ofstr);
	//ofstr.close();

	//cout<<"* MCsampling finished *"<<endl<<flush;
	//cout<<"**********************"<<endl<<flush;
};
//
///*void diagram::diagMC_S(){
//	int ox;
//       	bool valid;
//
//	//cout<<"thermalization"<<endl;
//	thermalization();
//
//	// initialize quantities
//	//long Sumsign = 0, Nmeasure = 0, Nvalid = 0, Nproposed = 0, Naccepted = 0;
//	long Nmeasure = 0, Nvalid = 0, Nproposed = 0, Naccepted = 0;
//
//	normalization[0].resize(this->Nordermax+1);
//	for(int i=0; i<normalization[0].size(); i++) normalization[0][i] = 0;
//
//	Owx.resize(Nw);
//	setSwx();
//
//	for(int iorder=0; iorder<Nordermax; iorder++)
//	for(int iflavor=0; iflavor<Nflavor; iflavor++) for(int jflavor=0; jflavor<Nflavor; jflavor++)
//		Swxacc[iorder][iflavor][jflavor].init();
//
//	//cout<<"* MCsampling started *"<<endl<<flush;
//	//printconfiguration(cout);
//	for(int imc=0; imc<Nmc; imc++){
//		//cout<<"*** "<<imc<<"th iterator ***"<<endl<<"update: ";
//		switch( static_cast<int>(Nupdate*unidist(mt)) ){
//			case 0:
//				//cout<<"addvertice"<<endl;
//				ox = addvertices();
//				break;
//			case 1:
//				//cout<<"removevertice"<<endl;
//				ox = removevertices();
//				break;
//			case 2:
//				//cout<<"reconnect"<<endl;
//				ox = reconnect();
//				break;
//			case 3:
//				//cout<<"intraswap"<<endl;
//				ox = intraswap();
//				break;
//			case 4:
//				//cout<<"absorbbubble"<<endl;
//				ox = absorbbubble();
//				break;
//			default:
//				//cout<<"ejectbubble"<<endl;
//				ox = ejectbubble();
//				break;
//		}
//		//printconfiguration(cout);
//		if(ox>1){
//			Nproposed++; Naccepted++;
//			if(Norder>1 && ox==3) setSwx();
//		}
//		else if(ox==1) Nproposed++;
//
//		if(imc%Nmeasureinterval==0){
//			Nmeasure++;
//			valid = (checkirreducibility()&&checkcompactness(true))&&checkcompactness(false);
//			if(valid){
//				Nvalid++;
//				//Sumsign += sign;
//				measure_S();
//			}
//		}
//	}
//	//printconfiguration(cout);
//	//Avesign = static_cast<double>(Sumsign)/Nvalid;
//	Rvalid = static_cast<double>(Nvalid)/Nmeasure;
//	Rproposal = static_cast<double>(Nproposed)/Nmc;
//	Racceptance = static_cast<double>(Naccepted)/Nproposed;
//	//cout<<"* MCsampling finished *"<<endl<<flush;
//};*/
///*void diagram::diagMC_S(int iorder){
//	//cout<<endl<<"* diagMC module *"<<endl;
//	int ox;
//	bool valid;
//
//	//cout<<"thermalization"<<endl;
//	thermalization();
//
//	// initialize quantities
//	//long Sumsign = 0, Nmeasure = 0, Nvalid = 0, Nproposed = 0, Naccepted = 0;
//	long Nmeasure = 0, Nvalid = 0, Nproposed = 0, Naccepted = 0;
//
//	normalization[0].resize(this->Nordermax+1);
//	for(int i=0; i<normalization[0].size(); i++) normalization[0][i] = 0;
//
//	Owx.resize(Nw);
//	setSwx();
//
//	for(int iflavor=0; iflavor<Nflavor; iflavor++) for(int jflavor=0; jflavor<Nflavor; jflavor++)
//		Swxacc[iorder-1][iflavor][jflavor].init();
//
//	//cout<<"* MCsampling started *"<<endl<<flush;
//	//printconfiguration(cout);
//	for(long long imc=0; imc<Nmc; imc++){
//		//cout<<"*** "<<imc<<"th iterator ***"<<endl<<"update: ";
//		switch( static_cast<int>(Nupdate*unidist(mt)) ){
//			case 0:
//				//cout<<"addvertice"<<endl;
//				ox = addvertices();
//				break;
//			case 1:
//				//cout<<"removevertice"<<endl;
//				ox = removevertices();
//				break;
//			case 2:
//				//cout<<"reconnect"<<endl;
//				ox = reconnect();
//				break;
//			case 3:
//				//cout<<"intraswap"<<endl;
//				ox = intraswap();
//				break;
//			case 4:
//				//cout<<"absorbbubble"<<endl;
//				ox = absorbbubble();
//				break;
//			default:
//				//cout<<"ejectbubble"<<endl;
//				ox = ejectbubble();
//				break;
//		}
//		//printconfiguration(cout);
//
//		if(ox>1){
//			Nproposed++; Naccepted++;
//			if(Norder==iorder && ox==3) setSwx();
//		}
//		else if(ox==1) Nproposed++;
//
//		if(imc%Nmeasureinterval==0){
//			Nmeasure++;
//			valid = (checkirreducibility()&&checkcompactness(true))&&checkcompactness(false);
//			if(valid){
//				Nvalid++;
//				//Sumsign += sign;
//				measure_S(iorder);
//			}
//		}
//	}
//	//printconfiguration(cout);
//
//	//Avesign = static_cast<double>(Sumsign)/Nvalid;
//	Rvalid = static_cast<double>(Nvalid)/Nmeasure;
//	Rproposal = static_cast<double>(Nproposed)/Nmc;
//	Racceptance = static_cast<double>(Naccepted)/Nproposed;
//};*/
///*void diagram::diagMC_P(){
//	//cout<<endl<<"* diagMC_P module *"<<endl;
//	int ox;
//	bool valid;
//
//	//cout<<"thermalization"<<endl;
//	//thermalization();
//
//	// initialize quantities
//	//long Sumsign = 0, Nmeasure = 0, Nvalid = 0, Nproposed = 0, Naccepted = 0;
//	long Nmeasure = 0, Nvalid = 0, Nproposed = 0, Naccepted = 0;
//
//	normalization[1].resize(this->Nordermax+1);
//	for(int i=0; i<normalization[1].size(); i++) normalization[1][i] = 0;
//
//
//	Owx.resize(Nw);
//	setPwx();
//
//	for(int iorder=0; iorder<Nordermax+1; iorder++)
//	for(int iflavor=0; iflavor<Nflavor; iflavor++) for(int jflavor=0; jflavor<Nflavor; jflavor++)
//		Pwxacc[iorder][iflavor][jflavor].init();
//
//	//cout<<"* MCsampling started *"<<endl<<flush;
//	//printconfiguration(cout);
//	for(long long imc=0; imc<Nmc; imc++){
//		//cout<<"*** "<<imc<<"th iterator ***"<<endl<<"update: ";
//		switch( static_cast<int>(Nupdate*unidist(mt)) ){
//		//switch( static_cast<int>(4*unidist(mt)) ){
//			case 0:
//				//cout<<"addvertice"<<endl;
//				ox = addvertices();
//				break;
//			case 1:
//				//cout<<"removevertice"<<endl;
//				ox = removevertices();
//				break;
//			case 2:
//				//cout<<"intraswap"<<endl;
//				ox = intraswap();
//				break;
//			case 3:
//				//cout<<"reconnect"<<endl;
//				ox = reconnect();
//				break;
//			case 4:
//				//cout<<"absorbbubble"<<endl;
//				ox = absorbbubble();
//				break;
//			default:
//				//cout<<"ejectbubble"<<endl;
//				ox = ejectbubble();
//				break;
//		}
//
//		//if(Norder==0) cout<<imc<<" ";
//		//if(imc==211) printconfiguration(cout);
//		//else if(imc==10000) printconfiguration(cout);
//
//		if(ox>1){
//			//printconfiguration(cout);
//			Nproposed++; Naccepted++;
//			//if(Norder==iorder && ox==3) setPwx();
//			setPwx();
//		}
//		else if(ox==1) Nproposed++;
//		if(imc%Nmeasureinterval==0){
//			Nmeasure++;
//			valid = (checkirreducibility(true)&&checkirreducibility(false))&&(checkcompactness(true)&&checkcompactness(false));
//			if(valid){
//				Nvalid++;
//				//Sumsign += sign;
//				measure_P();
//			}
//		}
//	}
//	//printconfiguration(cout);
//	//Avesign = static_cast<double>(Sumsign)/Nvalid;
//	Rvalid = static_cast<double>(Nvalid)/Nmeasure;
//	Rproposal = static_cast<double>(Nproposed)/Nmc;
//	Racceptance = static_cast<double>(Naccepted)/Nproposed;
//};*/
///*void diagram::diagMC_P(int iorder){
//	//cout<<endl<<"* diagMC_P module *"<<endl;
//	int ox;
//	bool valid;
//
//	//cout<<"thermalization"<<endl;
//	//thermalization();
//
//	// initialize quantities
//	long Nmeasure = 0, Nvalid = 0, Nproposed = 0, Naccepted = 0;
//
//	normalization[1].resize(this->Nordermax+1);
//	for(int i=0; i<normalization[1].size(); i++) normalization[1][i] = 0;
//
//	Owx.resize(Nw);
//	setPwx();
//
//	for(int iflavor=0; iflavor<Nflavor; iflavor++) for(int jflavor=0; jflavor<Nflavor; jflavor++)
//		Pwxacc[iorder][iflavor][jflavor].init();
//
//	//cout<<"* MCsampling started *"<<endl<<flush;
//	//printconfiguration(cout);
//	int count1 = 0, count2 = 0;
//	for(long long imc=0; imc<Nmc; imc++){
//		//cout<<"*** "<<imc<<"th iterator ***"<<endl<<"update: ";
//		switch( static_cast<int>(Nupdate*unidist(mt)) ){
//			case 0:
//				//cout<<"addvertice"<<endl;
//				ox = addvertices();
//				break;
//			case 1:
//				//cout<<"removevertice"<<endl;
//				ox = removevertices();
//				break;
//			case 2:
//				//cout<<"reconnect"<<endl;
//				ox = reconnect();
//				break;
//			case 3:
//				//cout<<"intraswap"<<endl;
//				ox = intraswap();
//				break;
//			case 4:
//				//cout<<"absorbbubble"<<endl;
//				ox = absorbbubble();
//				break;
//			default:
//				//cout<<"ejectbubble"<<endl;
//				ox = ejectbubble();
//				break;
//		}
//		if((checkirreducibility(true)&&checkirreducibility(false)) && (checkcompactness(true)&&checkcompactness(false))){ 
//			if(Norder==0) count1++;
//			else if(Norder==2 && Ws.size()==3 && measuringline->getiflavor()==measuringline->getjflavor()) count2++;
//		}
//		if(ox>1){
//			Nproposed++; Naccepted++;
//			if(Norder==iorder) setPwx();
//		}
//		else if(ox==1) Nproposed++;
//		if(imc%Nmeasureinterval==0){
//			Nmeasure++;
//			valid = (checkirreducibility(true)&&checkirreducibility(false))&&(checkcompactness(true)&&checkcompactness(false));
//			if(valid){
//				Nvalid++;
//				//Sumsign += sign;
//				measure_P(iorder);
//			}
//		}
//	}
//	//printconfiguration(cout);
//	//Avesign = static_cast<double>(Sumsign)/Nvalid;
//	Rvalid = static_cast<double>(Nvalid)/Nmeasure;
//	Rproposal = static_cast<double>(Nproposed)/Nmc;
//	Racceptance = static_cast<double>(Naccepted)/Nproposed;
//};*/
///*void diagram::diagMC(){
//	//cout<<endl<<"* diagMC for both S and P module *"<<endl;
//	int ox;
//	bool valid;
//
//	//cout<<"thermalization"<<endl;
//	//thermalization();
//
//	// initialize quantities
//	long Nmeasure = 0, Nvalid = 0, Nproposed = 0, Naccepted = 0;
//
//	for(int sp=0; sp<2; sp++){
//		normalization[sp].resize(this->Nordermax+1);
//		for(int i=0; i<normalization[sp].size(); i++) normalization[sp][i] = 0.;
//	}
//
//	Owx.resize(Nw);
//	if(measuringline->gettype()==0) setSwx();
//	else setPwx();
//
//	for(int iflavor=0; iflavor<Nflavor; iflavor++) for(int jflavor=0; jflavor<Nflavor; jflavor++){
//		for(int iorder=0; iorder<Nordermax; iorder++) Swxacc[iorder][iflavor][jflavor].init();
//		for(int iorder=0; iorder<Nordermax+1; iorder++) Pwxacc[iorder][iflavor][jflavor].init();
//	}
//
//	//cout<<"* MCsampling started *"<<endl<<flush;
//	//printconfiguration(cout);
//	for(long long imc=0; imc<Nmc; imc++){
//		//cout<<"*** "<<imc<<"th iterator ***"<<endl<<"update: ";
//		switch( static_cast<int>(Nupdate*unidist(mt)) ){
//			case 0:
//				//cout<<"addvertice"<<endl;
//				ox = addvertices();
//				break;
//			case 1:
//				//cout<<"removevertice"<<endl;
//				ox = removevertices();
//				break;
//			case 2:
//				//cout<<"interswap"<<endl;
//				ox = interswap();
//				//ox = intraswap();
//				break;
//			case 3:
//				//cout<<"reconnect"<<endl;
//				ox = reconnect();
//				break;
//			case 4:
//				//cout<<"absorbbubble"<<endl;
//				ox = absorbbubble();
//				break;
//			default:
//				//cout<<"ejectbubble"<<endl;
//				ox = ejectbubble();
//				break;
//		}
//
//		//printconfiguration(cout);
//
//		if(ox>1){
//			Nproposed++; Naccepted++;
//			if(measuringline->gettype()==0) setSwx();
//			else setPwx();
//		}
//		else if(ox==1) Nproposed++;
//
//		if(imc%Nmeasureinterval==0){
//			Nmeasure++;
//			valid = (checkirreducibility(true)&&checkirreducibility(false))&&(checkcompactness(true)&&checkcompactness(false));
//			if(valid){
//				Nvalid++;
//				//Sumsign += sign;
//				measure();
//			}
//		}
//	}
//	//printconfiguration(cout);
//
//	//Avesign = static_cast<double>(Sumsign)/Nvalid;
//	Rvalid = static_cast<double>(Nvalid)/Nmeasure;
//	Rproposal = static_cast<double>(Nproposed)/Nmc;
//	Racceptance = static_cast<double>(Naccepted)/Nproposed;
//
//	//cout<<"Avesign = "<<Avesign<<endl;
//	//cout<<"Rvalid = "<<Rvalid<<endl;
//	//cout<<"Rproposal = "<<Rproposal<<endl;
//	//cout<<"Racceptance = "<<Racceptance<<endl;
//};*/
///*void diagram::diagMC(int iorder){
//	//cout<<endl<<"* diagMC for both S and P module *"<<endl;
//	int ox;
//	bool valid;
//
//	//cout<<"thermalization"<<endl;
//	//thermalization();
//
//	// initialize quantities
//	long Nmeasure = 0, Nvalid = 0, Nproposed = 0, Naccepted = 0;
//
//	for(int sp=0; sp<2; sp++){
//		normalization[sp].resize(this->Nordermax+1);
//		for(int i=0; i<normalization[sp].size(); i++) normalization[sp][i] = 0.;
//	}
//
//	Owx.resize(Nw);
//	if(measuringline->gettype()==0) setSwx();
//	else setPwx();
//
//	for(int iflavor=0; iflavor<Nflavor; iflavor++) for(int jflavor=0; jflavor<Nflavor; jflavor++){
//		Swxacc[iorder-1][iflavor][jflavor].init();
//		Pwxacc[iorder][iflavor][jflavor].init();
//	}
//
//	//cout<<"* MCsampling started *"<<endl<<flush;
//	//printconfiguration(cout);
//	for(long long imc=0; imc<Nmc; imc++){
//		cout<<"*** "<<imc<<"th iterator ***"<<endl<<"update: ";
//		switch( static_cast<int>(Nupdate*unidist(mt)) ){
//			case 0:
//				cout<<"addvertice"<<endl;
//				ox = addvertices();
//				break;
//			case 1:
//				cout<<"removevertice"<<endl;
//				ox = removevertices();
//				break;
//			case 2:
//				cout<<"reconnect"<<endl;
//				ox = reconnect();
//				break;
//			case 3:
//				cout<<"interswap"<<endl;
//				ox = interswap();
//				break;
//			case 4:
//				cout<<"absorbbubble"<<endl;
//				ox = absorbbubble();
//				break;
//			default:
//				cout<<"ejectbubble"<<endl;
//				ox = ejectbubble();
//				break;
//		}
//
//		if(ox>1){
//			Nproposed++; Naccepted++;
//			//if(Norder==iorder && ox==3) setPwx();
//			if(Norder==iorder){
//				if(measuringline->gettype()==0) setSwx();
//				else setPwx();
//			}
//		}
//		else if(ox==1) Nproposed++;
//
//		if(imc%Nmeasureinterval==0){
//			Nmeasure++;
//			valid = (checkirreducibility()&&checkcompactness(true))&&checkcompactness(false);
//			if(valid){
//				Nvalid++;
//				//Sumsign += sign;
//				if(measuringline->gettype()==0) measure_S(iorder);
//				else measure_P(iorder);
//			}
//		}
//	}
//	//printconfiguration(cout);
//	//Avesign = static_cast<double>(Sumsign)/Nvalid;
//	Rvalid = static_cast<double>(Nvalid)/Nmeasure;
//	Rproposal = static_cast<double>(Nproposed)/Nmc;
//	Racceptance = static_cast<double>(Naccepted)/Nproposed;
//}*/
///*void diagram::bold_diagMC_G(){
//	//cout<<endl<<"* diagMC for S module *"<<endl;
//	int ox;
//	bool valid;
//
//	//cout<<"thermalization"<<endl;
//	//thermalization();
//
//	// initialize quantities
//	long Nmeasure = 0, Nvalid = 0, Nproposed = 0, Naccepted = 0;
//
//	normalization[0].resize(this->Nordermax+1);
//	for(int i=0; i<normalization[0].size(); i++) normalization[0][i] = 0.;
//
//	Owx.resize(Nw);
//	setSwx();
//
//	for(int iflavor=0; iflavor<Nflavor; iflavor++) for(int jflavor=0; jflavor<Nflavor; jflavor++){
//		for(int iorder=0; iorder<Nordermax; iorder++) Swxacc[iorder][iflavor][jflavor].init();
//	}
//
//	//cout<<"* MCsampling started *"<<endl<<flush;
//	//printconfiguration(cout);
//	for(long long imc=0; imc<Nmc; imc++){
//		//cout<<"*** "<<imc<<"th iterator ***"<<endl<<"update: ";
//		switch( static_cast<int>(Nupdate*unidist(mt)) ){
//			case 0:
//				//cout<<"addvertice"<<endl;
//				ox = addvertices();
//				break;
//			case 1:
//				//cout<<"removevertice"<<endl;
//				ox = removevertices();
//				break;
//			case 2:
//				//cout<<"interswap"<<endl;
//				ox = intraswap();
//				break;
//			default:
//				//cout<<"reconnect"<<endl;
//				ox = reconnect();
//				break;
//		}
//
//		if(ox>1){
//			Nproposed++; Naccepted++;
//			if(measuringline->gettype()==0) setSwx();
//			else setPwx();
//		}
//		else if(ox==1) Nproposed++;
//
//		if(imc%Nmeasureinterval==0){
//			Nmeasure++;
//			valid = (checkirreducibility(true)&&checkcompactness(true));
//			if(valid){
//				Nvalid++;
//				//Sumsign += sign;
//				measure();
//			}
//		}
//	}
//
//	//Avesign = static_cast<double>(Sumsign)/Nvalid;
//	Rvalid = static_cast<double>(Nvalid)/Nmeasure;
//	Rproposal = static_cast<double>(Nproposed)/Nmc;
//	Racceptance = static_cast<double>(Naccepted)/Nproposed;
//};*/
///*void diagram::bold_diagMC_GW(){
//	// cout<<endl<<"* bold diagMC for both S and P module *"<<endl;
//	int ox, updatetype;
//	bool valid, Gcompact, Wcompact;
//
//	// cout<<"thermalization"<<endl;
//	// thermalization();
//
//	// initialize quantities
//	long Nmeasure = 0, Nvalid = 0, Nproposed = 0, Naccepted = 0;
//	long Nirreducible = 0, NGcompact = 0, NWcompact = 0;
//	for(int sp=0; sp<2; sp++){
//		normalization[sp].resize(this->Nordermax+1);
//		for(int i=0; i<normalization[sp].size(); i++) normalization[sp][i] = 0;
//	}
//
//	Plocacc.resize(this->Nordermax);
//	for(int i=0; i<Plocacc.size(); i++){ Plocacc[i] = 0.; }
//
//	Owx.resize(Nw);
//	if(measuringline->gettype()==0) setSwx();
//	else if(measuringline->gettype()==2) setPwx();
//
//	for(int iorder=0; iorder<Nordermax; iorder++)
//	for(int iflavor=0; iflavor<Nflavor; iflavor++) for(int jflavor=0; jflavor<Nflavor; jflavor++){
//		Swxacc[iorder][iflavor][jflavor].init();
//		Pwxacc[iorder][iflavor][jflavor].init();
//	}
//
//	//cout<<"* MCsampling starts *"<<endl<<flush;
//	//printconfiguration(cout);
//	for(long long imc=0; imc<Nmc; imc++){
//		//cout<<"*** "<<imc<<"th iterator ***"<<endl<<"update: ";
//		//updatetype = static_cast<int>(Nupdate*unidist(mt));
//		updatetype = static_cast<int>((2*Nupdate-1)*unidist(mt));
//		switch( updatetype ){
//			case 0:
//				// cout<<"addvertice"<<endl<<flush;
//				ox = addvertices();
//				break;
//			case 1:
//				// cout<<"removevertice"<<endl<<flush;
//				ox = removevertices();
//				break;
//			case 2:
//				// cout<<"reconnect"<<endl<<flush;
//				ox = reconnect();
//				break;
//			case 3:
//				// cout<<"interswap"<<endl<<flush;
//				ox = interswap();
//				break;
//			default:
//				// cout<<"transformUW"<<endl<<flush;
//				ox = transformUW();
//				break;
//		};
// 
//		if(ox>1){
//			Nproposed++; Naccepted++;
//			if(measuringline->gettype()==0) setSwx();
//			else if(measuringline->gettype()==2) setPwx();
//		}
//		else if(ox==1) Nproposed++;
//
//		if(imc%Nmeasureinterval==0){
//			Nmeasure++; Nirreducible++; NGcompact++; NWcompact++; Nvalid++;
//			measure();
//		}
//	}
//	//cout<<"* MCsampling ends *"<<endl<<flush;
//
//	//Avesign = static_cast<double>(Sumsign)/Nvalid;
//	Rproposal = static_cast<double>(Nproposed)/Nmc;
//	Racceptance = static_cast<double>(Naccepted)/Nproposed;
//	Rvalid = static_cast<double>(Nvalid)/Nmeasure;
//	Rirreducible = static_cast<double>(Nirreducible)/Nmeasure;
//	RGcompact = static_cast<double>(NGcompact)/Nmeasure;
//	RWcompact = static_cast<double>(NWcompact)/Nmeasure;
//};
//void diagram::bold_diagMC_GW_S(){
//	// cout<<endl<<"* bold diagMC for both S and P module *"<<endl;
//	int ox, updatetype;
//	bool valid, Gcompact, Wcompact;
//
//	// cout<<"thermalization"<<endl;
//	//thermalization();
//
//	// initialize quantities
//	long Nmeasure = 0, Nvalid = 0, Nproposed = 0, Naccepted = 0;
//	long Nirreducible = 0, NGcompact = 0, NWcompact = 0;
//	for(int sp=0; sp<2; sp++){
//		normalization[sp].resize(this->Nordermax+1);
//		for(int i=0; i<normalization[sp].size(); i++) normalization[sp][i] = 0;
//	}
//
//	Plocacc.resize(this->Nordermax);
//	for(int i=0; i<Plocacc.size(); i++){ Plocacc[i] = 0.; }
//
//	Owx.resize(Nw);
//	if(measuringline->gettype()==0) setSwx();
//	else if(measuringline->gettype()==2) setPwx();
//
//	for(int iflavor=0; iflavor<Nflavor; iflavor++) for(int jflavor=0; jflavor<Nflavor; jflavor++){
//		for(int iorder=0; iorder<Nordermax; iorder++) Swxacc[iorder][iflavor][jflavor].init();
//		for(int iorder=0; iorder<Nordermax; iorder++) Pwxacc[iorder][iflavor][jflavor].init();
//	}
//
//	 //cout<<"* MCsampling started *"<<endl<<flush;
//	 //printconfiguration(cout);
//	for(long long imc=0; imc<Nmc; imc++){
//		//cout<<"*** "<<imc<<"th iterator ***"<<endl<<"update: ";
//		//updatetype = static_cast<int>(Nupdate*unidist(mt));
//		updatetype = static_cast<int>((2*Nupdate-1)*unidist(mt));
//		switch( updatetype ){
//			case 0:
//				// cout<<"addvertice"<<endl;
//				ox = addvertices();
//				break;
//			case 1:
//				// cout<<"removevertice"<<endl;
//				ox = removevertices();
//				break;
//			case 2:
//				// cout<<"reconnect"<<endl;
//				ox = reconnect();
//				break;
//			case 3:
//				// cout<<"interswap"<<endl;
//				ox = intraswap();
//				break;
//			default:
//				// cout<<"interswap"<<endl;
//				ox = transformUW();
//				break;
//		};
// 
//		if(ox>1){
//			Nproposed++; Naccepted++;
//			if(measuringline->gettype()==0) setSwx();
//			else if(measuringline->gettype()==2) setPwx();
//		}
//		else if(ox==1) Nproposed++;
//
//		if(imc%Nmeasureinterval==0){
//			Nmeasure++; Nirreducible++; NGcompact++; NWcompact++; Nvalid++;
//			measure();
//		}
//	}
//
//	//Avesign = static_cast<double>(Sumsign)/Nvalid;
//	Rproposal = static_cast<double>(Nproposed)/Nmc;
//	Racceptance = static_cast<double>(Naccepted)/Nproposed;
//	Rvalid = static_cast<double>(Nvalid)/Nmeasure;
//	Rirreducible = static_cast<double>(Nirreducible)/Nmeasure;
//	RGcompact = static_cast<double>(NGcompact)/Nmeasure;
//	RWcompact = static_cast<double>(NWcompact)/Nmeasure;
//};*/
//
void diagram::normalize(){
	//cout<<"before normalization: sampled factors"<<endl;
	//for(int i=0; i<normalization.size(); i++)
	//	cout<<normalization[i]<<" ";
	//cout<<endl;

	for(int i=1; i<normalization.size(); i++)
		normalization[i] /= normalization[0];

	for(int i=0; i<Nordermax; i++)
		Qacc[i] /= normalization[0];

	normalization[0] = 1.0;

	//cout<<"after normalization: sampled factors"<<endl;
	//for(int i=0; i<normalization.size(); i++)
	//	cout<<normalization[i]<<" ";
	//cout<<endl;
	//cout<<"*****************"<<endl;
};
void diagram::partialSum(){
	Qsum.initZero(NtG);
	for(int i=1; i<Nordermax; i++)
		Qsum += Qacc[i];
};

///*void diagram::normalize_S(const int& reftype, const int& reforder, const double& norm_i){
//	if(normalization[reftype][reforder]==0){
//		normalization[reftype][reforder]++;
//		cout<<endl<<"WARNING! ZERO COUNT ON "<<reforder<<"TH ORDER DIAGRAM, ARTIFICIAL COUNT IS ADDED!"<<endl<<flush;
//	}
//	double ratio = static_cast<double>(normalization[reftype][reforder])*Ri_N[reforder]/norm_i;
//
//	for(int iorder=0; iorder<(this->Nordermax); iorder++)
//	for(int iflavor=0; iflavor<Nflavor; iflavor++) for(int jflavor=0; jflavor<Nflavor; jflavor++){
//		Swxacc[iorder][iflavor][jflavor].fourier_xtok();
//		Swxacc[iorder][iflavor][jflavor] /= ratio;
//	}
//	//cout<<"*****************"<<endl;
//};
//void diagram::normalize_S(const int& reftype, const int& reforder, const double& norm_i, const int& iorder){
//	if(normalization[reftype][reforder]==0){
//		normalization[reftype][reforder]++;
//		cout<<endl<<"WARNING! ZERO COUNT ON "<<reforder<<"TH ORDER DIAGRAM, ARTIFICIAL COUNT IS ADDED!"<<endl<<flush;
//	}
//	double ratio = static_cast<double>(normalization[reftype][reforder])*Ri_N[reforder]/norm_i;
//	for(int iflavor=0; iflavor<Nflavor; iflavor++) for(int jflavor=0; jflavor<Nflavor; jflavor++){
//		Swxacc[iorder-1][iflavor][jflavor].fourier_xtok();
//		Swxacc[iorder-1][iflavor][jflavor] /= ratio;
//	}
//	//cout<<"*****************"<<endl;
//};
//void diagram::normalize_P(const int& reftype, const int& reforder, const double& norm_i){
//	if(normalization[reftype][reforder]==0){
//		normalization[reftype][reforder]++;
//		cout<<endl<<"WARNING! ZERO COUNT ON "<<reforder<<"TH ORDER DIAGRAM, ARTIFICIAL COUNT IS ADDED!"<<endl<<flush;
//	}
//	double ratio = static_cast<double>(normalization[reftype][reforder]*Ri_N[reforder])/norm_i;
//
//	for(int iorder=0; iorder<(this->Nordermax); iorder++){
//		Plocacc[iorder] /= ratio;
//		for(int iflavor=0; iflavor<Nflavor; iflavor++) for(int jflavor=0; jflavor<Nflavor; jflavor++){
//			Pwxacc[iorder][iflavor][jflavor].fourier_xtok();
//			Pwxacc[iorder][iflavor][jflavor] /= ratio;
//		}
//	}
//	//cout<<"*****************"<<endl;
//};
//void diagram::normalize_P(const int& reftype, const int& reforder, const double& norm_i, const int& iorder){
//	if(normalization[reftype][reforder]==0){
//		normalization[reftype][reforder]++;
//		cout<<endl<<"WARNING! ZERO COUNT ON "<<reforder<<"TH ORDER DIAGRAM, ARTIFICIAL COUNT IS ADDED!"<<endl<<flush;
//	}
//	double ratio = static_cast<double>(normalization[reftype][reforder])*Ri_N[reforder]/norm_i;
//
//	Plocacc[iorder] /= ratio;
//	for(int iflavor=0; iflavor<Nflavor; iflavor++) for(int jflavor=0; jflavor<Nflavor; jflavor++){
//		Pwxacc[iorder][iflavor][jflavor].fourier_xtok();
//		Pwxacc[iorder][iflavor][jflavor] /= ratio;
//	}
//	//cout<<"*****************"<<endl;
//};
//
////double diagram::getAvesign(){ return Avesign; };
//double diagram::getRproposal(){ return Rproposal; };
//double diagram::getRacceptance(){ return Racceptance; };
//double diagram::getRvalid(){ return Rvalid; };
//double diagram::getRirreducible(){ return Rirreducible; };
//double diagram::getRGcompact(){ return RGcompact; };
//double diagram::getRWcompact(){ return RWcompact; };
//long diagram::getnorm(const int& type_index, const int& order_index){ return normalization[type_index][order_index]; };
//double diagram::getPloc(const int& order_index){ return Plocacc[order_index]; };
//double diagram::getPlocSum(const int& order_index){
//	double sum = Plocacc[0];
//	for(int i=1; i<order_index+1; i++) sum += Plocacc[i];
//	return sum;
//};
//selfenergy_w diagram::getSw(const int& order_index, const int& iflavor_i, const int& jflavor_i){ return Swxacc[order_index][iflavor_i][jflavor_i]; };
//selfenergy_w diagram::getSwSum(const int& order_index, const int& iflavor_i, const int& jflavor_i){ 
//	selfenergy_w tmp = Swxacc[0][iflavor_i][jflavor_i];
//	for(int i=1; i<order_index; i++) tmp += Swxacc[i][iflavor_i][jflavor_i];
//	return tmp;
//};
//polarization_w diagram::getPw(const int& order_index, const int& iflavor_i, const int& jflavor_i){ return Pwxacc[order_index][iflavor_i][jflavor_i]; };
//polarization_w diagram::getPwSum(const int& order_index, const int& iflavor_i, const int& jflavor_i){ 
//	polarization_w tmp = Pwxacc[0][iflavor_i][jflavor_i];
//	for(int i=1; i<order_index+1; i++) tmp += Pwxacc[i][iflavor_i][jflavor_i];
//	return tmp;
//};*/

double diagram::trapezoidalWeight(const int& i, const int& N){
    if( (i==0) || (i==N) ) return 0.5;
    else return 1.0;
}
vector<double> diagram::getQseries(const int& iout_i, const int& iin_i){
	//cout<<"* getQseries *"<<endl<<flush;
	vector<double> Qseries_tmp(Nordermax+1,0.0);
    double wt, wt2, wt1;
	for(int it=0; it<NtG; it++){
        wt = trapezoidalWeight(it,NtG-1);
        for(int it2=0; it2<it+1; it2++){
            wt2 = trapezoidalWeight(it2,it);
            for(int it1=0; it1<it2+1; it1++){
                wt1 = trapezoidalWeight(it1,it2);
                for(int iorder=2; iorder<Nordermax+1; iorder++)
                    Qseries_tmp[iorder] += wt*wt2*wt1*Qacc[iorder-1].getQt(it,it2,it1,iout_i,iin_i);
                    //Qseries_tmp[iorder] += Qacc[iorder-1].getQt(it,it2,it1,iout_i,iin_i);
            }
        }
	}
	//for(int iorder=2; iorder<Nordermax+1; iorder++) Qseries_tmp[iorder] /= static_cast<double>(NtG*(NtG+1)*(NtG+2)/6);
	for(int iorder=2; iorder<Nordermax+1; iorder++) Qseries_tmp[iorder] *= pow(parameters::beta/(NtG-1),3);

	return Qseries_tmp;
};
vertex_4pt diagram::getQt(){
	return Qsum;
};
////
////Eigen::MatrixXd diagram::getQt(const int& n){
////	return Qsum.getQt(n);
////};
////
////void diagram::printlog(ostream& ostr){
////	//printconfiguration(ostr);
////	ostr<<endl<<"-----------------------------------------------------------" << endl;
////	ostr<<"# Norder average = "<<NorderAv<<endl;
////	ostr<<"# of updates / second = "<<static_cast<double>(Nmc)/time_span.count()<<endl;
////};
void diagram::printconfiguration(ostream& ostr){
	ostr<<endl;
	ostr<<"* diagram information *"<<endl;
	ostr<<"Norder : "<<Norder<<endl;
	//ostr<<"sign : "<<sign<<endl;
	////ostr<<"Nselfloop : "<<Nselfloop<<endl;
	//ostr<<"vertex list : "<<endl;
	//for(auto& it:vertices) ostr<<it<<endl;
	////ostr<<"fermionicloop list : "<<endl;
	////for(auto& it:fermionicloops) ostr<<it<<endl;
	//ostr<<endl;
	//ostr<<"W-line list : "<<endl;
	//for(auto& it:Ws) ostr<<it<<endl;
	//ostr<<endl;
	//ostr<<"G-line list : "<<endl;
	//for(auto& it:Gs) ostr<<it<<endl;
	//ostr<<endl;
	////ostr<<"SetUWInter : "<<endl;
	////ostr<<SetUWInter<<endl;
	////ostr<<endl<<"measuringline (miflavor:"<<miflavor<<",mjflavor:"<<mjflavor<<"): "<<(*measuringline)<<endl;
	//ostr<<endl<<"ref_line: "<<(*ref_line);
	//ostr<<endl<<"measuring_line: "<<(*measuring_line)<<endl;
	Conf.print(cout);
	ostr<<"***********************"<<endl<<endl<<flush;
};
////void diagram::print_normalization(ostream& ostr){
////	for(int i=0; i<normalization.size(); i++){
////		ostr<<setw(30)<<left<<i
////			<<setw(30)<<right<<normalization[i]<<endl;
////	}
////};
////void diagram::printQ(ostream& ostr, const int& iorder, const int& it){
////	Qacc[iorder-1].print(ostr,it);
////};
//
///*void diagram::printSwx(ostream& ostr, const int& iorder, const int& iflavor, const int& jflavor){
//	stringstream ss;
//	ostr<<setw(25)<<left<<"# w_n";
//	for(int ix=0; ix<Nx; ix++){
//		ss.str(""); ss<<"("<<Lattice.getxgrid(ix)(0)<<","<<Lattice.getxgrid(ix)(1)<<")";
//		ostr<<setw(50)<<right<<ss.str();
//	}
//	ostr<<endl;
//
//	ostr<<setprecision(15)<<scientific;
//	for(int n=0; n<Nw; n++){
//		ostr<<setw(25)<<left<<w(n);
//		for(int ix=0; ix<Nx; ix++){
//			ostr<<setw(25)<<right<<Swxacc[iorder-1][iflavor][jflavor].getOwx(n,ix).real()
//				<<setw(25)<<right<<Swxacc[iorder-1][iflavor][jflavor].getOwx(n,ix).imag();
//		}
//		ostr<<endl;
//	}
//};
//void diagram::printSwk(ostream& ostr, const int& iorder, const int& iflavor, const int& jflavor){
//	stringstream ss;
//	ostr<<setw(25)<<left<<"# w_n";
//	for(int ik=0; ik<Nx; ik++){
//		ss.str(""); ss<<"("<<Lattice.getkgrid(ik)(0)<<","<<Lattice.getkgrid(ik)(1)<<")";
//		ostr<<setw(50)<<right<<ss.str();
//	}
//	ostr<<endl;
//
//	ostr<<setprecision(15)<<scientific;
//	for(int n=0; n<Nw; n++){
//		ostr<<setw(25)<<left<<w(n);
//		for(int ik=0; ik<Nx; ik++){
//			ostr<<setw(25)<<right<<Swxacc[iorder-1][iflavor][jflavor].getOwk(n,ik).real()
//				<<setw(25)<<right<<Swxacc[iorder-1][iflavor][jflavor].getOwk(n,ik).imag();
//		}
//		ostr<<endl;
//	}
//};
//void diagram::printPwx(ostream& ostr, const int& iorder, const int& iflavor, const int& jflavor){
//	stringstream ss;
//	ostr<<setw(25)<<left<<"# w_n";
//	for(int ix=0; ix<Nx; ix++){
//		ss.str(""); ss<<"("<<Lattice.getxgrid(ix)(0)<<","<<Lattice.getxgrid(ix)(1)<<")";
//		ostr<<setw(50)<<right<<ss.str();
//	}
//	ostr<<endl;
//
//	ostr<<setprecision(15)<<scientific;
//	for(int n=0; n<Nw; n++){
//		ostr<<setw(25)<<left<<wB(n);
//		for(int ix=0; ix<Nx; ix++){
//			ostr<<setw(25)<<right<<Pwxacc[iorder][iflavor][jflavor].getOwx(n,ix).real()
//				<<setw(25)<<right<<Pwxacc[iorder][iflavor][jflavor].getOwx(n,ix).imag();
//		}
//		ostr<<endl;
//	}
//};
//void diagram::printPwk(ostream& ostr, const int& iorder, const int& iflavor, const int& jflavor){
//	stringstream ss;
//	ostr<<setw(25)<<left<<"# w_n";
//	for(int ik=0; ik<Nx; ik++){
//		ss.str(""); ss<<"("<<Lattice.getkgrid(ik)(0)<<","<<Lattice.getkgrid(ik)(1)<<")";
//		ostr<<setw(50)<<right<<ss.str();
//	}
//	ostr<<endl;
//
//	ostr<<setprecision(15)<<scientific;
//	for(int n=0; n<Nw; n++){
//		ostr<<setw(25)<<left<<wB(n);
//		for(int ik=0; ik<Nx; ik++){
//			ostr<<setw(25)<<right<<Pwxacc[iorder][iflavor][jflavor].getOwk(n,ik).real()
//				<<setw(25)<<right<<Pwxacc[iorder][iflavor][jflavor].getOwk(n,ik).imag();
//		}
//		ostr<<endl;
//	}
//};
//void diagram::printoutput(std::string foldername){
//	ofstream ofstr;
//	stringstream ss;
//	for(int i=1; i<Nordermax+1; i++){
//		ss.str(""); ss<<foldername<<"/Swx_"<<i<<"th_up.dat";
//		ofstr.open(ss.str().c_str());
//		printSwx(ofstr,i,1,1);
//		ofstr.close();
//
//		ss.str(""); ss<<foldername<<"/Swx_"<<i<<"th_dn.dat";
//		ofstr.open(ss.str().c_str());
//		printSwx(ofstr,i,0,0);
//		ofstr.close();
//
//
//		//ss.str(""); ss<<foldername<<"/Swk_"<<i<<"th_up.dat";
//		//ofstr.open(ss.str().c_str());
//		//printSwk(ofstr,i,1,1);
//		//ofstr.close();
//
//		//ss.str(""); ss<<foldername<<"/Swk_"<<i<<"th_dn.dat";
//		//ofstr.open(ss.str().c_str());
//		//printSwk(ofstr,i,0,0);
//		//ofstr.close();
//	}
//	for(int i=0; i<Nordermax+1; i++){
//		ss.str(""); ss<<foldername<<"/Pwx_"<<i<<"th_up_up.dat";
//		ofstr.open(ss.str().c_str());
//		printPwx(ofstr,i,1,1);
//		ofstr.close();
//
//		ss.str(""); ss<<foldername<<"/Pwx_"<<i<<"th_up_dn.dat";
//		ofstr.open(ss.str().c_str());
//		printPwx(ofstr,i,1,0);
//		ofstr.close();
//
//		ss.str(""); ss<<foldername<<"/Pwx_"<<i<<"th_dn_dn.dat";
//		ofstr.open(ss.str().c_str());
//		printPwx(ofstr,i,0,0);
//		ofstr.close();
//
//		ss.str(""); ss<<foldername<<"/Pwx_"<<i<<"th_dn_up.dat";
//		ofstr.open(ss.str().c_str());
//		printPwx(ofstr,i,0,1);
//		ofstr.close();
//
//		//ss.str(""); ss<<foldername<<"/Swk_"<<i<<"th_up.dat";
//		//ofstr.open(ss.str().c_str());
//		//printSwk(ofstr,i,1,1);
//		//ofstr.close();
//
//		//ss.str(""); ss<<foldername<<"/Swk_"<<i<<"th_dn.dat";
//		//ofstr.open(ss.str().c_str());
//		//printSwk(ofstr,i,0,0);
//		//ofstr.close();
//	}
//};
//void diagram::printoutput_S(std::string foldername){
//	ofstream ofstr;
//	stringstream ss;
//	for(int i=1; i<Nordermax+1; i++){
//		ss.str(""); ss<<foldername<<"/Swx_"<<i<<"th_up.dat";
//		ofstr.open(ss.str().c_str());
//		printSwx(ofstr,i,1,1);
//		ofstr.close();
//
//		ss.str(""); ss<<foldername<<"/Swx_"<<i<<"th_dn.dat";
//		ofstr.open(ss.str().c_str());
//		printSwx(ofstr,i,0,0);
//		ofstr.close();
//
//
//		//ss.str(""); ss<<foldername<<"/Swk_"<<i<<"th_up.dat";
//		//ofstr.open(ss.str().c_str());
//		//printSwk(ofstr,i,1,1);
//		//ofstr.close();
//
//		//ss.str(""); ss<<foldername<<"/Swk_"<<i<<"th_dn.dat";
//		//ofstr.open(ss.str().c_str());
//		//printSwk(ofstr,i,0,0);
//		//ofstr.close();
//	}
//};
//void diagram::printoutput_P(std::string foldername){
//	ofstream ofstr;
//	stringstream ss;
//	for(int i=0; i<Nordermax+1; i++){
//		ss.str(""); ss<<foldername<<"/Pwx_"<<i<<"th_up.dat";
//		ofstr.open(ss.str().c_str());
//		printPwx(ofstr,i,1,1);
//		ofstr.close();
//
//		ss.str(""); ss<<foldername<<"/Pwx_"<<i<<"th_dn.dat";
//		ofstr.open(ss.str().c_str());
//		printPwx(ofstr,i,0,0);
//		ofstr.close();
//
//
//		//ss.str(""); ss<<foldername<<"/Swk_"<<i<<"th_up.dat";
//		//ofstr.open(ss.str().c_str());
//		//printSwk(ofstr,i,1,1);
//		//ofstr.close();
//
//		//ss.str(""); ss<<foldername<<"/Swk_"<<i<<"th_dn.dat";
//		//ofstr.open(ss.str().c_str());
//		//printSwk(ofstr,i,0,0);
//		//ofstr.close();
//	}
//};*/
///*void diagram::updatebubbles(){
//	bubbles.clear();
//	
//	vector<bool> discovered(vertices.size(),false);
//	
//	list<vertex>::iterator v_ref = vertices.begin();
//	list<vertex>::iterator vi_iter;
//	list<vertex>::iterator vj_iter;
//
//	list<line>::iterator Gi_iter, Gj_iter;
//	list<line>::iterator UWi_iter, UWj_iter;
//
//	int count;
//	bool physicalpart;
//
//	for(int i=0; i<discovered.size(); i++){
//		if(!discovered[i]){
//			count = 0;
//			physicalpart = true;
//
//			vi_iter = next(vertices.begin(),i);
//			vj_iter = vi_iter;
//			discovered[i] = true;
//			Gi_iter = vi_iter->getconnected_Gout();
//			do{
//				Gj_iter = vj_iter->getconnected_Gout();
//				vj_iter = Gj_iter->getconnected_vj();
//				discovered[distance(v_ref,vj_iter)] = true;
//				count++;
//				if(!(Gj_iter->getphysical())) physicalpart = false;
//			} while(vi_iter!=vj_iter);
//			if(count==2 && physicalpart==true){
//				//vi_iter = Gi_iter->getconnected_vi();
//				//vj_iter = Gi_iter->getconnected_vj();
//				//UWi_iter = vi_iter->getconnected_W();
//				//UWj_iter = vj_iter->getconnected_W();
//				//if(UWi_iter!=UWj_iter && (UWi_iter->gettype()==1 && UWj_iter->gettype()==1))
//				//if(UWi_iter->gettype()==1 && UWj_iter->gettype()==1)
//				bubbles.emplace_front(Gi_iter,Gj_iter);
//			}
//		}
//	}
//};*/
///*void diagram::updatefermionicloops(){
//	//cout<<"** updatefermionicloops module starts **"<<endl;
//	fermionicloops.clear();
//	int Nvertices = vertices.size();
//	vector<bool> discovered(Nvertices,false);
//
//	list<list<vertex>::iterator> vs_tmp;
//
//	list<vertex>::iterator v_ref = vertices.begin();
//	list<vertex>::iterator vi_iter;
//	list<vertex>::iterator vj_iter;
//	list<line>::iterator UW_iter, Gi_out;
//	list<fermionicloop>::iterator loop_iter;
//
//	for(int i=0; i<Nvertices; i++){
//		if(!discovered[i]){
//			discovered[i] = true;
//			vs_tmp.clear();
//			vi_iter = next(v_ref,i);
//			vj_iter = vi_iter;
//			do{
//				vs_tmp.push_back(vj_iter);
//				Gi_out = vj_iter->getconnected_Gout();
//				vj_iter = Gi_out->getconnected_vj();
//				discovered[distance(v_ref,vj_iter)] = true;
//			} while(vi_iter!=vj_iter);
//
//			fermionicloops.emplace_front(vs_tmp);
//
//			loop_iter = fermionicloops.begin();
//			for(list<list<vertex>::iterator>::iterator iter=vs_tmp.begin(); iter!=vs_tmp.end(); ++iter)
//				(*iter)->setconnected_loop(loop_iter);
//		}
//	}
//
//	for(UW_iter=Ws.begin(); UW_iter!=Ws.end(); ++UW_iter){
//		vi_iter = UW_iter->getconnected_vi();
//		vj_iter = UW_iter->getconnected_vj();
//	}
//	//cout<<"** updatefermionicloops module ends **"<<endl;
//};
//void diagram::updateSetUWInter(){
//	SetUWInter.clear();
//	SetUWInter.setNfermionicloops( fermionicloops.size() );
//	SetUWInter.releaseloopflavor();
//
//	int iloop, jloop;
//	list<line>::iterator UWiter;
//	list<vertex>::iterator vi_iter, vj_iter;
//	for(UWiter=Ws.begin(); UWiter!=Ws.end(); ++UWiter){
//		vi_iter = UWiter->getconnected_vi();
//		vj_iter = UWiter->getconnected_vj();
//		if( vi_iter->getconnected_loop()!=vj_iter->getconnected_loop() ){
//			iloop = distance(fermionicloops.begin(),vi_iter->getconnected_loop());
//			jloop = distance(fermionicloops.begin(),vj_iter->getconnected_loop());
//			SetUWInter.push_back(*UWiter,iloop,jloop);
//		}
//	}
//
//	vi_iter = measuringline->getconnected_vi();
//	vj_iter = measuringline->getconnected_vj();
//	iloop = distance(fermionicloops.begin(),vi_iter->getconnected_loop());
//	jloop = distance(fermionicloops.begin(),vj_iter->getconnected_loop());
//	SetUWInter.fixloopflavor(iloop,miflavor);
//	SetUWInter.fixloopflavor(jloop,mjflavor);
//
//};
//bool diagram::checkconnection(){
//	list<vertex>::iterator v_iter;
//	list<line>::iterator G_iter, UW_iter;
//
//	list<fermionicloop>::iterator loop_iter, loop_iter_return;
//	list<list<vertex>::iterator> viters;
//	list<list<vertex>::iterator>::iterator viiter;
//
//	bool tmp = true;
//	for(v_iter=vertices.begin(); v_iter!=vertices.end(); ++v_iter){
//		G_iter = v_iter->getconnected_Gin();
//		if((v_iter->getflavor()!=G_iter->getjflavor()) || (v_iter->getx()!=G_iter->getxj())){ tmp = false; }
//		G_iter = v_iter->getconnected_Gout();
//		if((v_iter->getflavor()!=G_iter->getiflavor()) || (v_iter->getx()!=G_iter->getxi())){ tmp = false; }
//
//		UW_iter = v_iter->getconnected_W();
//		if( (!((v_iter->getflavor()==UW_iter->getiflavor()) && (v_iter->getx()==UW_iter->getxi())))
//				&& (!((v_iter->getflavor()==UW_iter->getjflavor()) && (v_iter->getx()==UW_iter->getxj()))))
//		{ tmp = false; }
//	}
//	for(G_iter=Gs.begin(); G_iter!=Gs.end(); ++G_iter){
//		v_iter = G_iter->getconnected_vi();
//		if((G_iter->getiflavor()!=v_iter->getflavor()) || (G_iter->getxi()!=v_iter->getx())){ tmp = false; }
//		v_iter = G_iter->getconnected_vj();
//		if((G_iter->getjflavor()!=v_iter->getflavor()) || (G_iter->getxj()!=v_iter->getx())){ tmp = false; }
//	}
//	for(UW_iter=Ws.begin(); UW_iter!=Ws.end(); ++UW_iter){
//		v_iter = UW_iter->getconnected_vi();
//		if((UW_iter->getiflavor()!=v_iter->getflavor()) || (UW_iter->getxi()!=v_iter->getx())){ tmp = false; }
//		v_iter = UW_iter->getconnected_vj();
//		if((UW_iter->getjflavor()!=v_iter->getflavor()) || (UW_iter->getxj()!=v_iter->getx())){ tmp = false; }
//	}
//	for(loop_iter=fermionicloops.begin(); loop_iter!=fermionicloops.end(); ++loop_iter){
//		viters = loop_iter->getvs();
//		for(viiter=viters.begin(); viiter!=viters.end(); ++viiter){
//			loop_iter_return = (*viiter)->getconnected_loop();
//			if(loop_iter!=loop_iter_return){ tmp = false; }
//		}
//	}
//	if(!tmp){
//		cout<<"***** ALERT! WRONG CONNECTION *****"<<endl;
//		printconfiguration(cout);
//		return false;
//	}
//
//	return true;
//};
//// When constructing adjLtmp max number of adjacent vetice will be 3/2*Nvertices
//bool diagram::checkconnectivity(){
//	//cout<<"** checkconnectivity module **"<<endl;
//	int Nvertices = vertices.size();
//	//if(Nvertices<3) return true;
//	//else{
//		graph topology;
//
//		vector<vector<int> > adjLtmp(Nvertices);
//		for(int i=0; i<Nvertices; i++) adjLtmp[i].reserve(3*Nvertices/2);
//
//		list<vertex>::iterator ref0 = vertices.begin();
//		list<line>::iterator iter;
//
//		int physical;
//		int iv, jv;
//		for(iter=Gs.begin(); iter!=Gs.end(); ++iter){
//			physical = iter->getphysical();
//			iv = distance(ref0,iter->getconnected_vi());
//			jv = distance(ref0,iter->getconnected_vj());
//			if((iv!=jv) && physical){
//				adjLtmp[iv].push_back(jv);
//				adjLtmp[jv].push_back(iv);
//			}
//		}
//		for(iter=Ws.begin(); iter!=Ws.end(); ++iter){
//			physical = iter->getphysical();
//			iv = distance(ref0,iter->getconnected_vi());
//			jv = distance(ref0,iter->getconnected_vj());
//			if((iv!=jv) && physical){
//				adjLtmp[iv].push_back(jv);
//				adjLtmp[jv].push_back(iv);
//			}
//		}
//
//		//cout<<"adjL info"<<endl;
//		//vector<int>::iterator viter;
//		//for(int i=0; i<Nvertices; i++){
//		//	for(viter=adjLtmp[i].begin(); viter!=adjLtmp[i].end(); viter++)
//		//		cout<<(*viter);
//		//	cout<<endl;
//		//}
//
//		//topology.set(Nvertices,adjLtmp,adjLtype);
//		topology.set(Nvertices,adjLtmp);
//
//		return topology.edge_connectivity();
//	//}
//};
//bool diagram::checkconnectivity(const vector<vector<int> >& GadjL_i, const vector<vector<int> >& WadjL_i){
//	//cout<<"** checkconnectivity module **"<<endl;
//	int Nvertices = GadjL_i.size();
//	//if(Nvertices<3) return true;
//	//else{
//		vector<vector<int> > adjLtmp(Nvertices);
//		for(int i=0; i<Nvertices; i++){
//			adjLtmp[i].reserve(3*Nvertices/2);
//			//adjLtype[i].reserve(3*Nvertices/2);
//			for(int j=0; j<GadjL_i[i].size(); j++)
//				adjLtmp[i].push_back(GadjL_i[i][j]);
//			for(int j=0; j<WadjL_i[i].size(); j++)
//				adjLtmp[i].push_back(WadjL_i[i][j]);
//		}
//
//		//vector<vector<int> > adjLtmp(Nvertices);
//		//for(int i=0; i<Nvertices; i++) adjLtmp[i].reserve(3*Nvertices/2);
//
//		//list<vertex>::iterator ref0 = vertices.begin();
//		//list<line>::iterator iter;
//
//		//int physical;
//		//int iv, jv;
//		//for(iter=Gs.begin(); iter!=Gs.end(); ++iter){
//		//	physical = iter->getphysical();
//		//	iv = distance(ref0,iter->getconnected_vi());
//		//	jv = distance(ref0,iter->getconnected_vj());
//		//	if((iv!=jv) && physical){
//		//		adjLtmp[iv].push_back(jv);
//		//		adjLtmp[jv].push_back(iv);
//		//	}
//		//}
//		//for(iter=Ws.begin(); iter!=Ws.end(); ++iter){
//		//	physical = iter->getphysical();
//		//	iv = distance(ref0,iter->getconnected_vi());
//		//	jv = distance(ref0,iter->getconnected_vj());
//		//	if((iv!=jv) && physical){
//		//		adjLtmp[iv].push_back(jv);
//		//		adjLtmp[jv].push_back(iv);
//		//	}
//		//}
//
//		//cout<<"adjL info"<<endl;
//		//vector<int>::iterator viter;
//		//for(int i=0; i<Nvertices; i++){
//		//	for(viter=adjLtmp[i].begin(); viter!=adjLtmp[i].end(); viter++)
//		//		cout<<(*viter);
//		//	cout<<endl;
//		//}
//
//		//topology.set(Nvertices,adjLtmp,adjLtype);
//		graph topology;
//		topology.set(Nvertices,adjLtmp);
//
//		return topology.edge_connectivity();
//	//}
//};
//
//bool diagram::checkconnectivityforremove(list<line>::iterator UWremove){
//	//cout<<"** checkconnectivity module **"<<endl;
//	int Nvertices = vertices.size();
//	graph topology;
//
//	vector<vector<int> > adjLtmp(Nvertices);
//	for(int i=0; i<Nvertices; i++) adjLtmp[i].reserve(3*Nvertices/2);
//
//	list<vertex>::iterator ref0 = vertices.begin();
//	list<line>::iterator iter;
//
//	int physical;
//	int iv, jv;
//	for(iter=Gs.begin(); iter!=Gs.end(); ++iter){
//		physical = iter->getphysical();
//		iv = distance(ref0,iter->getconnected_vi());
//		jv = distance(ref0,iter->getconnected_vj());
//		if((iv!=jv) && physical){
//			adjLtmp[iv].push_back(jv);
//			adjLtmp[jv].push_back(iv);
//		}
//	}
//	for(iter=Ws.begin(); iter!=Ws.end(); ++iter){
//		if(iter!=UWremove){
//			physical = iter->getphysical();
//			iv = distance(ref0,iter->getconnected_vi());
//			jv = distance(ref0,iter->getconnected_vj());
//			if((iv!=jv) && physical){
//				adjLtmp[iv].push_back(jv);
//				adjLtmp[jv].push_back(iv);
//			}
//		}
//	}
//	//for(iter=Ws.begin(); iter!=Ws.end(); ++iter){
//	//	if(iter!=UWremove){
//	//		physical = iter->getphysical();
//	//		iv = distance(ref0,iter->getconnected_vi());
//	//		jv = distance(ref0,iter->getconnected_vj());
//	//		if((iv!=jv) && physical){
//	//			adjLtmp[iv].push_back(jv);
//	//			adjLtmp[jv].push_back(iv);
//	//		}
//	//	}
//	//}
//
//	//cout<<"adjL info"<<endl;
//	//vector<int>::iterator viter;
//	//for(int i=0; i<Nvertices; i++){
//	//	for(viter=adjLtmp[i].begin(); viter!=adjLtmp[i].end(); viter++)
//	//		cout<<(*viter);
//	//	cout<<endl;
//	//}
//
//	topology.set(Nvertices,adjLtmp);
//
//	return topology.edge_connectivity();
//};*/
////bool diagram::check_irreducibility_for_remove(list<line>::iterator Wremove){
////	//cout<<"** check irreducibility module for remove update **"<<endl;
////	int Nvertices = vertices.size();
////	if(Nvertices<5) return true;
////	else{
////		//cout<<"real check module"<<endl<<flush;
////		//cout<<"Nvertices = "<<Nvertices<<endl<<flush;
////		//printconfiguration(cout);
////
////		graph topology;
////
////		vector<vector<int> > adjLtmp(Nvertices);
////		for(int i=0; i<Nvertices; i++) adjLtmp[i].reserve(3*Nvertices/2);
////
////		list<vertex>::iterator ref0 = vertices.begin();
////		list<line>::iterator iter;
////
////		list<vertex>::iterator vi_iter = Wremove->getconnected_vi();
////		list<vertex>::iterator vj_iter = Wremove->getconnected_vj();
////
////		list<line>::iterator Gi_in = vi_iter->getconnected_Gin();
////		list<line>::iterator Gj_in = vj_iter->getconnected_Gin();
////		list<line>::iterator Gi_out = vi_iter->getconnected_Gout();
////		list<line>::iterator Gj_out = vj_iter->getconnected_Gout();
////
////		int iv_remove = distance(ref0,vi_iter);
////		int jv_remove = distance(ref0,vj_iter);
////
////		//cout<<"iv_remove: "<<iv_remove<<endl
////		//	<<"jv_remove: "<<jv_remove<<endl<<flush;
////
////		int iv_prev = distance(ref0,Gi_in->getconnected_vi());
////		int jv_prev = distance(ref0,Gj_in->getconnected_vi());
////		int iv_next = distance(ref0,Gi_out->getconnected_vj());
////		int jv_next = distance(ref0,Gj_out->getconnected_vj());
////
////		//cout<<"iv_prev: "<<iv_prev<<endl
////		//	<<"jv_prev: "<<jv_prev<<endl
////		//	<<"iv_next: "<<iv_next<<endl
////		//	<<"jv_next: "<<jv_next<<endl<<flush;
////
////		//cout<<"G push_back"<<endl<<flush;
////		if(Gi_in==Gj_out){
////			//cout<<"case 1"<<endl;
////			if(Gi_out->getphysical()&&Gj_in->getphysical()){
////				adjLtmp[iv_next].push_back(jv_prev);
////				adjLtmp[jv_prev].push_back(iv_next);
////			}
////		}
////		else if(Gj_in==Gi_out){
////			//cout<<"case 2"<<endl;
////			if(Gi_in->getphysical()&&Gj_out->getphysical()){
////				adjLtmp[iv_prev].push_back(jv_next);
////				adjLtmp[jv_next].push_back(iv_prev);
////			}
////		}
////		else{
////			//cout<<"case 3"<<endl;
////			if(Gi_in->getphysical()&&Gi_out->getphysical()){
////				//cout<<"case 3, a"<<endl;
////				adjLtmp[iv_prev].push_back(iv_next);
////				adjLtmp[iv_next].push_back(iv_prev);
////			}
////			if(Gj_in->getphysical()&&Gj_out->getphysical()){
////				//cout<<"case 3, b"<<endl;
////				adjLtmp[jv_prev].push_back(jv_next);
////				adjLtmp[jv_next].push_back(jv_prev);
////			}
////		}
////
////		//cout<<"adjL info (Nvertices: "<<Nvertices<<")"<<endl;
////		//vector<int>::iterator viter;
////		//for(int i=0; i<Nvertices; i++){
////		//	for(viter=adjLtmp[i].begin(); viter!=adjLtmp[i].end(); viter++)
////		//		cout<<(*viter);
////		//	cout<<endl;
////		//}
////
////
////		int physical;
////		int iv, jv;
////		for(iter=Gs.begin(); iter!=Gs.end(); ++iter){
////			physical = iter->getphysical();
////			iv = distance(ref0,iter->getconnected_vi());
////			jv = distance(ref0,iter->getconnected_vj());
////			if((iv!=iv_remove) && (iv!=jv_remove) && (jv!=iv_remove) && (jv!=jv_remove) && (iv!=jv) && physical){
////				adjLtmp[iv].push_back(jv);
////				adjLtmp[jv].push_back(iv);
////			}
////		}
////
////		//cout<<"W push_back"<<endl<<flush;
////		for(iter=Ws.begin(); iter!=Ws.end(); ++iter){
////			if(iter!=Wremove){
////				physical = iter->getphysical();
////				iv = distance(ref0,iter->getconnected_vi());
////				jv = distance(ref0,iter->getconnected_vj());
////				if((iv!=jv) && physical){
////					adjLtmp[iv].push_back(jv);
////					adjLtmp[jv].push_back(iv);
////				}
////			}
////		}
////
////
////		//cout<<"adjL info (Nvertices: "<<Nvertices<<")"<<endl;
////		//vector<int>::iterator viter;
////		//for(int i=0; i<adjLtmp.size(); i++){
////		//	for(viter=adjLtmp[i].begin(); viter!=adjLtmp[i].end(); viter++)
////		//		cout<<(*viter);
////		//	cout<<endl;
////		//}
////
////		vector<int>::iterator viter;
////		int tmp = max(iv_remove,jv_remove);
////		adjLtmp.erase(adjLtmp.begin()+tmp);
////		for(int i=0; i<adjLtmp.size(); i++){
////			for(viter=adjLtmp[i].begin(); viter!=adjLtmp[i].end(); viter++)
////				if(*viter>tmp) (*viter)--;
////		}
////
////		tmp = min(iv_remove,jv_remove);
////		adjLtmp.erase(adjLtmp.begin()+tmp);
////		for(int i=0; i<adjLtmp.size(); i++){
////			for(viter=adjLtmp[i].begin(); viter!=adjLtmp[i].end(); viter++)
////				if(*viter>tmp) (*viter)--;
////		}
////
////		//cout<<"adjL info (Nvertices: "<<Nvertices<<")"<<endl;
////		//vector<int>::iterator viter;
////		//for(int i=0; i<adjLtmp.size(); i++){
////		//	for(viter=adjLtmp[i].begin(); viter!=adjLtmp[i].end(); viter++)
////		//		cout<<(*viter);
////		//	cout<<endl;
////		//}
////		//cout<<"at this point, adjLtmp is inserted to topology"<<endl<<flush;
////
////		//exit(EXIT_FAILURE);
////
////		//topology.set(Nvertices,adjLtmp);
////		topology.set(adjLtmp.size(),adjLtmp);
////		return topology.two_edge_connectivity();
////	}
////};
////bool diagram::check_irreducibility_for_swap(list<line>::iterator target){
////	//cout<<"** checkconnectivity module **"<<endl;
////	int Nvertices = vertices.size();
////
////	graph topology;
////
////	vector<vector<int> > adjLtmp(Nvertices);
////	for(int i=0; i<Nvertices; i++) adjLtmp[i].reserve(3*Nvertices/2);
////
////	list<vertex>::iterator ref0 = vertices.begin();
////	list<line>::iterator iter;
////
////	int isG, physical;
////	int iv, jv;
////	for(iter=Gs.begin(); iter!=Gs.end(); ++iter){
////		if(iter!=target && iter!=ref_line){
////			iv = distance(ref0,iter->getconnected_vi());
////			jv = distance(ref0,iter->getconnected_vj());
////			if(iv!=jv){
////				adjLtmp[iv].push_back(jv);
////				adjLtmp[jv].push_back(iv);
////			}
////		}
////	}
////	for(iter=Ws.begin(); iter!=Ws.end(); ++iter){
////		iv = distance(ref0,iter->getconnected_vi());
////		jv = distance(ref0,iter->getconnected_vj());
////		if(iv!=jv){
////			adjLtmp[iv].push_back(jv);
////			adjLtmp[jv].push_back(iv);
////		}
////	}
////
////	//cout<<"adjL info"<<endl;
////	//vector<int>::iterator viter;
////	//for(int i=0; i<Nvertices; i++){
////	//	for(viter=adjLtmp[i].begin(); viter!=adjLtmp[i].end(); viter++)
////	//		cout<<(*viter);
////	//	cout<<endl;
////	//}
////
////	topology.set(Nvertices,adjLtmp);
////	if(topology.edge_connectivity()){
////		if(Nvertices<3) return true;
////		else return topology.two_edge_connectivity();
////	}
////	else return false;
////
////	//topology.set(Nvertices,adjLtmp);
////	//return topology.two_edge_connectivity();
////};
//
///*bool diagram::checkconnectivityforreconnect(list<line>::iterator G_i, list<line>::iterator G_j){
//	//cout<<"** checkconnectivity module **"<<endl;
//	int Nvertices = vertices.size();
//	graph topology;
//
//	vector<vector<int> > adjLtmp(Nvertices);
//	for(int i=0; i<Nvertices; i++) adjLtmp[i].reserve(3*Nvertices/2);
//
//	list<vertex>::iterator ref0 = vertices.begin();
//	list<line>::iterator iter;
//
//	int physical;
//	int iv, jv;
//	for(iter=Gs.begin(); iter!=Gs.end(); ++iter){
//		physical = iter->getphysical();
//
//		if(iter==G_i){
//			//cout<<"1"<<flush;
//			iv = distance(ref0,iter->getconnected_vi());
//			jv = distance(ref0,G_j->getconnected_vj());
//		}
//		else if(iter==G_j){
//			//cout<<"2"<<flush;
//			iv = distance(ref0,iter->getconnected_vi());
//			jv = distance(ref0,G_i->getconnected_vj());
//		}
//		else{
//			//cout<<"0"<<flush;
//			iv = distance(ref0,iter->getconnected_vi());
//			jv = distance(ref0,iter->getconnected_vj());
//		}
//		if((iv!=jv) && physical){
//			adjLtmp[iv].push_back(jv);
//			adjLtmp[jv].push_back(iv);
//		}
//	}
//	//cout<<" "<<flush;
//	for(iter=Ws.begin(); iter!=Ws.end(); ++iter){
//		physical = iter->getphysical();
//
//		iv = distance(ref0,iter->getconnected_vi());
//		jv = distance(ref0,iter->getconnected_vj());
//		if((iv!=jv) && physical){
//			adjLtmp[iv].push_back(jv);
//			adjLtmp[jv].push_back(iv);
//		}
//	}
//	//for(iter=Ws.begin(); iter!=Ws.end(); ++iter){
//	//	physical = iter->getphysical();
//
//	//	iv = distance(ref0,iter->getconnected_vi());
//	//	jv = distance(ref0,iter->getconnected_vj());
//	//	if((iv!=jv) && physical){
//	//		adjLtmp[iv].push_back(jv);
//	//		adjLtmp[jv].push_back(iv);
//	//	}
//	//}
//
//	topology.set(Nvertices,adjLtmp);
//
//	return topology.edge_connectivity();
//};*/
///*bool diagram::checkconnectivityforinterswap(list<line>::iterator target){
//	//cout<<"** checkconnectivity module **"<<endl;
//	int Nvertices = vertices.size();
//	graph topology;
//
//	vector<vector<int> > adjLtmp(Nvertices);
//	for(int i=0; i<Nvertices; i++) adjLtmp[i].reserve(3*Nvertices/2);
//
//	list<vertex>::iterator ref0 = vertices.begin();
//	list<line>::iterator iter;
//
//	int isG, physical;
//	int iv, jv;
//	for(iter=Gs.begin(); iter!=Gs.end(); ++iter){
//		if(iter!=target){
//			iv = distance(ref0,iter->getconnected_vi());
//			jv = distance(ref0,iter->getconnected_vj());
//			if(iv!=jv){
//				adjLtmp[iv].push_back(jv);
//				adjLtmp[jv].push_back(iv);
//			}
//		}
//	}
//	for(iter=Us.begin(); iter!=Us.end(); ++iter){
//		if(iter!=target){
//			iv = distance(ref0,iter->getconnected_vi());
//			jv = distance(ref0,iter->getconnected_vj());
//			if(iv!=jv){
//				adjLtmp[iv].push_back(jv);
//				adjLtmp[jv].push_back(iv);
//			}
//		}
//	}
//	for(iter=Ws.begin(); iter!=Ws.end(); ++iter){
//		if(iter!=target){
//			iv = distance(ref0,iter->getconnected_vi());
//			jv = distance(ref0,iter->getconnected_vj());
//			if(iv!=jv){
//				adjLtmp[iv].push_back(jv);
//				adjLtmp[jv].push_back(iv);
//			}
//		}
//	}
//
//	//cout<<"adjL info"<<endl;
//	//vector<int>::iterator viter;
//	//for(int i=0; i<Nvertices; i++){
//	//	for(viter=adjLtmp[i].begin(); viter!=adjLtmp[i].end(); viter++)
//	//		cout<<(*viter);
//	//	cout<<endl;
//	//}
//
//	topology.set(Nvertices,adjLtmp);
//
//	return topology.edge_connectivity();
//};*/
//
////bool diagram::checkirreducibility(){
////	//cout<<"** checkirreducibility module **"<<endl;
////	int Nvertices = vertices.size();
////	if(Nvertices<3) return true;
////	else{
////		graph topology;
////
////		vector<vector<int> > adjLtmp(Nvertices);
////		for(int i=0; i<Nvertices; i++){
////			adjLtmp[i].reserve(3*Nvertices/2);
////		}
////
////		list<vertex>::iterator ref0 = vertices.begin();
////		list<line>::iterator iter;
////
////		int physical;
////		int iv, jv;
////		for(iter=Gs.begin(); iter!=Gs.end(); ++iter){
////			physical = iter->getphysical();
////
////			iv = distance(ref0,iter->getconnected_vi());
////			jv = distance(ref0,iter->getconnected_vj());
////			if((iv!=jv) && physical){
////				adjLtmp[iv].push_back(jv);
////				adjLtmp[jv].push_back(iv);
////			}
////		}
////		for(iter=Ws.begin(); iter!=Ws.end(); ++iter){
////			physical = iter->getphysical();
////
////			iv = distance(ref0,iter->getconnected_vi());
////			jv = distance(ref0,iter->getconnected_vj());
////			if((iv!=jv) && physical){
////				adjLtmp[iv].push_back(jv);
////				adjLtmp[jv].push_back(iv);
////			}
////		}
////
////		//cout<<"adjL info"<<endl;
////		//vector<int>::iterator viter;
////		//for(int i=0; i<Nvertices; i++){
////		//	for(viter=adjLtmp[i].begin(); viter!=adjLtmp[i].end(); viter++)
////		//		cout<<(*viter);
////		//	cout<<endl;
////		//}
////		//cout<<"adjL type info"<<endl;
////		//for(int i=0; i<Nvertices; i++){
////		//	for(viter=adjLtype[i].begin(); viter!=adjLtype[i].end(); viter++)
////		//		cout<<(*viter);
////		//	cout<<endl;
////		//}
////
////		topology.set(Nvertices,adjLtmp);
////		return topology.two_edge_connectivity();
////	}
////};
///*bool diagram::checkirreducibility(vector<vector<int> > GadjL_i, vector<vector<int> > WadjL_i){
//	//cout<<"** checkirreducibility module **"<<endl;
//	int Nvertices = GadjL_i.size();
//	if(Nvertices<3) return true;
//	else{
//		vector<vector<int> > adjLtmp(Nvertices);
//		for(int i=0; i<Nvertices; i++){
//			adjLtmp[i].reserve(3*Nvertices/2);
//			//adjLtype[i].reserve(3*Nvertices/2);
//			for(int j=0; j<GadjL_i[i].size(); j++)
//				adjLtmp[i].push_back(GadjL_i[i][j]);
//			for(int j=0; j<WadjL_i[i].size(); j++)
//				adjLtmp[i].push_back(WadjL_i[i][j]);
//		}
//
//		//cout<<"adjL info"<<endl;
//		//vector<int>::iterator viter;
//		//for(int i=0; i<Nvertices; i++){
//		//	for(viter=adjLtmp[i].begin(); viter!=adjLtmp[i].end(); viter++)
//		//		cout<<(*viter);
//		//	cout<<endl;
//		//}
//		//cout<<"adjL type info"<<endl;
//		//for(int i=0; i<Nvertices; i++){
//		//	for(viter=adjLtype[i].begin(); viter!=adjLtype[i].end(); viter++)
//		//		cout<<(*viter);
//		//	cout<<endl;
//		//}
//
//		graph topology;
//		//topology.set(Nvertices,adjLtmp,adjLtype);
//		topology.set(Nvertices,adjLtmp);
//		return topology.two_edge_connectivity();
//	}
//}*/
///*bool diagram::checkirreducibility(const bool& isG_i){
//	//cout<<"** checkirreducibility module **"<<endl;
//	int Nvertices = vertices.size();
//	if(Nvertices<3) return true;
//	else{
//		graph topology;
//
//		vector<vector<int> > adjLtmp(Nvertices);
//		vector<vector<int> > adjLtype(Nvertices);
//		for(int i=0; i<Nvertices; i++){
//			adjLtmp[i].reserve(3*Nvertices/2);
//			adjLtype[i].reserve(3*Nvertices/2);
//		}
//
//		list<vertex>::iterator ref0 = vertices.begin();
//		list<line>::iterator iter;
//
//		int physical;
//		int iv, jv;
//		for(iter=Gs.begin(); iter!=Gs.end(); ++iter){
//			physical = iter->getphysical();
//
//			iv = distance(ref0,iter->getconnected_vi());
//			jv = distance(ref0,iter->getconnected_vj());
//			if((iv!=jv) && physical){
//				adjLtype[iv].push_back(0);
//				adjLtype[jv].push_back(0);
//
//				adjLtmp[iv].push_back(jv);
//				adjLtmp[jv].push_back(iv);
//			}
//		}
//		for(iter=Ws.begin(); iter!=Ws.end(); ++iter){
//			physical = iter->getphysical();
//
//			iv = distance(ref0,iter->getconnected_vi());
//			jv = distance(ref0,iter->getconnected_vj());
//			if((iv!=jv) && physical){
//				adjLtype[iv].push_back(1);
//				adjLtype[jv].push_back(1);
//
//				adjLtmp[iv].push_back(jv);
//				adjLtmp[jv].push_back(iv);
//			}
//		}
//
//		//cout<<"adjL info"<<endl;
//		//vector<int>::iterator viter;
//		//for(int i=0; i<Nvertices; i++){
//		//	for(viter=adjLtmp[i].begin(); viter!=adjLtmp[i].end(); viter++)
//		//		cout<<(*viter);
//		//	cout<<endl;
//		//}
//		//cout<<"adjL type info"<<endl;
//		//for(int i=0; i<Nvertices; i++){
//		//	for(viter=adjLtype[i].begin(); viter!=adjLtype[i].end(); viter++)
//		//		cout<<(*viter);
//		//	cout<<endl;
//		//}
//
//		topology.set(Nvertices,adjLtmp,adjLtype);
//		return topology.two_edge_connectivity(2+(isG_i?1:0));
//	}
//};
//bool diagram::brutecheckirreducibility(){
//	//cout<<"** checkirreducibility module **"<<endl;
//	int Nvertices = vertices.size();
//	if(Nvertices<3) return true;
//	else{
//		graph topology;
//
//		vector<vector<int> > adjLtmp(Nvertices);
//		vector<vector<int> > adjLtype(Nvertices);
//
//		list<vertex>::iterator ref0 = vertices.begin();
//		list<line>::iterator iter;
//		list<line>::iterator iter_tmp;
//
//		bool connectivity;
//
//		int physical;
//		int iv, jv;
//		for(iter_tmp = Gs.begin(); iter_tmp!=Gs.end(); ++iter_tmp){
//			if(iter_tmp->getphysical()){
//				for(int i=0; i<Nvertices; i++){
//					adjLtmp[i].clear();
//					adjLtype[i].clear();
//				}
//				for(iter=Gs.begin(); iter!=Gs.end(); ++iter){
//					if(iter!=iter_tmp){
//						physical = iter->getphysical();
//
//						iv = distance(ref0,iter->getconnected_vi());
//						jv = distance(ref0,iter->getconnected_vj());
//						if((iv!=jv) && physical){
//							adjLtype[iv].push_back(0);
//							adjLtype[jv].push_back(0);
//
//							adjLtmp[iv].push_back(jv);
//							adjLtmp[jv].push_back(iv);
//						}
//					}
//				}
//				for(iter=Ws.begin(); iter!=Ws.end(); ++iter){
//					physical = iter->getphysical();
//
//					iv = distance(ref0,iter->getconnected_vi());
//					jv = distance(ref0,iter->getconnected_vj());
//					if((iv!=jv) && physical){
//						adjLtype[iv].push_back(1);
//						adjLtype[jv].push_back(1);
//
//						adjLtmp[iv].push_back(jv);
//						adjLtmp[jv].push_back(iv);
//					}
//				}
//				topology.set(Nvertices,adjLtmp,adjLtype);
//				connectivity = topology.edge_connectivity();
//				if( !connectivity ) return false;
//			}
//		}
//		for(iter_tmp = Ws.begin(); iter_tmp!=Ws.end(); ++iter_tmp){
//			if(iter_tmp->getphysical()){
//				for(int i=0; i<Nvertices; i++){
//					adjLtmp[i].clear();
//					adjLtype[i].clear();
//				}
//				for(iter=Gs.begin(); iter!=Gs.end(); ++iter){
//					physical = iter->getphysical();
//
//					iv = distance(ref0,iter->getconnected_vi());
//					jv = distance(ref0,iter->getconnected_vj());
//					if((iv!=jv) && physical){
//						adjLtype[iv].push_back(0);
//						adjLtype[jv].push_back(0);
//
//						adjLtmp[iv].push_back(jv);
//						adjLtmp[jv].push_back(iv);
//					}
//				}
//				for(iter=Ws.begin(); iter!=Ws.end(); ++iter){
//					if(iter!=iter_tmp){
//						physical = iter->getphysical();
//
//						iv = distance(ref0,iter->getconnected_vi());
//						jv = distance(ref0,iter->getconnected_vj());
//						if((iv!=jv) && physical){
//							adjLtype[iv].push_back(1);
//							adjLtype[jv].push_back(1);
//
//							adjLtmp[iv].push_back(jv);
//							adjLtmp[jv].push_back(iv);
//						}
//					}
//				}
//				topology.set(Nvertices,adjLtmp,adjLtype);
//				connectivity = topology.edge_connectivity();
//				if( !connectivity ) return false;
//			}
//		}
//
//		return true;
//
//		//cout<<"adjL info"<<endl;
//		//vector<int>::iterator viter;
//		//for(int i=0; i<Nvertices; i++){
//		//	for(viter=adjLtmp[i].begin(); viter!=adjLtmp[i].end(); viter++)
//		//		cout<<(*viter);
//		//	cout<<endl;
//		//}
//		//cout<<"adjL type info"<<endl;
//		//for(int i=0; i<Nvertices; i++){
//		//	for(viter=adjLtype[i].begin(); viter!=adjLtype[i].end(); viter++)
//		//		cout<<(*viter);
//		//	cout<<endl;
//		//}
//
//		//topology.set(Nvertices,adjLtmp,adjLtype);
//		//return topology.two_edge_connectivity(3);
//	}
//};*/
///*bool diagram::check_irreducibility_for_remove(list<line>::iterator Wremove){
//	//cout<<"** check irreducibility module for remove update **"<<endl;
//	int Nvertices = vertices.size();
//	if(Nvertices<5) return true;
//	else{
//		//cout<<"real check module"<<endl<<flush;
//		//cout<<"Nvertices = "<<Nvertices<<endl<<flush;
//		//printconfiguration(cout);
//
//		graph topology;
//
//		vector<vector<int> > adjLtmp(Nvertices);
//		for(int i=0; i<Nvertices; i++) adjLtmp[i].reserve(3*Nvertices/2);
//
//		list<vertex>::iterator ref0 = vertices.begin();
//		list<line>::iterator iter;
//
//		list<vertex>::iterator vi_iter = Wremove->getconnected_vi();
//		list<vertex>::iterator vj_iter = Wremove->getconnected_vj();
//
//		list<line>::iterator Gi_in = vi_iter->getconnected_Gin();
//		list<line>::iterator Gj_in = vj_iter->getconnected_Gin();
//		list<line>::iterator Gi_out = vi_iter->getconnected_Gout();
//		list<line>::iterator Gj_out = vj_iter->getconnected_Gout();
//
//		int iv_remove = distance(ref0,vi_iter);
//		int jv_remove = distance(ref0,vj_iter);
//
//		int iv_prev = distance(ref0,Gi_in->getconnected_vi());
//		int jv_prev = distance(ref0,Gj_in->getconnected_vi());
//		int iv_next = distance(ref0,Gi_out->getconnected_vj());
//		int jv_next = distance(ref0,Gj_out->getconnected_vj());
//
//		//cout<<"iv_prev: "<<iv_prev<<endl
//		//	<<"jv_prev: "<<jv_prev<<endl
//		//	<<"iv_next: "<<iv_next<<endl
//		//	<<"jv_next: "<<jv_next<<endl<<flush;
//
//		//cout<<"G push_back"<<endl<<flush;
//		if(Gi_in==Gj_out){
//			//cout<<"case 1"<<endl;
//			if(Gi_out->getphysical()&&Gj_in->getphysical()){
//				adjLtmp[iv_next].push_back(jv_prev);
//				adjLtmp[jv_prev].push_back(iv_next);
//			}
//		}
//		else if(Gj_in==Gi_out){
//			//cout<<"case 2"<<endl;
//			if(Gi_in->getphysical()&&Gj_out->getphysical()){
//				adjLtmp[iv_prev].push_back(jv_next);
//				adjLtmp[jv_next].push_back(iv_prev);
//			}
//		}
//		else{
//			//cout<<"case 3"<<endl;
//			if(Gi_in->getphysical()&&Gi_out->getphysical()){
//				//cout<<"case 3, a"<<endl;
//				adjLtmp[iv_prev].push_back(iv_next);
//				adjLtmp[iv_next].push_back(iv_prev);
//			}
//			if(Gj_in->getphysical()&&Gj_out->getphysical()){
//				//cout<<"case 3, b"<<endl;
//				adjLtmp[jv_prev].push_back(jv_next);
//				adjLtmp[jv_next].push_back(jv_prev);
//			}
//		}
//
//		//cout<<"adjL info (Nvertices: "<<Nvertices<<")"<<endl;
//		//vector<int>::iterator viter;
//		//for(int i=0; i<Nvertices; i++){
//		//	for(viter=adjLtmp[i].begin(); viter!=adjLtmp[i].end(); viter++)
//		//		cout<<(*viter);
//		//	cout<<endl;
//		//}
//
//		int physical;
//		int iv, jv;
//		for(iter=Gs.begin(); iter!=Gs.end(); ++iter){
//			physical = iter->getphysical();
//			iv = distance(ref0,iter->getconnected_vi());
//			jv = distance(ref0,iter->getconnected_vj());
//			if((iv!=iv_remove) && (iv!=jv_remove) && (jv!=iv_remove) && (jv!=jv_remove) && (iv!=jv) && physical){
//				adjLtmp[iv].push_back(jv);
//				adjLtmp[jv].push_back(iv);
//			}
//		}
//
//		//cout<<"U push_back"<<endl<<flush;
//		for(iter=Ws.begin(); iter!=Ws.end(); ++iter){
//			if(iter!=Wremove){
//				physical = iter->getphysical();
//				iv = distance(ref0,iter->getconnected_vi());
//				jv = distance(ref0,iter->getconnected_vj());
//				if((iv!=jv) && physical){
//					adjLtmp[iv].push_back(jv);
//					adjLtmp[jv].push_back(iv);
//				}
//			}
//		}
//
//
//		//cout<<"adjL info (Nvertices: "<<Nvertices<<")"<<endl;
//		//vector<int>::iterator viter;
//		//for(int i=0; i<adjLtmp.size(); i++){
//		//	for(viter=adjLtmp[i].begin(); viter!=adjLtmp[i].end(); viter++)
//		//		cout<<(*viter);
//		//	cout<<endl;
//		//}
//
//		vector<int>::iterator viter;
//		int tmp = max(iv_remove,jv_remove);
//		adjLtmp.erase(adjLtmp.begin()+tmp);
//		for(int i=0; i<adjLtmp.size(); i++){
//			for(viter=adjLtmp[i].begin(); viter!=adjLtmp[i].end(); viter++)
//				if(*viter>tmp) (*viter)--;
//		}
//
//		tmp = min(iv_remove,jv_remove);
//		adjLtmp.erase(adjLtmp.begin()+tmp);
//		for(int i=0; i<adjLtmp.size(); i++){
//			for(viter=adjLtmp[i].begin(); viter!=adjLtmp[i].end(); viter++)
//				if(*viter>tmp) (*viter)--;
//		}
//
//		//cout<<"adjL info (Nvertices: "<<Nvertices<<")"<<endl;
//		//vector<int>::iterator viter;
//		//for(int i=0; i<adjLtmp.size(); i++){
//		//	for(viter=adjLtmp[i].begin(); viter!=adjLtmp[i].end(); viter++)
//		//		cout<<(*viter);
//		//	cout<<endl;
//		//}
//		//cout<<"at this point, adjLtmp is inserted to topology"<<endl<<flush;
//
//		//exit(EXIT_FAILURE);
//
//		//topology.set(Nvertices,adjLtmp);
//		topology.set(adjLtmp.size(),adjLtmp);
//		return topology.two_edge_connectivity();
//	}
//};*/
///*bool diagram::checkirreducibilityforreconnect(list<line>::iterator G_i, list<line>::iterator G_j){
//	//cout<<"** checkconnectivity module **"<<endl;
//	int Nvertices = vertices.size();
//	graph topology;
//
//	vector<vector<int> > adjLtmp(Nvertices);
//	for(int i=0; i<Nvertices; i++) adjLtmp[i].reserve(3*Nvertices/2);
//
//	list<vertex>::iterator ref0 = vertices.begin();
//	list<line>::iterator iter;
//
//	int physical;
//	int iv, jv;
//	for(iter=Gs.begin(); iter!=Gs.end(); ++iter){
//		physical = iter->getphysical();
//
//		if(iter==G_i){
//			iv = distance(ref0,iter->getconnected_vi());
//			jv = distance(ref0,G_j->getconnected_vj());
//		}
//		else if(iter==G_j){
//			iv = distance(ref0,iter->getconnected_vi());
//			jv = distance(ref0,G_i->getconnected_vj());
//		}
//		else{
//			iv = distance(ref0,iter->getconnected_vi());
//			jv = distance(ref0,iter->getconnected_vj());
//		}
//		if((iv!=jv) && physical){
//			adjLtmp[iv].push_back(jv);
//			adjLtmp[jv].push_back(iv);
//		}
//	}
//	for(iter=Ws.begin(); iter!=Ws.end(); ++iter){
//		physical = iter->getphysical();
//
//		iv = distance(ref0,iter->getconnected_vi());
//		jv = distance(ref0,iter->getconnected_vj());
//		if((iv!=jv) && physical){
//			adjLtmp[iv].push_back(jv);
//			adjLtmp[jv].push_back(iv);
//		}
//	}
//	//cout<<"adjL info"<<endl;
//	//vector<int>::iterator viter;
//	//for(int i=0; i<Nvertices; i++){
//	//	for(viter=adjLtmp[i].begin(); viter!=adjLtmp[i].end(); viter++)
//	//		cout<<(*viter);
//	//	cout<<endl;
//	//}
//
//	topology.set(Nvertices,adjLtmp);
//
//	int OneTwoConnectivity = topology.one_two_edge_connectivity();
//	if(OneTwoConnectivity==0) return false;
//	else if(OneTwoConnectivity==1){
//		if(Nvertices<3) return true;
//		else return false;
//	}
//	else return true;
//
//	//if(topology.edge_connectivity()){
//	//	if(Nvertices<3) return true;
//	//	else return topology.two_edge_connectivity();
//	//}
//	//else return false;
//};
//bool diagram::checkirreducibilityforinterswap(list<line>::iterator target){
//	//cout<<"** checkconnectivity module **"<<endl;
//	int Nvertices = vertices.size();
//
//	graph topology;
//
//	vector<vector<int> > adjLtmp(Nvertices);
//	for(int i=0; i<Nvertices; i++) adjLtmp[i].reserve(3*Nvertices/2);
//
//	list<vertex>::iterator ref0 = vertices.begin();
//	list<line>::iterator iter;
//
//	int isG, physical;
//	int iv, jv;
//	for(iter=Gs.begin(); iter!=Gs.end(); ++iter){
//		if(iter!=target){
//			iv = distance(ref0,iter->getconnected_vi());
//			jv = distance(ref0,iter->getconnected_vj());
//			if(iv!=jv){
//				adjLtmp[iv].push_back(jv);
//				adjLtmp[jv].push_back(iv);
//			}
//		}
//	}
//	for(iter=Ws.begin(); iter!=Ws.end(); ++iter){
//		if(iter!=target){
//			iv = distance(ref0,iter->getconnected_vi());
//			jv = distance(ref0,iter->getconnected_vj());
//			if(iv!=jv){
//				adjLtmp[iv].push_back(jv);
//				adjLtmp[jv].push_back(iv);
//			}
//		}
//	}
//
//	//cout<<"adjL info"<<endl;
//	//vector<int>::iterator viter;
//	//for(int i=0; i<Nvertices; i++){
//	//	for(viter=adjLtmp[i].begin(); viter!=adjLtmp[i].end(); viter++)
//	//		cout<<(*viter);
//	//	cout<<endl;
//	//}
//
//	topology.set(Nvertices,adjLtmp);
//	if(topology.edge_connectivity()){
//		if(Nvertices<3) return true;
//		else return topology.two_edge_connectivity();
//	}
//	else return false;
//
//	//topology.set(Nvertices,adjLtmp);
//	//return topology.two_edge_connectivity();
//};*/
//
////bool diagram::checkcompactness(const bool& isG_i){
////	three_edge topology;
////	int Nvertices = vertices.size();
////	if(Nvertices<3) return true;
////
////	vector<vector<int> > adjLtmp;
////
////	vector<int> supervertices(Nvertices);
////	vector<bool> discovered(Nvertices,false);
////
////	list<vertex>::iterator ref0 = vertices.begin();
////	list<line>::iterator iter;
////
////	bool physical;
////	int iv, jv;
////	int Nsupervertice;
////	if(isG_i){
////		//cout<<"** check G-compactness module **"<<endl;
////		int isv_run = 0;
////		for(iter=Ws.begin(); iter!=Ws.end(); ++iter){
////			physical = iter->getphysical();
////
////			iv = distance(ref0,iter->getconnected_vi());
////			jv = distance(ref0,iter->getconnected_vj());
////
////			supervertices[iv] = isv_run;
////			supervertices[jv] = isv_run;
////			isv_run++;
////		}
////
////		Nsupervertice = isv_run;
////		adjLtmp.resize(Nsupervertice);
////		for(int i=0; i<Nsupervertice; i++) adjLtmp[i].reserve(vertices.size());
////
////		for(iter=Gs.begin(); iter!=Gs.end(); ++iter){
////			physical = iter->getphysical();
////
////			iv = supervertices[distance(ref0,iter->getconnected_vi())];
////			jv = supervertices[distance(ref0,iter->getconnected_vj())];
////
////			if(iv!=jv){
////				adjLtmp[iv].push_back(jv);
////				adjLtmp[jv].push_back(iv);
////			}
////		}
////	}
////	else{
////		//cout<<"** check W-compactness module **"<<endl;
////		int isv_run = 0;
////		list<vertex>::iterator outer;
////		list<vertex>::iterator inner;
////		for(outer=ref0; outer!=vertices.end(); ++outer){
////			//cout<<*outer<<endl;
////			if(!discovered[distance(ref0,outer)]){
////				inner = outer;
////				do{
////					iv = distance(ref0,inner);
////					discovered[iv] = true;
////					supervertices[iv] = isv_run;
////					iter = inner->getconnected_Gout();
////					inner = iter->getconnected_vj();
////				} while(inner!=outer);
////				isv_run++;
////			}
////		}
////
////		Nsupervertice = isv_run;
////		adjLtmp.resize(Nsupervertice);
////		for(int i=0; i<Nsupervertice; i++) adjLtmp[i].reserve(vertices.size()/2);
////
////		for(iter=Ws.begin(); iter!=Ws.end(); ++iter){
////			//physical = iter->getphysical();
////
////			iv = supervertices[distance(ref0,iter->getconnected_vi())];
////			jv = supervertices[distance(ref0,iter->getconnected_vj())];
////
////			if(iv!=jv){
////				adjLtmp[iv].push_back(jv);
////				adjLtmp[jv].push_back(iv);
////			}
////		}
////	}
////
////	topology.set(Nsupervertice,adjLtmp);
////
////	//if(Nsupervertice!=0){ if(!(topology.three_edge_connectivity())){
////	//printconfiguration(cout);
////	//cout<<"Nsupervertice = "<<Nsupervertice<<endl<<flush;
////	//for(auto i: adjLtmp){
////	//	for(auto j: i) cout<<j;
////	//	cout<<endl;
////	//}
////	//} }
////
////	if(Nsupervertice==0) return true;
////	else return topology.three_edge_connectivity();
////};
/////*bool diagram::checkGcompactness(vector<vector<int> > GadjL_i, vector<int> WmadjL_i,vector<vector<int> > WadjL_i){
////	//cout<<"** check G-compactness module **"<<endl;
////	int Nvertices = GadjL_i.size();
////	//if(GadjL_i.size()!=WadjL_i.size()){
////	//	cout<<" size of GadjL and WadjL is different! "<<endl;
////	//}
////	//cout<<"Nvertices = "<<Nvertices<<endl;
////	if(Nvertices<3) return true;
////
////	vector<vector<int> > adjLtmp;
////
////	vector<int> supervertices(Nvertices);
////	vector<bool> discovered(Nvertices,false);
////
////	bool physical;
////	int iv, jv;
////	int Nsupervertice;
////	int isv_run = 0;
////	for(int i=0; i<WmadjL_i.size(); i++){
////		discovered[WmadjL_i[i]] = true;
////		supervertice[WmadjL_i[i]] = isv_run;
////	}
////	if(WmadjL_i.size()>0) isv_run++;
////	for(int i=0; i<WadjL_i.size(); i++){
////		if(!discovered[i] ){
////			discovered[i] = true;
////			supervertice[i] = isv_run;
////			for(int j=0; j<WadjL_i[i].size(); j++){
////				discovered[WadjL_i[i][j]] = true;
////				supervertice[WadjL_i[i][j]] = isv_run;
////			}
////			if(WadjL_i[i].size()>0) isv_run++;
////		}
////	}
////
////	Nsupervertice = isv_run;
////	adjLtmp.resize(Nsupervertice);
////
////	//cout<<"supervertice"<<endl;
////	//for(int i=0; i<Nvertices; i++)
////	//	cout<<supervertice[i]<<" ";
////	//cout<<endl;
////	//for(int i=0; i<Nsupervertice; i++){
////	//	adjLtmp[i].reserve(Nvertices);
////	//}
////
////	//cout<<"Nsupervertice = "<<Nsupervertice<<endl;
////
////	for(int i=0; i<GadjL_i.size(); i++){
////		iv = supervertice[i];
////		for(int j=0; j<GadjL_i[i].size(); j++){
////			jv = supervertice[GadjL_i[i][j]];
////			if(iv!=jv)
////				adjLtmp[iv].push_back(jv);
////		}
////	}
////
////	//cout<<"adjL info"<<endl;
////	//vector<int>::iterator viter;
////	//for(int i=0; i<Nsupervertice; i++){
////	//	for(viter=adjLtmp[i].begin(); viter!=adjLtmp[i].end(); viter++)
////	//		cout<<(*viter);
////	//	cout<<endl;
////	//}
////
////	//for(auto i: adjLtmp){
////	//	for(auto j: i) cout<<j;
////	//	cout<<endl;
////	//}
////
////	//for(iter=Gs.begin(); iter!=Gs.end(); ++iter){
////	//	physical = iter->getphysical();
////
////	//	iv = supervertice[distance(ref0,iter->getconnected_vi())];
////	//	jv = supervertice[distance(ref0,iter->getconnected_vj())];
////
////	//	if((iv!=jv) && physical){
////	//		adjLtmp[iv].push_back(jv);
////	//		adjLtmp[jv].push_back(iv);
////	//	}
////	//}
////
////	three_edge topology;
////	topology.set(Nsupervertice,adjLtmp);
////
////	//if(Nsupervertice!=0){ if(!(topology.three_edge_connectivity())){
////	//printconfiguration(cout);
////	//cout<<"Nsupervertice = "<<Nsupervertice<<endl<<flush;
////	//for(auto i: adjLtmp){
////	//	for(auto j: i) cout<<j;
////	//	cout<<endl;
////	//}
////	//} }
////
////	if(Nsupervertice==0) return true;
////	else return topology.three_edge_connectivity();
////};
////bool diagram::checkWcompactness(vector<vector<int> > GadjL_i, vector<int> WmadjL_i, vector<vector<int> > WadjL_i){
////	//cout<<"** check W-compactness module **"<<endl;
////	int Nvertices = WadjL_i.size();
////	if(Nvertices<3) return true;
////
////	vector<vector<int> > adjLtmp;
////
////	vector<int> supervertices(Nvertices);
////	vector<bool> discovered(Nvertices,false);
////
////	bool physical;
////	int iv, jv;
////	int Nsupervertice;
////	int isv_run = 0;
////
////	//cout<<" GadjL_i "<<endl;
////	//for(auto i: GadjL_i){
////	//	for(auto j: i) cout<<j;
////	//	cout<<endl;
////	//}
////
////	graph Ggraph;
////	Ggraph.set(Nvertices,GadjL_i);
////	supervertice = Ggraph.assignsupervertices(Nsupervertice);
////	//cout<<"Nsupervertice = "<<Nsupervertice<<endl<<flush;
////	//cout<<"supervertices index"<<endl;
////	//for(int i=0; i<supervertices.size(); i++)
////	//	cout<<supervertice[i]<<" ";
////	//cout<<endl;
////
////
////	adjLtmp.resize(Nsupervertice);
////	for(int i=0; i<Nsupervertice; i++) adjLtmp[i].reserve(vertices.size()/2);
////
////	//cout<<"add measuring line to the adjacent vertices list"<<endl<<flush;
////	//cout<<" WmadjL_i "<<endl;
////	//for(auto i: WmadjL_i) cout<<i;
////	//cout<<endl;
////	if(WmadjL_i.size()!=0){
////		iv = supervertice[WmadjL_i[0]];
////		jv = supervertice[WmadjL_i[1]];
////		if(iv!=jv){
////			adjLtmp[iv].push_back(jv);
////			adjLtmp[jv].push_back(iv);
////		}
////	}
////
////	//cout<<" WadjL_i "<<endl;
////	//for(auto i: WadjL_i){
////	//	for(auto j: i) cout<<j;
////	//	cout<<endl;
////	//}
////	//cout<<endl;
////	//cout<<"add physical line to the adjacent vertices list"<<endl<<flush;
////	for(int i=0; i<WadjL_i.size(); i++){
////		iv = supervertice[i];
////		for(int j=0; j<WadjL_i[i].size(); j++){
////			jv = supervertice[WadjL_i[i][j]];
////
////			if(iv!=jv){
////				adjLtmp[iv].push_back(jv);
////				//adjLtmp[jv].push_back(iv);
////			}
////		}
////	}
////
////	//for(auto i: adjLtmp){
////	//	for(auto j: i) cout<<j;
////	//	cout<<endl;
////	//}
////
////	three_edge topology;
////	topology.set(Nsupervertice,adjLtmp);
////
////	//if(Nsupervertice!=0){ if(!(topology.three_edge_connectivity())){
////	//printconfiguration(cout);
////	//} }
////
////	if(Nsupervertice==0) return true;
////	else return topology.three_edge_connectivity();
////};*/
////
/////*bool diagram::brutecheckcompactness(const bool& isG){
////	//cout<<"** bruteforce check compactness module **"<<endl;
////	int Nvertices = vertices.size();
////	if(Nvertices<3) return true;
////	if(isG){
////		//cout<<" * G-part compactness module *"<<endl;
////		graph topology;
////
////		vector<vector<int> > adjLtmp(Nvertices);
////		vector<vector<int> > adjLtype(Nvertices);
////
////		list<vertex>::iterator ref0 = vertices.begin();
////		list<line>::iterator iter;
////		list<line>::iterator iter_cut1;
////		list<line>::iterator iter_cut2;
////
////		bool connectivity;
////
////		int physical;
////		int iv, jv;
////		for(iter_cut1 = Gs.begin(); iter_cut1!=Gs.end(); ++iter_cut1)
////		for(iter_cut2 = next(iter_cut1); iter_cut2!=Gs.end(); ++iter_cut2){
////			if( iter_cut1->getphysical() && iter_cut2->getphysical()){
////				for(int i=0; i<Nvertices; i++){
////					adjLtmp[i].clear();
////					adjLtype[i].clear();
////				}
////				for(iter=Gs.begin(); iter!=Gs.end(); ++iter){
////					if(iter!=iter_cut1 && iter!=iter_cut2){
////						physical = iter->getphysical();
////						//isG = iter->getisG();
////
////						iv = distance(ref0,iter->getconnected_vi());
////						jv = distance(ref0,iter->getconnected_vj());
////						if(iv!=jv){
////							adjLtype[iv].push_back(0);
////							adjLtype[jv].push_back(0);
////
////							adjLtmp[iv].push_back(jv);
////							adjLtmp[jv].push_back(iv);
////						}
////					}
////				}
////				for(iter=Ws.begin(); iter!=Ws.end(); ++iter){
////					//physical = iter->getphysical();
////
////					iv = distance(ref0,iter->getconnected_vi());
////					jv = distance(ref0,iter->getconnected_vj());
////					if(iv!=jv){
////						adjLtype[iv].push_back(1);
////						adjLtype[jv].push_back(1);
////
////						adjLtmp[iv].push_back(jv);
////						adjLtmp[jv].push_back(iv);
////					}
////				}
////				topology.set(Nvertices,adjLtmp,adjLtype);
////				connectivity =topology.edge_connectivity();
////				if( !connectivity ) return false;
////				//return topology.two_edge_connectivity(3);
////			}
////		}
////
////		return true;
////
////		//cout<<"adjL info"<<endl;
////		//vector<int>::iterator viter;
////		//for(int i=0; i<Nvertices; i++){
////		//	for(viter=adjLtmp[i].begin(); viter!=adjLtmp[i].end(); viter++)
////		//		cout<<(*viter);
////		//	cout<<endl;
////		//}
////		//cout<<"adjL type info"<<endl;
////		//for(int i=0; i<Nvertices; i++){
////		//	for(viter=adjLtype[i].begin(); viter!=adjLtype[i].end(); viter++)
////		//		cout<<(*viter);
////		//	cout<<endl;
////		//}
////
////		//topology.set(Nvertices,adjLtmp,adjLtype);
////		//return topology.two_edge_connectivity(3);
////	}
////	else{
////		//cout<<" * UW-part compactness module *"<<endl;
////		graph topology;
////
////		vector<vector<int> > adjLtmp(Nvertices);
////		vector<vector<int> > adjLtype(Nvertices);
////
////		list<vertex>::iterator ref0 = vertices.begin();
////		list<line>::iterator iter;
////		list<line>::iterator iter_cut1;
////		list<line>::iterator iter_cut2;
////
////		int NWs = Ws.size();
////		//int NWs = Us.size()+Ws.size();
////		//for(iter=Us.begin(); iter!=Us.end(); ++iter) Ws.push_back(*iter);
////		//for(iter=Ws.begin(); iter!=Ws.end(); ++iter) Ws.push_back(*iter);
////
////		bool connectivity;
////
////		int physical;
////		int iv, jv;
////		for(int iUW=0; iUW<NWs; iUW++) for(int jUW=iUW+1; jUW<NWs; jUW++){
////			iter_cut1 = next(Ws.begin(),iUW);
////			iter_cut2 = next(Ws.begin(),jUW);
////			if( iter_cut1->getphysical() && iter_cut2->getphysical()){
////				for(int i=0; i<Nvertices; i++){
////					adjLtmp[i].clear();
////					adjLtype[i].clear();
////				}
////				for(iter=Gs.begin(); iter!=Gs.end(); ++iter){
////					physical = iter->getphysical();
////
////					iv = distance(ref0,iter->getconnected_vi());
////					jv = distance(ref0,iter->getconnected_vj());
////					if((iv!=jv) && physical){
////						adjLtype[iv].push_back(0);
////						adjLtype[jv].push_back(0);
////
////						adjLtmp[iv].push_back(jv);
////						adjLtmp[jv].push_back(iv);
////					}
////				}
////				for(iter=Ws.begin(); iter!=Ws.end(); ++iter){
////					if(iter!=iter_cut1 && iter!=iter_cut2){
////						iv = distance(ref0,iter->getconnected_vi());
////						jv = distance(ref0,iter->getconnected_vj());
////						if(iv!=jv){
////							adjLtype[iv].push_back(1);
////							adjLtype[jv].push_back(1);
////
////							adjLtmp[iv].push_back(jv);
////							adjLtmp[jv].push_back(iv);
////						}
////					}
////				}
////
////				topology.set(Nvertices,adjLtmp,adjLtype);
////				connectivity = topology.edge_connectivity();
////				if( !connectivity ) return false;
////				//return topology.two_edge_connectivity(3);
////			}
////		}
////
////		return true;
////
////		//cout<<"adjL info"<<endl;
////		//vector<int>::iterator viter;
////		//for(int i=0; i<Nvertices; i++){
////		//	for(viter=adjLtmp[i].begin(); viter!=adjLtmp[i].end(); viter++)
////		//		cout<<(*viter);
////		//	cout<<endl;
////		//}
////		//cout<<"adjL type info"<<endl;
////		//for(int i=0; i<Nvertices; i++){
////		//	for(viter=adjLtype[i].begin(); viter!=adjLtype[i].end(); viter++)
////		//		cout<<(*viter);
////		//	cout<<endl;
////		//}
////
////		//topology.set(Nvertices,adjLtmp,adjLtype);
////		//return topology.two_edge_connectivity(3);
////	}
////};*/
////
////bool diagram::check_compactness_for_remove(list<line>::iterator Wremove){
////	//cout<<"** check compactness module for remove update **"<<endl;
////	//printconfiguration(cout);
////	//cout<<"removing W-line : "<<*Wremove<<endl;
////	int Nvertices = vertices.size();
////	if(Nvertices<5) return true;
////	else{
////		vector<vector<int> > GadjLtmp, WadjLtmp;
////
////		vector<int> supervertices(Nvertices);
////
////		list<vertex>::iterator ref0 = vertices.begin();
////
////		list<line>::iterator iter;
////
////		list<vertex>::iterator vi_iter = Wremove->getconnected_vi();
////		list<vertex>::iterator vj_iter = Wremove->getconnected_vj();
////
////		list<line>::iterator Gi_in = vi_iter->getconnected_Gin();
////		list<line>::iterator Gj_in = vj_iter->getconnected_Gin();
////		list<line>::iterator Gi_out = vi_iter->getconnected_Gout();
////		list<line>::iterator Gj_out = vj_iter->getconnected_Gout();
////
////		int iv_remove = distance(ref0,vi_iter);
////		int jv_remove = distance(ref0,vj_iter);
////
////		int iv_prev = distance(ref0,Gi_in->getconnected_vi());
////		int jv_prev = distance(ref0,Gj_in->getconnected_vi());
////		int iv_next = distance(ref0,Gi_out->getconnected_vj());
////		int jv_next = distance(ref0,Gj_out->getconnected_vj());
////
////
////		bool physical;
////		int iv, jv;
////		int Nsupervertices;
////		int isv_run;
////
////		//cout<<"G topology construct"<<endl;
////		three_edge Gtopology;
////
////		isv_run = 0;
////		for(iter=Ws.begin(); iter!=Ws.end(); ++iter){
////			if(iter!=Wremove){
////				iv = distance(ref0,iter->getconnected_vi());
////				jv = distance(ref0,iter->getconnected_vj());
////
////				supervertices[iv] = isv_run;
////				supervertices[jv] = isv_run;
////				isv_run++;
////			}
////		}
////
////		Nsupervertices = isv_run;
////		GadjLtmp.resize(Nsupervertices);
////		for(int i=0; i<Nsupervertices; i++) GadjLtmp[i].reserve(Nvertices);
////
////		//cout<<"G push_back"<<endl<<flush;
////		iv = supervertices[iv_prev]; jv = supervertices[iv_next];
////		if(iv!=jv){
////			GadjLtmp[iv].push_back(jv);
////			GadjLtmp[jv].push_back(iv);
////		}
////		iv = supervertices[jv_prev]; jv = supervertices[jv_next];
////		if(iv!=jv){
////			GadjLtmp[iv].push_back(jv);
////			GadjLtmp[jv].push_back(iv);
////		}
////		for(iter=Gs.begin(); iter!=Gs.end(); ++iter){
////			iv = distance(ref0,iter->getconnected_vi());
////			jv = distance(ref0,iter->getconnected_vj());
////			if((iv!=iv_remove) && (iv!=jv_remove) && (jv!=iv_remove) && (jv!=jv_remove) && (iv!=jv)){
////				GadjLtmp[supervertices[iv]].push_back(supervertices[jv]);
////				GadjLtmp[supervertices[jv]].push_back(supervertices[iv]);
////			}
////		}
////
////		//cout<<endl;
////		//for(auto i: GadjLtmp){
////		//	for(auto j: i) cout<<j;
////		//	cout<<endl;
////		//}
////
////		Gtopology.set(GadjLtmp.size(),GadjLtmp);
////		bool Gcompactness = Gtopology.three_edge_connectivity();
////
////		return Gcompactness;
////	}
////};
/////*bool diagram::checkcompactnessforreconnect(list<line>::iterator G_i, list<line>::iterator G_j){
////	int Nvertices = vertices.size();
////	if(Nvertices<3) return true;
////	else{
////		vector<vector<int> > GadjLtmp, WadjLtmp;;
////
////		vector<int> supervertices(Nvertices);
////
////		list<vertex>::iterator ref0 = vertices.begin();
////		list<line>::iterator iter;
////
////		bool physical;
////		int iv, jv;
////		int Nsupervertice;
////		int isv_run;
////
////		//cout<<"** check G-compactness module **"<<endl;
////		three_edge Gtopology;
////
////		isv_run = 0;
////		for(iter=Ws.begin(); iter!=Ws.end(); ++iter){
////			physical = iter->getphysical();
////
////			iv = distance(ref0,iter->getconnected_vi());
////			jv = distance(ref0,iter->getconnected_vj());
////
////			supervertices[iv] = isv_run;
////			supervertices[jv] = isv_run;
////			isv_run++;
////		}
////
////		Nsupervertice = isv_run;
////		GadjLtmp.resize(Nsupervertice);
////		for(int i=0; i<Nsupervertice; i++) GadjLtmp[i].reserve(Nvertices);
////
////		for(iter=Gs.begin(); iter!=Gs.end(); ++iter){
////			physical = iter->getphysical();
////
////			if(iter==G_i){
////				iv = supervertices[distance(ref0,iter->getconnected_vi())];
////				jv = supervertices[distance(ref0,G_j->getconnected_vj())];
////			}
////			else if(iter==G_j){
////				iv = supervertices[distance(ref0,iter->getconnected_vi())];
////				jv = supervertices[distance(ref0,G_i->getconnected_vj())];
////			}
////			else{
////				iv = supervertices[distance(ref0,iter->getconnected_vi())];
////				jv = supervertices[distance(ref0,iter->getconnected_vj())];
////			}
////			if((iv!=jv) && physical){
////				GadjLtmp[iv].push_back(jv);
////				GadjLtmp[jv].push_back(iv);
////			}
////		}
////
////		Gtopology.set(GadjLtmp.size(),GadjLtmp);
////		bool Gcompactness = Gtopology.three_edge_connectivity();
////
////		//cout<<"** check W-compactness module **"<<endl;
////		three_edge Wtopology;
////
////		isv_run = 0;
////		list<vertex>::iterator vi_iter, vj_iter;
////		vector<bool> discovered(vertices.size(),false);
////		for(vi_iter=ref0; vi_iter!=vertices.end(); ++vi_iter){
////			iv = distance(ref0,vi_iter);
////			if( !discovered[iv] ){
////				vj_iter = vi_iter;
////				do{
////					jv = distance(ref0,vj_iter);
////					discovered[jv] = true;
////					supervertices[jv] = isv_run;
////					iter = vj_iter->getconnected_Gout();
////					if(iter==G_i)
////						vj_iter = G_j->getconnected_vj();
////					else if(iter==G_j)
////						vj_iter = G_i->getconnected_vj();
////					else
////						vj_iter = iter->getconnected_vj();
////				} while(vj_iter!=vi_iter);
////				isv_run++;
////			}
////		}
////
////		Nsupervertice = isv_run;
////		WadjLtmp.resize(Nsupervertice);
////		for(int i=0; i<Nsupervertice; i++) WadjLtmp[i].reserve(Nvertices/2);
////
////		for(iter=Ws.begin(); iter!=Ws.end(); ++iter){
////			iv = supervertices[distance(ref0,iter->getconnected_vi())];
////			jv = supervertices[distance(ref0,iter->getconnected_vj())];
////
////			if(iv!=jv){
////				WadjLtmp[iv].push_back(jv);
////				WadjLtmp[jv].push_back(iv);
////			}
////		}
////
////		Wtopology.set(WadjLtmp.size(),WadjLtmp);
////		bool Wcompactness = Wtopology.three_edge_connectivity();
////
////		return Gcompactness && Wcompactness;
////
////		//if(Nsupervertice!=0){ if(!(topology.three_edge_connectivity())){
////		//printconfiguration(cout);
////		//cout<<"Nsupervertice = "<<Nsupervertice<<endl<<flush;
////		//for(auto i: adjLtmp){
////		//	for(auto j: i) cout<<j;
////		//	cout<<endl;
////		//}
////		//} }
////	}
////
////};*/
////bool diagram::check_compactness_for_swap(list<line>::iterator target){
////	//cout<<"** check compactness module for interswap update **"<<endl;
////	//printconfiguration(cout);
////	//cout<<"target measuringline : "<<*target<<endl<<flush;
////
////	int Nvertices = vertices.size();
////	if(Nvertices<3) return true;
////	else{
////		vector<vector<int> > GadjLtmp, WadjLtmp;
////
////		vector<int> supervertices(Nvertices);
////
////		list<vertex>::iterator ref0 = vertices.begin();
////		list<line>::iterator iter;
////
////		bool physical;
////		int iv, jv;
////		int Nsupervertice;
////		int isv_run;
////
////		//cout<<"** check G-compactness module **"<<endl;
////		isv_run = 0;
////		for(iter=Ws.begin(); iter!=Ws.end(); ++iter){
////			iv = distance(ref0,iter->getconnected_vi());
////			jv = distance(ref0,iter->getconnected_vj());
////
////			supervertices[iv] = isv_run;
////			supervertices[jv] = isv_run;
////			isv_run++;
////		}
////
////		Nsupervertice = isv_run;
////		GadjLtmp.resize(Nsupervertice);
////		for(int i=0; i<Nsupervertice; i++) GadjLtmp[i].reserve(Nvertices);
////
////		for(iter=Gs.begin(); iter!=Gs.end(); ++iter){
////			physical = !(iter==target);
////
////			iv = supervertices[distance(ref0,iter->getconnected_vi())];
////			jv = supervertices[distance(ref0,iter->getconnected_vj())];
////			if((iv!=jv) && physical){
////				GadjLtmp[iv].push_back(jv);
////				GadjLtmp[jv].push_back(iv);
////			}
////		}
////
////		//for(auto i: GadjLtmp){
////		//	for(auto j: i) cout<<j;
////		//	cout<<endl;
////		//}
////
////		three_edge Gtopology;
////		Gtopology.set(GadjLtmp.size(),GadjLtmp);
////		//bool Gcompactness = Gtopology.three_edge_connectivity();
////		return Gtopology.three_edge_connectivity();
////
////		//cout<<"** check W-compactness module **"<<endl;
////		/*isv_run = 0;
////		list<vertex>::iterator vi_iter, vj_iter;
////		vector<bool> discovered(vertices.size(),false);
////		for(vi_iter=ref0; vi_iter!=vertices.end(); ++vi_iter){
////			iv = distance(ref0,vi_iter);
////			if( !discovered[iv] ){
////				vj_iter = vi_iter;
////				do{
////					jv = distance(ref0,vj_iter);
////					discovered[jv] = true;
////					supervertice[jv] = isv_run;
////					iter = vj_iter->getconnected_Gout();
////					vj_iter = iter->getconnected_vj();
////				} while(vj_iter!=vi_iter);
////				isv_run++;
////			}
////		}
////
////		Nsupervertice = isv_run;
////		WadjLtmp.resize(Nsupervertice);
////		for(int i=0; i<Nsupervertice; i++) WadjLtmp[i].reserve(Nvertices/2);
////
////		for(iter=Ws.begin(); iter!=Ws.end(); ++iter){
////			iv = supervertice[distance(ref0,iter->getconnected_vi())];
////			jv = supervertice[distance(ref0,iter->getconnected_vj())];
////
////			if(iv!=jv){
////				WadjLtmp[iv].push_back(jv);
////				WadjLtmp[jv].push_back(iv);
////			}
////		}
////
////		//for(auto i: WadjLtmp){
////		//	for(auto j: i) cout<<j;
////		//	cout<<endl;
////		//}
////
////
////		three_edge Wtopology;
////		Wtopology.set(WadjLtmp.size(),WadjLtmp);
////		bool Wcompactness = Wtopology.three_edge_connectivity();*/
////
////		//return Gcompactness && Wcompactness;
////	}
////};
/////*bool diagram::checkbubble(){
////	if(vertices.size()<3) return false;
////	vector<bool> discovered(vertices.size(),false);
////	
////	list<vertex>::iterator v_ref = vertices.begin();
////	list<vertex>::iterator vi_iter;
////	list<vertex>::iterator vj_iter;
////
////	list<line>::iterator Gi_iter, Gj_iter;
////	list<line>::iterator UWi_iter, UWj_iter;
////
////	int count;
////	bool physicalpart;
////
////	for(int i=0; i<discovered.size(); i++){
////		if(!discovered[i]){
////			count = 0;
////			physicalpart = true;
////
////			vi_iter = next(vertices.begin(),i);
////			vj_iter = vi_iter;
////			discovered[i] = true;
////			Gi_iter = vi_iter->getconnected_Gout();
////			do{
////				Gj_iter = vj_iter->getconnected_Gout();
////				vj_iter = Gj_iter->getconnected_vj();
////				discovered[distance(v_ref,vj_iter)] = true;
////				count++;
////				if(!(Gj_iter->getphysical())) physicalpart = false;
////			} while(vi_iter!=vj_iter);
////			if(count==2 && physicalpart==true){
////				return true;
////			}
////		}
////	}
////	return false;
////};
////bool diagram::checkbubble(const vector<vector<int> >& GadjL_i){
////	int Nvertices = GadjL_i.size();
////	if(Nvertices<3) return false;
////
////	vector<bool> discovered(Nvertices);
////	for(int i=0; i<Nvertices; i++){
////		for(vector<bool>::iterator iter=discovered.begin(); iter!=discovered.end(); ++iter) *iter = false;
////		for(int j=0; j<GadjL_i[i].size(); j++){
////			if(!discovered[GadjL_i[i][j]]) discovered[GadjL_i[i][j]] = true;
////			else return true;
////		}
////	}
////	return false;
////};
////bool diagram::checkdiagramsign(){
////	//cout<<"** checkconnectivity module **"<<endl;
////	int Nfermionloop = fermionicloops.size();
////	if( measuringline->gettype()==0 ) Nfermionloop--;
////
////	double Gproduct = 1;
////	list<fermionicloop>::iterator iter;
////	for(iter=fermionicloops.begin(); iter!=fermionicloops.end(); ++iter){ 
////		if(
////			measuringline->gettype()==0
////			&&
////			measuringline->getconnected_vi()->getconnected_loop()==iter
////		){
////			Gproduct *= iter->getvaluerh();
////		}
////		else{
////			Gproduct *= (iter->getvaluerh() + iter->getvaluelh()); 
////		}
////	}
////	//list<line>::iterator iter;
////	//for(iter=Gs.begin(); iter!=Gs.end(); ++iter){ Gproduct *= (iter->getvalue()); }
////
////	double Wproduct = SetUWInter.getvalue();
////
////	int signdirect = (1-2*(Nfermionloop%2))
////			* static_cast<int>(Wproduct>0 ? 1 : -1)
////			* static_cast<int>(Gproduct>0 ? 1 : -1);
////
////	if(signdirect==sign) return true;
////	else{
////		cout<<"***** ALERT! WRONG DIAGRAM SIGN *****"<<endl;
////		printconfiguration(cout);
////		return false;
////	}
////}*/
/////*bool diagram::checkdiagramsign(){
////	//cout<<"** checkconnectivity module **"<<endl;
////	int Nvertices = vertices.size();
////	graph topology;
////
////	vector<vector<int> > adjLtmp(Nvertices);
////	vector<vector<int> > adjLtype(Nvertices);
////	//for(int i=0; i<Nvertices; i++){
////	//	adjLtmp[i].reserve(Nvertices);
////	//	adjLtype[i].reserve(Nvertices);
////	//}
////
////	list<vertex>::iterator ref0 = vertices.begin();
////	list<line>::iterator iter;
////
////	int physical;
////	int iv, jv;
////	for(iter=Gs.begin(); iter!=Gs.end(); ++iter){
////		physical = iter->getphysical();
////
////		iv = distance(ref0,iter->getconnected_vi());
////		jv = distance(ref0,iter->getconnected_vj());
////		if((iv!=jv) && physical){
////			adjLtype[iv].push_back(0);
////			adjLtype[jv].push_back(0);
////
////			adjLtmp[iv].push_back(jv);
////			adjLtmp[jv].push_back(iv);
////		}
////	}
////
////	//cout<<"adjL info"<<endl;
////	//vector<int>::iterator viter;
////	//for(int i=0; i<Nvertices; i++){
////	//	for(viter=adjLtmp[i].begin(); viter!=adjLtmp[i].end(); viter++)
////	//		cout<<(*viter);
////	//	cout<<endl;
////	//}
////
////	topology.set(Nvertices,adjLtmp,adjLtype);
////	int Nfermionloop = topology.getNloop() + Nselfloop;
////	if(measuringline->getconnected_vi()==measuringline->getconnected_vj()) Nfermionloop--;
////
////	double Gproduct = 1;
////	for(iter=Gs.begin(); iter!=Gs.end(); ++iter){ Gproduct *= (iter->getvalue()); }
////	double Wproduct = 1;
////	for(iter=Ws.begin(); iter!=Ws.end(); ++iter){ Wproduct *= (iter->getvalue()); }
////	//for(iter=Us.begin(); iter!=Us.end(); ++iter){ Wproduct *= (iter->getvalue()); }
////	//for(iter=Ws.begin(); iter!=Ws.end(); ++iter){ Wproduct *= (iter->getvalue()); }
////
////	int signdirect = (1-2*(Nfermionloop%2))
////			// *(1-2*(Norder%2))
////			* static_cast<int>(Wproduct>0 ? 1 : -1)
////			* static_cast<int>(Gproduct>0 ? 1 : -1);
////
////	//cout<<"Nfermionloop: "<<Nfermionloop<<endl;
////	//cout<<"Norder: "<<Norder<<endl;
////	//cout<<"Gproduct: "<<Gproduct<<endl;
////	//cout<<"signdirect: "<<signdirect<<endl;
////
////	if(signdirect==sign) return true;
////	else{
////		cout<<"***** ALERT! WRONG DIAGRAM SIGN *****"<<endl;
////		printconfiguration(cout);
////		return false;
////	}
////};*/
/////*bool diagram::checkNorder(){
////	int Norder_check = Ws.size();
////	if(measuringline->gettype()!=0) Norder_check--;
////
////	if(Norder==Norder_check && !(Norder>(this->Nordermax))) return true;
////	else{
////		cout<<"***** ALERT! WRONG DIAGRAM ORDER *****"<<endl;
////		printconfiguration(cout);
////		return false;
////	}
////};
////bool diagram::checkUWlines(){
////	int UW_type;
////	list<line>::iterator iter;
////	list<vertex>::iterator vi_iter, vj_iter;
////	for(iter=Ws.begin(); iter!=Ws.end(); ++iter){
////		UW_type = iter->gettype();
////		vi_iter = iter->getconnected_vi();
////		vj_iter = iter->getconnected_vj();
////		if( (UW_type==1) && !(vi_iter->getx()==vj_iter->getx()) )  return false;
////		else if( (UW_type==2) && (vi_iter->getx()==vj_iter->getx()) ) return false;
////	}
////	return true;
////	//if( 
////	//	(UW_type==1) 
////	//	&& !((vi_iter->getflavor()^vj_iter->getflavor()) 
////	//	&& (vi_iter->getx()==vj_iter->getx())) 
////	//) return false;
////};
////bool diagram::checksetUWinters(){
////	setUWinter trial;
////
////	trial.clear();
////	trial.setNfermionicloops( fermionicloops.size() );
////	trial.releaseloopflavor();
////
////	int iloop, jloop;
////	list<line>::iterator UWiter;
////	list<vertex>::iterator vi_iter, vj_iter;
////	for(UWiter=Ws.begin(); UWiter!=Ws.end(); ++UWiter){
////		vi_iter = UWiter->getconnected_vi();
////		vj_iter = UWiter->getconnected_vj();
////		if( vi_iter->getconnected_loop()!=vj_iter->getconnected_loop() ){
////			iloop = distance(fermionicloops.begin(),vi_iter->getconnected_loop());
////			jloop = distance(fermionicloops.begin(),vj_iter->getconnected_loop());
////			trial.push_back(*UWiter,iloop,jloop);
////		}
////	}
////
////	vi_iter = measuringline->getconnected_vi();
////	vj_iter = measuringline->getconnected_vj();
////	iloop = distance(fermionicloops.begin(),vi_iter->getconnected_loop());
////	jloop = distance(fermionicloops.begin(),vj_iter->getconnected_loop());
////	trial.fixloopflavor(iloop,miflavor);
////	trial.fixloopflavor(jloop,mjflavor);
////
////	trial.calvalue();
////
////	if( abs(trial.getvalue()-SetUWInter.getvalue())<1e-10 ){ return true; }
////	else{
////		cout<<"***** ALERT! WRONG SetUWinters value *****"<<endl;
////		cout<<"trial info."<<endl
////			<<trial<<endl;
////		cout<<"SetUWInter info."<<endl
////			<<SetUWInter<<endl;
////		printconfiguration(cout);
////		return false;
////	}
////};
////bool diagram::checkdress(){
////	list<line>::iterator iter;
////	for(iter=Gs.begin(); iter!=Gs.end(); ++iter)
////		if(iter->getdress()!=0) return false;
////	for(iter=Ws.begin(); iter!=Ws.end(); ++iter)
////		if(iter->getdress()!=0) return false;
////
////	return true;
////}
////bool diagram::checkNselfloop(){
////	//cout<<"** checkNselfloop module **"<<endl;
////	list<line>::iterator iter;
////	int Nselfloop_running = 0;
////	for(iter=Gs.begin(); iter!=Gs.end(); ++iter){
////		if(iter->getconnected_vi()==iter->getconnected_vj()) Nselfloop_running++;
////	}
////
////	if(Nselfloop==Nselfloop_running){
////		if((Nselfloop<NselfloopMax+1) || (Ws.size()==1)) return true;
////		else{
////			cout<<"***** ALERT! WRONG Nselfloop *****"<<endl;
////			printconfiguration(cout);
////			return false;
////		}
////	}
////	else{
////		cout<<"***** ALERT! WRONG Nselfloop *****"<<endl;
////		printconfiguration(cout);
////		return false;
////	}
////};
////bool diagram::checkfermionicloop(){
////	//cout<<"** checkNselfloop module **"<<endl;
////	list<fermionicloop>::iterator loop_iter;
////	bool isNverticeodd;
////	double valuerh_tmp, valuelh_tmp;
////	for(loop_iter=fermionicloops.begin(); loop_iter!=fermionicloops.end(); ++loop_iter){
////		if( abs(abs(loop_iter->getvaluerh()) - abs(loop_iter->getvaluelh()))>1.0e-10 ){
////			cout<<"***** ALERT! WRONG fermionic loop values *****"<<endl;
////			printconfiguration(cout);
////			return false;
////		}
////
////		valuerh_tmp = loop_iter->getvaluerh();
////		valuelh_tmp = loop_iter->getvaluelh();
////
////		//isNverticeodd = static_cast<bool>(loop_iter->getvssize()%2);
////		//if(
////		//	(isNverticeodd && abs(valuerh_tmp+valuelh_tmp)>1.0e-10) ||
////		//	((!isNverticeodd) && abs(valuerh_tmp-valuelh_tmp)>1.0e-10)
////		//){
////		//	cout<<"***** ALERT! WRONG fermionic loop values *****"<<endl;
////		//	cout<<"Nvertices of the given loop = "<<loop_iter->getvssize()<<endl;
////		//	cout<<"assigned valuerh = "<<valuerh_tmp<<endl;
////		//	cout<<"assigned valuelh = "<<valuelh_tmp<<endl;
////		//	printconfiguration(cout);
////		//	return false;
////		//}
////
////		loop_iter->calvalue();
////
////		if( 
////			(abs(loop_iter->getvaluerh()-valuerh_tmp)>1.0e-10) || 
////			(abs(loop_iter->getvaluelh()-valuelh_tmp)>1.0e-10)
////		){
////			cout<<"***** ALERT! WRONG fermionic loop values *****"<<endl;
////			cout<<"assigned valuerh = "<<valuerh_tmp<<endl;
////			cout<<"calculated valuerh = "<<loop_iter->getvaluerh()<<endl;
////			cout<<"assigned valuelh = "<<valuelh_tmp<<endl;
////			cout<<"calculated valuelh = "<<loop_iter->getvaluelh()<<endl;
////			printconfiguration(cout);
////			return false;
////		}
////	}
////	return true;
////};*/
////bool diagram::check_measuring_line(){
////	//cout<<"** check measuringline module **"<<endl;
////	if( measuring_line->getphysical() )
////		return false;
////	else 
////		return true;
////
////	/*int mcount = 0;
////	list<line>::iterator liter;
////	for(liter=Gs.begin(); liter!=Gs.end(); ++liter) if( !liter->getphysical() ) mcount++;
////	for(liter=Ws.begin(); liter!=Ws.end(); ++liter) if( !liter->getphysical() ) mcount++;
////
////	bool running = (mcount==1 && !measuringline->getphysical());
////	//cout<<"running = "<<running<<endl<<flush;
////
////	list<fermionicloop>::iterator loopi_iter, loopj_iter;
////	loopi_iter = measuringline->getconnected_vi()->getconnected_loop();
////	loopj_iter = measuringline->getconnected_vj()->getconnected_loop();
////
////	int iloop = distance(fermionicloops.begin(),loopi_iter);
////	int jloop = distance(fermionicloops.begin(),loopj_iter);
////
////	running &= (miflavor==SetUWInter.getloopflavor(iloop));
////	running &= (mjflavor==SetUWInter.getloopflavor(jloop));
////
////	running &= SetUWInter.getisloopflavorfixed(iloop);
////	running &= SetUWInter.getisloopflavorfixed(jloop);
////
////	for(int i=0; i<fermionicloops.size(); i++)
////		if( (i!=iloop && i!=jloop) && SetUWInter.getisloopflavorfixed(i)) running = false;
////	//cout<<"running = "<<running<<endl<<flush;
////
////	if( (measuringline->gettype()==1) && (miflavor==mjflavor) ) running = false;
////	//cout<<"running = "<<running<<endl<<flush;
////
////	//if(mcount==1 && !measuringline->getphysical()) return true;
////	if(running) return true;
////	else{
////		cout<<"***** ALERT! WRONG measuring line *****"<<endl;
////		printconfiguration(cout);
////		return false;
////	}*/
////};
////
////bool diagram::check_ref_line(){
////	//cout<<"** check ref_line module **"<<endl;
////	bool is_ref_line_physical = ref_line->getphysical();
////	bool is_end_point_zero;
////	double tj = ref_line->getconnected_vj()->getx().time;
////	const double eps = 1.0e-10;
////	if( is_ref_line_physical || fabs(tj)>eps ){
////		cout<<"***** ALERT! WRONG ref_line *****"<<endl;
////		printconfiguration(cout);
////		return false;
////	}
////	else return true;
////
////	/*int mcount = 0;
////	list<line>::iterator liter;
////	for(liter=Gs.begin(); liter!=Gs.end(); ++liter) if( !liter->getphysical() ) mcount++;
////	for(liter=Ws.begin(); liter!=Ws.end(); ++liter) if( !liter->getphysical() ) mcount++;
////
////	bool running = (mcount==1 && !measuringline->getphysical());
////	//cout<<"running = "<<running<<endl<<flush;
////
////	list<fermionicloop>::iterator loopi_iter, loopj_iter;
////	loopi_iter = measuringline->getconnected_vi()->getconnected_loop();
////	loopj_iter = measuringline->getconnected_vj()->getconnected_loop();
////
////	int iloop = distance(fermionicloops.begin(),loopi_iter);
////	int jloop = distance(fermionicloops.begin(),loopj_iter);
////
////	running &= (miflavor==SetUWInter.getloopflavor(iloop));
////	running &= (mjflavor==SetUWInter.getloopflavor(jloop));
////
////	running &= SetUWInter.getisloopflavorfixed(iloop);
////	running &= SetUWInter.getisloopflavorfixed(jloop);
////
////	for(int i=0; i<fermionicloops.size(); i++)
////		if( (i!=iloop && i!=jloop) && SetUWInter.getisloopflavorfixed(i)) running = false;
////	//cout<<"running = "<<running<<endl<<flush;
////
////	if( (measuringline->gettype()==1) && (miflavor==mjflavor) ) running = false;
////	//cout<<"running = "<<running<<endl<<flush;
////
////	//if(mcount==1 && !measuringline->getphysical()) return true;
////	if(running) return true;
////	else{
////		cout<<"***** ALERT! WRONG measuring line *****"<<endl;
////		printconfiguration(cout);
////		return false;
////	}*/
////};
////
////
/////*bool diagram::checkdetailedbalance(){
////	bool running = true;
////	running &= checkdetailedbalance_add_remove_vertices();
////	running &= checkdetailedbalance_remove_add_vertices();
////	running &= checkdetailedbalance_reconnect();
////	running &= checkdetailedbalance_interswap();
////	running &= checkdetailedbalance_transformUW();
////	
////	//running &= checkdetailedbalance_eject_absorb_bubble();
////	//running &= checkdetailedbalance_absorb_eject_bubble();
////	return running;
////};
////
////bool diagram::checkdetailedbalance_add_remove_vertices(){
////	//cout<<"** check detailed balance **"<<endl;
////	//cout<<" * addvertice / removevertice *"<<endl;
////
////
////	/////////////////
////	// add vertice //
////	/////////////////
////
////	double probability_add, probability_remove;
////	//if(Norder==(this->Nordermax)){ return true; }
////	if(Ws.size()==(this->Nordermax)){ return true; }
////	else{
////		int NdressChoice = 1;
////		int Ndress = 0;
////
////		int NGs = Gs.size();
////		int iG,jG;
////		select_two_index(NGs,iG,jG);
////
////		int m_type = measuringline->gettype();
////		int NWs_selection;
////		if(m_type==0){ NWs_selection = Ws.size()+1; }
////		else{ NWs_selection = Ws.size(); }
////
////		coordinate xi_n, xj_n;
////
////		int iflavor_n, jflavor_n;
////
////		list<line>::iterator Gi_iter = next(Gs.begin(),iG);
////		list<line>::iterator Gj_iter = next(Gs.begin(),jG);
////
////		int dNselfloop = 0;
////		if(Gi_iter->getconnected_vi()==Gi_iter->getconnected_vj()) dNselfloop--;
////		if(Gj_iter->getconnected_vi()==Gj_iter->getconnected_vj()) dNselfloop--;
////
////		list<vertex>::iterator vi_prev = Gi_iter->getconnected_vi();
////		list<vertex>::iterator vj_prev = Gj_iter->getconnected_vi();
////
////		list<fermionicloop>::iterator loopi_iter = vi_prev->getconnected_loop();
////		list<fermionicloop>::iterator loopj_iter = vj_prev->getconnected_loop();
////
////		bool sameloop = (loopi_iter==loopj_iter);
////
////		iflavor_n = Gi_iter->getjflavor();
////		jflavor_n = Gj_iter->getjflavor();
////
////		//compactness check
////		if(Norder==1 && Nselfloop!=0) return true;
////		double proposal = static_cast<double>(NdressChoice*NGs*(NGs-1))/2./NWs_selection;
////		int Wselection;
////
////		if(!sameloop){
////			Wselection = 1 + static_cast<int>( 2*unidist(mt) );
////			proposal *= 2.;
////			if(Wselection==1){
////				coordinate x_c = mean(
////						Gi_iter->getxi(),Gi_iter->getxj(),
////						Gj_iter->getxi(),Gj_iter->getxj()
////						);
////				xi_n.genrand(x_c);
////				xj_n = xi_n;
////				proposal *= beta
////					/xi_n.getsitegenprobability(x_c);
////			}
////			else{
////				coordinate xc_i = mean(Gi_iter->getxi(),Gi_iter->getxj());
////				coordinate xc_j = mean(Gj_iter->getxi(),Gj_iter->getxj());
////				xi_n.genrand(xc_i);
////				xj_n.genrand(xc_j);
////				proposal *= beta*beta
////					/xi_n.getsitegenprobability(xc_i)
////					/xj_n.getsitegenprobability(xc_j);
////			}
////		}
////		else{
////			Wselection = 2;
////
////			coordinate xc_i = mean(Gi_iter->getxi(),Gi_iter->getxj());
////			coordinate xc_j = mean(Gj_iter->getxi(),Gj_iter->getxj());
////			xi_n.genrand(xc_i);
////			xj_n.genrand(xc_j);
////			proposal *= beta*beta
////				/xi_n.getsitegenprobability(xc_i)
////				/xj_n.getsitegenprobability(xc_j);
////		}
////
////		// for probability check
////		line UWij(true,Wselection,Ndress,iflavor_n,jflavor_n,xi_n,xj_n);
////		line Gin_m = (*Gi_iter); Gin_m.setj(iflavor_n,xi_n);
////		line Gjn_m = (*Gj_iter); Gjn_m.setj(jflavor_n,xj_n);
////		line Gin_c(true,0,0,iflavor_n,Gi_iter->getjflavor(),xi_n,Gi_iter->getxj());
////		line Gjn_c(true,0,0,jflavor_n,Gj_iter->getjflavor(),xj_n,Gj_iter->getxj());
////
////		// measuring line selection
////		bool measuring_line_involved = false;
////		bool relocate_measuring_line = false;
////		if( Gi_iter->getphysical()^Gj_iter->getphysical() ){
////			measuring_line_involved = true;
////			proposal *= 2.;
////			if(unidist(mt)<0.5){
////				relocate_measuring_line = true;
////				Gin_m.setphysical(true);
////				Gjn_m.setphysical(true);
////				Gin_c.setphysical(Gi_iter->getphysical());
////				Gjn_c.setphysical(Gj_iter->getphysical());
////			}
////		}
////
////		probability_add 
////			= proposal
////			* Ri_N[Norder] / Ri_N[Norder+1];
////		
////			// * Gin_m.getvalue()*Gjn_m.getvalue()
////			// * Gin_c.getvalue()*Gjn_c.getvalue()
////			// / Gi_iter->getvalue() / Gj_iter->getvalue();
////
////		double loopivaluerh_n, loopivaluelh_n, loopjvaluerh_n, loopjvaluelh_n;
////
////		bool ismeasuringlineiniloop;
////		bool ismeasuringlineinjloop;
////
////		setUWinter trial;
////
////		if(sameloop){
////			ismeasuringlineiniloop = (measuringline->getconnected_vi()->getconnected_loop()==loopi_iter);
////
////			loopivaluerh_n = Gin_m.getvalueij()*Gjn_m.getvalueij()
////					* Gin_c.getvalueij()*Gjn_c.getvalueij()
////					/ Gi_iter->getvalueij() / Gj_iter->getvalueij()
////					* loopi_iter->getvaluerh();
////			loopivaluelh_n = Gin_m.getvalueji()*Gjn_m.getvalueji()
////					* Gin_c.getvalueji()*Gjn_c.getvalueji()
////					/ Gi_iter->getvalueji() / Gj_iter->getvalueji()
////					*loopi_iter->getvaluelh();
////
////			if(ismeasuringlineiniloop){
////				probability_add *= loopivaluerh_n / loopi_iter->getvaluerh()
////						* UWij.getvalueij(); 
////			}
////			else{
////				probability_add *= (loopivaluerh_n + loopivaluelh_n) 
////						/ ( loopi_iter->getvaluerh() + loopi_iter->getvaluelh() )
////						* UWij.getvalueij(); 
////			}
////		}
////		else{
////			ismeasuringlineiniloop = (measuringline->getconnected_vi()->getconnected_loop()==loopi_iter);
////			ismeasuringlineinjloop = (measuringline->getconnected_vi()->getconnected_loop()==loopj_iter);
////
////			loopivaluerh_n = Gin_m.getvalueij() * Gin_c.getvalueij()
////					/ Gi_iter->getvalueij() 
////					* loopi_iter->getvaluerh();
////			loopivaluelh_n = Gin_m.getvalueji() * Gin_c.getvalueji()
////					/ Gi_iter->getvalueji() 
////					*loopi_iter->getvaluelh();
////
////			loopjvaluerh_n = Gjn_m.getvalueij() * Gjn_c.getvalueij()
////					/ Gj_iter->getvalueij()
////					* loopj_iter->getvaluerh();
////			loopjvaluelh_n = Gjn_m.getvalueji() *Gjn_c.getvalueji()
////					/ Gj_iter->getvalueji()
////					*loopj_iter->getvaluelh();
////
////
////			//cout<<"different loop case"<<endl;
////			int iloop = distance(fermionicloops.begin(),loopi_iter);
////			int jloop = distance(fermionicloops.begin(),loopj_iter);
////			trial = SetUWInter;
////			trial.push_back(UWij,iloop,jloop);
////			trial.calvalue();
////					
////			if(ismeasuringlineiniloop){
////				probability_add *= loopivaluerh_n / loopi_iter->getvaluerh()
////						* (loopjvaluerh_n+loopjvaluelh_n) 
////						/ (loopj_iter->getvaluerh()+loopj_iter->getvaluelh())
////						* trial.getvalue() / SetUWInter.getvalue();
////			}
////			else if(ismeasuringlineinjloop){
////				probability_add *= (loopivaluerh_n+loopivaluelh_n) 
////						/ (loopi_iter->getvaluerh()+loopi_iter->getvaluelh())
////						* loopjvaluerh_n / loopj_iter->getvaluerh()
////						* trial.getvalue() / SetUWInter.getvalue();
////			}
////			else{
////				probability_add *= (loopivaluerh_n+loopivaluelh_n) 
////						/ (loopi_iter->getvaluerh()+loopi_iter->getvaluelh())
////						* (loopjvaluerh_n+loopjvaluelh_n) 
////						/ (loopj_iter->getvaluerh()+loopj_iter->getvaluelh())
////						* trial.getvalue() / SetUWInter.getvalue();
////			}
////		}
////
////		//printconfiguration(cout);
////		//cout<<"* addvertice info *"<<endl;
////		//cout<<"  Ndress = "<<Ndress<<endl;
////		//cout<<"  probability = "<<probability<<endl;
////		//cout<<"  proposal = "<<proposal<<endl;
////		//if(sameloop){
////		//	cout<<"  same loop case: "<<endl;
////		//	cout<<"    loopivaluerh_n = "<<loopivaluerh_n<<endl;
////		//	cout<<"    loopivaluelh_n = "<<loopivaluelh_n<<endl;
////		//}
////		//else{
////		//	cout<<"  different loop case: "<<endl;
////		//	cout<<"    loopivaluerh_n = "<<loopivaluerh_n<<endl;
////		//	cout<<"    loopivaluelh_n = "<<loopivaluelh_n<<endl;
////		//	cout<<"    loopjvaluerh_n = "<<loopjvaluerh_n<<endl;
////		//	cout<<"    loopjvaluelh_n = "<<loopjvaluelh_n<<endl;
////		//}
////		//exit(EXIT_FAILURE);
////		//cout<<"  UWij = "<<UWij<<endl;
////		//cout<<"  G.G.G.G = "<<Gin_m.getvalue()*Gjn_m.getvalue()*Gin_c.getvalue()*Gjn_c.getvalue()<<endl;
////		//cout<<"  Gi_iter->getvalue() * Gj_iter->getvalue() = "<<Gi_iter->getvalue() * Gj_iter->getvalue()<<endl;
////		//cout<<"  coordinates"<<endl;
////		//cout<<xi_n<<xj_n<<endl;
////		//cout<<trial<<endl;
////		//cout<<" trial.getvalue() = "<<trial.getvalue()<<endl;
////		//cout<<"*******************"<<endl;
////
////		int signfactor;
////		if(probability_add<0.) signfactor = -1;
////		else signfactor = 1;
////
////		if( abs(probability_add)>1.0e-10 ){
////			vertices.emplace_front(vertex(iflavor_n,xi_n));
////			list<vertex>::iterator vi_iter = vertices.begin();
////			vertices.emplace_front(vertex(jflavor_n,xj_n));
////			list<vertex>::iterator vj_iter = vertices.begin();
////
////			Ws.push_front(UWij); list<line>::iterator UW_iter = Ws.begin();
////			UW_iter->setconnected_vertices(vi_iter,vj_iter);
////
////			Gs.push_front(Gin_c); list<line>::iterator Gi_next = Gs.begin();
////			Gs.push_front(Gjn_c); list<line>::iterator Gj_next = Gs.begin();
////
////			Gi_iter->setj(iflavor_n,xi_n); Gi_iter->setphysical(Gin_m.getphysical());
////			Gj_iter->setj(jflavor_n,xj_n); Gj_iter->setphysical(Gjn_m.getphysical());
////
////			list<vertex>::iterator vi_next = Gi_iter->getconnected_vj();
////			list<vertex>::iterator vj_next = Gj_iter->getconnected_vj();
////			Gi_next->setconnected_vertices(vi_iter,vi_next);
////			Gj_next->setconnected_vertices(vj_iter,vj_next);
////			vi_next->setconnected_Gin(Gi_next);
////			vj_next->setconnected_Gin(Gj_next);
////
////			Gi_iter->setconnected_vj(vi_iter);
////			Gj_iter->setconnected_vj(vj_iter);
////
////			vi_iter->setconnected_lines(UW_iter,Gi_iter,Gi_next);
////			vj_iter->setconnected_lines(UW_iter,Gj_iter,Gj_next);
////
////			if(relocate_measuring_line){
////				for(measuringline=Gs.begin(); measuringline!=Gs.end(); ++measuringline)
////					if(!(measuringline->getphysical())) break;
////
////			}
////
////			// loop SetUWInter update
////			loopi_iter->insert(vi_prev,vi_iter);
////			loopj_iter->insert(vj_prev,vj_iter);
////			vi_iter->setconnected_loop(loopi_iter);
////			vj_iter->setconnected_loop(loopj_iter);
////
////			if(sameloop){
////				loopi_iter->setvalue(loopivaluerh_n,loopivaluelh_n);
////			}
////			else{ 
////				loopi_iter->setvalue(loopivaluerh_n,loopivaluelh_n);
////				loopj_iter->setvalue(loopjvaluerh_n,loopjvaluelh_n);
////				SetUWInter = trial; 
////			}
////
////			Norder += 1 + Ndress;
////			sign *= signfactor;
////			Nselfloop += dNselfloop;
////		}
////		else{ return true; }
////	}
////
////
////	///////////////////
////	// remove update //
////	///////////////////
////
////	int NWs= Ws.size();
////	if(NWs<2) return false;
////
////	int iUW, NWs_selection;
////	int m_type = measuringline->gettype();
////	list<line>::iterator UW_iter;
////
////	if(m_type==0){
////		NWs_selection = NWs;
////		//iUW = NWs*unidist(mt);
////		iUW = 0;
////		UW_iter = next(Ws.begin(),iUW);
////	}
////	else{
////		int iUW_measuringline;
////		iUW_measuringline = distance(Ws.begin(),measuringline);
////		NWs_selection = NWs - 1;
////		//select_one_excluding_one(NWs,iUW_measuringline,iUW);
////		iUW = 0;
////		UW_iter = next(Ws.begin(),iUW);
////	}
////
////	int Wselection = UW_iter->gettype();
////
////	int NdressChoice = 1;
////	int Ndress = UW_iter->getdress();
////
////	list<vertex>::iterator vi_iter = UW_iter->getconnected_vi();
////	list<vertex>::iterator vj_iter = UW_iter->getconnected_vj();
////	list<fermionicloop>::iterator loopi_iter = vi_iter->getconnected_loop();
////	list<fermionicloop>::iterator loopj_iter = vj_iter->getconnected_loop();
////
////	bool sameloop = (loopi_iter==loopj_iter);
////
////	list<line>::iterator Giin_iter = vi_iter->getconnected_Gin();
////	list<line>::iterator Giout_iter = vi_iter->getconnected_Gout();
////
////	list<line>::iterator Gjin_iter = vj_iter->getconnected_Gin();
////	list<line>::iterator Gjout_iter = vj_iter->getconnected_Gout();
////
////	list<vertex>::iterator vi_next = Giout_iter->getconnected_vj();
////	list<vertex>::iterator vj_next = Gjout_iter->getconnected_vj();
////
////	bool PhyGiin = Giin_iter->getphysical(), PhyGiout = Giout_iter->getphysical();
////	bool PhyGjin = Gjin_iter->getphysical(), PhyGjout = Gjout_iter->getphysical();
////	bool measuring_line_involved = false;
////	bool relocate_measuring_line = false;
////	if( (Giin_iter==Giout_iter) || (Gjin_iter==Gjout_iter) ){ return false; }
////	else if((!(PhyGiin&&PhyGiout))||(!(PhyGjin&&PhyGjout))){
////		measuring_line_involved = true;
////		if(PhyGiout^PhyGjout){ relocate_measuring_line = true; }
////	}
////
////	// count the change of the number of self-loop and fermionic-loop
////	int dNselfloop = 0;
////	if(vi_next->getconnected_Gout()==Giin_iter){ dNselfloop++; }
////	if(vj_next->getconnected_Gout()==Gjin_iter){ dNselfloop++; }
////	if(Nselfloop+dNselfloop>NselfloopMax){ if(Ws.size()!=2) return false; }
////
////	line Gi_n((Giin_iter->getphysical()&&Giout_iter->getphysical()),0,Giin_iter->getdress()+Giout_iter->getdress(),
////			Giin_iter->getiflavor(),Giout_iter->getjflavor(),
////			Giin_iter->getxi(),Giout_iter->getxj());
////	line Gj_n((Gjin_iter->getphysical()&&Gjout_iter->getphysical()),0,Gjin_iter->getdress()+Gjout_iter->getdress(),
////			Gjin_iter->getiflavor(),Gjout_iter->getjflavor(),
////			Gjin_iter->getxi(),Gjout_iter->getxj());
////
////	// check whether the update has relevance with the measuring line
////	// Since only Gout segment will be removed, 
////	double proposal = 2.*static_cast<double>(NWs_selection)/NdressChoice/(Gs.size()-2)/(Gs.size()-3);
////	if(!sameloop){
////		proposal /= 2.;
////		if(Wselection==1){
////			coordinate x_c = mean(
////					Giin_iter->getxi(),Giout_iter->getxj(),
////					Gjin_iter->getxi(),Gjout_iter->getxj()
////					);
////			proposal *= vi_iter->getx().getsitegenprobability(x_c)
////				/ beta;
////		}
////		else{
////			coordinate xc_i = mean(Giin_iter->getxi(),Giout_iter->getxj());
////			coordinate xc_j = mean(Gjin_iter->getxi(),Gjout_iter->getxj());
////			proposal *= vi_iter->getx().getsitegenprobability(xc_i)
////				* vj_iter->getx().getsitegenprobability(xc_j)
////				/ beta / beta;
////		}
////	}
////	else{
////		coordinate xc_i = mean(Giin_iter->getxi(),Giout_iter->getxj());
////		coordinate xc_j = mean(Gjin_iter->getxi(),Gjout_iter->getxj());
////		proposal *= vi_iter->getx().getsitegenprobability(xc_i)
////			* vj_iter->getx().getsitegenprobability(xc_j)
////			/ beta / beta;
////	}
////	if(measuring_line_involved) proposal /= 2.;
////
////	probability_remove 
////		= Ri_N[Norder] / Ri_N[Norder-1];
////
////		// * Gi_n.getvalue()*Gj_n.getvalue()
////		// / Giin_iter->getvalue() / Giout_iter->getvalue()
////		// / Gjin_iter->getvalue() / Gjout_iter->getvalue();
////
////	//cout<<"trial setUWinter update"<<endl<<flush;
////	//cout<<"removing line info."<<endl;
////	//cout<<*UW_iter<<endl;
////	//cout<<"previous setUWinter info."<<endl;
////	//cout<<SetUWInter<<endl;
////	//cout<<"trial setUWinter info."<<endl;
////	//cout<<trial<<endl;
////
////	double loopivaluerh_n, loopivaluelh_n, loopjvaluerh_n, loopjvaluelh_n;
////
////	bool ismeasuringlineiniloop;
////	bool ismeasuringlineinjloop;
////
////	setUWinter trial = SetUWInter;
////
////	if(sameloop){
////		ismeasuringlineiniloop = (measuringline->getconnected_vi()->getconnected_loop()==loopi_iter);
////
////		loopivaluerh_n = Gi_n.getvalueij()*Gj_n.getvalueij()
////				/ Giin_iter->getvalueij() / Giout_iter->getvalueij()
////				/ Gjin_iter->getvalueij() / Gjout_iter->getvalueij()
////				* loopi_iter->getvaluerh();
////
////		loopivaluelh_n = Gi_n.getvalueji()*Gj_n.getvalueji()
////				/ Giin_iter->getvalueji() / Giout_iter->getvalueji()
////				/ Gjin_iter->getvalueji() / Gjout_iter->getvalueji()
////				* loopi_iter->getvaluelh();
////
////		if(ismeasuringlineiniloop){
////			probability_remove *= loopivaluerh_n / loopi_iter->getvaluerh()
////					/ UW_iter->getvalueij(); 
////		}
////		else{
////			probability_remove *= (loopivaluerh_n+loopivaluelh_n)
////					/ (loopi_iter->getvaluerh()+loopi_iter->getvaluelh())
////					/ UW_iter->getvalueij(); 
////		}
////	}
////	else{
////		ismeasuringlineiniloop = (measuringline->getconnected_vi()->getconnected_loop()==loopi_iter);
////		ismeasuringlineinjloop = (measuringline->getconnected_vi()->getconnected_loop()==loopj_iter);
////
////		loopivaluerh_n = Gi_n.getvalueij()
////				/ Giin_iter->getvalueij() / Giout_iter->getvalueij()
////				* loopi_iter->getvaluerh();
////		loopivaluelh_n = Gi_n.getvalueji()
////				/ Giin_iter->getvalueji() / Giout_iter->getvalueji()
////				* loopi_iter->getvaluelh();
////
////		loopjvaluerh_n = Gj_n.getvalueij()
////				/ Gjin_iter->getvalueij() / Gjout_iter->getvalueij()
////				* loopj_iter->getvaluerh();
////		loopjvaluelh_n = Gj_n.getvalueji()
////				/ Gjin_iter->getvalueji() / Gjout_iter->getvalueji()
////				* loopj_iter->getvaluelh();
////		
////		trial.erase(*UW_iter);
////		trial.calvalue();
////
////		if(ismeasuringlineiniloop){
////			probability_remove *= loopivaluerh_n / loopi_iter->getvaluerh()
////					* (loopjvaluerh_n+loopjvaluelh_n)
////					/ (loopj_iter->getvaluerh() + loopj_iter->getvaluelh())
////					* trial.getvalue() / SetUWInter.getvalue();
////		}
////		else if(ismeasuringlineinjloop){
////			probability_remove *= (loopivaluerh_n+loopivaluelh_n)
////					/ (loopi_iter->getvaluerh() + loopi_iter->getvaluelh())
////					* loopjvaluerh_n
////					/ loopj_iter->getvaluerh() 
////					* trial.getvalue() / SetUWInter.getvalue();
////		}
////		else{
////			probability_remove *= (loopivaluerh_n+loopivaluelh_n)
////					/ (loopi_iter->getvaluerh() + loopi_iter->getvaluelh())
////					* (loopjvaluerh_n+loopjvaluelh_n)
////					/ (loopj_iter->getvaluerh() + loopj_iter->getvaluelh())
////					* trial.getvalue() / SetUWInter.getvalue();
////		}
////	}
////	probability_remove *= proposal;
////
////	int signfactor;
////	if(probability_remove<0.) signfactor = -1;
////	else signfactor = 1;
////
////	//printconfiguration(cout);
////	//cout<<"* removevertice info *"<<endl;
////	//cout<<" NWs = "<<NWs<<endl;
////	//cout<<" Ndress = "<<Ndress<<endl;
////	//cout<<" probability = "<<probability<<endl;
////	//cout<<" proposal = "<<proposal<<endl;
////	//if(sameloop){
////	//	cout<<"same loop case"<<endl;
////	//	cout<<"loopivaluerh_n = "<<loopivaluerh_n<<endl;
////	//	cout<<"loopivaluelh_n = "<<loopivaluelh_n<<endl;
////	//}
////	//else{
////	//	cout<<"different loop case"<<endl;
////	//	cout<<"loopivaluerh_n = "<<loopivaluerh_n<<endl;
////	//	cout<<"loopivaluelh_n = "<<loopivaluelh_n<<endl;
////	//	cout<<"loopjvaluerh_n = "<<loopjvaluerh_n<<endl;
////	//	cout<<"loopjvaluelh_n = "<<loopjvaluelh_n<<endl;
////	//}
////	//cout<<" UW_iter->getvalue() = "<<UW_iter->getvalue()<<endl;
////	//cout<<" Gi_n.getvalue() = "<<Gi_n.getvalue()<<endl;
////	//cout<<" Gj_n.getvalue() = "<<Gj_n.getvalue()<<endl;
////	//cout<<" Giin_iter->getvalue() = "<<Giin_iter->getvalue()<<endl;
////	//cout<<" Giout_iter->getvalue() = "<<Giout_iter->getvalue()<<endl;
////	//cout<<" Gjin_iter->getvalue() = "<<Gjin_iter->getvalue()<<endl;
////	//cout<<" Gjout_iter->getvalue() = "<<Gjout_iter->getvalue()<<endl;
////	//cout<<"*******************"<<endl;
////
////	{
////		//if(!(checkconnectivityforremove(UW_iter))){ return 0; }
////		//if(!(checkirreducibilityforremove(UW_iter))){ return 0; }
////		if(m_type==0){
////			if(checkirreducibilityforremove(UW_iter)){ 
////				if(!(checkcompactnessforremove(UW_iter))){ return false; }
////			}
////			else return false;
////		}
////		else{
////			if(!(checkcompactnessforremove(UW_iter))){ return false; }
////		}
////
////		vi_next->setconnected_Gin(Giin_iter);
////		Giin_iter->setphysical(Giin_iter->getphysical()&&Giout_iter->getphysical());
////		Giin_iter->setj(Giout_iter->getjflavor(),Giout_iter->getxj());
////		Giin_iter->setconnected_vj(Giout_iter->getconnected_vj());
////
////		vj_next = Gjout_iter->getconnected_vj();
////		vj_next->setconnected_Gin(Gjin_iter);
////		Gjin_iter->setphysical(Gjin_iter->getphysical()&&Gjout_iter->getphysical());
////		Gjin_iter->setj(Gjout_iter->getjflavor(),Gjout_iter->getxj());
////		Gjin_iter->setconnected_vj(Gjout_iter->getconnected_vj());
////
////		loopi_iter->erase(vi_iter);
////		loopj_iter->erase(vj_iter);
////		if(sameloop){
////			loopi_iter->setvalue(loopivaluerh_n,loopivaluelh_n);
////		}
////		else{ 
////			loopi_iter->setvalue(loopivaluerh_n,loopivaluelh_n);
////			loopj_iter->setvalue(loopjvaluerh_n,loopjvaluelh_n);
////			SetUWInter = trial; 
////		}
////
////		vertices.erase(vi_iter); vertices.erase(vj_iter);
////		Ws.erase(UW_iter);
////		Gs.erase(Giout_iter); Gs.erase(Gjout_iter);
////
////		if(relocate_measuring_line){
////			for(measuringline=Gs.begin(); measuringline!=Gs.end(); ++measuringline)
////				if(!(measuringline->getphysical())) break;
////		}
////
////		Norder -= 1 + Ndress;
////		sign *= signfactor;
////		Nselfloop += dNselfloop;
////	}
////
////	if( abs(probability_add*probability_remove - 1.)<1.0e-10 ){
////		cout<<"probability_addvertice = "<<probability_add<<endl;
////		cout<<"probability_removevertice = "<<probability_remove<<endl;
////		cout<<"probability_addvertice*probability_removevertice = "<<probability_add*probability_remove<<endl;
////		return true;
////	}
////	else{
////		cout<<"***** ALERT! DETAILED BALANCE CONDITION FAILS! *****"<<endl;
////		cout<<"probability_addvertice = "<<probability_add<<endl;
////		cout<<"probability_removevertice = "<<probability_remove<<endl;
////		cout<<"probability_addvertice*probability_removevertice = "<<probability_add*probability_remove<<endl;
////		printconfiguration(cout);
////		return false;
////	}
////};
////bool diagram::checkdetailedbalance_remove_add_vertices(){
////
////	////////////////////
////	// remove vertice //
////	////////////////////
////
////	int NWs= Ws.size();
////	if(NWs<2) return true;
////
////	int iUW, NWs_selection;
////	int m_type = measuringline->gettype();
////	list<line>::iterator UW_iter;
////
////	if(m_type==0){
////		NWs_selection = NWs;
////		iUW = NWs*unidist(mt);
////		UW_iter = next(Ws.begin(),iUW);
////	}
////	else{
////		int iUW_measuringline;
////		iUW_measuringline = distance(Ws.begin(),measuringline);
////		NWs_selection = NWs - 1;
////		select_one_excluding_one(NWs,iUW_measuringline,iUW);
////		UW_iter = next(Ws.begin(),iUW);
////	}
////
////	int Wselection = UW_iter->gettype();
////
////	int NdressChoice = 1;
////	int Ndress = UW_iter->getdress();
////
////	list<vertex>::iterator vi_iter = UW_iter->getconnected_vi();
////	list<vertex>::iterator vj_iter = UW_iter->getconnected_vj();
////	list<fermionicloop>::iterator loopi_iter = vi_iter->getconnected_loop();
////	list<fermionicloop>::iterator loopj_iter = vj_iter->getconnected_loop();
////
////	bool sameloop = (loopi_iter==loopj_iter);
////
////	list<line>::iterator Giin_iter = vi_iter->getconnected_Gin();
////	list<line>::iterator Giout_iter = vi_iter->getconnected_Gout();
////
////	list<line>::iterator Gjin_iter = vj_iter->getconnected_Gin();
////	list<line>::iterator Gjout_iter = vj_iter->getconnected_Gout();
////
////	list<vertex>::iterator vi_next = Giout_iter->getconnected_vj();
////	list<vertex>::iterator vj_next = Gjout_iter->getconnected_vj();
////
////	bool PhyGiin = Giin_iter->getphysical(), PhyGiout = Giout_iter->getphysical();
////	bool PhyGjin = Gjin_iter->getphysical(), PhyGjout = Gjout_iter->getphysical();
////	bool measuring_line_involved = false;
////	bool relocate_measuring_line = false;
////	if( (Giin_iter==Giout_iter) || (Gjin_iter==Gjout_iter) ){ return true; }
////	else if((!(PhyGiin&&PhyGiout))||(!(PhyGjin&&PhyGjout))){
////		measuring_line_involved = true;
////		if(PhyGiout^PhyGjout){ relocate_measuring_line = true; }
////	}
////
////	// count the change of the number of self-loop and fermionic-loop
////	int dNselfloop = 0;
////	if(vi_next->getconnected_Gout()==Giin_iter){ dNselfloop++; }
////	if(vj_next->getconnected_Gout()==Gjin_iter){ dNselfloop++; }
////	if(Nselfloop+dNselfloop>NselfloopMax){ if(Ws.size()!=2) return true; }
////
////	line Gi_n((Giin_iter->getphysical()&&Giout_iter->getphysical()),0,Giin_iter->getdress()+Giout_iter->getdress(),
////			Giin_iter->getiflavor(),Giout_iter->getjflavor(),
////			Giin_iter->getxi(),Giout_iter->getxj());
////	line Gj_n((Gjin_iter->getphysical()&&Gjout_iter->getphysical()),0,Gjin_iter->getdress()+Gjout_iter->getdress(),
////			Gjin_iter->getiflavor(),Gjout_iter->getjflavor(),
////			Gjin_iter->getxi(),Gjout_iter->getxj());
////
////	// check whether the update has relevance with the measuring line
////	// Since only Gout segment will be removed, 
////	double proposal = 2.*static_cast<double>(NWs_selection)/NdressChoice/(Gs.size()-2)/(Gs.size()-3);
////	if(!sameloop){
////		proposal /= 2.;
////		if(Wselection==1){
////			coordinate x_c = mean(
////					Giin_iter->getxi(),Giout_iter->getxj(),
////					Gjin_iter->getxi(),Gjout_iter->getxj()
////					);
////			proposal *= vi_iter->getx().getsitegenprobability(x_c)
////				/ beta;
////		}
////		else{
////			coordinate xc_i = mean(Giin_iter->getxi(),Giout_iter->getxj());
////			coordinate xc_j = mean(Gjin_iter->getxi(),Gjout_iter->getxj());
////			proposal *= vi_iter->getx().getsitegenprobability(xc_i)
////				* vj_iter->getx().getsitegenprobability(xc_j)
////				/ beta / beta;
////		}
////	}
////	else{
////		coordinate xc_i = mean(Giin_iter->getxi(),Giout_iter->getxj());
////		coordinate xc_j = mean(Gjin_iter->getxi(),Gjout_iter->getxj());
////		proposal *= vi_iter->getx().getsitegenprobability(xc_i)
////			* vj_iter->getx().getsitegenprobability(xc_j)
////			/ beta / beta;
////	}
////	if(measuring_line_involved) proposal /= 2.;
////
////	double probability_remove
////		= Ri_N[Norder] / Ri_N[Norder-1];
////
////		// * Gi_n.getvalue()*Gj_n.getvalue()
////		// / Giin_iter->getvalue() / Giout_iter->getvalue()
////		// / Gjin_iter->getvalue() / Gjout_iter->getvalue();
////
////	//cout<<"trial setUWinter update"<<endl<<flush;
////	//cout<<"removing line info."<<endl;
////	//cout<<*UW_iter<<endl;
////	//cout<<"previous setUWinter info."<<endl;
////	//cout<<SetUWInter<<endl;
////	//cout<<"trial setUWinter info."<<endl;
////	//cout<<trial<<endl;
////
////	double loopivaluerh_n, loopivaluelh_n, loopjvaluerh_n, loopjvaluelh_n;
////
////	bool ismeasuringlineiniloop;
////	bool ismeasuringlineinjloop;
////
////	setUWinter trial = SetUWInter;
////
////	if(sameloop){
////		ismeasuringlineiniloop = (measuringline->getconnected_vi()->getconnected_loop()==loopi_iter);
////
////		loopivaluerh_n = Gi_n.getvalueij()*Gj_n.getvalueij()
////				/ Giin_iter->getvalueij() / Giout_iter->getvalueij()
////				/ Gjin_iter->getvalueij() / Gjout_iter->getvalueij()
////				* loopi_iter->getvaluerh();
////
////		loopivaluelh_n = Gi_n.getvalueji()*Gj_n.getvalueji()
////				/ Giin_iter->getvalueji() / Giout_iter->getvalueji()
////				/ Gjin_iter->getvalueji() / Gjout_iter->getvalueji()
////				* loopi_iter->getvaluelh();
////
////		if(ismeasuringlineiniloop){
////			probability_remove *= loopivaluerh_n / loopi_iter->getvaluerh()
////					/ UW_iter->getvalueij(); 
////		}
////		else{
////			probability_remove *= (loopivaluerh_n+loopivaluelh_n)
////					/ (loopi_iter->getvaluerh()+loopi_iter->getvaluelh())
////					/ UW_iter->getvalueij(); 
////		}
////	}
////	else{
////		ismeasuringlineiniloop = (measuringline->getconnected_vi()->getconnected_loop()==loopi_iter);
////		ismeasuringlineinjloop = (measuringline->getconnected_vi()->getconnected_loop()==loopj_iter);
////
////		loopivaluerh_n = Gi_n.getvalueij()
////				/ Giin_iter->getvalueij() / Giout_iter->getvalueij()
////				* loopi_iter->getvaluerh();
////		loopivaluelh_n = Gi_n.getvalueji()
////				/ Giin_iter->getvalueji() / Giout_iter->getvalueji()
////				* loopi_iter->getvaluelh();
////
////		loopjvaluerh_n = Gj_n.getvalueij()
////				/ Gjin_iter->getvalueij() / Gjout_iter->getvalueij()
////				* loopj_iter->getvaluerh();
////		loopjvaluelh_n = Gj_n.getvalueji()
////				/ Gjin_iter->getvalueji() / Gjout_iter->getvalueji()
////				* loopj_iter->getvaluelh();
////		
////		trial.erase(*UW_iter);
////		trial.calvalue();
////
////		if(ismeasuringlineiniloop){
////			probability_remove *= loopivaluerh_n / loopi_iter->getvaluerh()
////					* (loopjvaluerh_n+loopjvaluelh_n)
////					/ (loopj_iter->getvaluerh() + loopj_iter->getvaluelh())
////					* trial.getvalue() / SetUWInter.getvalue();
////		}
////		else if(ismeasuringlineinjloop){
////			probability_remove *= (loopivaluerh_n+loopivaluelh_n)
////					/ (loopi_iter->getvaluerh() + loopi_iter->getvaluelh())
////					* loopjvaluerh_n
////					/ loopj_iter->getvaluerh() 
////					* trial.getvalue() / SetUWInter.getvalue();
////		}
////		else{
////			probability_remove *= (loopivaluerh_n+loopivaluelh_n)
////					/ (loopi_iter->getvaluerh() + loopi_iter->getvaluelh())
////					* (loopjvaluerh_n+loopjvaluelh_n)
////					/ (loopj_iter->getvaluerh() + loopj_iter->getvaluelh())
////					* trial.getvalue() / SetUWInter.getvalue();
////		}
////	}
////	probability_remove *= proposal;
////
////	int signfactor;
////	if(probability_remove<0.) signfactor = -1;
////	else signfactor = 1;
////
////	coordinate xi_n, xj_n;
////	xi_n = vi_iter->getx();
////	xj_n = vj_iter->getx();
////
////	//printconfiguration(cout);
////	//cout<<"* removevertice info *"<<endl;
////
////	//cout<<" NWs = "<<NWs<<endl;
////	//cout<<" Ndress = "<<Ndress<<endl;
////	//cout<<" probability = "<<probability<<endl;
////	//cout<<" proposal = "<<proposal<<endl;
////	//if(sameloop){
////	//	cout<<"same loop case"<<endl;
////	//	cout<<"loopivaluerh_o = "<<loopi_iter->getvaluerh()<<endl;
////	//	cout<<"loopivaluelh_o = "<<loopi_iter->getvaluelh()<<endl;
////	//	cout<<"loopivaluerh_n = "<<loopivaluerh_n<<endl;
////	//	cout<<"loopivaluelh_n = "<<loopivaluelh_n<<endl;
////	//}
////	//else{
////	//	cout<<"different loop case"<<endl;
////	//	cout<<"loopivaluerh_n = "<<loopivaluerh_n<<endl;
////	//	cout<<"loopivaluelh_n = "<<loopivaluelh_n<<endl;
////	//	cout<<"loopjvaluerh_n = "<<loopjvaluerh_n<<endl;
////	//	cout<<"loopjvaluelh_n = "<<loopjvaluelh_n<<endl;
////	//}
////	//cout<<" *UW_iter: "<<*UW_iter<<endl;
////	//cout<<" trial.getvalue() / SetUWInter.getvalue() = "<<trial.getvalue() / SetUWInter.getvalue()<<endl;
////	//cout<<" Gi_n.getvalue() = "<<Gi_n.getvalue()<<endl;
////	//cout<<" Gj_n.getvalue() = "<<Gj_n.getvalue()<<endl;
////	//cout<<" *Giin_iter: "<<*Giin_iter<<endl;
////	//cout<<" *Giout_ite: "<<*Giout_iter<<endl;
////	//cout<<" *Gjin_iter: "<<*Gjin_iter<<endl;
////	//cout<<" *Gjout_iter:  "<<*Gjout_iter<<endl;
////	//cout<<"*******************"<<endl;
////
////	if( abs(probability_remove)>1.0e-10 ){
////		//if(!(checkconnectivityforremove(UW_iter))){ return 0; }
////		//if(!(checkirreducibilityforremove(UW_iter))){ return 0; }
////		if(m_type==0){
////			if(checkirreducibilityforremove(UW_iter)){ 
////				if(!(checkcompactnessforremove(UW_iter))){ return true; }
////			}
////			else return true;
////		}
////		else{
////			if(!(checkcompactnessforremove(UW_iter))){ return true; }
////		}
////
////		vi_next->setconnected_Gin(Giin_iter);
////		Giin_iter->setphysical(Giin_iter->getphysical()&&Giout_iter->getphysical());
////		Giin_iter->setj(Giout_iter->getjflavor(),Giout_iter->getxj());
////		Giin_iter->setconnected_vj(Giout_iter->getconnected_vj());
////
////		vj_next = Gjout_iter->getconnected_vj();
////		vj_next->setconnected_Gin(Gjin_iter);
////		Gjin_iter->setphysical(Gjin_iter->getphysical()&&Gjout_iter->getphysical());
////		Gjin_iter->setj(Gjout_iter->getjflavor(),Gjout_iter->getxj());
////		Gjin_iter->setconnected_vj(Gjout_iter->getconnected_vj());
////
////		loopi_iter->erase(vi_iter);
////		loopj_iter->erase(vj_iter);
////		if(sameloop){
////			loopi_iter->setvalue(loopivaluerh_n,loopivaluelh_n);
////		}
////		else{ 
////			loopi_iter->setvalue(loopivaluerh_n,loopivaluelh_n);
////			loopj_iter->setvalue(loopjvaluerh_n,loopjvaluelh_n);
////			SetUWInter = trial; 
////		}
////
////		vertices.erase(vi_iter); vertices.erase(vj_iter);
////		Ws.erase(UW_iter);
////		Gs.erase(Giout_iter); Gs.erase(Gjout_iter);
////
////		if(relocate_measuring_line){
////			for(measuringline=Gs.begin(); measuringline!=Gs.end(); ++measuringline)
////				if(!(measuringline->getphysical())) break;
////		}
////
////		Norder -= 1 + Ndress;
////		sign *= signfactor;
////		Nselfloop += dNselfloop;
////
////		//cout<<"-->update accepted"<<endl;
////	}
////	else{ return true; }
////
////
////	/////////////////
////	// add vertice //
////	/////////////////
////
////	double probability_add;
////	//if(Norder==(this->Nordermax)){ return false; }
////	if(Ws.size()==(this->Nordermax)){ return false; }
////	else{
////		int NdressChoice = 1;
////		int Ndress = 0;
////
////		int NGs = Gs.size();
////		int iG,jG;
////		select_two_index(NGs,iG,jG);
////
////		int m_type = measuringline->gettype();
////		int NWs_selection;
////		if(m_type==0){ NWs_selection = Ws.size()+1; }
////		else{ NWs_selection = Ws.size(); }
////
////		//coordinate xi_n, xj_n;
////
////		int iflavor_n, jflavor_n;
////
////		//list<line>::iterator Gi_iter = next(Gs.begin(),iG);
////		//list<line>::iterator Gj_iter = next(Gs.begin(),jG);
////		list<line>::iterator Gi_iter = Giin_iter;
////		list<line>::iterator Gj_iter = Gjin_iter;
////
////		int dNselfloop = 0;
////		if(Gi_iter->getconnected_vi()==Gi_iter->getconnected_vj()) dNselfloop--;
////		if(Gj_iter->getconnected_vi()==Gj_iter->getconnected_vj()) dNselfloop--;
////
////		list<vertex>::iterator vi_prev = Gi_iter->getconnected_vi();
////		list<vertex>::iterator vj_prev = Gj_iter->getconnected_vi();
////
////		list<fermionicloop>::iterator loopi_iter = vi_prev->getconnected_loop();
////		list<fermionicloop>::iterator loopj_iter = vj_prev->getconnected_loop();
////
////		bool sameloop = (loopi_iter==loopj_iter);
////
////		iflavor_n = Gi_iter->getjflavor();
////		jflavor_n = Gj_iter->getjflavor();
////
////		//compactness check
////		if(Norder==1 && Nselfloop!=0) return false;
////		double proposal = static_cast<double>(NdressChoice*NGs*(NGs-1))/2./NWs_selection;
////		//int Wselection;
////
////		if(!sameloop){
////			//Wselection = 1 + static_cast<int>( 2*unidist(mt) );
////			proposal *= 2.;
////			if(Wselection==1){
////				coordinate x_c = mean(
////						Gi_iter->getxi(),Gi_iter->getxj(),
////						Gj_iter->getxi(),Gj_iter->getxj()
////						);
////				//xi_n.genrand(x_c);
////				//xj_n = xi_n;
////				proposal *= beta
////					/xi_n.getsitegenprobability(x_c);
////			}
////			else{
////				coordinate xc_i = mean(Gi_iter->getxi(),Gi_iter->getxj());
////				coordinate xc_j = mean(Gj_iter->getxi(),Gj_iter->getxj());
////				//xi_n.genrand(xc_i);
////				//xj_n.genrand(xc_j);
////				proposal *= beta*beta
////					/xi_n.getsitegenprobability(xc_i)
////					/xj_n.getsitegenprobability(xc_j);
////			}
////		}
////		else{
////			Wselection = 2;
////
////			coordinate xc_i = mean(Gi_iter->getxi(),Gi_iter->getxj());
////			coordinate xc_j = mean(Gj_iter->getxi(),Gj_iter->getxj());
////			//xi_n.genrand(xc_i);
////			//xj_n.genrand(xc_j);
////			proposal *= beta*beta
////				/xi_n.getsitegenprobability(xc_i)
////				/xj_n.getsitegenprobability(xc_j);
////		}
////
////		// for probability check
////		line UWij(true,Wselection,Ndress,iflavor_n,jflavor_n,xi_n,xj_n);
////		line Gin_m = (*Gi_iter); Gin_m.setj(iflavor_n,xi_n);
////		line Gjn_m = (*Gj_iter); Gjn_m.setj(jflavor_n,xj_n);
////		line Gin_c(true,0,0,iflavor_n,Gi_iter->getjflavor(),xi_n,Gi_iter->getxj());
////		line Gjn_c(true,0,0,jflavor_n,Gj_iter->getjflavor(),xj_n,Gj_iter->getxj());
////
////		// measuring line selection
////		bool measuring_line_involved = false;
////		//bool relocate_measuring_line = false;
////		if( Gi_iter->getphysical()^Gj_iter->getphysical() ){
////			measuring_line_involved = true;
////			proposal *= 2.;
////			//if(unidist(mt)<0.5){
////			if(relocate_measuring_line){
////				//relocate_measuring_line = true;
////				Gin_m.setphysical(true);
////				Gjn_m.setphysical(true);
////				Gin_c.setphysical(Gi_iter->getphysical());
////				Gjn_c.setphysical(Gj_iter->getphysical());
////			}
////		}
////
////		probability_add 
////			= proposal
////			* Ri_N[Norder] / Ri_N[Norder+1];
////		
////			// * Gin_m.getvalue()*Gjn_m.getvalue()
////			// * Gin_c.getvalue()*Gjn_c.getvalue()
////			// / Gi_iter->getvalue() / Gj_iter->getvalue();
////
////		double loopivaluerh_n, loopivaluelh_n, loopjvaluerh_n, loopjvaluelh_n;
////
////		bool ismeasuringlineiniloop;
////		bool ismeasuringlineinjloop;
////
////		setUWinter trial;
////
////		if(sameloop){
////			ismeasuringlineiniloop = (measuringline->getconnected_vi()->getconnected_loop()==loopi_iter);
////
////			loopivaluerh_n = Gin_m.getvalueij()*Gjn_m.getvalueij()
////					* Gin_c.getvalueij()*Gjn_c.getvalueij()
////					/ Gi_iter->getvalueij() / Gj_iter->getvalueij()
////					* loopi_iter->getvaluerh();
////			loopivaluelh_n = Gin_m.getvalueji()*Gjn_m.getvalueji()
////					* Gin_c.getvalueji()*Gjn_c.getvalueji()
////					/ Gi_iter->getvalueji() / Gj_iter->getvalueji()
////					*loopi_iter->getvaluelh();
////
////			if(ismeasuringlineiniloop){
////				probability_add *= loopivaluerh_n / loopi_iter->getvaluerh()
////						* UWij.getvalueij(); 
////			}
////			else{
////				probability_add *= (loopivaluerh_n + loopivaluelh_n) 
////						/ ( loopi_iter->getvaluerh() + loopi_iter->getvaluelh() )
////						* UWij.getvalueij(); 
////			}
////		}
////		else{
////			ismeasuringlineiniloop = (measuringline->getconnected_vi()->getconnected_loop()==loopi_iter);
////			ismeasuringlineinjloop = (measuringline->getconnected_vi()->getconnected_loop()==loopj_iter);
////
////			loopivaluerh_n = Gin_m.getvalueij() * Gin_c.getvalueij()
////					/ Gi_iter->getvalueij() 
////					* loopi_iter->getvaluerh();
////			loopivaluelh_n = Gin_m.getvalueji() * Gin_c.getvalueji()
////					/ Gi_iter->getvalueji() 
////					*loopi_iter->getvaluelh();
////
////			loopjvaluerh_n = Gjn_m.getvalueij() * Gjn_c.getvalueij()
////					/ Gj_iter->getvalueij()
////					* loopj_iter->getvaluerh();
////			loopjvaluelh_n = Gjn_m.getvalueji() *Gjn_c.getvalueji()
////					/ Gj_iter->getvalueji()
////					*loopj_iter->getvaluelh();
////
////
////			//cout<<"different loop case"<<endl;
////			int iloop = distance(fermionicloops.begin(),loopi_iter);
////			int jloop = distance(fermionicloops.begin(),loopj_iter);
////			trial = SetUWInter;
////			trial.push_back(UWij,iloop,jloop);
////			trial.calvalue();
////					
////			if(ismeasuringlineiniloop){
////				probability_add *= loopivaluerh_n / loopi_iter->getvaluerh()
////						* (loopjvaluerh_n+loopjvaluelh_n) 
////						/ (loopj_iter->getvaluerh()+loopj_iter->getvaluelh())
////						* trial.getvalue() / SetUWInter.getvalue();
////			}
////			else if(ismeasuringlineinjloop){
////				probability_add *= (loopivaluerh_n+loopivaluelh_n) 
////						/ (loopi_iter->getvaluerh()+loopi_iter->getvaluelh())
////						* loopjvaluerh_n / loopj_iter->getvaluerh()
////						* trial.getvalue() / SetUWInter.getvalue();
////			}
////			else{
////				probability_add *= (loopivaluerh_n+loopivaluelh_n) 
////						/ (loopi_iter->getvaluerh()+loopi_iter->getvaluelh())
////						* (loopjvaluerh_n+loopjvaluelh_n) 
////						/ (loopj_iter->getvaluerh()+loopj_iter->getvaluelh())
////						* trial.getvalue() / SetUWInter.getvalue();
////			}
////		}
////
////		//printconfiguration(cout);
////		//cout<<"* addvertice info *"<<endl;
////		//cout<<"  Ndress = "<<Ndress<<endl;
////		//cout<<"  probability = "<<probability<<endl;
////		//cout<<"  proposal = "<<proposal<<endl;
////		//if(sameloop){
////		//	cout<<"  same loop case: "<<endl;
////		//	cout<<"    loopivaluerh_o = "<<loopi_iter->getvaluerh()<<endl;
////		//	cout<<"    loopivaluelh_o = "<<loopi_iter->getvaluelh()<<endl;
////		//	cout<<"    loopivaluerh_n = "<<loopivaluerh_n<<endl;
////		//	cout<<"    loopivaluelh_n = "<<loopivaluelh_n<<endl;
////		//}
////		//else{
////		//	cout<<"  different loop case: "<<endl;
////		//	cout<<"    loopivaluerh_n = "<<loopivaluerh_n<<endl;
////		//	cout<<"    loopivaluelh_n = "<<loopivaluelh_n<<endl;
////		//	cout<<"    loopjvaluerh_n = "<<loopjvaluerh_n<<endl;
////		//	cout<<"    loopjvaluelh_n = "<<loopjvaluelh_n<<endl;
////		//}
////		//exit(EXIT_FAILURE);
////		//cout<<"  UWij = "<<UWij<<endl;
////		//cout<<"  trial.getvalue() / SetUWInter.getvalue() = "<<trial.getvalue() / SetUWInter.getvalue()<<endl;
////
////		//cout<<"Gin_m: "<<Gin_m<<endl;
////		//cout<<"Gjn_m: "<<Gjn_m<<endl;
////		//cout<<"Gin_c: "<<Gin_c<<endl;
////		//cout<<"Gjn_c: "<<Gjn_c<<endl;
////
////		//cout<<"  G.G.G.G = "<<Gin_m.getvalue()*Gjn_m.getvalue()*Gin_c.getvalue()*Gjn_c.getvalue()<<endl;
////		//cout<<"  Gi_iter->getvalue() * Gj_iter->getvalue() = "<<Gi_iter->getvalue() * Gj_iter->getvalue()<<endl;
////		//cout<<"  coordinates"<<endl;
////		//cout<<xi_n<<xj_n<<endl;
////		//cout<<trial<<endl;
////		//cout<<" trial.getvalue() = "<<trial.getvalue()<<endl;
////		//cout<<"*******************"<<endl;
////
////		int signfactor;
////		if(probability_add<0.) signfactor = -1;
////		else signfactor = 1;
////
////		{
////			vertices.emplace_front(vertex(iflavor_n,xi_n));
////			list<vertex>::iterator vi_iter = vertices.begin();
////			vertices.emplace_front(vertex(jflavor_n,xj_n));
////			list<vertex>::iterator vj_iter = vertices.begin();
////
////			Ws.push_front(UWij); list<line>::iterator UW_iter = Ws.begin();
////			UW_iter->setconnected_vertices(vi_iter,vj_iter);
////
////			Gs.push_front(Gin_c); list<line>::iterator Gi_next = Gs.begin();
////			Gs.push_front(Gjn_c); list<line>::iterator Gj_next = Gs.begin();
////
////			Gi_iter->setj(iflavor_n,xi_n); Gi_iter->setphysical(Gin_m.getphysical());
////			Gj_iter->setj(jflavor_n,xj_n); Gj_iter->setphysical(Gjn_m.getphysical());
////
////			list<vertex>::iterator vi_next = Gi_iter->getconnected_vj();
////			list<vertex>::iterator vj_next = Gj_iter->getconnected_vj();
////			Gi_next->setconnected_vertices(vi_iter,vi_next);
////			Gj_next->setconnected_vertices(vj_iter,vj_next);
////			vi_next->setconnected_Gin(Gi_next);
////			vj_next->setconnected_Gin(Gj_next);
////
////			Gi_iter->setconnected_vj(vi_iter);
////			Gj_iter->setconnected_vj(vj_iter);
////
////			vi_iter->setconnected_lines(UW_iter,Gi_iter,Gi_next);
////			vj_iter->setconnected_lines(UW_iter,Gj_iter,Gj_next);
////
////			if(relocate_measuring_line){
////				for(measuringline=Gs.begin(); measuringline!=Gs.end(); ++measuringline)
////					if(!(measuringline->getphysical())) break;
////
////			}
////
////			// loop SetUWInter update
////			loopi_iter->insert(vi_prev,vi_iter);
////			loopj_iter->insert(vj_prev,vj_iter);
////			vi_iter->setconnected_loop(loopi_iter);
////			vj_iter->setconnected_loop(loopj_iter);
////
////			if(sameloop){
////				loopi_iter->setvalue(loopivaluerh_n,loopivaluelh_n);
////			}
////			else{ 
////				loopi_iter->setvalue(loopivaluerh_n,loopivaluelh_n);
////				loopj_iter->setvalue(loopjvaluerh_n,loopjvaluelh_n);
////				SetUWInter = trial; 
////			}
////
////			Norder += 1 + Ndress;
////			sign *= signfactor;
////			Nselfloop += dNselfloop;
////		}
////	}
////
////	int NWs= Ws.size();
////	if(NWs<2) return true;
////
////	int iUW, NWs_selection;
////	int m_type = measuringline->gettype();
////	list<line>::iterator UW_iter;
////
////	if(m_type==0){
////		NWs_selection = NWs;
////		iUW = NWs*unidist(mt);
////		UW_iter = next(Ws.begin(),iUW);
////	}
////	else{
////		int iUW_measuringline;
////		iUW_measuringline = distance(Ws.begin(),measuringline);
////		NWs_selection = NWs - 1;
////		select_one_excluding_one(NWs,iUW_measuringline,iUW);
////		UW_iter = next(Ws.begin(),iUW);
////	}
////
////	int Wselection = UW_iter->gettype();
////
////	int NdressChoice = 1;
////	int Ndress = UW_iter->getdress();
////
////	list<vertex>::iterator vi_iter = UW_iter->getconnected_vi();
////	list<vertex>::iterator vj_iter = UW_iter->getconnected_vj();
////	list<fermionicloop>::iterator loopi_iter = vi_iter->getconnected_loop();
////	list<fermionicloop>::iterator loopj_iter = vj_iter->getconnected_loop();
////
////	bool sameloop = (loopi_iter==loopj_iter);
////
////	list<line>::iterator Giin_iter = vi_iter->getconnected_Gin();
////	list<line>::iterator Giout_iter = vi_iter->getconnected_Gout();
////
////	list<line>::iterator Gjin_iter = vj_iter->getconnected_Gin();
////	list<line>::iterator Gjout_iter = vj_iter->getconnected_Gout();
////
////	list<vertex>::iterator vi_next = Giout_iter->getconnected_vj();
////	list<vertex>::iterator vj_next = Gjout_iter->getconnected_vj();
////
////	bool PhyGiin = Giin_iter->getphysical(), PhyGiout = Giout_iter->getphysical();
////	bool PhyGjin = Gjin_iter->getphysical(), PhyGjout = Gjout_iter->getphysical();
////	bool measuring_line_involved = false;
////	bool relocate_measuring_line = false;
////	if( (Giin_iter==Giout_iter) || (Gjin_iter==Gjout_iter) ){ return true; }
////	else if((!(PhyGiin&&PhyGiout))||(!(PhyGjin&&PhyGjout))){
////		measuring_line_involved = true;
////		if(PhyGiout^PhyGjout){ relocate_measuring_line = true; }
////	}
////
////	// count the change of the number of self-loop and fermionic-loop
////	int dNselfloop = 0;
////	if(vi_next->getconnected_Gout()==Giin_iter){ dNselfloop++; }
////	if(vj_next->getconnected_Gout()==Gjin_iter){ dNselfloop++; }
////	if(Nselfloop+dNselfloop>NselfloopMax){ if(Ws.size()!=2) return true; }
////
////	line Gi_n((Giin_iter->getphysical()&&Giout_iter->getphysical()),0,Giin_iter->getdress()+Giout_iter->getdress(),
////			Giin_iter->getiflavor(),Giout_iter->getjflavor(),
////			Giin_iter->getxi(),Giout_iter->getxj());
////	line Gj_n((Gjin_iter->getphysical()&&Gjout_iter->getphysical()),0,Gjin_iter->getdress()+Gjout_iter->getdress(),
////			Gjin_iter->getiflavor(),Gjout_iter->getjflavor(),
////			Gjin_iter->getxi(),Gjout_iter->getxj());
////
////	// check whether the update has relevance with the measuring line
////	// Since only Gout segment will be removed, 
////	double proposal = 2.*static_cast<double>(NWs_selection)/NdressChoice/(Gs.size()-2)/(Gs.size()-3);
////	if(!sameloop){
////		proposal /= 2.;
////		if(Wselection==1){
////			coordinate x_c = mean(
////					Giin_iter->getxi(),Giout_iter->getxj(),
////					Gjin_iter->getxi(),Gjout_iter->getxj()
////					);
////			proposal *= vi_iter->getx().getsitegenprobability(x_c)
////				/ beta;
////		}
////		else{
////			coordinate xc_i = mean(Giin_iter->getxi(),Giout_iter->getxj());
////			coordinate xc_j = mean(Gjin_iter->getxi(),Gjout_iter->getxj());
////			proposal *= vi_iter->getx().getsitegenprobability(xc_i)
////				* vj_iter->getx().getsitegenprobability(xc_j)
////				/ beta / beta;
////		}
////	}
////	else{
////		coordinate xc_i = mean(Giin_iter->getxi(),Giout_iter->getxj());
////		coordinate xc_j = mean(Gjin_iter->getxi(),Gjout_iter->getxj());
////		proposal *= vi_iter->getx().getsitegenprobability(xc_i)
////			* vj_iter->getx().getsitegenprobability(xc_j)
////			/ beta / beta;
////	}
////	if(measuring_line_involved) proposal /= 2.;
////
////	double probability_remove 
////		= Ri_N[Norder] / Ri_N[Norder-1]
////		* Gi_n.getvalue()*Gj_n.getvalue()
////		/ Giin_iter->getvalue() / Giout_iter->getvalue()
////		/ Gjin_iter->getvalue() / Gjout_iter->getvalue();
////
////	cout<<"loop part factor for probability_remove: "
////		<<Gi_n.getvalue()*Gj_n.getvalue()
////		/ Giin_iter->getvalue() / Giout_iter->getvalue()
////		/ Gjin_iter->getvalue() / Gjout_iter->getvalue()
////		<<endl;
////
////
////	setUWinter trial = SetUWInter;
////	if(sameloop){ probability_remove /= UW_iter->getvalue(); }
////	else{
////		trial.erase(*UW_iter);
////		trial.calvalue();
////		probability_remove *= trial.getvalue()/SetUWInter.getvalue();
////	}
////	probability_remove *= proposal;
////
////	int signfactor;
////	if(probability_remove<0.) signfactor = -1;
////	else signfactor = 1;
////
////	coordinate xi_n, xj_n;
////	xi_n = vi_iter->getx();
////	xj_n = vj_iter->getx();
////
////	if( abs(probability_remove)>1.0e-10 ){
////		cout<<"original SetUWInter"<<endl<<SetUWInter<<endl;
////		if(m_type==0){
////			if(checkirreducibilityforremove(UW_iter)){ 
////				if(!(checkcompactnessforremove(UW_iter))){ return true; }
////			}
////			else return true;
////		}
////		else{
////			if(!(checkcompactnessforremove(UW_iter))){ return true; }
////		}
////
////		vi_next->setconnected_Gin(Giin_iter);
////		Giin_iter->setphysical(Giin_iter->getphysical()&&Giout_iter->getphysical());
////		Giin_iter->setj(Giout_iter->getjflavor(),Giout_iter->getxj());
////		Giin_iter->setconnected_vj(Giout_iter->getconnected_vj());
////
////		vj_next = Gjout_iter->getconnected_vj();
////		vj_next->setconnected_Gin(Gjin_iter);
////		Gjin_iter->setphysical(Gjin_iter->getphysical()&&Gjout_iter->getphysical());
////		Gjin_iter->setj(Gjout_iter->getjflavor(),Gjout_iter->getxj());
////		Gjin_iter->setconnected_vj(Gjout_iter->getconnected_vj());
////
////		loopi_iter->erase(vi_iter);
////		loopj_iter->erase(vj_iter);
////		if(!sameloop){ SetUWInter = trial; }
////
////		cout<<"removed UWline: "<<*UW_iter<<endl;
////		cout<<"intermediate SetUWInter"<<endl<<SetUWInter<<endl;
////
////		vertices.erase(vi_iter); vertices.erase(vj_iter);
////		Ws.erase(UW_iter);
////		Gs.erase(Giout_iter); Gs.erase(Gjout_iter);
////
////		if(relocate_measuring_line){
////			for(measuringline=Gs.begin(); measuringline!=Gs.end(); ++measuringline)
////				if(!(measuringline->getphysical())) break;
////		}
////
////		Norder -= 1 + Ndress;
////		sign *= signfactor;
////		Nselfloop += dNselfloop;
////	}
////	else{ return false; }
////
////	double probability_add;
////	if(Norder==(this->Nordermax)){ return false; }
////	else{
////		int NdressChoice = 1;
////		int Ndress = 0;
////
////		int NGs = Gs.size();
////		int iG,jG;
////		select_two_index(NGs,iG,jG);
////
////		int m_type = measuringline->gettype();
////		int NWs_selection;
////		if(m_type==0){ NWs_selection = Ws.size()+1; }
////		else{ NWs_selection = Ws.size(); }
////
////		//coordinate xi_n, xj_n;
////
////		int iflavor_n, jflavor_n;
////
////		//list<line>::iterator Gi_iter = next(Gs.begin(),iG);
////		//list<line>::iterator Gj_iter = next(Gs.begin(),jG);
////		list<line>::iterator Gi_iter = Giin_iter;
////		list<line>::iterator Gj_iter = Gjin_iter;
////
////		int dNselfloop = 0;
////		if(Gi_iter->getconnected_vi()==Gi_iter->getconnected_vj()) dNselfloop--;
////		if(Gj_iter->getconnected_vi()==Gj_iter->getconnected_vj()) dNselfloop--;
////
////		list<vertex>::iterator vi_prev = Gi_iter->getconnected_vi();
////		list<vertex>::iterator vj_prev = Gj_iter->getconnected_vi();
////
////		list<fermionicloop>::iterator loopi_iter = vi_prev->getconnected_loop();
////		list<fermionicloop>::iterator loopj_iter = vj_prev->getconnected_loop();
////
////		bool sameloop = (loopi_iter==loopj_iter);
////
////		iflavor_n = Gi_iter->getjflavor();
////		jflavor_n = Gj_iter->getjflavor();
////
////		//compactness check
////		if(Norder==1 && Nselfloop!=0) return false;
////		double proposal = static_cast<double>(NdressChoice*NGs*(NGs-1))/2./NWs_selection;
////		//int Wselection;
////
////		if(!sameloop){
////			//Wselection = 1 + static_cast<int>( 2*unidist(mt) );
////			proposal *= 2.;
////			if(Wselection==1){
////				coordinate x_c = mean(
////						Gi_iter->getxi(),Gi_iter->getxj(),
////						Gj_iter->getxi(),Gj_iter->getxj()
////						);
////				//xi_n.genrand(x_c);
////				//xj_n = xi_n;
////				proposal *= beta
////					/xi_n.getsitegenprobability(x_c);
////			}
////			else{
////				coordinate xc_i = mean(Gi_iter->getxi(),Gi_iter->getxj());
////				coordinate xc_j = mean(Gj_iter->getxi(),Gj_iter->getxj());
////				//xi_n.genrand(xc_i);
////				//xj_n.genrand(xc_j);
////				proposal *= beta*beta
////					/xi_n.getsitegenprobability(xc_i)
////					/xj_n.getsitegenprobability(xc_j);
////			}
////		}
////		else{
////			Wselection = 2;
////
////			coordinate xc_i = mean(Gi_iter->getxi(),Gi_iter->getxj());
////			coordinate xc_j = mean(Gj_iter->getxi(),Gj_iter->getxj());
////			//xi_n.genrand(xc_i);
////			//xj_n.genrand(xc_j);
////			proposal *= beta*beta
////				/xi_n.getsitegenprobability(xc_i)
////				/xj_n.getsitegenprobability(xc_j);
////		}
////
////		// for probability check
////		line UWij(true,Wselection,Ndress,iflavor_n,jflavor_n,xi_n,xj_n);
////		line Gin_m = (*Gi_iter); Gin_m.setj(iflavor_n,xi_n);
////		line Gjn_m = (*Gj_iter); Gjn_m.setj(jflavor_n,xj_n);
////		line Gin_c(true,0,0,iflavor_n,Gi_iter->getjflavor(),xi_n,Gi_iter->getxj());
////		line Gjn_c(true,0,0,jflavor_n,Gj_iter->getjflavor(),xj_n,Gj_iter->getxj());
////
////		// measuring line selection
////		bool measuring_line_involved = false;
////		//bool relocate_measuring_line = false;
////		if( Gi_iter->getphysical()^Gj_iter->getphysical() ){
////			measuring_line_involved = true;
////			proposal *= 2.;
////			//if(unidist(mt)<0.5){
////			if(relocate_measuring_line){
////				//relocate_measuring_line = true;
////				Gin_m.setphysical(true);
////				Gjn_m.setphysical(true);
////				Gin_c.setphysical(Gi_iter->getphysical());
////				Gjn_c.setphysical(Gj_iter->getphysical());
////			}
////		}
////
////		probability_add
////			= proposal
////			* Ri_N[Norder] / Ri_N[Norder+1]
////			* Gin_m.getvalue()*Gjn_m.getvalue()
////			* Gin_c.getvalue()*Gjn_c.getvalue()
////			/ Gi_iter->getvalue() / Gj_iter->getvalue();
////
////		cout<<"loop part factor for probability_add: "
////		<<Gin_m.getvalue()*Gjn_m.getvalue()
////		* Gin_c.getvalue()*Gjn_c.getvalue()
////		/ Gi_iter->getvalue() / Gj_iter->getvalue()
////		<<endl;
////		cout<<"1/(loop part factor for probability_add): "
////		<<1.0 / Gin_m.getvalue() / Gjn_m.getvalue()
////		/ Gin_c.getvalue() / Gjn_c.getvalue()
////		* Gi_iter->getvalue() * Gj_iter->getvalue()
////		<<endl;
////
////		setUWinter trial;
////		if(sameloop){ probability_add *= UWij.getvalue(); }
////		else{
////			int iloop = distance(fermionicloops.begin(),loopi_iter);
////			int jloop = distance(fermionicloops.begin(),loopj_iter);
////			trial = SetUWInter;
////			trial.push_back(UWij,iloop,jloop);
////			trial.calvalue();
////			probability_add *= trial.getvalue() / SetUWInter.getvalue();
////		}
////
////		int signfactor;
////		if(probability_add<0.) signfactor = -1;
////		else signfactor = 1;
////
////		{
////			vertices.emplace_front(vertex(iflavor_n,xi_n));
////			list<vertex>::iterator vi_iter = vertices.begin();
////			vertices.emplace_front(vertex(jflavor_n,xj_n));
////			list<vertex>::iterator vj_iter = vertices.begin();
////
////			Ws.push_front(UWij); list<line>::iterator UW_iter = Ws.begin();
////			UW_iter->setconnected_vertices(vi_iter,vj_iter);
////
////			cout<<"added UWline: "<<*UW_iter<<endl;
////
////			Gs.push_front(Gin_c); list<line>::iterator Gi_next = Gs.begin();
////			Gs.push_front(Gjn_c); list<line>::iterator Gj_next = Gs.begin();
////
////			Gi_iter->setj(iflavor_n,xi_n); Gi_iter->setphysical(Gin_m.getphysical());
////			Gj_iter->setj(jflavor_n,xj_n); Gj_iter->setphysical(Gjn_m.getphysical());
////
////			list<vertex>::iterator vi_next = Gi_iter->getconnected_vj();
////			list<vertex>::iterator vj_next = Gj_iter->getconnected_vj();
////			Gi_next->setconnected_vertices(vi_iter,vi_next);
////			Gj_next->setconnected_vertices(vj_iter,vj_next);
////			vi_next->setconnected_Gin(Gi_next);
////			vj_next->setconnected_Gin(Gj_next);
////
////			Gi_iter->setconnected_vj(vi_iter);
////			Gj_iter->setconnected_vj(vj_iter);
////
////			vi_iter->setconnected_lines(UW_iter,Gi_iter,Gi_next);
////			vj_iter->setconnected_lines(UW_iter,Gj_iter,Gj_next);
////
////			if(relocate_measuring_line){
////				for(measuringline=Gs.begin(); measuringline!=Gs.end(); ++measuringline)
////					if(!(measuringline->getphysical())) break;
////
////			}
////
////			// loop update
////			loopi_iter->insert(vi_prev,vi_iter);
////			loopj_iter->insert(vj_prev,vj_iter);
////			vi_iter->setconnected_loop(loopi_iter);
////			vj_iter->setconnected_loop(loopj_iter);
////
////			if(!sameloop){ SetUWInter = trial; }
////			cout<<"coming back SetUWInter"<<endl<<SetUWInter<<endl;
////
////			Norder += 1 + Ndress;
////			sign *= signfactor;
////			Nselfloop += dNselfloop;
////
////		}
////	}
////
////	if( abs(probability_add*probability_remove - 1.)<1.0e-10 ){
////		cout<<"probability_removevertice = "<<probability_remove<<endl;
////		cout<<"probability_addvertice = "<<probability_add<<endl;
////		cout<<"probability_remove*probability_add = "<<probability_add*probability_remove<<endl;
////		return true;
////	}
////	else{
////		cout<<"***** ALERT! DETAILED BALANCE CONDITION FAILS! *****"<<endl;
////		cout<<"probability_removevertice = "<<probability_remove<<endl;
////		cout<<"probability_addvertice = "<<probability_add<<endl;
////		cout<<"probability_remove*probability_add = "<<probability_add*probability_remove<<endl;
////		printconfiguration(cout);
////		return false;
////	}
////};*/
////
bool diagram::checkSign(){
	if( sign*Conf.getWeight() )
		return true;
	else{
		cout<<"***** ALERT! Sign is wrong! *****"<<endl;
		return false;
	}
};
bool diagram::checkDetailedBalanceInsertRemove(){
	//cout<<"** check detailed balance **"<<endl;
	//cout<<" * addvertice / removevertice *"<<endl;
	if(Norder==(this->Nordermax)){ return true; }
	else{
		//cout<<"* before insertion *"<<endl;
		//printconfiguration(cout);

		// ** ADDVERTICE **
		int Nvs = Conf.getVerticesSize();
		int NWs = Conf.getWsSize();
		int iv, jv;
		select_two_index(Nvs,iv,jv);

		coordinate xi_n, xj_n, xb;

		xb.settime(parameters::beta);

		double ti_interval, tj_interval;
		ti_interval = Conf.getVertex(iv+1).time - Conf.getVertex(iv).time;
		if(jv==Nvs-1){
			tj_interval = xb.time - Conf.getVertex(jv).time;
		}
		else{
			tj_interval = Conf.getVertex(jv+1).time - Conf.getVertex(jv).time;
		}

		double proposal = static_cast<double>(Nvs*(Nvs-1))/2./NWs
				* ti_interval * tj_interval;

		//cout<<"proposal = "<<proposal<<endl;
		//cout<<"ti_interval = "<<ti_interval<<endl;
		//cout<<"tj_interval = "<<tj_interval<<endl;

		if( iv==Conf.getim2()-1 || jv==Conf.getim2()-1 ){
			proposal *= 2.;
		}

		Conf.proposeInsertVertices(iv, jv);

		double probabilityInsert
			= proposal
			* Ri_N[Norder] / Ri_N[Norder+1]
			* Conf.getProposedWeight()
			/ Conf.getWeight();

		int signfactor;
		if(probabilityInsert<0.) signfactor = -1;
		else signfactor = 1;

		if( abs(probabilityInsert)>1.0e-10 ){
			Conf.updateInsertVertices();

			Norder += 1;
			sign *= signfactor;
		}
		else{
			//cout<<"Weight zero update: you need to check!"<<endl;
			//Conf.printProposed(cout);
			//return false;
			return true;
		}

		//cout<<"* after insertion *"<<endl;
		//printconfiguration(cout);
		//cout<<"  Wij: "<<Wij<<endl;
		//if( measuring_line->getphysical() ){
		//	exit(EXIT_FAILURE);
		//}

		// ** REMOVEVERTICE **
		NWs= Conf.getWsSize();
		if(NWs<3) return false;

		Nvs = Conf.getVerticesSize();
		int iW = NWs-1;
		int iv_new = Conf.getiv(iW);
		int jv_new = Conf.getjv(iW);
		int im2 = Conf.getim2();

		if( (iv_new==0 || iv_new==1) && im2==2) return false;
		if( (jv_new==Nvs-1 || jv_new==Nvs-2) && im2==Nvs-2) return false;

		Conf.proposeRemoveVertices(iW);

		proposal = 2.*static_cast<double>(NWs-1)/(Nvs-2)/(Nvs-3)
			/ Conf.get_tinterval(iv_new) / Conf.get_tinterval(jv_new);

		//cout<<"proposal = "<<proposal<<endl;
		//cout<<"ti_interval = "<<Conf.get_tinterval(iv_new)<<endl;
		//cout<<"tj_interval = "<<Conf.get_tinterval(jv_new)<<endl;

		if( iv_new==im2-1 || iv_new==im2 ) proposal /= 2.;
		else if( jv_new==im2-1 || jv_new==im2 ) proposal /= 2.;

		double probabilityRemove 
			= proposal
			* Ri_N[Norder] / Ri_N[Norder-1]
			* Conf.getProposedWeight()
			/ Conf.getWeight();

		if(probabilityRemove<0.) signfactor = -1;
		else signfactor = 1;
		{
			if( Conf.checkProposedIrreducibility() ){
				if( !(Conf.checkProposedCompactness()) ){ return false; }
			}
			else return false;

			Conf.updateRemoveVertices();

			Norder -= 1;
			sign *= signfactor;
		}
		//cout<<"* after removal *"<<endl;
		//printconfiguration(cout);

		if( abs(probabilityInsert*probabilityRemove - 1.)<1.0e-10 ){
			//cout<<"probability_insert = "<<probabilityInsert<<endl;
			//cout<<"probability_remove = "<<probabilityRemove<<endl;
			//cout<<"probability_insert*probability_remove = "<<probabilityInsert*probabilityRemove<<endl;
			return true;
		}
		else{
			cout<<"***** ALERT! DETAILED BALANCE CONDITION FAILS! *****"<<endl;
			cout<<"probabilityInsert = "<<probabilityInsert<<endl;
			cout<<"probabilityRemove = "<<probabilityRemove<<endl;
			cout<<"probabilityInsert*probabilityRemove = "<<probabilityInsert*probabilityRemove<<endl;
			printconfiguration(cout);
			return false;
		}
	}
};
bool diagram::checkDetailedBalanceRemoveInsert(){
	//cout<<"before remove"<<endl;
	//printconfiguration(cout);

	int NWs= Conf.getWsSize();
	if(NWs<3) return true;
	int Nvs = Conf.getVerticesSize();
	int iW = 1+(NWs-1)*unidist(mt);
	int iv_new = Conf.getiv(iW);
	int jv_new = Conf.getjv(iW);
	int im2 = Conf.getim2();

	//cout<<"(iv_new,jv_new) = ("<<iv_new<<","<<jv_new<<")"<<endl;

	if( (iv_new==0 || iv_new==1) && im2==2) return true;
	if( (jv_new==Nvs-1 || jv_new==Nvs-2) && im2==Nvs-2) return true;

	Conf.proposeRemoveVertices(iW);

	coordinate vi_new = Conf.getProposedVi();
	coordinate vj_new = Conf.getProposedVj();

	double proposal = 2.*static_cast<double>(NWs-1)/(Nvs-2)/(Nvs-3)
		/ Conf.get_tinterval(iv_new) / Conf.get_tinterval(jv_new);

	bool is_on_top;
	if( iv_new==im2-1 || iv_new==im2 ){
		proposal /= 2.;
		if( iv_new==im2-1 ) is_on_top = true;
		else is_on_top = false;
	}
	else if( jv_new==im2-1 || jv_new==im2 ){
		proposal /= 2.;

		if( jv_new==im2-1 ) is_on_top = true;
		else is_on_top = false;
	}

	double probabilityRemove 
		= proposal
		* Ri_N[Norder] / Ri_N[Norder-1]
		* Conf.getProposedWeight()
		/ Conf.getWeight();

	int signfactor;
	if(probabilityRemove<0.) signfactor = -1;
	else signfactor = 1;

	if( abs(probabilityRemove)>1.0e-15 ){
		if( Conf.checkProposedIrreducibility() ){
			if( !(Conf.checkProposedCompactness()) ){ return true; }
		}
		else return true;
		Conf.updateRemoveVertices();
		Norder -= 1;
		sign *= signfactor;
	}
	else{
		//cout<<"Warning: vanishing probabilityRemove"<<endl<<flush;
		//return false;
		return true;
	}

	//cout<<"after remove"<<endl;
	//printconfiguration(cout);

	iv_new -= 1;
	jv_new -= 2;

	//cout<<"(iv_new,jv_new) = ("<<iv_new<<","<<jv_new<<")"<<endl;
	//cout<<"(vi_new,vj_new) = ("<<vi_new<<","<<vj_new<<")"<<endl;

	Nvs = Conf.getVerticesSize();
	NWs = Conf.getWsSize();

	coordinate xi_n, xj_n, xb;

	xb.settime(parameters::beta);

	double ti_interval, tj_interval;
	ti_interval = Conf.getVertex(iv_new+1).time - Conf.getVertex(iv_new).time;
	if(jv_new==Nvs-1){
		tj_interval = xb.time - Conf.getVertex(jv_new).time;
	}
	else{
		tj_interval = Conf.getVertex(jv_new+1).time - Conf.getVertex(jv_new).time;
	}

	proposal = static_cast<double>(Nvs*(Nvs-1))/2./NWs
			* ti_interval * tj_interval;

	if( iv_new==Conf.getim2()-1 || jv_new==Conf.getim2()-1 ){
		proposal *= 2.;
	}

	Conf.proposeInsertVertices(iv_new, vi_new, jv_new, vj_new, is_on_top);

	double probabilityInsert
		= proposal
		* Ri_N[Norder] / Ri_N[Norder+1]
		* Conf.getProposedWeight()
		/ Conf.getWeight();

	if(probabilityInsert<0.) signfactor = -1;
	else signfactor = 1;
	{
		Conf.updateInsertVertices();

		Norder += 1;
		sign *= signfactor;

	}
	//cout<<"after insert"<<endl;
	//printconfiguration(cout);

	if( abs(probabilityInsert*probabilityRemove - 1.)<1.0e-10 ){
		//cout<<"probability_insert = "<<probabilityInsert<<endl;
		//cout<<"probability_remove = "<<probabilityRemove<<endl;
		//cout<<"probability_insert*probability_remove = "<<probabilityInsert*probabilityRemove<<endl;
		return true;
	}
	else{
		cout<<"***** ALERT! DETAILED BALANCE CONDITION FAILS! *****"<<endl;
		cout<<"probabilityInsert = "<<probabilityInsert<<endl;
		cout<<"probabilityRemove = "<<probabilityRemove<<endl;
		cout<<"probabilityInsert*probabilityRemove = "<<probabilityInsert*probabilityRemove<<endl;
		printconfiguration(cout);
		return false;
	}
};

////bool diagram::check_detailed_balance_remove_insert(){
////	//cout<<"** check detailed balance **"<<endl;
////	//cout<<" * removevertice / addvertice *"<<endl;
////	
////	/////////////////////
////	// REMOVE VERTICES //
////	/////////////////////
////	int NWs= Ws.size();
////	if(NWs<2) return true;
////
////	int iW;
////	list<line>::iterator W_iter;
////
////	iW = (NWs-1)*unidist(mt);
////	W_iter = next(Ws.begin(),1+iW);
////
////	list<vertex>::iterator vi_iter = W_iter->getconnected_vi();
////	list<vertex>::iterator vj_iter = W_iter->getconnected_vj();
////
////	coordinate xi_n = vi_iter->getx();
////	coordinate xj_n = vj_iter->getx();
////
////	list<line>::iterator Giin_iter = vi_iter->getconnected_Gin();
////	list<line>::iterator Giout_iter = vi_iter->getconnected_Gout();
////
////	list<line>::iterator Gjin_iter = vj_iter->getconnected_Gin();
////	list<line>::iterator Gjout_iter = vj_iter->getconnected_Gout();
////
////	list<vertex>::iterator vi_next = Giout_iter->getconnected_vj();
////	list<vertex>::iterator vj_next = Gjout_iter->getconnected_vj();
////
////	bool PhyGiin = Giin_iter->getphysical(), PhyGiout = Giout_iter->getphysical();
////	bool PhyGjin = Gjin_iter->getphysical(), PhyGjout = Gjout_iter->getphysical();
////	bool unphysical_line_involved = !(PhyGiin && PhyGiout && PhyGjin && PhyGjout);
////
////	if( NWs!=2 
////		&& 
////		(
////		 (
////		  ( vj_iter==prev(vertices.end(),1) || vj_iter==prev(vertices.end(),2) ) 
////		  && !(prev(Gs.end(),3)->getphysical())
////		  )
////		 || 
////		 ( 
////		  ( vi_iter==vertices.begin() || vi_iter==next(vertices.begin()) ) 
////		  && !(next(Gs.begin())->getphysical())
////		  )
////		 )
////		) {
////		return true;
////	}
////
////	bool measuring_line_choice = false;
////	bool relocate_measuring_line = false;
////	bool relocate_ref_line = false;
////	if( Giin_iter==measuring_line || Giout_iter==measuring_line ){
////		if( Giin_iter!=Gs.begin() ){
////			measuring_line_choice = true;
////		}
////		else{
////			measuring_line_choice = false;
////		}
////		if( Giout_iter==measuring_line ){ relocate_measuring_line = true; }
////	}
////
////	if( Gjout_iter!=ref_line ){
////		if( Gjin_iter==measuring_line || Gjout_iter==measuring_line ){
////			measuring_line_choice = true;
////			if( Gjout_iter==measuring_line ){ relocate_measuring_line = true; }
////		}
////	}
////	else{
////		relocate_ref_line = true;
////	}
////
////	line Gi_n((Giin_iter->getphysical()&&Giout_iter->getphysical()), 0,
////			Giin_iter->getiflavor(),Giout_iter->getjflavor(),
////			Giin_iter->getxi(),Giout_iter->getxj());
////	line Gj_n((Gjin_iter->getphysical()&&Gjout_iter->getphysical()), 0,
////			Gjin_iter->getiflavor(),Gjout_iter->getjflavor(),
////			Gjin_iter->getxi(),Gjout_iter->getxj());
////
////	double proposal = 2.*static_cast<double>(NWs-1)/(Gs.size()-2)/(Gs.size()-3)
////		/ (Gi_n.getxj().gettime() - Gi_n.getxi().gettime())
////		/ fmod(Gj_n.getxj().gettime() - Gj_n.getxi().gettime() + beta, beta);
////
////	if(measuring_line_choice) proposal /= 2.;
////
////	double probability_remove
////		= proposal
////		* Ri_N[Norder] / Ri_N[Norder-1]
////		* Gi_n.getvalueij()*Gj_n.getvalueij()
////		/ Giin_iter->getvalueij() / Giout_iter->getvalueij()
////		/ Gjin_iter->getvalueij() / Gjout_iter->getvalueij()
////		/ W_iter->getvalueij();
////
////	int signfactor;
////	if(probability_remove<0.) signfactor = -1;
////	else signfactor = 1;
////
////	//cout<<"* removeVertices info *"<<endl;
////	//cout<<"  proposal = "<<proposal<<endl;
////	//cout<<"  probability = "<<probability_remove<<endl;
////	//cout<<"***********************"<<endl;
////	//printconfiguration(cout);
////
////	int iG, jG;
////	if( abs(probability_remove)>1.0e-10 ){
////		//cout<<"Wremove: "<<*W_iter<<endl<<flush;
////		//cout<<"relocate_measuring_line: "<<relocate_measuring_line<<endl<<flush;
////		if(check_irreducibility_for_remove(W_iter)){ 
////			if(!(check_compactness_for_remove(W_iter))){ return true; }
////		}
////		else return true;
////
////		vi_next->setconnected_Gin(Giin_iter);
////		Giin_iter->setphysical(Giin_iter->getphysical()&&Giout_iter->getphysical());
////		Giin_iter->setj(Giout_iter->getjflavor(),Giout_iter->getxj());
////		Giin_iter->setconnected_vj(Giout_iter->getconnected_vj());
////
////		vj_next = Gjout_iter->getconnected_vj();
////		vj_next->setconnected_Gin(Gjin_iter);
////		Gjin_iter->setphysical(Gjin_iter->getphysical()&&Gjout_iter->getphysical());
////		Gjin_iter->setj(Gjout_iter->getjflavor(),Gjout_iter->getxj());
////		Gjin_iter->setconnected_vj(Gjout_iter->getconnected_vj());
////
////		vertices.erase(vi_iter); vertices.erase(vj_iter);
////		Ws.erase(W_iter);
////		Gs.erase(Giout_iter); Gs.erase(Gjout_iter);
////
////		if( relocate_measuring_line ){
////			if( !Giin_iter->getphysical() )
////				measuring_line = Giin_iter;
////			else
////				measuring_line = Gjin_iter;
////		}
////
////		if(relocate_ref_line){
////			ref_line = Gjin_iter;
////		}
////
////		Norder -= 1;
////		sign *= signfactor;
////
////		iG = distance(Gs.begin(),Giin_iter);
////		jG = distance(Gs.begin(),Gjin_iter);
////	}
////
////	//printconfiguration(cout);
////	//cout<<"iG: "<<iG<<endl
////	//	<<"jG: "<<jG<<endl;
////
////	////////////////////
////	// INSERT VERTICE //
////	////////////////////
////
////	double probability_insert;
////	if(Ws.size()==(this->Nordermax)){ return false; }
////	else{
////		int NGs = Gs.size();
////		NWs = Ws.size();
////
////
////		int iflavor_n, jflavor_n;
////
////		list<line>::iterator Gi_iter = next(Gs.begin(),iG);
////		list<line>::iterator Gj_iter = next(Gs.begin(),jG);
////
////		list<vertex>::iterator vi_Gi = Gi_iter->getconnected_vi();
////		list<vertex>::iterator vj_Gi = Gi_iter->getconnected_vj();
////
////		list<vertex>::iterator vi_Gj = Gj_iter->getconnected_vi();
////		list<vertex>::iterator vj_Gj = Gj_iter->getconnected_vj();
////
////		iflavor_n = Gi_iter->getjflavor();
////		jflavor_n = Gj_iter->getjflavor();
////
////		//xi_n.genrand(vi_Gi->getx(),vj_Gi->getx());
////		//xj_n.genrand(vi_Gj->getx(),vj_Gj->getx());
////
////		proposal = static_cast<double>(NGs*(NGs-1))/2./NWs
////				* (vj_Gi->getx().gettime() - vi_Gi->getx().gettime()) 
////			 	* fmod(vj_Gj->getx().gettime() - vi_Gj->getx().gettime() + beta,beta);
////
////		// for probability check
////		line Wij(true,1,iflavor_n,jflavor_n,xi_n,xj_n);
////		line Gin_m = (*Gi_iter); Gin_m.setj(iflavor_n,xi_n);
////		line Gjn_m = (*Gj_iter); Gjn_m.setj(jflavor_n,xj_n);
////		line Gin_c(true,0,iflavor_n,Gi_iter->getjflavor(),xi_n,Gi_iter->getxj());
////		line Gjn_c(true,0,jflavor_n,Gj_iter->getjflavor(),xj_n,Gj_iter->getxj());
////
////		// measuring line selection
////		bool measuring_line_involved = false;
////		bool relocate_measuring_line = false;
////
////		if( Gi_iter==measuring_line ){
////			measuring_line_involved = true;
////			if( Gi_iter!=Gs.begin() ){
////				proposal *= 2.;
////				if(unidist(mt)<0.5){
////					relocate_measuring_line = true;
////					Gin_m.setphysical(true);
////					Gin_c.setphysical(false);
////				}
////			}
////			else{
////				relocate_measuring_line = true;
////				Gin_m.setphysical(true);
////				Gin_c.setphysical(false);
////			}
////		}
////
////		bool ref_line_involved = false;
////		bool relocate_ref_line = false;
////		if( Gj_iter!=ref_line ){
////			if( Gj_iter == measuring_line ){
////				measuring_line_involved = true;
////				proposal *= 2.;
////				if(unidist(mt)<0.5){
////					relocate_measuring_line = true;
////					Gjn_m.setphysical(true);
////					Gjn_c.setphysical(false);
////				}
////			}
////		}
////		else{
////			ref_line_involved = true;
////			relocate_ref_line = true;
////			Gjn_m.setphysical(true);
////			Gjn_c.setphysical(false);
////		}
////
////		probability_insert
////			= proposal
////			* Ri_N[Norder] / Ri_N[Norder+1]
////			* Wij.getvalueij()
////			* Gin_m.getvalueij()*Gjn_m.getvalueij()
////			* Gin_c.getvalueij()*Gjn_c.getvalueij()
////			/ Gi_iter->getvalueij() / Gj_iter->getvalueij();
////
////		if(probability_insert<0.) signfactor = -1;
////		else signfactor = 1;
////
////
////		{
////			vertices.emplace(next(vi_Gi),vertex(iflavor_n,xi_n));
////			vertices.emplace(next(vi_Gj),vertex(jflavor_n,xj_n));
////			list<vertex>::iterator vi_iter = next(vi_Gi);
////			list<vertex>::iterator vj_iter = next(vi_Gj);
////
////			Ws.push_back(Wij); list<line>::iterator W_iter = prev(Ws.end());
////			W_iter->setconnected_vertices(vi_iter,vj_iter);
////			//cout<<*W_iter<<endl;
////
////			Gs.insert(next(Gi_iter),Gin_c); list<line>::iterator Gi_next = next(Gi_iter);
////			Gs.insert(next(Gj_iter),Gjn_c); list<line>::iterator Gj_next = next(Gj_iter);
////
////			Gi_iter->setj(iflavor_n,xi_n); Gi_iter->setphysical(Gin_m.getphysical());
////			Gj_iter->setj(jflavor_n,xj_n); Gj_iter->setphysical(Gjn_m.getphysical());
////
////			list<vertex>::iterator vi_next = Gi_iter->getconnected_vj();
////			list<vertex>::iterator vj_next = Gj_iter->getconnected_vj();
////
////			Gi_next->setconnected_vertices(vi_iter,vi_next);
////			Gj_next->setconnected_vertices(vj_iter,vj_next);
////			vi_next->setconnected_Gin(Gi_next);
////			vj_next->setconnected_Gin(Gj_next);
////
////			Gi_iter->setconnected_vj(vi_iter);
////			Gj_iter->setconnected_vj(vj_iter);
////
////			vi_iter->setconnected_lines(W_iter,Gi_iter,Gi_next);
////			vj_iter->setconnected_lines(W_iter,Gj_iter,Gj_next);
////
////			if( relocate_measuring_line ){
////				if( !Gi_next->getphysical() )
////					measuring_line = Gi_next;
////				else
////					measuring_line = Gj_next;
////			}
////			if( relocate_ref_line ){
////				ref_line = prev(Gs.end());
////			}
////
////			Norder += 1;
////		}
////	}
////	if( abs(probability_insert*probability_remove - 1.)<1.0e-10 ){
////		cout<<"probability_remove = "<<probability_remove<<endl;
////		cout<<"probability_insert = "<<probability_insert<<endl;
////		cout<<"probability_insert*probability_remove = "<<probability_insert*probability_remove<<endl;
////		return true;
////	}
////	else{
////		cout<<"***** ALERT! DETAILED BALANCE CONDITION FAILS! *****"<<endl;
////		cout<<"probability_remove = "<<probability_remove<<endl;
////		cout<<"probability_insert = "<<probability_insert<<endl;
////		cout<<"probability_insert*probability_remove = "<<probability_insert*probability_remove<<endl;
////		printconfiguration(cout);
////		return false;
////	}
////};
bool diagram::checkDetailedBalanceShiftIm2(){
	//cout<<"* shiftIm2 update *"<<endl;
	//cout<<"before the first shiftIm2"<<endl;
	//printconfiguration(cout);

	int Nvs = Conf.getVerticesSize();

	if(Nvs<5) return true;

	int im2 = Conf.getim2();

	int dim2_new;
	select_one_excluding_one(Nvs-3,im2-2,dim2_new);
	int im2_new = 2 + dim2_new;

	Conf.proposeShiftIm2(im2_new);

	double probability1 = Conf.getProposedWeight() / Conf.getWeight();

	int signfactor;
	if(probability1<0.) signfactor = -1;
	else signfactor = 1;

	if( abs(probability1)>1.0e-15 ){
		if( Conf.checkProposedIrreducibility() ){
			if( !(Conf.checkProposedCompactness()) ){ return true; }
		}
		else return true;

		Conf.updateShiftIm2();
		sign *= signfactor;
	}
	else{
		//cout<<"Vanishing probability1"<<endl;
		//printconfiguration(cout);
		//return false;
		return true;
	}

	//cout<<"after the first shiftIm2"<<endl;
	//printconfiguration(cout);

	im2_new = im2;

	Conf.proposeShiftIm2(im2_new);

	double probability2 = Conf.getProposedWeight() / Conf.getWeight();

	if(probability2<0.) signfactor = -1;
	else signfactor = 1;

	if( abs(probability2)>1.0e-15 ){
		if( Conf.checkProposedIrreducibility() ){
			if( !(Conf.checkProposedCompactness()) ){ return true; }
		}
		else return true;

		Conf.updateShiftIm2();
		sign *= signfactor;
	}
	else{
		//cout<<"Vanishing probability2"<<endl;
		//printconfiguration(cout);
		//return false;
		return true;
	}

	//cout<<"after the second shiftIm2"<<endl;
	//printconfiguration(cout);

	if( abs(probability1*probability2 - 1.)<1.0e-10 ){
		//cout<<"probability_insert = "<<probabilityInsert<<endl;
		//cout<<"probability_remove = "<<probabilityRemove<<endl;
		//cout<<"probability_insert*probability_remove = "<<probabilityInsert*probabilityRemove<<endl;
		return true;
	}
	else{
		cout<<"***** ALERT! DETAILED BALANCE CONDITION FAILS! *****"<<endl;
		cout<<"probability1 = "<<probability1<<endl;
		cout<<"probability2 = "<<probability2<<endl;
		cout<<"probability1*probability2 = "<<probability1*probability2<<endl;
		printconfiguration(cout);
		return false;
	}
};
bool diagram::checkDetailedBalanceShiftVertex(){
	int Nvs = Conf.getVerticesSize();
	int iv = 1 + (Nvs-1)*unidist(mt);
	double t_old = Conf.getVertex(iv).time;

	Conf.proposeShift(iv);

	double t_new = Conf.getVertex(iv).time;

	double proposedWeight = Conf.getProposedWeight();
	double currentWeight = Conf.getWeight();

	double probability1 = proposedWeight / currentWeight;

	int signfactor;
	if(probability1<0.) signfactor = -1;
	else signfactor = 1;

	if( abs(probability1)>1.0e-15 ){
		Conf.updateShift();
		sign *= signfactor;
	}
	else{
		//cout<<"Vanishing probability1"<<endl;
		//printconfiguration(cout);
		//return false;
		return true;
	}

	Conf.proposeShift(iv,t_old);

	proposedWeight = Conf.getProposedWeight();
	currentWeight = Conf.getWeight();

	double probability2 = proposedWeight / currentWeight;
	if(probability2<0.) signfactor = -1;
	else signfactor = 1;

	if( abs(probability2)>1.0e-15 ){
		Conf.updateShift();
		sign *= signfactor;

	}
	else{
		cout<<"Vanishing probability2"<<endl;
		printconfiguration(cout);
		return false;
	}

	if( abs(probability1*probability2 - 1.)<1.0e-10 ){
		return true;
	}
	else{
		cout<<"***** ALERT! DETAILED BALANCE CONDITION FAILS! *****"<<endl;
		cout<<"probability1 = "<<probability1<<endl;
		cout<<"probability2 = "<<probability2<<endl;
		cout<<"probability1*probability2 = "<<probability1*probability2<<endl;
		printconfiguration(cout);
		return false;
	}
};
bool diagram::checkDetailedBalanceFlipExtFlavor(){
	int NExtV = 4;
	int iv = NExtV*unidist(mt);

	Conf.proposeFlip(iv);

	double probability1 = Conf.getProposedWeight() / Conf.getWeight();

	int signfactor;
	if(probability1<0.) signfactor = -1;
	else signfactor = 1;

	if( abs(probability1)>1.0e-15 ){
		Conf.updateFlip();
		sign *= signfactor;
	}
	else{
		//cout<<"Vanishing probability1"<<endl;
		//printconfiguration(cout);
		//return false;
		return true;
	}

	Conf.proposeFlip(iv);

	double probability2 = Conf.getProposedWeight() / Conf.getWeight();
	if(probability2<0.) signfactor = -1;
	else signfactor = 1;

	if( abs(probability2)>1.0e-15 ){
		Conf.updateFlip();
		sign *= signfactor;
	}
	else{
		//cout<<"Vanishing probability2"<<endl;
		//printconfiguration(cout);
		//return false;
		return true;
	}

	if( abs(probability1*probability2 - 1.)<1.0e-10 ){
		return true;
	}
	else{
		cout<<"***** ALERT! DETAILED BALANCE CONDITION FAILS! *****"<<endl;
		cout<<"probability1 = "<<probability1<<endl;
		cout<<"probability2 = "<<probability2<<endl;
		cout<<"probability1*probability2 = "<<probability1*probability2<<endl;
		printconfiguration(cout);
		return false;
	}
};

////bool diagram::check_detailed_balance_shift(){
////	int NGs = Gs.size();
////	int iG = (NGs-1)*unidist(mt);
////
////	list<line>::iterator G_iter = next(Gs.begin(),iG);
////	list<line>::iterator G_next = next(G_iter);
////	list<line>::iterator W_iter = G_iter->getconnected_vj()->getconnected_W();
////
////	list<vertex>::iterator v_iter = G_iter->getconnected_vj();
////
////	coordinate x_i = G_iter->getconnected_vi()->getx();
////	coordinate x_j = G_next->getconnected_vj()->getx();
////
////	coordinate x_o = G_iter->getconnected_vj()->getx();
////	coordinate x_n;
////	if( G_iter->getconnected_vi()==G_next->getconnected_vj() ){
////		x_n.genrand();
////	}
////	else{
////		x_n.genrand(x_i,x_j);
////	}
////
////	line Gn = *G_iter; 
////	line Gn_next  = *G_next;
////	line Wn  = *W_iter;
////
////	Gn.setxj(x_n);
////	Gn_next.setxi(x_n);
////	if( W_iter->getconnected_vi()==G_iter->getconnected_vj() )
////		Wn.setxi(x_n);
////	else
////		Wn.setxj(x_n);
////
////	double probability_1 = Gn.getvalueij() * Gn_next.getvalueij() * Wn.getvalueij()
////			   / G_iter->getvalueij() / G_next->getvalueij() / W_iter->getvalueij();
////
////	int signfactor;
////	if(probability_1<0.) signfactor = -1;
////	else signfactor = 1;
////
////	//cout<<"** shiftVertex info **"<<endl;
////	//cout<<"   probability = "<<probability<<endl;
////	//cout<<"   Gn.getvalueij() = "<<Gn.getvalueij()<<endl;
////	//cout<<"   Gn_next.getvalueij() = "<<Gn_next.getvalueij()<<endl;
////	//cout<<"   Wn.getvalueij() = "<<Wn.getvalueij()<<endl;
////	//cout<<"   G_iter->getvalueij() = "<<G_iter->getvalueij()<<endl;
////	//cout<<"   G_next->getvalueij() = "<<G_next->getvalueij()<<endl;
////	//cout<<"   W_iter->getvalueij() = "<<W_iter->getvalueij()<<endl;
////	//cout<<"**********************"<<endl;
////
////	if( abs(probability_1)>1.0e-15 ){
////		G_iter->setxj(x_n);
////		G_next->setxi(x_n);
////		if( W_iter->getconnected_vi()==G_iter->getconnected_vj() )
////			W_iter->setxi(x_n);
////		else
////			W_iter->setxj(x_n);
////
////		v_iter->setx(x_n);
////		
////		sign *= signfactor;
////	}
////
////	Gn = *G_iter; 
////	Gn_next  = *G_next;
////	Wn  = *W_iter;
////
////	Gn.setxj(x_o);
////	Gn_next.setxi(x_o);
////	if( W_iter->getconnected_vi()==G_iter->getconnected_vj() )
////		Wn.setxi(x_o);
////	else
////		Wn.setxj(x_o);
////
////	double probability_2 = Gn.getvalueij() * Gn_next.getvalueij() * Wn.getvalueij()
////			   / G_iter->getvalueij() / G_next->getvalueij() / W_iter->getvalueij();
////
////	if(probability_2<0.) signfactor = -1;
////	else signfactor = 1;
////
////	{
////		G_iter->setxj(x_o);
////		G_next->setxi(x_o);
////		if( W_iter->getconnected_vi()==G_iter->getconnected_vj() )
////			W_iter->setxi(x_o);
////		else
////			W_iter->setxj(x_o);
////
////		v_iter->setx(x_o);
////		
////		sign *= signfactor;
////	}
////
////	if( abs(probability_1*probability_2 - 1.)<1.0e-10 ){
////		cout<<"probability_shift_1 = "<<probability_1<<endl;
////		cout<<"probability_shift_2 = "<<probability_2<<endl;
////		cout<<"probability_1*probability_2 = "<<probability_1*probability_2<<endl;
////		return true;
////	}
////	else{
////		cout<<"***** ALERT! DETAILED BALANCE CONDITION FAILS! *****"<<endl;
////		cout<<"probability_shift_1 = "<<probability_1<<endl;
////		cout<<"probability_shift_2 = "<<probability_2<<endl;
////		cout<<"probability_1*probability_2 = "<<probability_1*probability_2<<endl;
////		return false;
////	}
////};
////
/////*bool diagram::checkdetailedbalance_reconnect(){
////	/////////////////////
////	// first reconnect //
////	/////////////////////
////
////	int Nvs = vertices.size();
////
////	int iv, jv;
////	select_two_index(Nvs, iv, jv);
////
////	list<vertex>::iterator vi_iter = next(vertices.begin(),iv);
////	list<vertex>::iterator vj_iter = next(vertices.begin(),jv);
////
////	list<line>::iterator Giin_iter = vi_iter->getconnected_Gin();
////	list<line>::iterator Gjin_iter = vj_iter->getconnected_Gin();
////
////	line Giin_new(*Gjin_iter);
////	line Gjin_new(*Giin_iter);
////	Giin_new.setj(Giin_iter->getjflavor(),Giin_iter->getxj());
////	Gjin_new.setj(Gjin_iter->getjflavor(),Gjin_iter->getxj());
////	Giin_new.setconnected_vertices(Gjin_iter->getconnected_vi(),vi_iter);
////	Gjin_new.setconnected_vertices(Giin_iter->getconnected_vi(),vj_iter);
////
////	int dNselfloop = 0;
////	if(Giin_iter==vi_iter->getconnected_Gout()) dNselfloop--;
////	if(Gjin_iter==vj_iter->getconnected_Gout()) dNselfloop--;
////	if(Giin_new.getconnected_vi()==Giin_new.getconnected_vj()) dNselfloop++;
////	if(Gjin_new.getconnected_vi()==Gjin_new.getconnected_vj()) dNselfloop++;
////	if(Ws.size()!=1 && Nselfloop+dNselfloop>NselfloopMax) return true;
////
////	// probability check
////	double probability_first
////		= -1.; // Fermionic sign factor
////
////	list<fermionicloop>::iterator loopi_iter = vi_iter->getconnected_loop();
////	list<fermionicloop>::iterator loopj_iter = vj_iter->getconnected_loop();
////
////	bool sameloop = (loopi_iter==loopj_iter);
////
////	setUWinter trial;
////	list<line>::iterator UW_iter;
////	list<fermionicloop>::iterator loop1_iter;
////	list<fermionicloop>::iterator loop2_iter;
////
////	list<list<vertex>::iterator> vset_split1, vset_split2, vset_combined;
////
////	double loop1_valuerh_n, loop1_valuelh_n, loop2_valuerh_n, loop2_valuelh_n;
////	double loop_valuerh_n, loop_valuelh_n;
////
////	if( sameloop ){
////		//cout<<"loop split case"<<endl<<flush;
////		vset_split1 = loopi_iter->getpath(vi_iter,vj_iter);
////		vset_split2 = loopi_iter->getpath(vj_iter,vi_iter);
////
////		int iloop = distance(fermionicloops.begin(),loopi_iter);
////
////		//cout<<"loop to be split index : "<<iloop<<endl;
////		//cout<<*loopi_iter;
////
////		fermionicloop loop1_pseudo(vset_split1);
////		fermionicloop loop2_pseudo(vset_split2);
////
////		loop1_pseudo.calvalue();
////		loop2_pseudo.calvalue();
////
////		cout<<"loop1_pseudo"<<endl<<loop1_pseudo<<endl;
////		cout<<"loop2_pseudo"<<endl<<loop2_pseudo<<endl;
////		cout<<"Giin_new.getvalueij() = "<<Giin_new.getvalueij()<<endl<<"Gjin_iter->getvalueij() = "<<Gjin_iter->getvalueij()<<endl;
////		cout<<"Giin_new.getvalueji() = "<<Giin_new.getvalueji()<<endl<<"Gjin_iter->getvalueji() = "<<Gjin_iter->getvalueji()<<endl;
////		cout<<"Gjin_new.getvalueij() = "<<Gjin_new.getvalueij()<<endl<<"Giin_iter->getvalueij() = "<<Giin_iter->getvalueij()<<endl;
////		cout<<"Gjin_new.getvalueji() = "<<Gjin_new.getvalueji()<<endl<<"Giin_iter->getvalueji() = "<<Giin_iter->getvalueji()<<endl;
////
////		loop1_valuerh_n = loop1_pseudo.getvaluerh() * Giin_new.getvalueij() / Gjin_iter->getvalueij();
////		loop1_valuelh_n = loop1_pseudo.getvaluelh() * Giin_new.getvalueji() / Gjin_iter->getvalueji();
////		loop2_valuerh_n = loop2_pseudo.getvaluerh() * Gjin_new.getvalueij() / Giin_iter->getvalueij();
////		loop2_valuelh_n = loop2_pseudo.getvaluelh() * Gjin_new.getvalueji() / Giin_iter->getvalueji();
////
////		trial.setNfermionicloops(fermionicloops.size()+1);
////		bool v1_in, v2_in;
////		list<list<vertex>::iterator>::iterator v_iiter;
////		list<vertex>::iterator v1_iter;
////		list<vertex>::iterator v2_iter;
////		for(UW_iter=Ws.begin(); UW_iter!=Ws.end(); ++UW_iter){
////			v1_iter = UW_iter->getconnected_vi();
////			v2_iter = UW_iter->getconnected_vj();
////			loop1_iter = v1_iter->getconnected_loop();
////			loop2_iter = v2_iter->getconnected_loop();
////
////			if( (loop1_iter==loopi_iter) && (loop2_iter==loopi_iter) ){
////				v1_in = false; v2_in = false;
////				for(v_iiter=vset_split1.begin(); v_iiter!=vset_split1.end(); ++v_iiter){
////					if(v1_iter==*v_iiter) v1_in = true;
////					else if(v2_iter==*v_iiter) v2_in = true;
////				}
////				if(v1_in^v2_in){
////					probability_first /= UW_iter->getvalueij();
////					if(v1_in) trial.push_back(*UW_iter,iloop,iloop+1);
////					else trial.push_back(*UW_iter,iloop+1,iloop);
////				}
////			}
////			if(loop1_iter!=loop2_iter){
////				if(loop1_iter==loopi_iter){
////					v1_in = false;
////					for(v_iiter=vset_split1.begin(); v_iiter!=vset_split1.end(); ++v_iiter){
////						if(v1_iter==*v_iiter){
////							v1_in = true;
////							break;
////						}
////					}
////					if(v1_in){
////						trial.push_back(
////							*UW_iter,
////							iloop,
////							gen_index_split(iloop,distance(fermionicloops.begin(),loop2_iter))
////						);
////					}
////					else{
////						trial.push_back(
////							*UW_iter,
////							iloop+1,
////							gen_index_split(iloop,distance(fermionicloops.begin(),loop2_iter))
////						);
////					}
////				}
////				else if(loop2_iter==loopi_iter){
////					v1_in = false;
////					for(v_iiter=vset_split1.begin(); v_iiter!=vset_split1.end(); ++v_iiter){
////						if(v2_iter==*v_iiter){
////							v1_in = true;
////							break;
////						}
////					}
////					if(v1_in){
////						trial.push_back(
////							*UW_iter,
////							gen_index_split(iloop,distance(fermionicloops.begin(),loop1_iter)),
////							iloop
////						);
////					}
////					else{
////						trial.push_back(
////							*UW_iter,
////							gen_index_split(iloop,distance(fermionicloops.begin(),loop1_iter)),
////							iloop+1
////						);
////					}
////				}
////				else{
////					trial.push_back(
////						*UW_iter,
////						gen_index_split(iloop,distance(fermionicloops.begin(),loop1_iter)),
////						gen_index_split(iloop,distance(fermionicloops.begin(),loop2_iter))
////					);
////				}
////			}
////		}
////
////		v1_iter = measuringline->getconnected_vi();
////		v2_iter = measuringline->getconnected_vj();
////		loop1_iter = v1_iter->getconnected_loop();
////		loop2_iter = v2_iter->getconnected_loop();
////		bool vi_1in = false; 
////		bool vj_1in = false;
////		for(v_iiter=vset_split1.begin(); v_iiter!=vset_split1.end(); ++v_iiter){
////			if(v1_iter==*v_iiter){ vi_1in = true; }
////			if(v2_iter==*v_iiter){ vj_1in = true; }
////		}
////		bool vi_2in = false; 
////		bool vj_2in = false;
////		for(v_iiter=vset_split2.begin(); v_iiter!=vset_split2.end(); ++v_iiter){
////			if(v1_iter==*v_iiter){ vi_2in = true; }
////			if(v2_iter==*v_iiter){ vj_2in = true; }
////		}
////
////		if( measuringline->gettype()==0 ){
////			if(vi_1in){
////				trial.fixloopflavor(iloop, miflavor);
////				probability_first *= loop1_valuerh_n * 0.5*(loop2_valuerh_n+loop2_valuelh_n)
////						/ loopi_iter->getvaluerh();
////			}
////			else if(vi_2in){
////				trial.fixloopflavor(iloop+1, miflavor);
////				probability_first *= 0.5*(loop1_valuerh_n+loop1_valuelh_n) * loop2_valuerh_n
////						/ loopi_iter->getvaluerh();
////			}
////			else{
////				trial.fixloopflavor(
////					gen_index_split(iloop,distance(fermionicloops.begin(),loop1_iter)),
////					miflavor
////				);
////				probability_first *= 0.5
////						* (loop1_valuerh_n+loop1_valuelh_n) * (loop2_valuerh_n+loop2_valuelh_n)
////						/ (loopi_iter->getvaluerh() + loopi_iter->getvaluelh());
////			}
////		}
////		else{
////			if(vi_1in) trial.fixloopflavor(iloop, miflavor);
////			else if(vi_2in) trial.fixloopflavor(iloop+1, miflavor);
////			else{
////				trial.fixloopflavor(
////					gen_index_split(iloop,distance(fermionicloops.begin(),loop1_iter)),
////					miflavor
////				);
////			}
////			if(vj_1in) trial.fixloopflavor(iloop, mjflavor);
////			else if(vj_2in) trial.fixloopflavor(iloop+1, mjflavor);
////			else{
////				trial.fixloopflavor(
////					gen_index_split(iloop,distance(fermionicloops.begin(),loop2_iter)),
////					mjflavor
////				);
////			}
////			probability_first *= 0.5
////					* (loop1_valuerh_n+loop1_valuelh_n) * (loop2_valuerh_n+loop2_valuelh_n)
////					/ (loopi_iter->getvaluerh() + loopi_iter->getvaluelh());
////		}
////
////		trial.calvalue();
////		probability_first *= trial.getvalue() / SetUWInter.getvalue();
////
////		//printconfiguration(cout);
////		//cout<<"trial.calvalue()"<<endl<<flush;
////		//cout<<trial<<endl;
////		//cout<<"symmetrized UWinter value: "<<trial.getvalue()<<endl;
////	}
////	else{
////		//cout<<"loop merge case"<<endl<<flush;
////		//cout<<*loopi_iter<<*loopj_iter<<endl;
////		// * Giin_new.getvalue() * Gjin_new.getvalue() 
////		// / Giin_iter->getvalue() / Gjin_iter->getvalue();
////		loop_valuerh_n = loopi_iter->getvaluerh()*loopj_iter->getvaluerh()
////					* Giin_new.getvalueij() * Gjin_new.getvalueij() 
////					/ Giin_iter->getvalueij() / Gjin_iter->getvalueij();
////		loop_valuelh_n = loopi_iter->getvaluelh()*loopj_iter->getvaluelh()
////					* Giin_new.getvalueji() * Gjin_new.getvalueji() 
////					/ Giin_iter->getvalueji() / Gjin_iter->getvalueji();
////
////		int iloop = distance(fermionicloops.begin(),loopi_iter);
////		int jloop = distance(fermionicloops.begin(),loopj_iter);
////
////		if( 
////			(SetUWInter.getisloopflavorfixed(iloop) && SetUWInter.getisloopflavorfixed(jloop))
////			&& 
////			(SetUWInter.getloopflavor(iloop)!=SetUWInter.getloopflavor(jloop))
////		){ return true; }
////
////		//cout<<"construct trial"<<endl;
////		trial.setNfermionicloops(fermionicloops.size()-1);
////		int i1loop, i2loop;
////		list<line>::iterator UW_iter;
////		for(UW_iter=Ws.begin(); UW_iter!=Ws.end(); ++UW_iter){
////			loop1_iter = UW_iter->getconnected_vi()->getconnected_loop();
////			loop2_iter = UW_iter->getconnected_vj()->getconnected_loop();
////
////			i1loop = gen_index_merge(iloop,jloop,distance(fermionicloops.begin(),loop1_iter));
////			i2loop = gen_index_merge(iloop,jloop,distance(fermionicloops.begin(),loop2_iter));
////
////			//cout<<"("<<i1loop<<";"<<distance(fermionicloops.begin(),loop1_iter)<<")"
////			//	<<"("<<i2loop<<";"<<distance(fermionicloops.begin(),loop2_iter)<<")"<<endl;
////
////			if(loop1_iter!=loop2_iter){
////				if(i1loop==i2loop){ probability_first *= UW_iter->getvalueij(); }
////				else{ trial.push_back(*UW_iter,i1loop,i2loop); }
////			}
////		}
////		//cout<<"mark measuringline"<<endl;
////		loop1_iter = measuringline->getconnected_vi()->getconnected_loop();
////		loop2_iter = measuringline->getconnected_vj()->getconnected_loop();
////		i1loop = gen_index_merge(iloop,jloop,distance(fermionicloops.begin(),loop1_iter));
////		i2loop = gen_index_merge(iloop,jloop,distance(fermionicloops.begin(),loop2_iter));
////		trial.fixloopflavor( i1loop, miflavor );
////		trial.fixloopflavor( i2loop, mjflavor );
////
////		if(measuringline->gettype()==0){
////			if(loop1_iter==loopi_iter){
////				probability_first *= loop_valuerh_n 
////						/ loopi_iter->getvaluerh() 
////						/ 0.5/(loopj_iter->getvaluerh()+loopj_iter->getvaluelh());
////			}
////			else if(loop1_iter==loopj_iter){
////				probability_first *= loop_valuerh_n 
////						/ 0.5/(loopi_iter->getvaluerh()+loopi_iter->getvaluelh())
////						/ loopj_iter->getvaluerh();
////			}
////			else{
////				probability_first *= (loop_valuerh_n+loop_valuelh_n) 
////						/ 0.5
////						/ (loopi_iter->getvaluerh()+loopi_iter->getvaluelh())
////						/ (loopj_iter->getvaluerh()+loopj_iter->getvaluelh());
////			}
////		}
////		else{
////			probability_first *= (loop_valuerh_n+loop_valuelh_n) 
////					/ 0.5
////					/ (loopi_iter->getvaluerh()+loopi_iter->getvaluelh())
////					/ (loopj_iter->getvaluerh()+loopj_iter->getvaluelh());
////		}
////		//cout<<"calculate trial value"<<endl;
////		//cout<<"symmetrized UWinter value: "<<trial.getvalue()<<endl;
////		trial.calvalue();
////
////		probability_first *= trial.getvalue() / SetUWInter.getvalue();
////
////	}
////
////	printconfiguration(cout);
////	cout<<"* reconnect info *"<<endl;
////	if(sameloop){
////		cout<<"  loop split case"<<endl;
////		cout<<"    loop1_valuerh_n = "<<loop1_valuerh_n<<endl;
////		cout<<"    loop1_valuelh_n = "<<loop1_valuelh_n<<endl;
////		cout<<"    loop2_valuerh_n = "<<loop2_valuerh_n<<endl;
////		cout<<"    loop2_valuelh_n = "<<loop2_valuelh_n<<endl;
////	}
////	else{
////		cout<<"  loop merge case"<<endl;
////		cout<<"    loop_valuerh_n = "<<loop_valuerh_n<<endl;
////		cout<<"    loop_valuelh_n = "<<loop_valuelh_n<<endl;
////	}
////
////	//cout<<"  probability = "<<probability<<endl;
////	//cout<<"  trial info. "<<endl<<trial<<endl;
////	//cout<<"  SetUWInter.getvalue() = "<<SetUWInter.getvalue()<<endl;
////	//cout<<"*******************"<<endl;
////
////	int signfactor;
////	if(probability_first<0.) signfactor = -1;
////	else signfactor = 1;
////
////	if( abs(probability_first)>1.0e-10 ){
////		//cout<<"compactness check"<<endl;
////		//if(!(checkconnectivityforreconnect(Giin_iter,Gjin_iter))) return 0;
////		if(!(checkirreducibilityforreconnect(Giin_iter,Gjin_iter))) return true;
////		if(!(checkcompactnessforreconnect(Giin_iter,Gjin_iter))) return true;
////
////		//cout<<"start to update quantities"<<endl;
////		bool measuring_line_involved = false;
////		if( Giin_iter->getphysical()^Gjin_iter->getphysical() )
////			measuring_line_involved = true;
////
////		vi_iter->setconnected_Gin(Gjin_iter);
////		vj_iter->setconnected_Gin(Giin_iter);
////		Giin_iter->setxj(vj_iter->getx());
////		Gjin_iter->setxj(vi_iter->getx());
////		Giin_iter->setconnected_vj(vj_iter);
////		Gjin_iter->setconnected_vj(vi_iter);
////
////		//cout<<"loop quantities update"<<endl;
////		list<list<vertex>::iterator>::iterator v_iiter;
////		if( sameloop ){
////			//cout<<"loop split case"<<endl;
////			loopi_iter->set(vset_split1);
////			fermionicloops.insert(next(loopi_iter),fermionicloop(vset_split2));
////			loopj_iter = next(loopi_iter);
////
////			loopi_iter->setvalue(loop1_valuerh_n,loop1_valuelh_n);
////			loopj_iter->setvalue(loop2_valuerh_n,loop2_valuelh_n);
////
////			for(v_iiter=vset_split1.begin(); v_iiter!=vset_split1.end(); ++v_iiter)
////				(*v_iiter)->setconnected_loop(loopi_iter);
////			for(v_iiter=vset_split2.begin(); v_iiter!=vset_split2.end(); ++v_iiter)
////				(*v_iiter)->setconnected_loop(loopj_iter);
////		}
////		else{
////			//cout<<"loop merge case"<<endl;
////			//cout<<*loopi_iter<<*loopj_iter<<endl;
////			vset_combined = loopi_iter->getpath(vi_iter,vi_iter);
////			vset_combined.splice(vset_combined.end(),loopj_iter->getpath(vj_iter,vj_iter));
////			int iloop = distance(fermionicloops.begin(),loopi_iter);
////			int jloop = distance(fermionicloops.begin(),loopj_iter);
////			//cout<<"vset_combined"<<endl;
////			//for(v_iiter=vset_combined.begin(); v_iiter!=vset_combined.end(); ++v_iiter)
////			//	cout<<(**v_iiter);
////			//cout<<endl;
////			if(iloop<jloop){
////				loopi_iter->set(vset_combined);
////				fermionicloops.erase(loopj_iter);
////
////				loopi_iter->setvalue(loop_valuerh_n,loop_valuelh_n);
////
////				for(v_iiter=vset_combined.begin(); v_iiter!=vset_combined.end(); ++v_iiter)
////					(*v_iiter)->setconnected_loop(loopi_iter);
////			}
////			else{
////				loopj_iter->set(vset_combined);
////				fermionicloops.erase(loopi_iter);
////
////				loopj_iter->setvalue(loop_valuerh_n,loop_valuelh_n);
////
////				for(v_iiter=vset_combined.begin(); v_iiter!=vset_combined.end(); ++v_iiter)
////					(*v_iiter)->setconnected_loop(loopj_iter);
////			}
////		}
////		SetUWInter = trial;
////		//updateSetUWInter();
////		//SetUWInter.setvalue(trial.getvalue());
////
////		sign *= signfactor;
////		Nselfloop += dNselfloop;
////
////		//printconfiguration(cout);
////		//exit(EXIT_FAILURE);
////	}
////	else{ return true; }
////
////
////	//////////////////////
////	// second reconnect //
////	//////////////////////
////
////	vi_iter = next(vertices.begin(),iv);
////	vj_iter = next(vertices.begin(),jv);
////
////	Giin_iter = vi_iter->getconnected_Gin();
////	Gjin_iter = vj_iter->getconnected_Gin();
////
////	Giin_new = (*Gjin_iter);
////	Gjin_new = (*Giin_iter);
////	Giin_new.setj(Giin_iter->getjflavor(),Giin_iter->getxj());
////	Gjin_new.setj(Gjin_iter->getjflavor(),Gjin_iter->getxj());
////	Giin_new.setconnected_vertices(Gjin_iter->getconnected_vi(),vi_iter);
////	Gjin_new.setconnected_vertices(Giin_iter->getconnected_vi(),vj_iter);
////
////	dNselfloop = 0;
////	if(Giin_iter==vi_iter->getconnected_Gout()) dNselfloop--;
////	if(Gjin_iter==vj_iter->getconnected_Gout()) dNselfloop--;
////	if(Giin_new.getconnected_vi()==Giin_new.getconnected_vj()) dNselfloop++;
////	if(Gjin_new.getconnected_vi()==Gjin_new.getconnected_vj()) dNselfloop++;
////	if(Ws.size()!=1 && Nselfloop+dNselfloop>NselfloopMax) return false;
////
////	// probability check
////	double probability_second
////		= -1.; // Fermionic sign factor
////
////	loopi_iter = vi_iter->getconnected_loop();
////	loopj_iter = vj_iter->getconnected_loop();
////
////	sameloop = (loopi_iter==loopj_iter);
////
////	trial.clear();
////
////	//setUWinter trial;
////	//list<line>::iterator UW_iter;
////	//list<fermionicloop>::iterator loop1_iter;
////	//list<fermionicloop>::iterator loop2_iter;
////
////	//list<list<vertex>::iterator> vset_split1, vset_split2, vset_combined;
////
////	//double loop1_valuerh_n, loop1_valuelh_n, loop2_valuerh_n, loop2_valuelh_n;
////	//double loop_valuerh_n, loop_valuelh_n;
////
////	if( sameloop ){
////		//cout<<"loop split case"<<endl<<flush;
////
////		vset_split1 = loopi_iter->getpath(vi_iter,vj_iter);
////		vset_split2 = loopi_iter->getpath(vj_iter,vi_iter);
////
////		int iloop = distance(fermionicloops.begin(),loopi_iter);
////
////		//cout<<"loop to be split index : "<<iloop<<endl;
////		//cout<<*loopi_iter;
////
////		fermionicloop loop1_pseudo(vset_split1);
////		fermionicloop loop2_pseudo(vset_split2);
////
////		loop1_pseudo.calvalue();
////		loop2_pseudo.calvalue();
////
////		cout<<"loop1_pseudo"<<endl<<loop1_pseudo<<endl;
////		cout<<"loop2_pseudo"<<endl<<loop2_pseudo<<endl;
////		cout<<"Giin_new.getvalueij() = "<<Giin_new.getvalueij()<<endl<<"Gjin_iter->getvalueij() = "<<Gjin_iter->getvalueij()<<endl;
////		cout<<"Giin_new.getvalueji() = "<<Giin_new.getvalueji()<<endl<<"Gjin_iter->getvalueji() = "<<Gjin_iter->getvalueji()<<endl;
////		cout<<"Gjin_new.getvalueij() = "<<Gjin_new.getvalueij()<<endl<<"Giin_iter->getvalueij() = "<<Giin_iter->getvalueij()<<endl;
////		cout<<"Gjin_new.getvalueji() = "<<Gjin_new.getvalueji()<<endl<<"Giin_iter->getvalueji() = "<<Giin_iter->getvalueji()<<endl;
////
////		loop1_valuerh_n = loop1_pseudo.getvaluerh() * Giin_new.getvalueij() / Gjin_iter->getvalueij();
////		loop1_valuelh_n = loop1_pseudo.getvaluelh() * Giin_new.getvalueji() / Gjin_iter->getvalueji();
////		loop2_valuerh_n = loop2_pseudo.getvaluerh() * Gjin_new.getvalueij() / Giin_iter->getvalueij();
////		loop2_valuelh_n = loop2_pseudo.getvaluelh() * Gjin_new.getvalueji() / Giin_iter->getvalueji();
////
////		trial.setNfermionicloops(fermionicloops.size()+1);
////		bool v1_in, v2_in;
////		list<list<vertex>::iterator>::iterator v_iiter;
////		list<vertex>::iterator v1_iter;
////		list<vertex>::iterator v2_iter;
////		for(UW_iter=Ws.begin(); UW_iter!=Ws.end(); ++UW_iter){
////			v1_iter = UW_iter->getconnected_vi();
////			v2_iter = UW_iter->getconnected_vj();
////			loop1_iter = v1_iter->getconnected_loop();
////			loop2_iter = v2_iter->getconnected_loop();
////
////			if( (loop1_iter==loopi_iter) && (loop2_iter==loopi_iter) ){
////				v1_in = false; v2_in = false;
////				for(v_iiter=vset_split1.begin(); v_iiter!=vset_split1.end(); ++v_iiter){
////					if(v1_iter==*v_iiter) v1_in = true;
////					else if(v2_iter==*v_iiter) v2_in = true;
////				}
////				if(v1_in^v2_in){
////					probability_second /= UW_iter->getvalueij();
////					if(v1_in) trial.push_back(*UW_iter,iloop,iloop+1);
////					else trial.push_back(*UW_iter,iloop+1,iloop);
////				}
////			}
////			if(loop1_iter!=loop2_iter){
////				if(loop1_iter==loopi_iter){
////					v1_in = false;
////					for(v_iiter=vset_split1.begin(); v_iiter!=vset_split1.end(); ++v_iiter){
////						if(v1_iter==*v_iiter){
////							v1_in = true;
////							break;
////						}
////					}
////					if(v1_in){
////						trial.push_back(
////							*UW_iter,
////							iloop,
////							gen_index_split(iloop,distance(fermionicloops.begin(),loop2_iter))
////						);
////					}
////					else{
////						trial.push_back(
////							*UW_iter,
////							iloop+1,
////							gen_index_split(iloop,distance(fermionicloops.begin(),loop2_iter))
////						);
////					}
////				}
////				else if(loop2_iter==loopi_iter){
////					v1_in = false;
////					for(v_iiter=vset_split1.begin(); v_iiter!=vset_split1.end(); ++v_iiter){
////						if(v2_iter==*v_iiter){
////							v1_in = true;
////							break;
////						}
////					}
////					if(v1_in){
////						trial.push_back(
////							*UW_iter,
////							gen_index_split(iloop,distance(fermionicloops.begin(),loop1_iter)),
////							iloop
////						);
////					}
////					else{
////						trial.push_back(
////							*UW_iter,
////							gen_index_split(iloop,distance(fermionicloops.begin(),loop1_iter)),
////							iloop+1
////						);
////					}
////				}
////				else{
////					trial.push_back(
////						*UW_iter,
////						gen_index_split(iloop,distance(fermionicloops.begin(),loop1_iter)),
////						gen_index_split(iloop,distance(fermionicloops.begin(),loop2_iter))
////					);
////				}
////			}
////		}
////
////		v1_iter = measuringline->getconnected_vi();
////		v2_iter = measuringline->getconnected_vj();
////		loop1_iter = v1_iter->getconnected_loop();
////		loop2_iter = v2_iter->getconnected_loop();
////		bool vi_1in = false; 
////		bool vj_1in = false;
////		for(v_iiter=vset_split1.begin(); v_iiter!=vset_split1.end(); ++v_iiter){
////			if(v1_iter==*v_iiter){ vi_1in = true; }
////			if(v2_iter==*v_iiter){ vj_1in = true; }
////		}
////		bool vi_2in = false; 
////		bool vj_2in = false;
////		for(v_iiter=vset_split2.begin(); v_iiter!=vset_split2.end(); ++v_iiter){
////			if(v1_iter==*v_iiter){ vi_2in = true; }
////			if(v2_iter==*v_iiter){ vj_2in = true; }
////		}
////
////		if( measuringline->gettype()==0 ){
////			if(vi_1in){
////				trial.fixloopflavor(iloop, miflavor);
////				probability_second *= loop1_valuerh_n * 0.5*(loop2_valuerh_n+loop2_valuelh_n)
////						/ loopi_iter->getvaluerh();
////			}
////			else if(vi_2in){
////				trial.fixloopflavor(iloop+1, miflavor);
////				probability_second *= 0.5*(loop1_valuerh_n+loop1_valuelh_n) * loop2_valuerh_n
////						/ loopi_iter->getvaluerh();
////			}
////			else{
////				trial.fixloopflavor(
////					gen_index_split(iloop,distance(fermionicloops.begin(),loop1_iter)),
////					miflavor
////				);
////				probability_second *= 0.5
////						* (loop1_valuerh_n+loop1_valuelh_n) * (loop2_valuerh_n+loop2_valuelh_n)
////						/ (loopi_iter->getvaluerh() + loopi_iter->getvaluelh());
////			}
////		}
////		else{
////			if(vi_1in) trial.fixloopflavor(iloop, miflavor);
////			else if(vi_2in) trial.fixloopflavor(iloop+1, miflavor);
////			else{
////				trial.fixloopflavor(
////					gen_index_split(iloop,distance(fermionicloops.begin(),loop1_iter)),
////					miflavor
////				);
////			}
////			if(vj_1in) trial.fixloopflavor(iloop, mjflavor);
////			else if(vj_2in) trial.fixloopflavor(iloop+1, mjflavor);
////			else{
////				trial.fixloopflavor(
////					gen_index_split(iloop,distance(fermionicloops.begin(),loop2_iter)),
////					mjflavor
////				);
////			}
////			probability_second *= 0.5
////					* (loop1_valuerh_n+loop1_valuelh_n) * (loop2_valuerh_n+loop2_valuelh_n)
////					/ (loopi_iter->getvaluerh() + loopi_iter->getvaluelh());
////		}
////
////		trial.calvalue();
////		probability_second *= trial.getvalue() / SetUWInter.getvalue();
////
////		//printconfiguration(cout);
////		//cout<<"trial.calvalue()"<<endl<<flush;
////		//cout<<trial<<endl;
////		//cout<<"symmetrized UWinter value: "<<trial.getvalue()<<endl;
////	}
////	else{
////		//cout<<"loop merge case"<<endl<<flush;
////		//cout<<*loopi_iter<<*loopj_iter<<endl;
////		// * Giin_new.getvalue() * Gjin_new.getvalue() 
////		// / Giin_iter->getvalue() / Gjin_iter->getvalue();
////		loop_valuerh_n = loopi_iter->getvaluerh()*loopj_iter->getvaluerh()
////					* Giin_new.getvalueij() * Gjin_new.getvalueij() 
////					/ Giin_iter->getvalueij() / Gjin_iter->getvalueij();
////		loop_valuelh_n = loopi_iter->getvaluelh()*loopj_iter->getvaluelh()
////					* Giin_new.getvalueji() * Gjin_new.getvalueji() 
////					/ Giin_iter->getvalueji() / Gjin_iter->getvalueji();
////
////		int iloop = distance(fermionicloops.begin(),loopi_iter);
////		int jloop = distance(fermionicloops.begin(),loopj_iter);
////
////		if( 
////			(SetUWInter.getisloopflavorfixed(iloop) && SetUWInter.getisloopflavorfixed(jloop))
////			&& 
////			(SetUWInter.getloopflavor(iloop)!=SetUWInter.getloopflavor(jloop))
////		){ return false; }
////
////		//cout<<"construct trial"<<endl;
////		trial.setNfermionicloops(fermionicloops.size()-1);
////		int i1loop, i2loop;
////		list<line>::iterator UW_iter;
////		for(UW_iter=Ws.begin(); UW_iter!=Ws.end(); ++UW_iter){
////			loop1_iter = UW_iter->getconnected_vi()->getconnected_loop();
////			loop2_iter = UW_iter->getconnected_vj()->getconnected_loop();
////
////			i1loop = gen_index_merge(iloop,jloop,distance(fermionicloops.begin(),loop1_iter));
////			i2loop = gen_index_merge(iloop,jloop,distance(fermionicloops.begin(),loop2_iter));
////
////			//cout<<"("<<i1loop<<";"<<distance(fermionicloops.begin(),loop1_iter)<<")"
////			//	<<"("<<i2loop<<";"<<distance(fermionicloops.begin(),loop2_iter)<<")"<<endl;
////
////			if(loop1_iter!=loop2_iter){
////				if(i1loop==i2loop){ probability_second *= UW_iter->getvalueij(); }
////				else{ trial.push_back(*UW_iter,i1loop,i2loop); }
////			}
////		}
////		//cout<<"mark measuringline"<<endl;
////		loop1_iter = measuringline->getconnected_vi()->getconnected_loop();
////		loop2_iter = measuringline->getconnected_vj()->getconnected_loop();
////		i1loop = gen_index_merge(iloop,jloop,distance(fermionicloops.begin(),loop1_iter));
////		i2loop = gen_index_merge(iloop,jloop,distance(fermionicloops.begin(),loop2_iter));
////		trial.fixloopflavor( i1loop, miflavor );
////		trial.fixloopflavor( i2loop, mjflavor );
////
////		if(measuringline->gettype()==0){
////			if(loop1_iter==loopi_iter){
////				probability_second *= loop_valuerh_n 
////						/ loopi_iter->getvaluerh() 
////						/ 0.5/(loopj_iter->getvaluerh()+loopj_iter->getvaluelh());
////			}
////			else if(loop1_iter==loopj_iter){
////				probability_second *= loop_valuerh_n 
////						/ 0.5/(loopi_iter->getvaluerh()+loopi_iter->getvaluelh())
////						/ loopj_iter->getvaluerh();
////			}
////			else{
////				probability_second *= (loop_valuerh_n+loop_valuelh_n) 
////						/ 0.5
////						/ (loopi_iter->getvaluerh()+loopi_iter->getvaluelh())
////						/ (loopj_iter->getvaluerh()+loopj_iter->getvaluelh());
////			}
////		}
////		else{
////			probability_second *= (loop_valuerh_n+loop_valuelh_n) 
////					/ 0.5
////					/ (loopi_iter->getvaluerh()+loopi_iter->getvaluelh())
////					/ (loopj_iter->getvaluerh()+loopj_iter->getvaluelh());
////		}
////		//cout<<"calculate trial value"<<endl;
////		//cout<<"symmetrized UWinter value: "<<trial.getvalue()<<endl;
////		trial.calvalue();
////
////		probability_second *= trial.getvalue() / SetUWInter.getvalue();
////
////	}
////
////	printconfiguration(cout);
////	cout<<"* reconnect info *"<<endl;
////	if(sameloop){
////		cout<<"  loop split case"<<endl;
////		cout<<"    loop1_valuerh_n = "<<loop1_valuerh_n<<endl;
////		cout<<"    loop1_valuelh_n = "<<loop1_valuelh_n<<endl;
////		cout<<"    loop2_valuerh_n = "<<loop2_valuerh_n<<endl;
////		cout<<"    loop2_valuelh_n = "<<loop2_valuelh_n<<endl;
////	}
////	else{
////		cout<<"  loop merge case"<<endl;
////		cout<<"    loop_valuerh_n = "<<loop_valuerh_n<<endl;
////		cout<<"    loop_valuelh_n = "<<loop_valuelh_n<<endl;
////	}
////
////	//cout<<"  probability = "<<probability<<endl;
////	//cout<<"  trial info. "<<endl<<trial<<endl;
////	//cout<<"  SetUWInter.getvalue() = "<<SetUWInter.getvalue()<<endl;
////	//cout<<"*******************"<<endl;
////
////	//int signfactor;
////	if(probability_second<0.) signfactor = -1;
////	else signfactor = 1;
////
////	{
////		//cout<<"compactness check"<<endl;
////		//if(!(checkconnectivityforreconnect(Giin_iter,Gjin_iter))) return 0;
////		if(!(checkirreducibilityforreconnect(Giin_iter,Gjin_iter))) return false;
////		if(!(checkcompactnessforreconnect(Giin_iter,Gjin_iter))) return false;
////
////		//cout<<"start to update quantities"<<endl;
////		bool measuring_line_involved = false;
////		if( Giin_iter->getphysical()^Gjin_iter->getphysical() )
////			measuring_line_involved = true;
////
////		vi_iter->setconnected_Gin(Gjin_iter);
////		vj_iter->setconnected_Gin(Giin_iter);
////		Giin_iter->setxj(vj_iter->getx());
////		Gjin_iter->setxj(vi_iter->getx());
////		Giin_iter->setconnected_vj(vj_iter);
////		Gjin_iter->setconnected_vj(vi_iter);
////
////		//cout<<"loop quantities update"<<endl;
////		list<list<vertex>::iterator>::iterator v_iiter;
////		if( sameloop ){
////			//cout<<"loop split case"<<endl;
////			loopi_iter->set(vset_split1);
////			fermionicloops.insert(next(loopi_iter),fermionicloop(vset_split2));
////			loopj_iter = next(loopi_iter);
////
////			loopi_iter->setvalue(loop1_valuerh_n,loop1_valuelh_n);
////			loopj_iter->setvalue(loop2_valuerh_n,loop2_valuelh_n);
////
////			for(v_iiter=vset_split1.begin(); v_iiter!=vset_split1.end(); ++v_iiter)
////				(*v_iiter)->setconnected_loop(loopi_iter);
////			for(v_iiter=vset_split2.begin(); v_iiter!=vset_split2.end(); ++v_iiter)
////				(*v_iiter)->setconnected_loop(loopj_iter);
////		}
////		else{
////			//cout<<"loop merge case"<<endl;
////			//cout<<*loopi_iter<<*loopj_iter<<endl;
////			vset_combined = loopi_iter->getpath(vi_iter,vi_iter);
////			vset_combined.splice(vset_combined.end(),loopj_iter->getpath(vj_iter,vj_iter));
////			int iloop = distance(fermionicloops.begin(),loopi_iter);
////			int jloop = distance(fermionicloops.begin(),loopj_iter);
////			//cout<<"vset_combined"<<endl;
////			//for(v_iiter=vset_combined.begin(); v_iiter!=vset_combined.end(); ++v_iiter)
////			//	cout<<(**v_iiter);
////			//cout<<endl;
////			if(iloop<jloop){
////				loopi_iter->set(vset_combined);
////				fermionicloops.erase(loopj_iter);
////
////				loopi_iter->setvalue(loop_valuerh_n,loop_valuelh_n);
////
////				for(v_iiter=vset_combined.begin(); v_iiter!=vset_combined.end(); ++v_iiter)
////					(*v_iiter)->setconnected_loop(loopi_iter);
////			}
////			else{
////				loopj_iter->set(vset_combined);
////				fermionicloops.erase(loopi_iter);
////
////				loopj_iter->setvalue(loop_valuerh_n,loop_valuelh_n);
////
////				for(v_iiter=vset_combined.begin(); v_iiter!=vset_combined.end(); ++v_iiter)
////					(*v_iiter)->setconnected_loop(loopj_iter);
////			}
////		}
////		SetUWInter = trial;
////		//updateSetUWInter();
////		//SetUWInter.setvalue(trial.getvalue());
////
////		sign *= signfactor;
////		Nselfloop += dNselfloop;
////
////		//printconfiguration(cout);
////		//exit(EXIT_FAILURE);
////	}
////
////	if( abs(probability_first*probability_second - 1.)<1.0e-10 ){
////		cout<<"probability_first = "<<probability_first<<endl;
////		cout<<"probability_second = "<<probability_second<<endl;
////		cout<<"probability_first*probability_second = "<<probability_first*probability_second<<endl;
////		return true;
////	}
////	else{
////		cout<<"***** ALERT! DETAILED BALANCE CONDITION FAILS! *****"<<endl;
////		cout<<"probability_first = "<<probability_first<<endl;
////		cout<<"probability_second = "<<probability_second<<endl;
////		cout<<"probability_first*probability_second = "<<probability_first*probability_second<<endl;
////		printconfiguration(cout);
////		return false;
////	}
////
////};*/
////
////
/////*bool diagram::checkdetailedbalance_reconnect(){
////	int Nvs = vertices.size();
////
////	int iv, jv;
////	select_two_index(Nvs, iv, jv);
////
////	list<vertex>::iterator vi_iter = next(vertices.begin(),iv);
////	list<vertex>::iterator vj_iter = next(vertices.begin(),jv);
////
////	list<line>::iterator Giin_iter = vi_iter->getconnected_Gin();
////	list<line>::iterator Gjin_iter = vj_iter->getconnected_Gin();
////
////	line Giin_new(*Giin_iter);
////	line Gjin_new(*Gjin_iter);
////	Giin_new.setj(Gjin_iter->getjflavor(),Gjin_iter->getxj());
////	Gjin_new.setj(Giin_iter->getjflavor(),Giin_iter->getxj());
////	Giin_new.setconnected_vertices(Giin_iter->getconnected_vi(),vj_iter);
////	Gjin_new.setconnected_vertices(Gjin_iter->getconnected_vi(),vi_iter);
////
////	int dNselfloop = 0;
////	if(Giin_iter==vi_iter->getconnected_Gout()) dNselfloop--;
////	if(Gjin_iter==vj_iter->getconnected_Gout()) dNselfloop--;
////	if(Giin_new.getconnected_vi()==Giin_new.getconnected_vj()) dNselfloop++;
////	if(Gjin_new.getconnected_vi()==Gjin_new.getconnected_vj()) dNselfloop++;
////	if(Ws.size()!=1 && Nselfloop+dNselfloop>NselfloopMax) return true;
////
////	double probability_first
////		= -1. // Fermionic sign factor
////		* Giin_new.getvalue() * Gjin_new.getvalue() 
////		/ Giin_iter->getvalue() / Gjin_iter->getvalue();
////
////	list<fermionicloop>::iterator loopi_iter = vi_iter->getconnected_loop();
////	list<fermionicloop>::iterator loopj_iter = vj_iter->getconnected_loop();
////
////	bool sameloop = (loopi_iter==loopj_iter);
////	int iloop, jloop;
////
////	setUWinter trial;
////	list<line>::iterator UW_iter;
////	list<fermionicloop>::iterator loop1_iter, loop2_iter;
////
////	list<list<vertex>::iterator> vset_split1, vset_split2, vset_combined;
////	if( sameloop ){
////		vset_split1 = loopi_iter->getpath(vi_iter,vj_iter);
////		vset_split2 = loopi_iter->getpath(vj_iter,vi_iter);
////
////		iloop = distance(fermionicloops.begin(),loopi_iter);
////
////		bool v1_in, v2_in;
////		trial.setNfermionicloops(fermionicloops.size()+1);
////		list<list<vertex>::iterator>::iterator v_iiter;
////		list<vertex>::iterator v1_iter;
////		list<vertex>::iterator v2_iter;
////		for(UW_iter=Ws.begin(); UW_iter!=Ws.end(); ++UW_iter){
////			v1_iter = UW_iter->getconnected_vi();
////			v2_iter = UW_iter->getconnected_vj();
////			loop1_iter = v1_iter->getconnected_loop();
////			loop2_iter = v2_iter->getconnected_loop();
////
////			if( (loop1_iter==loopi_iter) && (loop2_iter==loopi_iter) ){
////				v1_in = false; v2_in = false;
////				for(v_iiter=vset_split1.begin(); v_iiter!=vset_split1.end(); ++v_iiter){
////					if(v1_iter==*v_iiter) v1_in = true;
////					else if(v2_iter==*v_iiter) v2_in = true;
////				}
////				if(v1_in^v2_in){
////					probability_first /= UW_iter->getvalue();
////					if(v1_in) trial.push_back(*UW_iter,iloop,iloop+1);
////					else trial.push_back(*UW_iter,iloop+1,iloop);
////				}
////			}
////			if(loop1_iter!=loop2_iter){
////				if(loop1_iter==loopi_iter){
////					v1_in = false;
////					for(v_iiter=vset_split1.begin(); v_iiter!=vset_split1.end(); ++v_iiter){
////						if(v1_iter==*v_iiter){
////							v1_in = true;
////							break;
////						}
////					}
////					if(v1_in){
////						trial.push_back(
////							*UW_iter,
////							iloop,
////							gen_index_split(iloop,distance(fermionicloops.begin(),loop2_iter))
////						);
////					}
////					else{
////						trial.push_back(
////							*UW_iter,
////							iloop+1,
////							gen_index_split(iloop,distance(fermionicloops.begin(),loop2_iter))
////						);
////					}
////				}
////				else if(loop2_iter==loopi_iter){
////					v1_in = false;
////					for(v_iiter=vset_split1.begin(); v_iiter!=vset_split1.end(); ++v_iiter){
////						if(v2_iter==*v_iiter){
////							v1_in = true;
////							break;
////						}
////					}
////					if(v1_in){
////						trial.push_back(
////							*UW_iter,
////							gen_index_split(iloop,distance(fermionicloops.begin(),loop1_iter)),
////							iloop
////						);
////					}
////					else{
////						trial.push_back(
////							*UW_iter,
////							gen_index_split(iloop,distance(fermionicloops.begin(),loop1_iter)),
////							iloop+1
////						);
////					}
////				}
////				else{
////					trial.push_back(
////						*UW_iter,
////						gen_index_split(iloop,distance(fermionicloops.begin(),loop1_iter)),
////						gen_index_split(iloop,distance(fermionicloops.begin(),loop2_iter))
////					);
////				}
////			}
////		}
////
////		v1_iter = measuringline->getconnected_vi();
////		v2_iter = measuringline->getconnected_vj();
////		loop1_iter = v1_iter->getconnected_loop();
////		loop2_iter = v2_iter->getconnected_loop();
////		bool vi_1in = false; 
////		bool vj_1in = false;
////		for(v_iiter=vset_split1.begin(); v_iiter!=vset_split1.end(); ++v_iiter){
////			if(v1_iter==*v_iiter){ vi_1in = true; }
////			if(v2_iter==*v_iiter){ vj_1in = true; }
////		}
////		bool vi_2in = false; 
////		bool vj_2in = false;
////		for(v_iiter=vset_split2.begin(); v_iiter!=vset_split2.end(); ++v_iiter){
////			if(v1_iter==*v_iiter){ vi_2in = true; }
////			if(v2_iter==*v_iiter){ vj_2in = true; }
////		}
////
////		if( measuringline->gettype()==0 ){
////			if(vi_1in) trial.fixloopflavor(iloop, miflavor);
////			else if(vi_2in) trial.fixloopflavor(iloop+1, miflavor);
////			else{
////				trial.fixloopflavor(
////					gen_index_split(iloop,distance(fermionicloops.begin(),loop1_iter)),
////					miflavor
////				);
////			}
////		}
////		else{
////			if(vi_1in) trial.fixloopflavor(iloop, miflavor);
////			else if(vi_2in) trial.fixloopflavor(iloop+1, miflavor);
////			else{
////				trial.fixloopflavor(
////					gen_index_split(iloop,distance(fermionicloops.begin(),loop1_iter)),
////					miflavor
////				);
////			}
////			if(vj_1in) trial.fixloopflavor(iloop, mjflavor);
////			else if(vj_2in) trial.fixloopflavor(iloop+1, mjflavor);
////			else{
////				trial.fixloopflavor(
////					gen_index_split(iloop,distance(fermionicloops.begin(),loop2_iter)),
////					mjflavor
////				);
////			}
////		}
////
////		trial.calvalue();
////		probability_first *= trial.getvalue()/SetUWInter.getvalue();
////	}
////	else{
////		int iloop = distance(fermionicloops.begin(),loopi_iter);
////		int jloop = distance(fermionicloops.begin(),loopj_iter);
////
////		trial.setNfermionicloops(fermionicloops.size()-1);
////		int i1loop, i2loop;
////		list<line>::iterator UW_iter;
////		for(UW_iter=Ws.begin(); UW_iter!=Ws.end(); ++UW_iter){
////			loop1_iter = UW_iter->getconnected_vi()->getconnected_loop();
////			loop2_iter = UW_iter->getconnected_vj()->getconnected_loop();
////
////			i1loop = gen_index_merge(iloop,jloop,distance(fermionicloops.begin(),loop1_iter));
////			i2loop = gen_index_merge(iloop,jloop,distance(fermionicloops.begin(),loop2_iter));
////
////			if(loop1_iter!=loop2_iter){
////				if(i1loop==i2loop){ probability_first *= UW_iter->getvalue(); }
////				else{ trial.push_back(*UW_iter,i1loop,i2loop); }
////			}
////		}
////		loop1_iter = measuringline->getconnected_vi()->getconnected_loop();
////		loop2_iter = measuringline->getconnected_vj()->getconnected_loop();
////		i1loop = gen_index_merge(iloop,jloop,distance(fermionicloops.begin(),loop1_iter));
////		i2loop = gen_index_merge(iloop,jloop,distance(fermionicloops.begin(),loop2_iter));
////		trial.fixloopflavor( i1loop, miflavor );
////		trial.fixloopflavor( i2loop, mjflavor );
////
////		trial.calvalue();
////		probability_first *= trial.getvalue()/SetUWInter.getvalue();
////
////	}
////
////	int signfactor;
////	if(probability_first<0.) signfactor = -1;
////	else signfactor = 1;
////
////	if(abs(probability_first)>1.0e-10){
////		//if(!(checkconnectivityforreconnect(Giin_iter,Gjin_iter))) return 0;
////		if(!(checkirreducibilityforreconnect(Giin_iter,Gjin_iter))) return true;
////		if(!(checkcompactnessforreconnect(Giin_iter,Gjin_iter))) return true;
////
////		bool measuring_line_involved = false;
////		if( Giin_iter->getphysical()^Gjin_iter->getphysical() )
////			measuring_line_involved = true;
////
////		vi_iter->setconnected_Gin(Gjin_iter);
////		vj_iter->setconnected_Gin(Giin_iter);
////		Giin_iter->setxj(vj_iter->getx());
////		Gjin_iter->setxj(vi_iter->getx());
////		Giin_iter->setconnected_vj(vj_iter);
////		Gjin_iter->setconnected_vj(vi_iter);
////
////		list<list<vertex>::iterator>::iterator v_iiter;
////		if( sameloop ){
////			loopi_iter->set(vset_split1);
////			fermionicloops.insert(next(loopi_iter),fermionicloop(vset_split2));
////			loopj_iter = next(loopi_iter);
////			for(v_iiter=vset_split1.begin(); v_iiter!=vset_split1.end(); ++v_iiter)
////				(*v_iiter)->setconnected_loop(loopi_iter);
////			for(v_iiter=vset_split2.begin(); v_iiter!=vset_split2.end(); ++v_iiter)
////				(*v_iiter)->setconnected_loop(loopj_iter);
////		}
////		else{
////			vset_combined = loopi_iter->getpath(vi_iter,vi_iter);
////			vset_combined.splice(vset_combined.end(),loopj_iter->getpath(vj_iter,vj_iter));
////			int iloop = distance(fermionicloops.begin(),loopi_iter);
////			int jloop = distance(fermionicloops.begin(),loopj_iter);
////			if(iloop<jloop){
////				loopi_iter->set(vset_combined);
////				fermionicloops.erase(loopj_iter);
////				for(v_iiter=vset_combined.begin(); v_iiter!=vset_combined.end(); ++v_iiter)
////					(*v_iiter)->setconnected_loop(loopi_iter);
////			}
////			else{
////				loopj_iter->set(vset_combined);
////				fermionicloops.erase(loopi_iter);
////				for(v_iiter=vset_combined.begin(); v_iiter!=vset_combined.end(); ++v_iiter)
////					(*v_iiter)->setconnected_loop(loopj_iter);
////			}
////		}
////		SetUWInter = trial;
////
////		sign *= signfactor;
////		Nselfloop += dNselfloop;
////	}
////	else{ return true; }
////
////	vi_iter = next(vertices.begin(),iv);
////	vj_iter = next(vertices.begin(),jv);
////
////	Giin_iter = vi_iter->getconnected_Gin();
////	Gjin_iter = vj_iter->getconnected_Gin();
////
////	Giin_new = *Giin_iter;
////	Gjin_new = *Gjin_iter;
////	Giin_new.setj(Gjin_iter->getjflavor(),Gjin_iter->getxj());
////	Gjin_new.setj(Giin_iter->getjflavor(),Giin_iter->getxj());
////	Giin_new.setconnected_vertices(Giin_iter->getconnected_vi(),vj_iter);
////	Gjin_new.setconnected_vertices(Gjin_iter->getconnected_vi(),vi_iter);
////
////	dNselfloop = 0;
////	if(Giin_iter==vi_iter->getconnected_Gout()) dNselfloop--;
////	if(Gjin_iter==vj_iter->getconnected_Gout()) dNselfloop--;
////	if(Giin_new.getconnected_vi()==Giin_new.getconnected_vj()) dNselfloop++;
////	if(Gjin_new.getconnected_vi()==Gjin_new.getconnected_vj()) dNselfloop++;
////	if(Ws.size()!=1 && Nselfloop+dNselfloop>NselfloopMax) return true;
////
////	double probability_second
////		= -1. // Fermionic sign factor
////		* Giin_new.getvalue() * Gjin_new.getvalue() 
////		/ Giin_iter->getvalue() / Gjin_iter->getvalue();
////
////	loopi_iter = vi_iter->getconnected_loop();
////	loopj_iter = vj_iter->getconnected_loop();
////
////	sameloop = (loopi_iter==loopj_iter);
////
////	trial.clear();
////	if( sameloop ){
////		vset_split1 = loopi_iter->getpath(vi_iter,vj_iter);
////		vset_split2 = loopi_iter->getpath(vj_iter,vi_iter);
////
////		iloop = distance(fermionicloops.begin(),loopi_iter);
////
////		bool v1_in, v2_in;
////		trial.setNfermionicloops(fermionicloops.size()+1);
////		list<list<vertex>::iterator>::iterator v_iiter;
////		list<vertex>::iterator v1_iter;
////		list<vertex>::iterator v2_iter;
////		for(UW_iter=Ws.begin(); UW_iter!=Ws.end(); ++UW_iter){
////			v1_iter = UW_iter->getconnected_vi();
////			v2_iter = UW_iter->getconnected_vj();
////			loop1_iter = v1_iter->getconnected_loop();
////			loop2_iter = v2_iter->getconnected_loop();
////
////			if( (loop1_iter==loopi_iter) && (loop2_iter==loopi_iter) ){
////				v1_in = false; v2_in = false;
////				for(v_iiter=vset_split1.begin(); v_iiter!=vset_split1.end(); ++v_iiter){
////					if(v1_iter==*v_iiter) v1_in = true;
////					else if(v2_iter==*v_iiter) v2_in = true;
////				}
////				if(v1_in^v2_in){
////					probability_second /= UW_iter->getvalue();
////					if(v1_in) trial.push_back(*UW_iter,iloop,iloop+1);
////					else trial.push_back(*UW_iter,iloop+1,iloop);
////				}
////			}
////			if(loop1_iter!=loop2_iter){
////				if(loop1_iter==loopi_iter){
////					v1_in = false;
////					for(v_iiter=vset_split1.begin(); v_iiter!=vset_split1.end(); ++v_iiter){
////						if(v1_iter==*v_iiter){
////							v1_in = true;
////							break;
////						}
////					}
////					if(v1_in){
////						trial.push_back(
////							*UW_iter,
////							iloop,
////							gen_index_split(iloop,distance(fermionicloops.begin(),loop2_iter))
////						);
////					}
////					else{
////						trial.push_back(
////							*UW_iter,
////							iloop+1,
////							gen_index_split(iloop,distance(fermionicloops.begin(),loop2_iter))
////						);
////					}
////				}
////				else if(loop2_iter==loopi_iter){
////					v1_in = false;
////					for(v_iiter=vset_split1.begin(); v_iiter!=vset_split1.end(); ++v_iiter){
////						if(v2_iter==*v_iiter){
////							v1_in = true;
////							break;
////						}
////					}
////					if(v1_in){
////						trial.push_back(
////							*UW_iter,
////							gen_index_split(iloop,distance(fermionicloops.begin(),loop1_iter)),
////							iloop
////						);
////					}
////					else{
////						trial.push_back(
////							*UW_iter,
////							gen_index_split(iloop,distance(fermionicloops.begin(),loop1_iter)),
////							iloop+1
////						);
////					}
////				}
////				else{
////					trial.push_back(
////						*UW_iter,
////						gen_index_split(iloop,distance(fermionicloops.begin(),loop1_iter)),
////						gen_index_split(iloop,distance(fermionicloops.begin(),loop2_iter))
////					);
////				}
////			}
////		}
////
////		v1_iter = measuringline->getconnected_vi();
////		v2_iter = measuringline->getconnected_vj();
////		loop1_iter = v1_iter->getconnected_loop();
////		loop2_iter = v2_iter->getconnected_loop();
////		bool vi_1in = false; 
////		bool vj_1in = false;
////		for(v_iiter=vset_split1.begin(); v_iiter!=vset_split1.end(); ++v_iiter){
////			if(v1_iter==*v_iiter){ vi_1in = true; }
////			if(v2_iter==*v_iiter){ vj_1in = true; }
////		}
////		bool vi_2in = false; 
////		bool vj_2in = false;
////		for(v_iiter=vset_split2.begin(); v_iiter!=vset_split2.end(); ++v_iiter){
////			if(v1_iter==*v_iiter){ vi_2in = true; }
////			if(v2_iter==*v_iiter){ vj_2in = true; }
////		}
////
////		if( measuringline->gettype()==0 ){
////			if(vi_1in) trial.fixloopflavor(iloop, miflavor);
////			else if(vi_2in) trial.fixloopflavor(iloop+1, miflavor);
////			else{
////				trial.fixloopflavor(
////					gen_index_split(iloop,distance(fermionicloops.begin(),loop1_iter)),
////					miflavor
////				);
////			}
////		}
////		else{
////			if(vi_1in) trial.fixloopflavor(iloop, miflavor);
////			else if(vi_2in) trial.fixloopflavor(iloop+1, miflavor);
////			else{
////				trial.fixloopflavor(
////					gen_index_split(iloop,distance(fermionicloops.begin(),loop1_iter)),
////					miflavor
////				);
////			}
////			if(vj_1in) trial.fixloopflavor(iloop, mjflavor);
////			else if(vj_2in) trial.fixloopflavor(iloop+1, mjflavor);
////			else{
////				trial.fixloopflavor(
////					gen_index_split(iloop,distance(fermionicloops.begin(),loop2_iter)),
////					mjflavor
////				);
////			}
////		}
////
////		trial.calvalue();
////		probability_second *= trial.getvalue()/SetUWInter.getvalue();
////	}
////	else{
////		int iloop = distance(fermionicloops.begin(),loopi_iter);
////		int jloop = distance(fermionicloops.begin(),loopj_iter);
////
////		trial.setNfermionicloops(fermionicloops.size()-1);
////		int i1loop, i2loop;
////		list<line>::iterator UW_iter;
////		for(UW_iter=Ws.begin(); UW_iter!=Ws.end(); ++UW_iter){
////			loop1_iter = UW_iter->getconnected_vi()->getconnected_loop();
////			loop2_iter = UW_iter->getconnected_vj()->getconnected_loop();
////
////			i1loop = gen_index_merge(iloop,jloop,distance(fermionicloops.begin(),loop1_iter));
////			i2loop = gen_index_merge(iloop,jloop,distance(fermionicloops.begin(),loop2_iter));
////
////			if(loop1_iter!=loop2_iter){
////				if(i1loop==i2loop){ probability_second *= UW_iter->getvalue(); }
////				else{ trial.push_back(*UW_iter,i1loop,i2loop); }
////			}
////		}
////		loop1_iter = measuringline->getconnected_vi()->getconnected_loop();
////		loop2_iter = measuringline->getconnected_vj()->getconnected_loop();
////		i1loop = gen_index_merge(iloop,jloop,distance(fermionicloops.begin(),loop1_iter));
////		i2loop = gen_index_merge(iloop,jloop,distance(fermionicloops.begin(),loop2_iter));
////		trial.fixloopflavor( i1loop, miflavor );
////		trial.fixloopflavor( i2loop, mjflavor );
////
////		trial.calvalue();
////		probability_second *= trial.getvalue()/SetUWInter.getvalue();
////
////	}
////
////	if(probability_second<0.) signfactor = -1;
////	else signfactor = 1;
////
////	{
////		if(!(checkirreducibilityforreconnect(Giin_iter,Gjin_iter))){
////			cout<<"***** ALERT! DETAILED BALANCE CONDITION FAILS! *****"<<endl;
////			cout<<"fail to pass the irreducibility check when going back to original conf."<<endl;
////			return false;
////		}
////		if(!(checkcompactnessforreconnect(Giin_iter,Gjin_iter))){
////			cout<<"***** ALERT! DETAILED BALANCE CONDITION FAILS! *****"<<endl;
////			cout<<"fail to pass the compactness check when going back to original conf."<<endl;
////			return false;
////		}
////
////		bool measuring_line_involved = false;
////		if( Giin_iter->getphysical()^Gjin_iter->getphysical() )
////			measuring_line_involved = true;
////
////		vi_iter->setconnected_Gin(Gjin_iter);
////		vj_iter->setconnected_Gin(Giin_iter);
////		Giin_iter->setxj(vj_iter->getx());
////		Gjin_iter->setxj(vi_iter->getx());
////		Giin_iter->setconnected_vj(vj_iter);
////		Gjin_iter->setconnected_vj(vi_iter);
////
////		list<list<vertex>::iterator>::iterator v_iiter;
////		if( sameloop ){
////			loopi_iter->set(vset_split1);
////			fermionicloops.insert(next(loopi_iter),fermionicloop(vset_split2));
////			loopj_iter = next(loopi_iter);
////			for(v_iiter=vset_split1.begin(); v_iiter!=vset_split1.end(); ++v_iiter)
////				(*v_iiter)->setconnected_loop(loopi_iter);
////			for(v_iiter=vset_split2.begin(); v_iiter!=vset_split2.end(); ++v_iiter)
////				(*v_iiter)->setconnected_loop(loopj_iter);
////		}
////		else{
////			vset_combined = loopi_iter->getpath(vi_iter,vi_iter);
////			vset_combined.splice(vset_combined.end(),loopj_iter->getpath(vj_iter,vj_iter));
////			int iloop = distance(fermionicloops.begin(),loopi_iter);
////			int jloop = distance(fermionicloops.begin(),loopj_iter);
////			if(iloop<jloop){
////				loopi_iter->set(vset_combined);
////				fermionicloops.erase(loopj_iter);
////				for(v_iiter=vset_combined.begin(); v_iiter!=vset_combined.end(); ++v_iiter)
////					(*v_iiter)->setconnected_loop(loopi_iter);
////			}
////			else{
////				loopj_iter->set(vset_combined);
////				fermionicloops.erase(loopi_iter);
////				for(v_iiter=vset_combined.begin(); v_iiter!=vset_combined.end(); ++v_iiter)
////					(*v_iiter)->setconnected_loop(loopj_iter);
////			}
////		}
////		SetUWInter = trial;
////
////		sign *= signfactor;
////		Nselfloop += dNselfloop;
////	}
////
////	if( abs(probability_first*probability_second - 1.)<1.0e-10 ){
////		cout<<"probability_first = "<<probability_first<<endl;
////		cout<<"probability_second = "<<probability_second<<endl;
////		cout<<"probability_first*probability_second = "<<probability_first*probability_second<<endl;
////		return true;
////	}
////	else{
////		cout<<"***** ALERT! DETAILED BALANCE CONDITION FAILS! *****"<<endl;
////		cout<<"probability_first = "<<probability_first<<endl;
////		cout<<"probability_second = "<<probability_second<<endl;
////		cout<<"probability_first*probability_second = "<<probability_first*probability_second<<endl;
////		printconfiguration(cout);
////		return false;
////	}
////};*/
////
/////*bool diagram::checkdetailedbalance_interswap(){
////	cout<<"** interswap detailed balance check **"<<endl;
////
////	/////////////////////
////	// first interswap //
////	/////////////////////
////
////	cout<<"*  first interswap  *"<<endl;
////
////	int index_j;
////	int baredNorder;
////	list<line>::iterator line_j_iter;
////	if(measuringline->gettype()==0){
////		index_j = static_cast<int>( Ws.size()*unidist(mt) );
////		line_j_iter = next(Ws.begin(),index_j);
////		baredNorder = -1;
////	}
////	else{
////		index_j = static_cast<int>( Gs.size()*unidist(mt) );
////		line_j_iter = next(Gs.begin(),index_j);
////		baredNorder = +1;
////	}
////
////	int dNorder = baredNorder;
////	if(Norder+dNorder>(this->Nordermax)){ return true; }
////
////	line line_i_new(*measuringline);
////	line_i_new.setphysical(true);
////
////	double proposal = 1.0;
////	double probability_first = 
////		-1 // Fermionic sign factor
////		* Ri_N[Norder] / Ri_N[Norder+dNorder];
////	
////	//int miflavor_n, mjflavor_n;
////	//int iloop, jloop;
////	//list<vertex>::iterator vi_iter, vj_iter;
////	//list<fermionicloop>::iterator loopi_iter, loopj_iter;
////
////	int miflavor_n, mjflavor_n;
////	int iloop, jloop;
////	list<vertex>::iterator vi_iter, vj_iter;
////	list<fermionicloop>::iterator loop_iter;
////	list<fermionicloop>::iterator loopi_iter, loopj_iter;
////
////	cout<<"measuring line: "<<*measuringline<<endl;
////	cout<<"target line: "<<*line_j_iter<<endl;
////	cout<<"line_i_new: "<<line_i_new<<endl;
////
////	cout<<"original SetUWInter"<<endl<<SetUWInter;
////
////	setUWinter trial = SetUWInter;
////	if(measuringline->gettype()==0){
////		vi_iter = line_j_iter->getconnected_vi();
////		vj_iter = line_j_iter->getconnected_vj();
////		loopi_iter = vi_iter->getconnected_loop();
////		loopj_iter = vj_iter->getconnected_loop();
////		if(loopi_iter==loopj_iter){
////			iloop  = distance(fermionicloops.begin(), loopi_iter);
////
////			miflavor_n = Nflavor*unidist(mt);
////			mjflavor_n = miflavor_n;
////
////			trial.releaseloopflavor();
////			trial.fixloopflavor(iloop,miflavor_n);
////
////			probability_first *= line_i_new.getvalue() / line_j_iter->getvalue();
////		}
////		else{
////			iloop  = distance(fermionicloops.begin(), loopi_iter);
////			jloop  = distance(fermionicloops.begin(), loopj_iter);
////
////			proposal *= 2.;
////			miflavor_n = Nflavor*unidist(mt);
////			mjflavor_n = Nflavor*unidist(mt);
////
////			trial.setphysical(*line_j_iter,false);
////			trial.releaseloopflavor();
////			trial.fixloopflavor(iloop,miflavor_n);
////			trial.fixloopflavor(jloop,mjflavor_n);
////			trial.calvalue();
////
////			probability_first *= line_i_new.getvalue() * trial.getvalue() / SetUWInter.getvalue();
////		}
////	}
////	else{
////		proposal /= 2.;
////		miflavor_n = Nflavor*unidist(mt);
////		mjflavor_n = miflavor_n;
////
////		vi_iter = line_j_iter->getconnected_vi();
////		loopi_iter = vi_iter->getconnected_loop();
////		iloop  = distance(fermionicloops.begin(), loopi_iter);
////
////		if( measuringline->getconnected_vi()->getconnected_loop() 
////			== measuringline->getconnected_vj()->getconnected_loop() )
////			probability_first *= line_i_new.getvalue();
////		else
////			trial.setphysical(*measuringline,true);
////
////		trial.releaseloopflavor();
////		trial.fixloopflavor(iloop,miflavor_n);
////		trial.calvalue();
////		probability_first *= 1./line_j_iter->getvalue() * trial.getvalue() / SetUWInter.getvalue();
////	}
////
////	double loop_valuerh_n, loop_valuelh_n;
////
////	bool ismeasuringlineGline = (measuringline->gettype()==0);
////	setUWinter trial = SetUWInter;
////	if(ismeasuringlineGline){
////		loop_iter = measuringline->getconnected_vi()->getconnected_loop();
////		loop_valuerh_n = loop_iter->getvaluerh()*line_i_new.getvalueij();
////		loop_valuelh_n = loop_iter->getvaluelh()*line_i_new.getvalueji();
////
////		vi_iter = line_j_iter->getconnected_vi();
////		vj_iter = line_j_iter->getconnected_vj();
////		loopi_iter = vi_iter->getconnected_loop();
////		loopj_iter = vj_iter->getconnected_loop();
////		if(loopi_iter==loopj_iter){
////			iloop  = distance(fermionicloops.begin(), loopi_iter);
////
////			miflavor_n = Nflavor*unidist(mt);
////			mjflavor_n = miflavor_n;
////
////			trial.releaseloopflavor();
////			trial.fixloopflavor(iloop,miflavor_n);
////			trial.calvalue();
////
////			probability_first *= 0.5*(loop_valuerh_n + loop_valuelh_n) / loop_iter->getvaluerh() 
////					/ line_j_iter->getvalueij()
////					* trial.getvalue() / SetUWInter.getvalue();
////		}
////		else{
////			iloop  = distance(fermionicloops.begin(), loopi_iter);
////			jloop  = distance(fermionicloops.begin(), loopj_iter);
////
////			if(line_j_iter->gettype()==1){
////				miflavor_n = Nflavor*unidist(mt);
////				mjflavor_n = (miflavor_n^1);
////			}
////			else{
////				proposal *= 2.;
////				miflavor_n = Nflavor*unidist(mt);
////				mjflavor_n = Nflavor*unidist(mt);
////			}
////
////			trial.setphysical(*line_j_iter,false);
////			trial.releaseloopflavor();
////			trial.fixloopflavor(iloop,miflavor_n);
////			trial.fixloopflavor(jloop,mjflavor_n);
////			trial.calvalue();
////
////			probability_first *= 0.5*(loop_valuerh_n + loop_valuelh_n) / loop_iter->getvaluerh() 
////				* trial.getvalue() / SetUWInter.getvalue();
////		}
////	}
////	else{
////		loop_iter = line_j_iter->getconnected_vi()->getconnected_loop();
////		loop_valuerh_n = loop_iter->getvaluerh() / line_j_iter->getvalueij();
////		loop_valuelh_n = loop_iter->getvaluelh() / line_j_iter->getvalueji();
////
////		miflavor_n = Nflavor*unidist(mt);
////		mjflavor_n = miflavor_n;
////
////		vi_iter = line_j_iter->getconnected_vi();
////		loopi_iter = vi_iter->getconnected_loop();
////		iloop  = distance(fermionicloops.begin(), loopi_iter);
////
////		if( measuringline->getconnected_vi()->getconnected_loop()
////			== measuringline->getconnected_vj()->getconnected_loop() ){
////			probability_first *= line_i_new.getvalueij();
////		}
////		else{
////			if(measuringline->gettype()==2){
////				proposal /= 2.;
////			}
////			trial.setphysical(*measuringline,true);
////		}
////
////		trial.releaseloopflavor();
////		trial.fixloopflavor(iloop,miflavor_n);
////		trial.calvalue();
////		probability_first *= loop_valuerh_n / 0.5/(loop_iter->getvaluerh()+loop_iter->getvaluelh())
////				* trial.getvalue() / SetUWInter.getvalue();
////	}
////	probability_first *= proposal;
////
////	cout<<"  probability_first = "<<probability_first<<endl;
////	cout<<"  proposal = "<<proposal<<endl;
////	cout<<"  trial.getvalue() / SetUWInter.getvalue() = "<<trial.getvalue() / SetUWInter.getvalue()<<endl;
////	cout<<"  loop_valuerh_n = "<<loop_valuerh_n<<endl;
////	cout<<"  loop_valuelh_n = "<<loop_valuelh_n<<endl;
////
////	int signfactor;
////	if(probability_first<0.) signfactor = -1;
////	else signfactor = 1;
////
////	list<line>::iterator line_j_prev;
////	int miflavor_o, mjflavor_o;
////	if(abs(probability_first)>1.0e-10){
////		if(!(checkirreducibilityforinterswap(line_j_iter))){ return true; }
////		if(!(checkcompactnessforinterswap(line_j_iter))){ return true; }
////
////		loop_iter->setvalue(loop_valuerh_n,loop_valuelh_n);
////
////		measuringline->setphysical(true);
////		line_j_iter->setphysical(false);
////		line_j_prev = measuringline;
////		measuringline = line_j_iter;
////
////		SetUWInter = trial;
////		miflavor_o = miflavor;
////		mjflavor_o = mjflavor;
////		miflavor = miflavor_n;
////		mjflavor = mjflavor_n;
////
////		Norder += dNorder;
////		sign *= signfactor;
////	}
////	else{ return true; }
////
////	cout<<"intermediate SetUWInter"<<endl<<SetUWInter;
////
////	////////////////////////////////////////////////
////	// update for returning to the original conf. //
////	////////////////////////////////////////////////
////
////	cout<<"*  second interswap  *"<<endl;
////
////	line_j_iter = line_j_prev;
////	if(measuringline->gettype()==0){
////		//index_j = static_cast<int>( Ws.size()*unidist(mt) );
////		//line_j_iter = next(Ws.begin(),index_j);
////		baredNorder = -1;
////	}
////	else{
////		//index_j = static_cast<int>( Gs.size()*unidist(mt) );
////		//line_j_iter = next(Gs.begin(),index_j);
////		baredNorder = +1;
////	}
////
////	cout<<"measuring line: "<<flush<<*measuringline<<endl;
////	cout<<"target line: "<<flush<<*line_j_iter<<endl;
////
////	dNorder = baredNorder;
////	if(Norder+dNorder>(this->Nordermax)){ 
////		cout<<"***** ALERT! DETAILED BALANCE CONDITION FAILS! *****"<<endl;
////		cout<<"fail to pass the Norder check when going back to original conf."<<endl;
////		return false; 
////	}
////
////	line_i_new = *measuringline;
////	line_i_new.setphysical(true);
////	cout<<"line_i_new : "<<line_i_new<<endl;
////
////	proposal = 1.0;
////	double probability_second = 
////		-1 // Fermionic sign factor
////		* Ri_N[Norder] / Ri_N[Norder+dNorder];
////		
////
////	ismeasuringlineGline = (measuringline->gettype()==0);
////	//trial.clear();
////	trial = SetUWInter;
////	if(ismeasuringlineGline){
////		loop_iter = measuringline->getconnected_vi()->getconnected_loop();
////		loop_valuerh_n = loop_iter->getvaluerh()*line_i_new.getvalueij();
////		loop_valuelh_n = loop_iter->getvaluelh()*line_i_new.getvalueji();
////
////		vi_iter = line_j_iter->getconnected_vi();
////		vj_iter = line_j_iter->getconnected_vj();
////		loopi_iter = vi_iter->getconnected_loop();
////		loopj_iter = vj_iter->getconnected_loop();
////		if(loopi_iter==loopj_iter){
////			iloop  = distance(fermionicloops.begin(), loopi_iter);
////
////			miflavor_n = miflavor_o;
////			mjflavor_n = mjflavor_o;
////
////			trial.releaseloopflavor();
////			trial.fixloopflavor(iloop,miflavor_n);
////			trial.calvalue();
////
////			probability_second *= 0.5*(loop_valuerh_n + loop_valuelh_n) / loop_iter->getvaluerh() 
////					/ line_j_iter->getvalueij()
////					* trial.getvalue() / SetUWInter.getvalue();
////		}
////		else{
////			iloop  = distance(fermionicloops.begin(), loopi_iter);
////			jloop  = distance(fermionicloops.begin(), loopj_iter);
////
////			if(line_j_iter->gettype()==1){
////				miflavor_n = miflavor_o;
////				mjflavor_n = mjflavor_o;
////			}
////			else{
////				proposal *= 2.;
////				miflavor_n = miflavor_o;
////				mjflavor_n = mjflavor_o;
////			}
////
////			trial.setphysical(*line_j_iter,false);
////			trial.releaseloopflavor();
////			trial.fixloopflavor(iloop,miflavor_n);
////			trial.fixloopflavor(jloop,mjflavor_n);
////			trial.calvalue();
////
////			probability_second *= 0.5*(loop_valuerh_n + loop_valuelh_n) / loop_iter->getvaluerh() 
////				* trial.getvalue() / SetUWInter.getvalue();
////		}
////	}
////	else{
////		loop_iter = line_j_iter->getconnected_vi()->getconnected_loop();
////		loop_valuerh_n = loop_iter->getvaluerh() / line_j_iter->getvalueij();
////		loop_valuelh_n = loop_iter->getvaluelh() / line_j_iter->getvalueji();
////
////		miflavor_n = miflavor_o;
////		mjflavor_n = mjflavor_o;
////
////		vi_iter = line_j_iter->getconnected_vi();
////		loopi_iter = vi_iter->getconnected_loop();
////		iloop  = distance(fermionicloops.begin(), loopi_iter);
////
////		if( measuringline->getconnected_vi()->getconnected_loop()
////			== measuringline->getconnected_vj()->getconnected_loop() ){
////			probability_second *= line_i_new.getvalueij();
////		}
////		else{
////			if(measuringline->gettype()==2){
////				proposal /= 2.;
////			}
////			trial.setphysical(*measuringline,true);
////		}
////
////		trial.releaseloopflavor();
////		trial.fixloopflavor(iloop,miflavor_n);
////		trial.calvalue();
////		probability_second *= loop_valuerh_n / 0.5/(loop_iter->getvaluerh()+loop_iter->getvaluelh())
////				* trial.getvalue() / SetUWInter.getvalue();
////	}
////	probability_second *= proposal;
////	trial = SetUWInter;
////	if(measuringline->gettype()==0){
////		vi_iter = line_j_iter->getconnected_vi();
////		vj_iter = line_j_iter->getconnected_vj();
////		loopi_iter = vi_iter->getconnected_loop();
////		loopj_iter = vj_iter->getconnected_loop();
////		if(loopi_iter==loopj_iter){
////			iloop  = distance(fermionicloops.begin(), loopi_iter);
////
////			miflavor_n = miflavor_o;
////			mjflavor_n = mjflavor_o;
////
////			trial.releaseloopflavor();
////			trial.fixloopflavor(iloop,miflavor_n);
////
////			probability_second *= line_i_new.getvalue() / line_j_iter->getvalue();
////		}
////		else{
////			iloop  = distance(fermionicloops.begin(), loopi_iter);
////			jloop  = distance(fermionicloops.begin(), loopj_iter);
////
////			proposal *= 2.;
////			miflavor_n = miflavor_o;
////			mjflavor_n = mjflavor_o;
////
////			trial.setphysical(*line_j_iter,false);
////			trial.releaseloopflavor();
////			trial.fixloopflavor(iloop,miflavor_n);
////			trial.fixloopflavor(jloop,mjflavor_n);
////			trial.calvalue();
////
////			probability_second *= line_i_new.getvalue() * trial.getvalue() / SetUWInter.getvalue();
////		}
////	}
////	else{
////		proposal /= 2.;
////		miflavor_n = miflavor_o;
////		mjflavor_n = mjflavor_o;
////
////		vi_iter = line_j_iter->getconnected_vi();
////		loopi_iter = vi_iter->getconnected_loop();
////		iloop  = distance(fermionicloops.begin(), loopi_iter);
////
////		if( measuringline->getconnected_vi()->getconnected_loop() 
////			== measuringline->getconnected_vj()->getconnected_loop() )
////			probability_second *= line_i_new.getvalue();
////		else
////			trial.setphysical(*measuringline,true);
////
////		trial.releaseloopflavor();
////		trial.fixloopflavor(iloop,miflavor_n);
////		trial.calvalue();
////		probability_second *= 1./line_j_iter->getvalue() * trial.getvalue() / SetUWInter.getvalue();
////	}
////
////	cout<<"  probability_second = "<<probability_second<<endl;
////	cout<<"  proposal = "<<proposal<<endl;
////	cout<<"  trial.getvalue() / SetUWInter.getvalue() = "<<trial.getvalue() / SetUWInter.getvalue()<<endl;
////	cout<<"  loop_valuerh_n = "<<loop_valuerh_n<<endl;
////	cout<<"  loop_valuelh_n = "<<loop_valuelh_n<<endl;
////
////
////	if(probability_second<0.) signfactor = -1;
////	else signfactor = 1;
////
////	{
////		if(!(checkirreducibilityforinterswap(line_j_iter))){ 
////			cout<<"***** ALERT! DETAILED BALANCE CONDITION FAILS! *****"<<endl;
////			cout<<"fail to pass the irreducibility check when going back to original conf."<<endl;
////			return false; 
////		}
////		if(!(checkcompactnessforinterswap(line_j_iter))){ 
////			cout<<"***** ALERT! DETAILED BALANCE CONDITION FAILS! *****"<<endl;
////			cout<<"fail to pass the compactness check when going back to original conf."<<endl;
////			return false; 
////		}
////
////		loop_iter->setvalue(loop_valuerh_n,loop_valuelh_n);
////
////		measuringline->setphysical(true);
////		line_j_iter->setphysical(false);
////		measuringline = line_j_iter;
////
////		SetUWInter = trial;
////		miflavor = miflavor_n;
////		mjflavor = mjflavor_n;
////
////		Norder += dNorder;
////		sign *= signfactor;
////	}
////
////	cout<<"returning SetUWInter"<<endl<<SetUWInter;
////
////	if( abs(probability_first*probability_second - 1.)<1.0e-10 ){
////		cout<<"probability_first = "<<probability_first<<endl;
////		cout<<"probability_second = "<<probability_second<<endl;
////		cout<<"probability_first*probability_second = "<<probability_first*probability_second<<endl;
////		return true;
////	}
////	else{
////		cout<<"***** ALERT! DETAILED BALANCE CONDITION FAILS! *****"<<endl;
////		cout<<"probability_first = "<<probability_first<<endl;
////		cout<<"probability_second = "<<probability_second<<endl;
////		cout<<"probability_first*probability_second = "<<probability_first*probability_second<<endl;
////		printconfiguration(cout);
////		return false;
////	}
////};
////bool diagram::checkdetailedbalance_intraswap(){
////	cout<<"** intraswap detailed balance check **"<<endl;
////
////	/////////////////////
////	// first interswap //
////	/////////////////////
////
////	cout<<"*  first intraswap  *"<<endl;
////	int index_i, index_j;
////	list<line>::iterator line_j_iter;
////	bool ismeasuringlineGline = (measuringline->gettype()==0);
////	if(ismeasuringlineGline){
////		index_i = distance(Gs.begin(),measuringline);
////		select_one_excluding_one(Gs.size(),index_i,index_j);
////		line_j_iter = next(Gs.begin(),index_j);
////	}
////	else{
////		index_i = distance(Ws.begin(),measuringline);
////		select_one_excluding_one(Ws.size(),index_i,index_j);
////		line_j_iter = next(Gs.begin(),index_j);
////	}
////
////	line line_i_new(*measuringline);
////	line_i_new.setphysical(true);
////
////	double proposal = 1.0;
////	double probability_first = 1.0;
////	
////	int miflavor_n, mjflavor_n;
////	int iloop, jloop;
////	list<fermionicloop>::iterator loopi_iter, loopj_iter;
////
////	cout<<"measuring line: "<<*measuringline<<endl;
////	cout<<"target line: "<<*line_j_iter<<endl;
////	cout<<"line_i_new: "<<line_i_new<<endl;
////
////	cout<<"original SetUWInter"<<endl<<SetUWInter;
////
////	double loopi_valuerh_n, loopi_valuelh_n;
////	double loopj_valuerh_n, loopj_valuelh_n;
////
////	vector<bool> isloopflavorfixed_n(fermionicloops.size(),false);
////	vector<int> loopflavor_n(fermionicloops.size(),0);
////	double valtrial;
////
////	bool islinei_inter, islinej_inter;
////	if(ismeasuringlineGline){
////		loopi_iter = measuringline->getconnected_vi()->getconnected_loop();
////		loopj_iter = line_j_iter->getconnected_vi()->getconnected_loop();
////		if(loopi_iter==loopj_iter){
////			loopi_valuerh_n = loopi_iter->getvaluerh() 
////					* line_i_new.getvalueij() 
////					/ line_j_iter->getvalueij();
////			loopi_valuelh_n = loopi_iter->getvaluelh() 
////					* line_i_new.getvalueji() 
////					/ line_j_iter->getvalueji();
////
////			miflavor_n = Nflavor*unidist(mt);
////			mjflavor_n = miflavor_n;
////
////			jloop = distance(fermionicloops.begin(),loopj_iter);
////			isloopflavorfixed_n[jloop] = true;
////			loopflavor_n[jloop] = miflavor_n;
////			valtrial = SetUWInter.getvalue();
////
////			probability_first *= loopi_valuerh_n / loopi_iter->getvaluerh();
////		}
////		else{
////			loopi_valuerh_n = loopi_iter->getvaluerh() 
////					* line_i_new.getvalueij();
////			loopi_valuelh_n = loopi_iter->getvaluelh() 
////					* line_i_new.getvalueji();
////			loopj_valuerh_n = loopj_iter->getvaluerh() 
////					/ line_j_iter->getvalueij();
////			loopj_valuelh_n = loopj_iter->getvaluelh() 
////					/ line_j_iter->getvalueji();
////
////			miflavor_n = Nflavor*unidist(mt);
////			mjflavor_n = miflavor_n;
////
////			jloop = distance(fermionicloops.begin(),loopj_iter);
////			isloopflavorfixed_n[jloop] = true;
////			loopflavor_n[jloop] = miflavor_n;
////			valtrial = SetUWInter.calvalueafterinterswap(isloopflavorfixed_n,loopflavor_n);
////
////			probability_first *= 0.5*(loopi_valuerh_n + loopi_valuelh_n) / loopi_iter->getvaluerh() 
////					* loopj_valuerh_n / 0.5/(loopj_iter->getvaluerh()+loopj_iter->getvaluelh())
////					* valtrial / SetUWInter.getvalue();
////		}
////	}
////	else{
////		islinei_inter = (measuringline->getconnected_vi()->getconnected_loop()
////					!= measuringline->getconnected_vj()->getconnected_loop());
////		islinej_inter = (line_j_iter->getconnected_vi()->getconnected_loop()
////					!= line_j_iter->getconnected_vj()->getconnected_loop());
////
////		if(!islinej_inter){
////			miflavor_n = Nflavor*unidist(mt);
////			mjflavor_n = miflavor_n;
////
////			if(!islinei_inter){
////				probability_first *= line_i_new.getvalueij() / line_j_iter->getvalueij();
////				valtrial = SetUWInter.getvalue();
////			}
////			else{
////				proposal /= 2.;
////
////				loopi_iter = line_j_iter->getconnected_vi()->getconnected_loop();
////
////				iloop = distance(fermionicloops.begin(),loopi_iter);
////				isloopflavorfixed_n[iloop] = true;
////				loopflavor_n[iloop] = miflavor_n;
////				valtrial = SetUWInter.calvalueafterinterswap(*measuringline,isloopflavorfixed_n,loopflavor_n);
////
////				probability_first *= 1. / line_j_iter->getvalueij()
////						* valtrial / SetUWInter.getvalue();
////			}
////		}
////		else{
////			miflavor_n = Nflavor*unidist(mt);
////			mjflavor_n = Nflavor*unidist(mt);
////
////			loopi_iter = line_j_iter->getconnected_vi()->getconnected_loop();
////			loopj_iter = line_j_iter->getconnected_vj()->getconnected_loop();
////
////			iloop = distance(fermionicloops.begin(),loopi_iter);
////			jloop = distance(fermionicloops.begin(),loopj_iter);
////			isloopflavorfixed_n[iloop] = true;
////			isloopflavorfixed_n[jloop] = true;
////			loopflavor_n[iloop] = miflavor_n;
////			loopflavor_n[jloop] = mjflavor_n;
////
////			if(!islinei_inter){
////				proposal *= 2.;
////				valtrial = SetUWInter.calvalueafterinterswap(*line_j_iter,isloopflavorfixed_n,loopflavor_n);
////			}
////			else{
////				valtrial = SetUWInter.calvalueafterintraswap(*measuringline,*line_j_iter,isloopflavorfixed_n,loopflavor_n);
////			}
////		}
////	}
////
////	probability_first *= proposal;
////
////	int signfactor;
////	if(probability_first<0.) signfactor = -1;
////	else signfactor = 1;
////
////	//if(NdressChoice==3 && measuringline->getisG() && Ndress==2){
////	//printconfiguration(cout);
////	//cout<<"* interswap info *"<<endl;
////	//cout<<"measuringline = "<<*measuringline<<endl;
////	//cout<<"line_i_new = "<<line_i_new<<endl;
////	//cout<<"line_j_iter = "<<*line_j_iter<<endl;
////	//cout<<"Norder = "<<Norder<<endl;
////	//cout<<"dNorder = "<<dNorder<<endl<<flush;
////	//cout<<"probability = "<<probability<<endl;
////	//cout<<"signfactor = "<<signfactor<<endl;
////	//cout<<"line_i_new.getvalue() = "<<line_i_new.getvalue()<<endl;
////	//cout<<"line_j_iter->getvalue() = "<<line_j_iter->getvalue()<<endl;
////	//cout<<endl<<endl<<flush;
////	//cout<<"type of measuringline = "<<measuringline->gettype()<<endl;
////	//cout<<"(index_i,index_j) = ("<<index_i<<","<<index_j<<")"<<endl;
////	//cout<<"NdressChoice = "<<NdressChoice<<endl;
////	//cout<<"NdressChoiceInverse = "<<(1 + Nordermax - Norder + line_j_iter->getdress())<<endl;
////	//cout<<"Ndress = "<<Ndress<<endl<<flush;
////	//}
////
////	list<line>::iterator line_j_prev;
////	int miflavor_o, mjflavor_o;
////	if(abs(probability_first)>1.0e-10){
////		//if(!(checkconnectivityforinterswap(line_j_iter))){ return 0; }
////		if(!(checkirreducibilityforinterswap(line_j_iter))){ return true; }
////		if(!(checkcompactnessforinterswap(line_j_iter))){ return true; }
////
////		SetUWInter.setisloopflavorfixed(isloopflavorfixed_n);
////		SetUWInter.setloopflavor(loopflavor_n);
////		if( ismeasuringlineGline ){
////			if(loopi_iter==loopj_iter)
////				loopi_iter->setvalue(loopi_valuerh_n,loopi_valuelh_n);
////			else{
////				loopi_iter->setvalue(loopi_valuerh_n,loopi_valuelh_n);
////				loopj_iter->setvalue(loopj_valuerh_n,loopj_valuelh_n);
////			}
////		}
////		else{
////			if(islinei_inter)
////				SetUWInter.setphysical(*measuringline,true);
////			if(islinej_inter)
////				SetUWInter.setphysical(*line_j_iter,false);
////		}
////		SetUWInter.setvalue(valtrial);
////
////		measuringline->setphysical(true);
////		line_j_iter->setphysical(false);
////		line_j_prev = measuringline;
////		measuringline = line_j_iter;
////
////		miflavor_o = miflavor;
////		mjflavor_o = mjflavor;
////		miflavor = miflavor_n;
////		mjflavor = mjflavor_n;
////
////		sign *= signfactor;
////
////	}
////	else{ return true; }
////
////	cout<<"intermediate SetUWInter"<<endl<<SetUWInter;
////
////	////////////////////////////////////////////////
////	// update for returning to the original conf. //
////	////////////////////////////////////////////////
////
////	cout<<"*  second interswap  *"<<endl;
////
////	line_j_iter = line_j_prev;
////
////	cout<<"measuring line: "<<flush<<*measuringline<<endl;
////	cout<<"target line: "<<flush<<*line_j_iter<<endl;
////
////	line_i_new = *measuringline;
////	line_i_new.setphysical(true);
////
////	proposal = 1.0;
////	double probability_second = 1.0;
////	
////	for(int i=0; i<fermionicloops.size(); i++){
////		isloopflavorfixed_n[i] = false;
////		loopflavor_n[i] = 0;
////	}
////
////	if(ismeasuringlineGline){
////		loopi_iter = measuringline->getconnected_vi()->getconnected_loop();
////		loopj_iter = line_j_iter->getconnected_vi()->getconnected_loop();
////		if(loopi_iter==loopj_iter){
////			loopi_valuerh_n = loopi_iter->getvaluerh() 
////					* line_i_new.getvalueij() 
////					/ line_j_iter->getvalueij();
////			loopi_valuelh_n = loopi_iter->getvaluelh() 
////					* line_i_new.getvalueji() 
////					/ line_j_iter->getvalueji();
////
////			miflavor_n = miflavor_o;
////			mjflavor_n = mjflavor_o;
////
////			jloop = distance(fermionicloops.begin(),loopj_iter);
////			isloopflavorfixed_n[jloop] = true;
////			loopflavor_n[jloop] = miflavor_n;
////			valtrial = SetUWInter.getvalue();
////
////			probability_second *= loopi_valuerh_n / loopi_iter->getvaluerh();
////		}
////		else{
////			loopi_valuerh_n = loopi_iter->getvaluerh() 
////					* line_i_new.getvalueij();
////			loopi_valuelh_n = loopi_iter->getvaluelh() 
////					* line_i_new.getvalueji();
////			loopj_valuerh_n = loopj_iter->getvaluerh() 
////					/ line_j_iter->getvalueij();
////			loopj_valuelh_n = loopj_iter->getvaluelh() 
////					/ line_j_iter->getvalueji();
////
////			miflavor_n = miflavor_o;
////			mjflavor_n = mjflavor_o;
////
////			jloop = distance(fermionicloops.begin(),loopj_iter);
////			isloopflavorfixed_n[jloop] = true;
////			loopflavor_n[jloop] = miflavor_n;
////			valtrial = SetUWInter.calvalueafterinterswap(isloopflavorfixed_n,loopflavor_n);
////
////			probability_second *= 0.5*(loopi_valuerh_n + loopi_valuelh_n) / loopi_iter->getvaluerh() 
////					* loopj_valuerh_n / 0.5/(loopj_iter->getvaluerh()+loopj_iter->getvaluelh())
////					* valtrial / SetUWInter.getvalue();
////		}
////	}
////	else{
////		islinei_inter = (measuringline->getconnected_vi()->getconnected_loop()
////					!= measuringline->getconnected_vj()->getconnected_loop());
////		islinej_inter = (line_j_iter->getconnected_vi()->getconnected_loop()
////					!= line_j_iter->getconnected_vj()->getconnected_loop());
////
////		if(!islinej_inter){
////			miflavor_n = miflavor_o;
////			mjflavor_n = mjflavor_o;
////
////			if(!islinei_inter){
////				probability_second *= line_i_new.getvalueij() / line_j_iter->getvalueij();
////				valtrial = SetUWInter.getvalue();
////			}
////			else{
////				proposal /= 2.;
////
////				loopi_iter = line_j_iter->getconnected_vi()->getconnected_loop();
////
////				iloop = distance(fermionicloops.begin(),loopi_iter);
////				isloopflavorfixed_n[iloop] = true;
////				loopflavor_n[iloop] = miflavor_n;
////				valtrial = SetUWInter.calvalueafterinterswap(*measuringline,isloopflavorfixed_n,loopflavor_n);
////
////				probability_second *= 1. / line_j_iter->getvalueij()
////						* valtrial / SetUWInter.getvalue();
////			}
////		}
////		else{
////			miflavor_n = miflavor_o;
////			mjflavor_n = mjflavor_o;
////
////			loopi_iter = line_j_iter->getconnected_vi()->getconnected_loop();
////			loopj_iter = line_j_iter->getconnected_vj()->getconnected_loop();
////
////			iloop = distance(fermionicloops.begin(),loopi_iter);
////			jloop = distance(fermionicloops.begin(),loopj_iter);
////			isloopflavorfixed_n[iloop] = true;
////			isloopflavorfixed_n[jloop] = true;
////			loopflavor_n[iloop] = miflavor_n;
////			loopflavor_n[jloop] = mjflavor_n;
////
////			if(!islinei_inter){
////				proposal *= 2.;
////				valtrial = SetUWInter.calvalueafterinterswap(*line_j_iter,isloopflavorfixed_n,loopflavor_n);
////			}
////			else{
////				valtrial = SetUWInter.calvalueafterintraswap(*measuringline,*line_j_iter,isloopflavorfixed_n,loopflavor_n);
////			}
////		}
////	}
////
////	probability_second *= proposal;
////
////	if(probability_second<0.) signfactor = -1;
////	else signfactor = 1;
////
////	//if(NdressChoice==3 && measuringline->getisG() && Ndress==2){
////	//printconfiguration(cout);
////	//cout<<"* interswap info *"<<endl;
////	//cout<<"measuringline = "<<*measuringline<<endl;
////	//cout<<"line_i_new = "<<line_i_new<<endl;
////	//cout<<"line_j_iter = "<<*line_j_iter<<endl;
////	//cout<<"Norder = "<<Norder<<endl;
////	//cout<<"dNorder = "<<dNorder<<endl<<flush;
////	//cout<<"probability = "<<probability<<endl;
////	//cout<<"signfactor = "<<signfactor<<endl;
////	//cout<<"line_i_new.getvalue() = "<<line_i_new.getvalue()<<endl;
////	//cout<<"line_j_iter->getvalue() = "<<line_j_iter->getvalue()<<endl;
////	//cout<<endl<<endl<<flush;
////	//cout<<"type of measuringline = "<<measuringline->gettype()<<endl;
////	//cout<<"(index_i,index_j) = ("<<index_i<<","<<index_j<<")"<<endl;
////	//cout<<"NdressChoice = "<<NdressChoice<<endl;
////	//cout<<"NdressChoiceInverse = "<<(1 + Nordermax - Norder + line_j_iter->getdress())<<endl;
////	//cout<<"Ndress = "<<Ndress<<endl<<flush;
////	//}
////
////	if(abs(probability_second)>1.0e-10){
////		//if(!(checkconnectivityforinterswap(line_j_iter))){ return 0; }
////		if(!(checkirreducibilityforinterswap(line_j_iter))){ return false; }
////		if(!(checkcompactnessforinterswap(line_j_iter))){ return false; }
////
////		SetUWInter.setisloopflavorfixed(isloopflavorfixed_n);
////		SetUWInter.setloopflavor(loopflavor_n);
////		if( ismeasuringlineGline ){
////			if(loopi_iter==loopj_iter)
////				loopi_iter->setvalue(loopi_valuerh_n,loopi_valuelh_n);
////			else{
////				loopi_iter->setvalue(loopi_valuerh_n,loopi_valuelh_n);
////				loopj_iter->setvalue(loopj_valuerh_n,loopj_valuelh_n);
////			}
////		}
////		else{
////			if(islinei_inter)
////				SetUWInter.setphysical(*measuringline,true);
////			if(islinej_inter)
////				SetUWInter.setphysical(*line_j_iter,false);
////		}
////		SetUWInter.setvalue(valtrial);
////
////		measuringline->setphysical(true);
////		line_j_iter->setphysical(false);
////		measuringline = line_j_iter;
////
////		miflavor = miflavor_n;
////		mjflavor = mjflavor_n;
////
////		sign *= signfactor;
////
////	}
////	else{ return false; }
////
////
////	cout<<"returning SetUWInter"<<endl<<SetUWInter;
////
////	if( abs(probability_first*probability_second - 1.)<1.0e-10 ){
////		cout<<"probability_first = "<<probability_first<<endl;
////		cout<<"probability_second = "<<probability_second<<endl;
////		cout<<"probability_first*probability_second = "<<probability_first*probability_second<<endl;
////		return true;
////	}
////	else{
////		cout<<"***** ALERT! DETAILED BALANCE CONDITION FAILS! *****"<<endl;
////		cout<<"probability_first = "<<probability_first<<endl;
////		cout<<"probability_second = "<<probability_second<<endl;
////		cout<<"probability_first*probability_second = "<<probability_first*probability_second<<endl;
////		printconfiguration(cout);
////		return false;
////	}
////};*/
////
////
/////*bool diagram::checkdetailedbalance_interswap(){
////	//cout<<"* interswap module *"<<endl<<flush;
////	int index_j;
////	int baredNorder;
////	list<line>::iterator line_j_iter;
////	if(measuringline->gettype()==0){
////		index_j = static_cast<int>( (Us.size()+Ws.size())*unidist(mt) );
////		if(index_j<Us.size())
////			line_j_iter = next(Us.begin(),index_j);
////		else
////			line_j_iter = next(Ws.begin(),index_j-Us.size());
////		baredNorder = -1;
////		}
////	else{
////		index_j = static_cast<int>( Gs.size()*unidist(mt) );
////		line_j_iter = next(Gs.begin(),index_j);
////		baredNorder = +1;
////	}
////
////	int dNorder = baredNorder;
////	if(Norder+dNorder>(this->Nordermax)){ return true; }
////	if(!(checkconnectivityforinterswap(line_j_iter))){ return true; }
////
////	if(Norder==1 && Nselfloop==2){
////		cout<<"there is problem"<<endl;
////		exit(EXIT_FAILURE);
////	}
////
////	line line_i_new(*measuringline);
////	line_i_new.setdress(measuringline->getdress());
////	line_i_new.setphysical(true);
////
////	double proposal = 1.0;
////
////	// Wow if I turn off below two lines, it works! Super strange... Orz...
////	if(measuringline->gettype()==0) proposal *= 0.5;
////	else proposal *= 2.0;
////
////	double probability_first = 
////		-1 // Fermionic sign factor
////		* proposal
////		* Ri_N[Norder] / Ri_N[Norder+dNorder]
////		* line_i_new.getvalue() / line_j_iter->getvalue();
////
////	int signfactor;
////	if(probability_first<0.) signfactor = -1;
////	else signfactor = 1;
////
////	//if(NdressChoice==3 && measuringline->getisG() && Ndress==2){
////	//printconfiguration(cout);
////	cout<<"* interswap info *"<<endl;
////	cout<<"measuringline = "<<*measuringline<<endl;
////	cout<<"line_i_new = "<<line_i_new<<endl;
////	cout<<"line_j_iter = "<<*line_j_iter<<endl;
////	cout<<"Norder = "<<Norder<<endl;
////	cout<<"dNorder = "<<dNorder<<endl<<flush;
////	cout<<"proposal = "<<proposal<<endl;
////	cout<<"probability_first = "<<probability_first<<endl;
////	cout<<"signfactor = "<<signfactor<<endl;
////	cout<<"line_i_new.getvalue() = "<<line_i_new.getvalue()<<endl;
////	cout<<"line_j_iter->getvalue() = "<<line_j_iter->getvalue()<<endl;
////	cout<<endl<<endl<<flush;
////	//cout<<"type of measuringline = "<<measuringline->gettype()<<endl;
////	//cout<<"(index_i,index_j) = ("<<index_i<<","<<index_j<<")"<<endl;
////	//cout<<"NdressChoice = "<<NdressChoice<<endl;
////	//cout<<"NdressChoiceInverse = "<<(1 + Nordermax - Norder + line_j_iter->getdress())<<endl;
////	//cout<<"Ndress = "<<Ndress<<endl<<flush;
////	//}
////
////	list<line>::iterator line_j_prev;
////	{
////		measuringline->setphysical(true);
////		line_j_iter->setphysical(false);
////		line_j_prev = measuringline;
////		measuringline = line_j_iter;
////
////		Norder += dNorder;
////		sign *= signfactor;
////	}
////
////	//cout<<"* interswap module *"<<endl<<flush;
////	line_j_iter = line_j_prev;
////	if(measuringline->gettype()==0){
////		//index_j = static_cast<int>( (Us.size()+Ws.size())*unidist(mt) );
////		//if(index_j<Us.size())
////		//	line_j_iter = next(Us.begin(),index_j);
////		//else
////		//	line_j_iter = next(Ws.begin(),index_j-Us.size());
////		baredNorder = -1;
////		}
////	else{
////		//index_j = static_cast<int>( Gs.size()*unidist(mt) );
////		//line_j_iter = next(Gs.begin(),index_j);
////		baredNorder = +1;
////	}
////
////	dNorder = baredNorder;
////	if(Norder+dNorder>(this->Nordermax)){ return false; }
////	if(!(checkconnectivityforinterswap(line_j_iter))){ return false; }
////
////	if(Norder==1 && Nselfloop==2){
////		cout<<"there is problem"<<endl;
////		exit(EXIT_FAILURE);
////	}
////
////	line_i_new = (*measuringline);
////	line_i_new.setdress(measuringline->getdress());
////	line_i_new.setphysical(true);
////
////	proposal = 1.0;
////
////	// Wow if I turn off below two lines, it works! Super strange... Orz...
////	if(measuringline->gettype()==0) proposal *= 0.5;
////	else proposal *= 2.0;
////
////	double probability_second = 
////		-1 // Fermionic sign factor
////		* proposal
////		* Ri_N[Norder] / Ri_N[Norder+dNorder]
////		* line_i_new.getvalue() / line_j_iter->getvalue();
////
////	if(probability_second<0.) signfactor = -1;
////	else signfactor = 1;
////
////	//if(NdressChoice==3 && measuringline->getisG() && Ndress==2){
////	//printconfiguration(cout);
////	cout<<"* interswap info *"<<endl;
////	cout<<"measuringline = "<<*measuringline<<endl;
////	cout<<"line_i_new = "<<line_i_new<<endl;
////	cout<<"line_j_iter = "<<*line_j_iter<<endl;
////	cout<<"Norder = "<<Norder<<endl;
////	cout<<"dNorder = "<<dNorder<<endl<<flush;
////	cout<<"proposal = "<<proposal<<endl;
////	cout<<"probability_second = "<<probability_second<<endl;
////	cout<<"signfactor = "<<signfactor<<endl;
////	cout<<"line_i_new.getvalue() = "<<line_i_new.getvalue()<<endl;
////	cout<<"line_j_iter->getvalue() = "<<line_j_iter->getvalue()<<endl;
////	cout<<endl<<endl<<flush;
////	//cout<<"type of measuringline = "<<measuringline->gettype()<<endl;
////	//cout<<"(index_i,index_j) = ("<<index_i<<","<<index_j<<")"<<endl;
////	//cout<<"NdressChoice = "<<NdressChoice<<endl;
////	//cout<<"NdressChoiceInverse = "<<(1 + Nordermax - Norder + line_j_iter->getdress())<<endl;
////	//cout<<"Ndress = "<<Ndress<<endl<<flush;
////	//}
////
////	{
////		measuringline->setphysical(true);
////		line_j_iter->setphysical(false);
////		measuringline = line_j_iter;
////
////		Norder += dNorder;
////		sign *= signfactor;
////	}
////
////	if( abs(probability_first*probability_second - 1.)<1.0e-10 ){
////		cout<<"probability_first = "<<probability_first<<endl;
////		cout<<"probability_second = "<<probability_second<<endl;
////		cout<<"probability_first*probability_second = "<<probability_first*probability_second<<endl;
////		return true;
////	}
////	else{
////		cout<<"***** ALERT! DETAILED BALANCE CONDITION FAILS! *****"<<endl;
////		cout<<"probability_first = "<<probability_first<<endl;
////		cout<<"probability_second = "<<probability_second<<endl;
////		cout<<"probability_first*probability_second = "<<probability_first*probability_second<<endl;
////		printconfiguration(cout);
////		return false;
////	}
////
////};*/
////
/////*bool diagram::checkdetailedbalance_transformUW(){
////	cout<<"** checkdetailedbalance_transformUW module **"<<endl;
////	list<list<line>::iterator> UWinters;
////	list<line>::iterator UW_iter;
////	list<fermionicloop>::iterator loopi_iter, loopj_iter;
////	for(UW_iter=Ws.begin(); UW_iter!=Ws.end(); ++UW_iter){
////		loopi_iter = UW_iter->getconnected_vi()->getconnected_loop();
////		loopj_iter = UW_iter->getconnected_vj()->getconnected_loop();
////		if( loopi_iter!=loopj_iter && UW_iter->getphysical() ) UWinters.push_back(UW_iter);
////	}
////	if(UWinters.size()==0) return true;
////	else{
////		/////////////////////
////		// first transform //
////		/////////////////////
////
////		int NUWinters = UWinters.size();
////		int iUWinters = static_cast<int>(NUWinters*unidist(mt));
////		int ijselection = static_cast<int>(2*unidist(mt));
////
////		UW_iter = (*next(UWinters.begin(),iUWinters));
////		int UWtype = UW_iter->gettype();
////		list<vertex>::iterator v_iter, v_oppo;
////		if(ijselection==0){
////			v_iter = UW_iter->getconnected_vi();
////			v_oppo = UW_iter->getconnected_vj();
////		}
////		else{
////			v_iter = UW_iter->getconnected_vj();
////			v_oppo = UW_iter->getconnected_vi();
////		}
////
////		double proposal;
////		double probability_first;
////
////		coordinate x_n, x_o;
////		line UW_n = *UW_iter;
////		if(UWtype==1){
////			x_n.genrand(v_iter->getx());
////
////			proposal = beta/(v_iter->getx().getsitegenprobability(x_n));
////			if(ijselection==0) 
////				UW_n = line(UW_iter->getphysical(),2,UW_iter->getdress(),UW_iter->getiflavor(),UW_iter->getjflavor(),x_n,UW_iter->getxj());
////			else 
////				UW_n = line(UW_iter->getphysical(),2,UW_iter->getdress(),UW_iter->getiflavor(),UW_iter->getjflavor(),UW_iter->getxi(),x_n);
////		}
////		else{
////			x_n = v_oppo->getx();
////
////			proposal = (v_iter->getx().getsitegenprobability(x_n))/beta;
////			if(ijselection==0)
////				UW_n = line(UW_iter->getphysical(),1,UW_iter->getdress(),UW_iter->getiflavor(),UW_iter->getjflavor(),x_n,UW_iter->getxj());
////			else 
////				UW_n = line(UW_iter->getphysical(),1,UW_iter->getdress(),UW_iter->getiflavor(),UW_iter->getjflavor(),UW_iter->getxi(),x_n);
////		}
////
////		list<fermionicloop>::iterator loop_iter = v_iter->getconnected_loop();
////		bool ismeasuringlineinloop = ( measuringline->getconnected_vi()->getconnected_loop()==loop_iter );
////
////		setUWinter trial = SetUWInter;
////		trial.transform(*UW_iter,UW_n);
////		trial.calvalue();
////
////		list<line>::iterator Gin_iter = v_iter->getconnected_Gin();
////		list<line>::iterator Gout_iter = v_iter->getconnected_Gout();
////		bool is_selfloop = (Gin_iter==Gout_iter);
////		line Gin_n = *Gin_iter;
////		line Gout_n = *Gout_iter;
////
////		double loop_valuerh_n, loop_valuelh_n;
////
////		if(is_selfloop){
////			Gin_n.setx(x_n,x_n);
////			probability_first = proposal 
////				* trial.getvalue() / SetUWInter.getvalue();
////		}
////		else{
////			Gin_n.setxj(x_n);
////			Gout_n.setxi(x_n);
////			loop_valuerh_n = Gin_n.getvalueij() * Gout_n.getvalueij() 
////					/ Gin_iter->getvalueij() / Gout_iter->getvalueij() 
////					* loop_iter->getvaluerh();
////			loop_valuelh_n = Gin_n.getvalueji() * Gout_n.getvalueji() 
////					/ Gin_iter->getvalueji() / Gout_iter->getvalueji() 
////					* loop_iter->getvaluelh();
////
////			if(ismeasuringlineinloop){
////				probability_first = proposal
////						* loop_valuerh_n / loop_iter->getvaluerh()
////						* trial.getvalue() / SetUWInter.getvalue();
////			}
////			else{
////				probability_first = proposal
////						* (loop_valuerh_n+loop_valuelh_n) 
////						/ (loop_iter->getvaluerh()+loop_iter->getvaluelh())
////						* trial.getvalue() / SetUWInter.getvalue();
////			}
////		}
////
////		if(is_selfloop){
////			Gin_n.setx(x_n,x_n);
////			probability_first = proposal
////				* Gin_n.getvalue() * trial.getvalue()
////				/ Gin_iter->getvalue() / SetUWInter.getvalue();
////		}
////		else{
////			Gin_n.setxj(x_n);
////			Gout_n.setxi(x_n);
////			probability_first = proposal
////				* Gin_n.getvalue() * Gout_n.getvalue() * trial.getvalue()
////				/ Gin_iter->getvalue() / Gout_iter->getvalue() / SetUWInter.getvalue();
////		}
////
////		//cout<<"* first update "<<endl;
////		//cout<<"proposal = "<<proposal<<endl;
////		//cout<<"Gin_n.getvalue()      	= "<<Gin_n.getvalue() 		<<endl;
////		//cout<<"Gout_n.getvalue()     	= "<<Gout_n.getvalue() 		<<endl;
////		//cout<<"trial.getvalue()      	= "<<trial.getvalue()		<<endl;
////		//cout<<"Gin_iter->getvalue()  	= "<<Gin_iter->getvalue() 	<<endl;
////		//cout<<"Gout_iter->getvalue()    = "<<Gout_iter->getvalue()	<<endl;
////		//cout<<"SetUWInter.getvalue()    = "<<SetUWInter.getvalue()	<<endl;
////		//cout<<"trial info."<<endl<<trial<<endl;
////		//cout<<"SetUWInter info."<<endl<<SetUWInter<<endl;
////
////		int signfactor;
////		if(probability_first<0.) signfactor = -1;
////		else signfactor = 1;
////
////		x_o = v_iter->getx();
////		if(abs(probability_first)>1.0e-10){
////			v_iter->setx(x_n);
////			UW_iter->assign(UW_n);
////			if(is_selfloop) 
////				Gin_iter->assign(Gin_n);
////			else{
////				Gin_iter->assign(Gin_n);
////				Gout_iter->assign(Gout_n);
////				loop_iter->setvalue(loop_valuerh_n,loop_valuelh_n);
////			}
////
////			SetUWInter = trial;
////
////			sign *= signfactor;
////		}
////		else{
////			return true;
////		}
////
////
////		//////////////////////
////		// second transform //
////		//////////////////////
////
////		double probability_second;
////
////		UWtype = UW_iter->gettype();
////		UW_n = *UW_iter;
////		if(UWtype==1){
////			x_n = x_o;
////
////			proposal = beta/(v_iter->getx().getsitegenprobability(x_n));
////			if(ijselection==0) 
////				UW_n = line(UW_iter->getphysical(),2,UW_iter->getdress(),UW_iter->getiflavor(),UW_iter->getjflavor(),x_n,UW_iter->getxj());
////			else 
////				UW_n = line(UW_iter->getphysical(),2,UW_iter->getdress(),UW_iter->getiflavor(),UW_iter->getjflavor(),UW_iter->getxi(),x_n);
////		}
////		else{
////			x_n = x_o;
////
////			proposal = (v_iter->getx().getsitegenprobability(x_n))/beta;
////			if(ijselection==0)
////				UW_n = line(UW_iter->getphysical(),1,UW_iter->getdress(),UW_iter->getiflavor(),UW_iter->getjflavor(),x_n,UW_iter->getxj());
////			else 
////				UW_n = line(UW_iter->getphysical(),1,UW_iter->getdress(),UW_iter->getiflavor(),UW_iter->getjflavor(),UW_iter->getxi(),x_n);
////		}
////
////		trial = SetUWInter;
////		trial.transform(*UW_iter,UW_n);
////		trial.calvalue();
////
////		Gin_iter = v_iter->getconnected_Gin();
////		Gout_iter = v_iter->getconnected_Gout();
////		is_selfloop = (Gin_iter==Gout_iter);
////		Gin_n = *Gin_iter;
////		Gout_n = *Gout_iter;
////
////		if(is_selfloop){
////			Gin_n.setx(x_n,x_n);
////			probability_second = proposal 
////				* trial.getvalue() / SetUWInter.getvalue();
////		}
////		else{
////			Gin_n.setxj(x_n);
////			Gout_n.setxi(x_n);
////			loop_valuerh_n = Gin_n.getvalueij() * Gout_n.getvalueij() 
////					/ Gin_iter->getvalueij() / Gout_iter->getvalueij() 
////					* loop_iter->getvaluerh();
////			loop_valuelh_n = Gin_n.getvalueji() * Gout_n.getvalueji() 
////					/ Gin_iter->getvalueji() / Gout_iter->getvalueji() 
////					* loop_iter->getvaluelh();
////
////			if(ismeasuringlineinloop){
////				probability_second = proposal
////						* loop_valuerh_n / loop_iter->getvaluerh()
////						* trial.getvalue() / SetUWInter.getvalue();
////			}
////			else{
////				probability_second = proposal
////						* (loop_valuerh_n+loop_valuelh_n) 
////						/ (loop_iter->getvaluerh()+loop_iter->getvaluelh())
////						* trial.getvalue() / SetUWInter.getvalue();
////			}
////		}
////
////		if(is_selfloop){
////			Gin_n.setx(x_n,x_n);
////			probability_second = proposal
////				* Gin_n.getvalue() * trial.getvalue()
////				/ Gin_iter->getvalue() / SetUWInter.getvalue();
////		}
////		else{
////			Gin_n.setxj(x_n);
////			Gout_n.setxi(x_n);
////			probability_second = proposal
////				* Gin_n.getvalue() * Gout_n.getvalue() * trial.getvalue()
////				/ Gin_iter->getvalue() / Gout_iter->getvalue() / SetUWInter.getvalue();
////		}
////		//cout<<"* second update "<<endl;
////		//cout<<"proposal = "<<proposal<<endl;
////		//cout<<"Gin_n.getvalue()      	= "<<Gin_n.getvalue() 		<<endl;
////		//cout<<"Gout_n.getvalue()     	= "<<Gout_n.getvalue() 		<<endl;
////		//cout<<"trial.getvalue()      	= "<<trial.getvalue()		<<endl;
////		//cout<<"Gin_iter->getvalue()  	= "<<Gin_iter->getvalue() 	<<endl;
////		//cout<<"Gout_iter->getvalue()    = "<<Gout_iter->getvalue()	<<endl;
////		//cout<<"SetUWInter.getvalue()    = "<<SetUWInter.getvalue()	<<endl;
////		//cout<<"trial info."<<endl<<trial<<endl;
////		//cout<<"SetUWInter info."<<endl<<SetUWInter<<endl;
////
////		if(probability_second<0.) signfactor = -1;
////		else signfactor = 1;
////
////		x_o = v_iter->getx();
////		if(abs(probability_second)>1.0e-10){
////			v_iter->setx(x_n);
////			UW_iter->assign(UW_n);
////			if(is_selfloop) 
////				Gin_iter->assign(Gin_n);
////			else{
////				Gin_iter->assign(Gin_n);
////				Gout_iter->assign(Gout_n);
////				loop_iter->setvalue(loop_valuerh_n,loop_valuelh_n);
////			}
////			SetUWInter = trial;
////			sign *= signfactor;
////		}
////
////		if( abs(probability_first*probability_second - 1.)<1.0e-10 ){
////			cout<<"probability_first = "<<probability_first<<endl;
////			cout<<"probability_second = "<<probability_second<<endl;
////			cout<<"probability_first*probability_second = "<<probability_first*probability_second<<endl;
////			return true;
////		}
////		else{
////			cout<<"***** ALERT! DETAILED BALANCE CONDITION FOR TRANSFORMUW UPDATE FAILS! *****"<<endl;
////			cout<<"probability_first = "<<probability_first<<endl;
////			cout<<"probability_second = "<<probability_second<<endl;
////			cout<<"probability_first*probability_second = "<<probability_first*probability_second<<endl;
////			printconfiguration(cout);
////			return false;
////		}
////	}
////};*/
////
/////*bool diagram::checkdetailedbalance_loopspinflip(){
////	cout<<"** check detailed balance for loop spin flip update **"<<endl<<flush;
////	cout<<"*  first update  *"<<endl;
////	cout<<"   updatefermionicloops"<<endl<<flush;
////	updatefermionicloops();
////
////	//for(auto fl:fermionicloops) cout<<fl;
////	
////	int Nloop = fermionicloops.size();
////	int iloop = static_cast<int>(Nloop*unidist(mt));
////
////	vector<coordinate> xs;
////	vector<int> UWtypes;
////
////	list<fermionicloop>::iterator fl_iter = next(fermionicloops.begin(),iloop);
////
////	xs = fl_iter->getxvector();
////	UWtypes = fl_iter->getUWtypes();
////
////	fl_iter->proposeloopspinflip();
////
////	double probability_first = fl_iter->get_loopspinflip_probability();
////	//cout<<"probability = "<<probability<<endl;
////
////	int signfactor;
////	if(probability_first<0.) signfactor = -1;
////	else signfactor = 1;
////
////	fl_iter->loopspinflip();
////	sign *= signfactor;
////
////	cout<<"*  first update  *"<<endl;
////
////	updatefermionicloops();
////	fl_iter = next(fermionicloops.begin(),iloop);
////	fl_iter->proposeloopspinflip(xs,UWtypes);
////	double probability_second = fl_iter->get_loopspinflip_probability();
////
////	if(probability_second<0.) signfactor = -1;
////	else signfactor = 1;
////
////	fl_iter->loopspinflip();
////	sign *= signfactor;
////
////	if( abs(probability_first*probability_second - 1.)<1.0e-10 ){
////		cout<<"probability_first = "<<probability_first<<endl;
////		cout<<"probability_second = "<<probability_second<<endl;
////		cout<<"probability_first*probability_second = "<<probability_first*probability_second<<endl;
////		return true;
////	}
////	else{
////		cout<<"***** ALERT! DETAILED BALANCE CONDITION FAILS! *****"<<endl;
////		cout<<"probability_first = "<<probability_first<<endl;
////		cout<<"probability_second = "<<probability_second<<endl;
////		cout<<"probability_first*probability_second = "<<probability_first*probability_second<<endl;
////		printconfiguration(cout);
////		return false;
////	}
////
////	//return 0;
////};*/
////
/////*int diagram::interswap(){
////	//cout<<"* interswap module *"<<endl<<flush;
////	int Nlines = Gs.size()+Us.size()+Ws.size();
////	int index_i, index_j;
////	//int baredNorder;
////	list<line>::iterator line_j_iter;
////
////	if(measuringline->gettype()==0){
////		index_i = distance(Gs.begin(),measuringline);
////		//baredNorder = -1;
////	}
////	else if(measuringline->gettype()==1){
////		index_i = Gs.size()+distance(Us.begin(),measuringline);
////		//baredNorder = +1;
////	}
////	else{
////		index_i = Gs.size()+Us.size()+distance(Ws.begin(),measuringline);
////		//baredNorder = +1;
////	}
////
////	select_one_excluding_one(Nlines,index_i,index_j);
////
////	if(index_j<Gs.size()) line_j_iter = next(Gs.begin(),index_j);
////	else if(index_j<Gs.size()+Us.size()) line_j_iter = next(Us.begin(),index_j-Gs.size());
////	else line_j_iter = next(Ws.begin(),index_j-Gs.size()-Us.size());
////
////
////	int dNorder;
////	if(measuringline->gettype()==0 && line_j_iter->gettype()!=0) dNorder = -1;
////	else if(measuringline->gettype()!=0 && line_j_iter->gettype()==0) dNorder = +1;
////	else dNorder = 0;
////
////	//if(Gs.size()!=2*(Us.size()+Ws.size())) exit(EXIT_FAILURE);
////	//int dNorder = baredNorder + measuringline->getdress() - line_j_iter->getdress();
////	//cout<<"dNorder = "<<dNorder<<endl;
////	if(Norder+dNorder>(this->Nordermax)){
////		//cout<<"Noder exceeds the limit!"<<endl<<flush;
////		//printconfiguration(cout);
////		return 0;
////	}
////	if(!(checkconnectivityforinterswap(line_j_iter))){
////		//cout<<"interswap connectivity check"<<endl;
////		//printconfiguration(cout);
////		return 0;
////	}
////
////	//if(Norder==1 && Nselfloop==2){
////	//	cout<<"there is problem"<<endl;
////	//	exit(EXIT_FAILURE);
////	//}
////
////	//int NdressChoice = 1 + Nordermax - Norder - baredNorder - line_j_iter->getdress();
////	//if(NdressChoice<1) return 0;
////	
////	//int Ndress = static_cast<int>(NdressChoice*unidist(mt));
////	//int dNorder = baredNorder + Ndress - line_j_iter->getdress();
////
////	//cout<<"isG: "<<measuringline->getisG()<<endl<<flush;
////	//cout<<"NdressChoice = "<<NdressChoice<<endl;
////	//cout<<"Ndress = "<<Ndress<<endl<<flush;
////	//cout<<"dNorder = "<<dNorder<<endl<<flush;
////	//cout<<"measuringline info."<<endl<<flush;
////	//cout<<(*measuringline)<<endl<<flush;
////
////	line line_i_new(*measuringline);
////	line_i_new.setdress(measuringline->getdress());
////	line_i_new.setphysical(true);
////
////	//double proposal = static_cast<double>(NdressChoice) / (1 + Nordermax - Norder + line_j_iter->getdress());
////	double proposal = 1.0;
////
////	// Wow if I turn off below two lines, it works! Super strange... Orz...
////	//if(measuringline->gettype()==0) proposal *= 0.5;
////	//else proposal *= 2.0;
////
////	// Fermionic sign factor
////	if(measuringline->gettype()*line_j_iter->gettype()==0 && measuringline->gettype()!=line_j_iter->gettype()) 
////		proposal *= -1;
////
////	double probability = proposal
////		* Ri_N[Norder] / Ri_N[Norder+dNorder]
////		* line_i_new.getvalue() / line_j_iter->getvalue();
////
////	int signfactor;
////	if(probability<0.) signfactor = -1;
////	else signfactor = 1;
////
////	//if(NdressChoice==3 && measuringline->getisG() && Ndress==2){
////	//printconfiguration(cout);
////	//cout<<"* interswap info *"<<endl;
////	//cout<<"measuringline = "<<*measuringline<<endl;
////	//cout<<"line_i_new = "<<line_i_new<<endl;
////	//cout<<"line_j_iter = "<<*line_j_iter<<endl;
////	//cout<<"Norder = "<<Norder<<endl;
////	//cout<<"dNorder = "<<dNorder<<endl<<flush;
////	//cout<<"proposal = "<<proposal<<endl;
////	//cout<<"probability = "<<probability<<endl;
////	//cout<<"signfactor = "<<signfactor<<endl;
////	//cout<<"line_i_new.getvalue() = "<<line_i_new.getvalue()<<endl;
////	//cout<<"line_j_iter->getvalue() = "<<line_j_iter->getvalue()<<endl;
////	//cout<<endl<<endl<<flush;
////	//cout<<"type of measuringline = "<<measuringline->gettype()<<endl;
////	//cout<<"(index_i,index_j) = ("<<index_i<<","<<index_j<<")"<<endl;
////	//cout<<"NdressChoice = "<<NdressChoice<<endl;
////	//cout<<"NdressChoiceInverse = "<<(1 + Nordermax - Norder + line_j_iter->getdress())<<endl;
////	//cout<<"Ndress = "<<Ndress<<endl<<flush;
////	//}
////
////	if(unidist(mt)<abs(probability)){
////		//printconfiguration(cout);
////		measuringline->setphysical(true);
////		line_j_iter->setphysical(false);
////		measuringline = line_j_iter;
////
////		Norder += dNorder;
////		sign *= signfactor;
////
////		//cout<<"-->update accepted"<<endl;
////		return 3;
////	}
////
////	//cout<<"-->update rejected"<<endl;
////	return 1;
////};*/
/////*bool diagram::checkUlines(){
////	list<line>::iterator iter;
////	list<vertex>::iterator vi_iter, vj_iter;
////	for(iter=Us.begin(); iter!=Us.end(); ++iter){
////		vi_iter = iter->getconnected_vi();
////		vj_iter = iter->getconnected_vj();
////		if( !((vi_iter->getflavor()^vj_iter->getflavor()) && (vi_iter->getx()==vj_iter->getx())) ) return false;
////	}
////	return true;
////};
////bool diagram::checkWlines(){
////	list<line>::iterator iter;
////	list<vertex>::iterator vi_iter, vj_iter;
////	for(iter=Ws.begin(); iter!=Ws.end(); ++iter){
////		vi_iter = iter->getconnected_vi();
////		vj_iter = iter->getconnected_vj();
////		if( (vi_iter->getflavor()^vj_iter->getflavor()) && (vi_iter->getx()==vj_iter->getx()))	return false;
////	}
////	return true;
////};*/
////
