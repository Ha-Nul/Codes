#include <iterator>
#include <algorithm>
#include <fstream>
//#include "green.hpp"
#include "configuration.hpp"
#include "parameters.hpp"
#include "graph.hpp"
using namespace std;
using namespace parameters;


///////////////////
// configuration //
///////////////////

const Eigen::Matrix2d configuration::Zero2x2 = Eigen::MatrixXd::Zero(2,2);
const Eigen::Matrix2d configuration::Identity2x2 = Eigen::MatrixXd::Identity(2,2);

configuration::configuration(){};
void configuration::init(const vector<Eigen::Matrix2d>& Gt_i){
	//mWt.initWt_spinBoson( parameters::NWtau, parameters::g_coupling, parameters::w_boson);
	//Gt.initOp( Eigen::MatrixXd::Identity(2,2) );
	mWt.initWtSingleModeWaveguide( parameters::NWtau, parameters::g_coupling, parameters::Omega11, parameters::omega_c, parameters::Nkx );
	Gt.init( Gt_i );

	//ofstream ofstr("WtDiagMC.in");
	//mWt.print(ofstr);
	//ofstr.close();

	//for(int iv=0; iv<NvExt; iv++){
	//	extFlavor(iv/Nflavor,iv%Nflavor) = static_cast<int>(Nflavor*unidist(mt));
	//}
	extFlavor(0,0) = 1;
	extFlavor(0,1) = 1;
	extFlavor(1,0) = 1;
	extFlavor(1,1) = 1;

	im2 = 2;
	v_beta.settime(parameters::beta);

	vector<coordinate> xi(NvExt);
	xi[0].setref();
	for(auto it=next(xi.begin()); it!=xi.end(); ++it) it->genrand();
	sort(xi.begin(), xi.end());
	for(int i=0; i<NvExt;  i++) vertices.push_back(xi[i]);

	for(int iv=1; iv<NvExt; iv++) indexExtVertices.push_back(iv);

	Ws.push_back( mWt.getWt(vertices[0],vertices[2]) );
	Ws.push_back( mWt.getWt(vertices[1],vertices[3]) );

	list_v_to_W.push_back(0);
	list_v_to_W.push_back(1);
	list_v_to_W.push_back(0);
	list_v_to_W.push_back(1);

	list_W_to_vi.push_back(0);
	list_W_to_vi.push_back(1);

	list_W_to_vj.push_back(2);
	list_W_to_vj.push_back(3);

	GtTop.push_back( Gt.getSzGtSz(vertices[0],vertices[1]) );
	GtBottom.push_back( Gt.getSzGtSz(vertices[2],vertices[3]) );

	calculate_WProduct();
	calculate_GtTopProduct();
	calculate_GtBottomProduct();
};
void configuration::init(const vector<double>& mWt_i,const vector<Eigen::Matrix2d>& Gt_i){
	mWt.init( mWt_i );
	//Gt.initOp( Eigen::MatrixXd::Identity(2,2) );
	Gt.init( Gt_i );

	//for(int iv=0; iv<NvExt; iv++){
	//	extFlavor(iv/Nflavor,iv%Nflavor) = static_cast<int>(Nflavor*unidist(mt));
	//}
	extFlavor(0,0) = 1;
	extFlavor(0,1) = 1;
	extFlavor(1,0) = 1;
	extFlavor(1,1) = 1;

	im2 = 2;
	v_beta.settime(parameters::beta);

	vector<coordinate> xi(NvExt);
	xi[0].setref();
	for(auto it=next(xi.begin()); it!=xi.end(); ++it) it->genrand();
	sort(xi.begin(), xi.end());
	for(int i=0; i<NvExt;  i++) vertices.push_back(xi[i]);

	for(int iv=1; iv<NvExt; iv++) indexExtVertices.push_back(iv);

	Ws.push_back( mWt.getWt(vertices[0],vertices[2]) );
	Ws.push_back( mWt.getWt(vertices[1],vertices[3]) );

	list_v_to_W.push_back(0);
	list_v_to_W.push_back(1);
	list_v_to_W.push_back(0);
	list_v_to_W.push_back(1);

	list_W_to_vi.push_back(0);
	list_W_to_vi.push_back(1);

	list_W_to_vj.push_back(2);
	list_W_to_vj.push_back(3);

	GtTop.push_back( Gt.getSzGtSz(vertices[0],vertices[1]) );
	GtBottom.push_back( Gt.getSzGtSz(vertices[2],vertices[3]) );

	calculate_WProduct();
	calculate_GtTopProduct();
	calculate_GtBottomProduct();
};
int configuration::getVerticesSize(){ return vertices.size(); };
int configuration::getWsSize(){ return Ws.size(); };
int configuration::getGtTopSize(){ return GtTop.size(); };
int configuration::getim2(){ return im2; };
int configuration::getiv(const int& iW){ return list_W_to_vi[iW]; };
int configuration::getjv(const int& iW){ return list_W_to_vj[iW]; };
int configuration::getiW(const int& iv){ return list_v_to_W[iv]; };
int configuration::getiv_partner(const int& iv_i){
	int iW = list_v_to_W[iv_i];
	return (list_W_to_vi[iW]==iv_i ? list_W_to_vj[iW] : list_W_to_vi[iW]);
};
double configuration::get_tinterval(const int& iv_i){ 
	if(iv_i==vertices.size()-1) return v_beta.time - vertices[iv_i-1].time;
	else return vertices[iv_i+1].time - vertices[iv_i-1].time;
}
double configuration::getWProduct(){ return WProduct; };
coordinate configuration::getVertex(const int& iv_i){ return vertices[iv_i]; };
coordinate configuration::getProposedVi(){ return vi_new; };
coordinate configuration::getProposedVj(){ return vj_new; };
bool configuration::getIsOnTop(){ return is_on_top; };
void configuration::calculate_WProduct(){
	WProduct = 1.0;
	for(auto& it: Ws) WProduct *= it;
};
void configuration::calculate_GtTopProduct(){
	GtTopProduct = Identity2x2;
	for(auto it=GtTop.rbegin(); it!=GtTop.rend(); ++it) GtTopProduct *= (*it);
};
void configuration::calculate_GtTopProductProposed(){
	GtTopProductProposed = Identity2x2;
	for(auto it=GtTopProposed.rbegin(); it!=GtTopProposed.rend(); ++it) GtTopProductProposed *= (*it);
};
void configuration::calculate_GtBottomProduct(){
	GtBottomProduct = Identity2x2;
	for(auto it=GtBottom.rbegin(); it!=GtBottom.rend(); ++it) GtBottomProduct *= (*it);
};
void configuration::calculate_GtBottomProductProposed(){
	GtBottomProductProposed = Identity2x2;
	for(auto it=GtBottomProposed.rbegin(); it!=GtBottomProposed.rend(); ++it) GtBottomProductProposed *= (*it);
};
double configuration::getWeight(){
	return WProduct
		* GtTopProduct(extFlavor(0,1), extFlavor(0,0))
		* GtBottomProduct(extFlavor(1,0), extFlavor(1,1));
};
double configuration::getProposedWeight(){
	return WProductProposed
		* GtTopProductProposed(extFlavorProposed(0,1), extFlavorProposed(0,0))
		* GtBottomProductProposed(extFlavorProposed(1,0), extFlavorProposed(1,1));
};
vector<coordinate> configuration::getExtVertices(){
	vector<coordinate> extVertices;

	//cout<<endl<<"indexExtVertices"<<endl;
	//for(auto& iterRef: indexExtVertices) cout<<iterRef;
	//cout<<endl;

	for(int i=0; i<NvExt-1; i++){
		extVertices.push_back( vertices[indexExtVertices[i]] );
	}

	return extVertices;
}
vector<int> configuration::getFlavorExt(){
	vector<int> VertexExt2pt;
	VertexExt2pt.push_back( extFlavor(0,0)*Nflavor + extFlavor(1,0) );
	VertexExt2pt.push_back( extFlavor(0,1)*Nflavor + extFlavor(1,1) );

	return VertexExt2pt;
};
//void configuration::initProposal(){
//	//WsProposed = Ws;
//	//GtTopProposed = GtTop;
//	//GtBottomProposed = GtBottom;
//};
void configuration::proposeInsertVertices(const int& iv, const int& jv){
	//cout<<endl<<"(iv,jv): ("<<iv<<","<<jv<<")"<<endl;
	iv_new = iv;
	jv_new = jv;

	int Nvs = vertices.size();
	vi_new.genrand( vertices[iv],vertices[iv+1] );
	if(jv==Nvs-1) vj_new.genrand(vertices[jv],v_beta);
	else vj_new.genrand(vertices[jv],vertices[jv+1]);

	if( iv==im2-1 || jv==im2-1 ){
		if(unidist(mt)<0.5) is_on_top = true;
		else is_on_top = false;
	}

	if(iv<im2){
		if(iv==im2-1 && is_on_top==false) is_vi_on_top = false;
		else is_vi_on_top = true;
	}
	else{ is_vi_on_top = false; }
	if(jv<im2){
		if(jv==im2-1 && is_on_top==false) is_vj_on_top = false;
		else is_vj_on_top = true;
	}
	else{ is_vj_on_top = false; }

	verticesProposed = vertices;
	verticesProposed.insert(verticesProposed.begin()+jv+1, vj_new);
	verticesProposed.insert(verticesProposed.begin()+iv+1, vi_new);

	//cout<<"Hmm, strange"<<endl<<flush;

	//int Nv_new = Nvs+2;
	if( is_vi_on_top&&is_vj_on_top ){
		//cout<<"both top"<<endl;
		im2_new = im2+2;

		GtTopProposed = GtTop;
		if(jv==im2-1){
			GtTopProposed.back() =  Gt.getGtSz(vertices[jv-1],vertices[jv]) ;
			GtTopProposed.insert( GtTopProposed.end(), Gt.getSzGtSz(vertices[jv],vj_new) );
		}
		else if(jv==im2-2){
			GtTopProposed.insert( GtTopProposed.begin()+jv, Gt.getGtSz(vertices[jv],vj_new) );
			GtTopProposed.back() =  Gt.getSzGtSz(vj_new,vertices[jv+1]) ;
		}
		else{
			GtTopProposed.insert( GtTopProposed.begin()+jv, Gt.getGtSz(vertices[jv],vj_new) );
			GtTopProposed[jv+1] =  Gt.getGtSz(vj_new,vertices[jv+1]) ;
		}

		GtTopProposed.insert( GtTopProposed.begin()+iv, Gt.getGtSz(vertices[iv],vi_new) );
		GtTopProposed[iv+1] = Gt.getGtSz(vi_new,vertices[iv+1]);

		calculate_GtTopProductProposed();
		GtBottomProductProposed = GtBottomProduct;
	}
	else if( !(is_vi_on_top||is_vj_on_top) ){
		//cout<<"both down"<<endl;
		im2_new = im2;

		GtBottomProposed = GtBottom;

		if(jv==Nvs-1){
			GtBottomProposed.back() =  Gt.getGtSz(vertices[jv-1],vertices[jv]) ;
			GtBottomProposed.insert( GtBottomProposed.end(), Gt.getSzGtSz(vertices[jv],vj_new) );
		}
		else if(jv==Nvs-2){
			GtBottomProposed.insert( GtBottomProposed.begin()+jv-im2, Gt.getGtSz(vertices[jv],vj_new) );
			GtBottomProposed.back() =  Gt.getSzGtSz(vj_new,vertices[jv+1]) ;
		}
		else{
			GtBottomProposed.insert( GtBottomProposed.begin()+jv-im2, Gt.getGtSz(vertices[jv],vj_new) );
			GtBottomProposed[jv-im2+1] =  Gt.getGtSz(vj_new,vertices[jv+1]) ;
		}

		if(iv==im2-1)
			GtBottomProposed.insert( GtBottomProposed.begin(), Gt.getGtSz(vi_new,vertices[im2]));
		else{
			GtBottomProposed.insert( GtBottomProposed.begin()+iv-im2, Gt.getGtSz(vertices[iv],vi_new) );
			GtBottomProposed[iv-im2+1] = Gt.getGtSz(vi_new,vertices[iv+1]);
		}

		calculate_GtBottomProductProposed();
		GtTopProductProposed = GtTopProduct;
	}
	else{
		//cout<<"one top one bottom"<<endl;
		im2_new = im2+1;
		GtTopProposed = GtTop;
		GtBottomProposed = GtBottom;

		//cout<<"vi_new: "<<vi_new<<endl;
		//cout<<"vj_new: "<<vj_new<<endl;

		if(jv==Nvs-1){
			GtBottomProposed.back() =  Gt.getGtSz(vertices[jv-1],vertices[jv]) ;
			GtBottomProposed.push_back( Gt.getSzGtSz(vertices[jv],vj_new) );
		}
		else if(jv==Nvs-2){
			GtBottomProposed.insert( GtBottomProposed.begin()+jv-im2, Gt.getGtSz(vertices[jv],vj_new) );
			GtBottomProposed.back() =  Gt.getSzGtSz(vj_new,vertices[jv+1]) ;
		}
		else if(jv==im2-1)
			GtBottomProposed.insert( GtBottomProposed.begin(), Gt.getGtSz(vj_new,vertices[im2]));
		else{
			GtBottomProposed.insert( GtBottomProposed.begin()+jv-im2, Gt.getGtSz(vertices[jv],vj_new) );
			GtBottomProposed[jv-im2+1] =  Gt.getGtSz(vj_new,vertices[jv+1]) ;
		}

		//cout<<"where?"<<endl<<flush;

		if(iv==im2-1){
			GtTopProposed.back() =  Gt.getGtSz(vertices[iv-1],vertices[iv]) ;
			GtTopProposed.push_back( Gt.getSzGtSz(vertices[iv],vi_new) );
		 }
		else if(iv==im2-2){
			GtTopProposed.insert( GtTopProposed.begin()+iv, Gt.getGtSz(vertices[iv],vi_new) );
			GtTopProposed.back() =  Gt.getSzGtSz(vi_new,vertices[iv+1]) ;
		}
		else{
			GtTopProposed.insert( GtTopProposed.begin()+iv, Gt.getGtSz(vertices[iv],vi_new) );
			GtTopProposed[iv+1] =  Gt.getGtSz(vi_new,vertices[iv+1]) ;
		}

		//cout<<"there?"<<endl<<flush;

		//if(jv==im2-1)
		//	GtBottomProposed.insert( GtBottomProposed.begin(), Gt.getGtSz(vi_new,vertices[im2]));
		//else
		//	GtBottomProposed.insert( GtBottom.begin()+jv-im2, Gt.getGtSz(vertices[jv],vj_new) );
		//if(jv!=im2+GtBottom.size()-1) GtBottomProposed[jv-im2+1] =  Gt.getGtSz(vj_new,vertices[jv-im2+1]);

		//GtTopProposed.insert( GtTopProposed.begin()+iv, Gt.getGtSz(vertices[iv],vj_new) );
		//if(iv!=GtTop.size()-1) GtTopProposed[iv+1] =  Gt.getGtSz(vj_new,vertices[iv+1]) ;

		//for(int i=0; i<GtTopProposed.size()-1; i++){
		//	GtTopProposed[i] = Gt.getGtSz(verticesProposed[i],verticesProposed[i+1]);
		//}
		//GtTopProposed.back() = Gt.getSzGtSz(verticesProposed[im2_new-2],verticesProposed[im2_new-1]);
		//
		//for(int i=im2_new; i<GtBottomProposed.size()-1; i++){
		//	GtBottomProposed[i] = Gt.getGtSz(verticesProposed[i],verticesProposed[i+1]);
		//}
		//GtBottomProposed.back() = Gt.getSzGtSz(verticesProposed[Nv_new-2],verticesProposed[Nv_new-1]);

		calculate_GtTopProductProposed();
		calculate_GtBottomProductProposed();
	}

	extFlavorProposed = extFlavor;

	Wt_new = mWt.getWt(vi_new, vj_new);
	WProductProposed = WProduct * Wt_new;
};
void configuration::proposeInsertVertices(const int& iv, const coordinate& vi_new_i, const int& jv, const coordinate& vj_new_i, const bool& is_on_top_i){
	//cout<<endl<<"(iv,jv): ("<<iv<<","<<jv<<")"<<endl;
	iv_new = iv;
	jv_new = jv;

	int Nvs = vertices.size();
	//vi_new.genrand( vertices[iv],vertices[iv+1] );
	//if(jv==Nvs-1) vj_new.genrand(vertices[jv],v_beta);
	//else vj_new.genrand(vertices[jv],vertices[jv+1]);
	vi_new = vi_new_i;
	vj_new = vj_new_i;

	//if( iv==im2-1 || jv==im2-1 ){
	//	if(unidist(mt)<0.5) is_on_top = true;
	//	else is_on_top = false;
	//}
	is_on_top = is_on_top_i;

	if(iv<im2){
		if(iv==im2-1 && is_on_top==false) is_vi_on_top = false;
		else is_vi_on_top = true;
	}
	else{ is_vi_on_top = false; }
	if(jv<im2){
		if(jv==im2-1 && is_on_top==false) is_vj_on_top = false;
		else is_vj_on_top = true;
	}
	else{ is_vj_on_top = false; }

	verticesProposed = vertices;
	verticesProposed.insert(verticesProposed.begin()+jv+1, vj_new);
	verticesProposed.insert(verticesProposed.begin()+iv+1, vi_new);

	//cout<<"Hmm, strange"<<endl<<flush;

	int Nv_new = Nvs+2;
	if( is_vi_on_top&&is_vj_on_top ){
		//cout<<"both top"<<endl;
		im2_new = im2+2;

		GtTopProposed = GtTop;
		if(jv==im2-1){
			GtTopProposed.back() =  Gt.getGtSz(vertices[jv-1],vertices[jv]) ;
			GtTopProposed.insert( GtTopProposed.end(), Gt.getSzGtSz(vertices[jv],vj_new) );
		}
		else if(jv==im2-2){
			GtTopProposed.insert( GtTopProposed.begin()+jv, Gt.getGtSz(vertices[jv],vj_new) );
			GtTopProposed.back() =  Gt.getSzGtSz(vj_new,vertices[jv+1]) ;
		}
		else{
			GtTopProposed.insert( GtTopProposed.begin()+jv, Gt.getGtSz(vertices[jv],vj_new) );
			GtTopProposed[jv+1] =  Gt.getGtSz(vj_new,vertices[jv+1]) ;
		}

		GtTopProposed.insert( GtTopProposed.begin()+iv, Gt.getGtSz(vertices[iv],vi_new) );
		GtTopProposed[iv+1] = Gt.getGtSz(vi_new,vertices[iv+1]);

		calculate_GtTopProductProposed();
		GtBottomProductProposed = GtBottomProduct;
	}
	else if( !(is_vi_on_top||is_vj_on_top) ){
		//cout<<"both down"<<endl;
		im2_new = im2;

		GtBottomProposed = GtBottom;

		if(jv==Nvs-1){
			GtBottomProposed.back() =  Gt.getGtSz(vertices[jv-1],vertices[jv]) ;
			GtBottomProposed.insert( GtBottomProposed.end(), Gt.getSzGtSz(vertices[jv],vj_new) );
		}
		else if(jv==Nvs-2){
			GtBottomProposed.insert( GtBottomProposed.begin()+jv-im2, Gt.getGtSz(vertices[jv],vj_new) );
			GtBottomProposed.back() =  Gt.getSzGtSz(vj_new,vertices[jv+1]) ;
		}
		else{
			GtBottomProposed.insert( GtBottomProposed.begin()+jv-im2, Gt.getGtSz(vertices[jv],vj_new) );
			GtBottomProposed[jv-im2+1] =  Gt.getGtSz(vj_new,vertices[jv+1]) ;
		}

		if(iv==im2-1)
			GtBottomProposed.insert( GtBottomProposed.begin(), Gt.getGtSz(vi_new,vertices[im2]));
		else{
			GtBottomProposed.insert( GtBottomProposed.begin()+iv-im2, Gt.getGtSz(vertices[iv],vi_new) );
			GtBottomProposed[iv-im2+1] = Gt.getGtSz(vi_new,vertices[iv+1]);
		}

		calculate_GtBottomProductProposed();
		GtTopProductProposed = GtTopProduct;
	}
	else{
		//cout<<"one top one bottom"<<endl;
		im2_new = im2+1;
		GtTopProposed = GtTop;
		GtBottomProposed = GtBottom;

		//cout<<"vi_new: "<<vi_new<<endl;
		//cout<<"vj_new: "<<vj_new<<endl;

		if(jv==Nvs-1){
			GtBottomProposed.back() =  Gt.getGtSz(vertices[jv-1],vertices[jv]) ;
			GtBottomProposed.push_back( Gt.getSzGtSz(vertices[jv],vj_new) );
		}
		else if(jv==Nvs-2){
			GtBottomProposed.insert( GtBottomProposed.begin()+jv-im2, Gt.getGtSz(vertices[jv],vj_new) );
			GtBottomProposed.back() =  Gt.getSzGtSz(vj_new,vertices[jv+1]) ;
		}
		else if(jv==im2-1)
			GtBottomProposed.insert( GtBottomProposed.begin(), Gt.getGtSz(vj_new,vertices[im2]));
		else{
			GtBottomProposed.insert( GtBottomProposed.begin()+jv-im2, Gt.getGtSz(vertices[jv],vj_new) );
			GtBottomProposed[jv-im2+1] =  Gt.getGtSz(vj_new,vertices[jv+1]) ;
		}

		//cout<<"where?"<<endl<<flush;

		if(iv==im2-1){
			GtTopProposed.back() =  Gt.getGtSz(vertices[iv-1],vertices[iv]) ;
			GtTopProposed.push_back( Gt.getSzGtSz(vertices[iv],vi_new) );
		 }
		else if(iv==im2-2){
			GtTopProposed.insert( GtTopProposed.begin()+iv, Gt.getGtSz(vertices[iv],vi_new) );
			GtTopProposed.back() =  Gt.getSzGtSz(vi_new,vertices[iv+1]) ;
		}
		else{
			GtTopProposed.insert( GtTopProposed.begin()+iv, Gt.getGtSz(vertices[iv],vi_new) );
			GtTopProposed[iv+1] =  Gt.getGtSz(vi_new,vertices[iv+1]) ;
		}

		//cout<<"there?"<<endl<<flush;

		//if(jv==im2-1)
		//	GtBottomProposed.insert( GtBottomProposed.begin(), Gt.getGtSz(vi_new,vertices[im2]));
		//else
		//	GtBottomProposed.insert( GtBottom.begin()+jv-im2, Gt.getGtSz(vertices[jv],vj_new) );
		//if(jv!=im2+GtBottom.size()-1) GtBottomProposed[jv-im2+1] =  Gt.getGtSz(vj_new,vertices[jv-im2+1]);

		//GtTopProposed.insert( GtTopProposed.begin()+iv, Gt.getGtSz(vertices[iv],vj_new) );
		//if(iv!=GtTop.size()-1) GtTopProposed[iv+1] =  Gt.getGtSz(vj_new,vertices[iv+1]) ;

		//for(int i=0; i<GtTopProposed.size()-1; i++){
		//	GtTopProposed[i] = Gt.getGtSz(verticesProposed[i],verticesProposed[i+1]);
		//}
		//GtTopProposed.back() = Gt.getSzGtSz(verticesProposed[im2_new-2],verticesProposed[im2_new-1]);
		//
		//for(int i=im2_new; i<GtBottomProposed.size()-1; i++){
		//	GtBottomProposed[i] = Gt.getGtSz(verticesProposed[i],verticesProposed[i+1]);
		//}
		//GtBottomProposed.back() = Gt.getSzGtSz(verticesProposed[Nv_new-2],verticesProposed[Nv_new-1]);

		calculate_GtTopProductProposed();
		calculate_GtBottomProductProposed();
	}

	extFlavorProposed = extFlavor;

	Wt_new = mWt.getWt(vi_new, vj_new);
	WProductProposed = WProduct * Wt_new;
};

void configuration::updateInsertVertices(){
	vertices = verticesProposed;

	Ws.push_back( Wt_new );
	WProduct = WProductProposed;

	if( is_vi_on_top&&is_vj_on_top ){
		GtTop = GtTopProposed;
		GtTopProduct = GtTopProductProposed;
	}
	else if( !(is_vi_on_top||is_vj_on_top) ){
		GtBottom = GtBottomProposed;
		GtBottomProduct = GtBottomProductProposed;
	}
	else{
		GtTop = GtTopProposed;
		GtTopProduct = GtTopProductProposed;

		GtBottom = GtBottomProposed;
		GtBottomProduct = GtBottomProductProposed;
	}

	list_v_to_W.insert(list_v_to_W.begin()+jv_new+1, Ws.size()-1);
	list_v_to_W.insert(list_v_to_W.begin()+iv_new+1, Ws.size()-1);

	for(auto& iterRef: list_W_to_vi){
		if( iterRef>jv_new ) iterRef++;
		if( iterRef>iv_new ) iterRef++;
	}
	for(auto& iterRef: list_W_to_vj){
		if( iterRef>jv_new ) iterRef++;
		if( iterRef>iv_new ) iterRef++;
	}

	list_W_to_vi.push_back(iv_new+1);
	list_W_to_vj.push_back(jv_new+2);

	//if( iv_new==im2-1 || jv_new==im2-1 ){
	//	if( is_on_top ) indexExtVertices[1] += 1;
	//	else indexExtVertices[2] += 1;
	//}

	im2 = im2_new;
	indexExtVertices[0] = im2-1;
	indexExtVertices[1] = im2;
	indexExtVertices[2] = vertices.size()-1;
};
void configuration::proposeRemoveVertices(const int& iW_i){
	iW_new = iW_i;
	iv_new = list_W_to_vi[iW_i];
	jv_new = list_W_to_vj[iW_i];
	vi_new = vertices[iv_new];
	vj_new = vertices[jv_new];

	//cout<<"(iv_new,jv_new): ("<<iv_new<<","<<jv_new<<")"<<endl;

	int Nvs = vertices.size();

	extFlavorProposed = extFlavor;
	
	if(iv_new<im2) is_vi_on_top = true;
	else is_vi_on_top = false;

	if(jv_new<im2) is_vj_on_top = true;
	else is_vj_on_top = false;

	//cout<<"Here?"<<endl<<flush;
	if( is_vi_on_top&&is_vj_on_top ){
		//cout<<"both on top"<<endl;
		im2_new = im2-2;
		GtTopProposed = GtTop;
		if(jv_new==im2-1){
			GtTopProposed.erase( GtTopProposed.end()-1 );
			GtTopProposed.back() = Gt.getSzGtSz(vertices[jv_new-2],vertices[jv_new-1]);
		}
		else if(jv_new==im2-2){
			GtTopProposed.erase(GtTopProposed.end()-2);
			GtTopProposed.back() = Gt.getSzGtSz(vertices[jv_new-1],vertices[jv_new+1]);
		}
		else{
			GtTopProposed.erase(GtTopProposed.begin()+jv_new);
			GtTopProposed[jv_new-1] = Gt.getGtSz(vertices[jv_new-1],vertices[jv_new+1]);
		}

		GtTopProposed.erase(GtTopProposed.begin()+iv_new);
		if( iv_new==GtTopProposed.size() )
			GtTopProposed[iv_new-1] = Gt.getSzGtSz(vertices[iv_new-1],vertices[iv_new+1]);
		else
			GtTopProposed[iv_new-1] = Gt.getGtSz(vertices[iv_new-1],vertices[iv_new+1]);

		calculate_GtTopProductProposed();
		GtBottomProductProposed = GtBottomProduct;
	}
	else if( !(is_vi_on_top||is_vj_on_top) ){
		//cout<<"both on bottom"<<endl;
		//cout<<"(iv_new,jv_new)=("<<iv_new<<","<<jv_new<<")"<<endl;
		im2_new = im2;
		GtBottomProposed = GtBottom;
		if(jv_new==Nvs-1){
			GtBottomProposed.erase( GtBottomProposed.end()-1 );
			GtBottomProposed.back() =  Gt.getSzGtSz(vertices[jv_new-2],vertices[jv_new-1]);
		}
		else if(jv_new==Nvs-2){
			GtBottomProposed.erase( GtBottomProposed.end()-2 );
			GtBottomProposed.back() =  Gt.getSzGtSz(vertices[jv_new-1],vertices[jv_new+1]);
		}
		else{
			GtBottomProposed.erase( GtBottomProposed.begin()+jv_new-im2 );
			GtBottomProposed[jv_new-im2-1] = Gt.getGtSz(vertices[jv_new-1],vertices[jv_new+1]);
		}
		//cout<<"Should be here!"<<endl<<flush;

		if(iv_new==im2)
			GtBottomProposed.erase( GtBottomProposed.begin() );
		else{
			GtBottomProposed.erase( GtBottomProposed.begin()+iv_new-im2 );
			if( iv_new-im2==GtBottomProposed.size() )
				GtBottomProposed[iv_new-im2-1] = Gt.getSzGtSz(vertices[iv_new-1],vertices[iv_new+1]);
			else
				GtBottomProposed[iv_new-im2-1] = Gt.getGtSz(vertices[iv_new-1],vertices[iv_new+1]);
		}

		calculate_GtBottomProductProposed();
		GtTopProductProposed = GtTopProduct;

		//cout<<"It's exhausting"<<endl<<flush;
	}
	else{
		//cout<<"one on top one on bottom"<<endl;
		//cout<<"iv_new = "<<iv_new<<endl;
		//cout<<"jv_new = "<<jv_new<<endl;

		im2_new = im2-1;
		GtTopProposed = GtTop;
		GtBottomProposed = GtBottom;

		if(jv_new==Nvs-1){
			//cout<<"jv_new==Nvs-1"<<endl<<flush;
			GtBottomProposed.erase( GtBottomProposed.end()-1 );
			GtBottomProposed.back() =  Gt.getSzGtSz(vertices[jv_new-2],vertices[jv_new-1]);
		}
		else if(jv_new==Nvs-2){
			GtBottomProposed.erase( GtBottomProposed.end()-2 );
			GtBottomProposed.back() =  Gt.getSzGtSz(vertices[jv_new-1],vertices[jv_new+1]);
		}
		else if(jv_new==im2){
			GtBottomProposed.erase( GtBottomProposed.begin() );
		}
		else{
			GtBottomProposed.erase(GtBottomProposed.begin()+jv_new-im2);
			GtBottomProposed[jv_new-im2-1] = Gt.getGtSz(vertices[jv_new-1],vertices[jv_new+1]);
		}


		if(iv_new==im2-1){
			//cout<<"iv_new==im2-1"<<endl<<flush;

			GtTopProposed.erase( GtTopProposed.end()-1 );

			//cout<<"Hmm"<<endl<<flush;
			//cout<<"vertices[iv_new-2] = "<<vertices[iv_new-2]<<endl<<flush;
			//cout<<"vertices[iv_new-1] = "<<vertices[iv_new-1]<<endl<<flush;
			//cout<<"Gt.getSzGtSz(vertices[iv_new-2],vertices[iv_new-1]) = "<<endl<<Gt.getSzGtSz(vertices[iv_new-2],vertices[iv_new-1])<<endl;
			//cout<<"GtBottomProposed list: "<<endl;
			//for(auto& it : GtBottomProposed) cout<<it<<endl<<"---"<<endl;
			//cout<<endl;
			//cout<<"GtTopProposed list: "<<endl;
			//for(auto& it : GtTopProposed) cout<<it<<endl<<"---"<<endl;
			//cout<<endl;
			//cout<<"Still fine, isn't it?"<<endl;

			GtTopProposed.back() = Gt.getSzGtSz(vertices[iv_new-2],vertices[iv_new-1]);
		}
		else if(iv_new==im2-2){
			GtTopProposed.erase(GtTopProposed.end()-2);
			GtTopProposed.back() = Gt.getSzGtSz(vertices[iv_new-1],vertices[iv_new+1]);
		}
		else{
			GtTopProposed.erase(GtTopProposed.begin()+iv_new);
			GtTopProposed[iv_new-1] = Gt.getGtSz(vertices[iv_new-1],vertices[iv_new+1]);
		}
		//cout<<"Seriously where?"<<endl<<flush;

		calculate_GtTopProductProposed();
		calculate_GtBottomProductProposed();
	}
	//cout<<"Where?"<<endl<<flush;

	//im2_new = GtTopProposed.size()+1;
	//cout<<"** im2_new = "<<im2_new<<" in proposeRemoveVertices **"<<endl<<flush;

	Wt_new = Ws[iW_i];
	WProductProposed = WProduct / Wt_new;
	//cout<<"Wt_new = "<<Wt_new<<endl;
	//cout<<"WProduct = "<<WProduct<<endl;
	//cout<<"WProductProposed = "<<WProductProposed<<endl;

	list_v_to_W_proposed = list_v_to_W;
	list_W_to_vi_proposed = list_W_to_vi;
	list_W_to_vj_proposed = list_W_to_vj;

	list_v_to_W_proposed.erase(list_v_to_W_proposed.begin()+jv_new);
	list_v_to_W_proposed.erase(list_v_to_W_proposed.begin()+iv_new);
	for(auto& iterRef: list_v_to_W_proposed){
		if(iterRef>iW_i) iterRef--;
	}

	list_W_to_vi_proposed.erase(list_W_to_vi_proposed.begin()+iW_i);
	list_W_to_vj_proposed.erase(list_W_to_vj_proposed.begin()+iW_i);

	for(auto& iterRef: list_W_to_vi_proposed){
		if( iterRef>jv_new ) iterRef--;
		if( iterRef>iv_new ) iterRef--;
	}
	for(auto& iterRef: list_W_to_vj_proposed){
		if( iterRef>jv_new ) iterRef--;
		if( iterRef>iv_new ) iterRef--;
	}
};
void configuration::updateRemoveVertices(){
	im2 = im2_new;
	list_v_to_W = list_v_to_W_proposed;
	list_W_to_vi = list_W_to_vi_proposed;
	list_W_to_vj = list_W_to_vj_proposed;

	vertices.erase(vertices.begin()+jv_new);
	vertices.erase(vertices.begin()+iv_new);

	Ws.erase(Ws.begin()+iW_new);
	WProduct = WProductProposed;

	if( is_vi_on_top&&is_vj_on_top ){
		GtTop = GtTopProposed;
		GtTopProduct = GtTopProductProposed;
	}
	else if( !(is_vi_on_top||is_vj_on_top) ){
		GtBottom = GtBottomProposed;
		GtBottomProduct = GtBottomProductProposed;
	}
	else{
		GtTop = GtTopProposed;
		GtTopProduct = GtTopProductProposed;

		GtBottom = GtBottomProposed;
		GtBottomProduct = GtBottomProductProposed;
	}
	indexExtVertices[0] = im2-1;
	indexExtVertices[1] = im2;
	indexExtVertices[2] = vertices.size()-1;

};
/*void configuration::proposeSwapWline(const int& iv_i, const int& jv_i){
	extFlavorProposed = extFlavor;

	iW_new= list_v_to_W[iv_i];
	iW_new2 = list_v_to_W[jv_i];

	int iv_partner = ( list_W_to_vi[iW_new]==iv_i ? list_W_to_vj[iW_new] : list_W_to_vi[iW_new] );
	int jv_partner = ( list_W_to_vi[iW_new2]==iv_i ? list_W_to_vj[iW_new2] : list_W_to_vi[iW_new2] );


	list_v_to_W_proposed = list_v_to_W;
	list_W_to_vi_proposed = list_W_to_vi;
	list_W_to_vj_proposed = list_W_to_vj;

	//list_v_to_W_proposed[iv_i] = iW_new2;
	//list_v_to_W_proposed[jv_i] = iW_new;

	if(iv_i<jv_partner){
		Wt_new = mWt.getWt(iv_i,jv_partner);

		list_W_to_vi_proposed[iW_new] = iv_i;
		list_W_to_vj_proposed[iW_new] = jv_partner;
	}
	else{
		Wt_new = mWt.getWt(jv_partner,iv_i);

		list_W_to_vi_proposed[iW_new] = jv_partner;
		list_W_to_vj_proposed[iW_new] = iv_i;
	}

	if(jv_i<iv_partner){
		Wt_new2 = mWt.getWt(jv_i,iv_partner);

		list_W_to_vi_proposed[iW_new2] = jv_i;
		list_W_to_vj_proposed[iW_new2] = iv_partner;
	}
	else{
		Wt_new2 = mWt.getWt(iv_partner,jv_i);

		list_W_to_vi_proposed[iW_new2] = iv_partner;
		list_W_to_vj_proposed[iW_new2] = jv_i;
	}

	WProductProposed = WProduct * Wt_new * Wt_new2 / Ws[iW_new] / Ws[iW_new2];
	GtTopProductProposed = GtTopProduct;
	GtBottomProductProposed = GtBottomProduct;
};
void configuration::updateSwapWline(){
	list_v_to_W = list_v_to_W_proposed;
	list_W_to_vi = list_W_to_vi_proposed;
	list_W_to_vj = list_W_to_vj_proposed;

	Ws[iW_new] = Wt_new;
	Ws[iW_new2] = Wt_new2;

	WProduct = WProductProposed;
};*/
void configuration::proposeShiftIm2(const int& im2_new_i){
	im2_new = im2_new_i;

	GtTopProposed = GtTop;
	GtBottomProposed = GtBottom;

	if(im2_new<im2){
		GtBottomProposed.insert( GtBottomProposed.begin(), Gt.getGtSz(vertices[im2-1],vertices[im2]) );
		GtTopProposed.back() = Gt.getGtSz(vertices[im2-2],vertices[im2-1]);
		for(int dim2=1; dim2<im2-im2_new; dim2++){
			GtBottomProposed.insert( GtBottomProposed.begin(), GtTopProposed.back() );
			GtTopProposed.erase( GtTopProposed.end()-1 );
		}
		GtTopProposed.erase( GtTopProposed.end()-1 );
		GtTopProposed.back() = Gt.getSzGtSz(vertices[im2_new-2],vertices[im2_new-1]);
	}
	else{
		GtTopProposed.back() = Gt.getGtSz(vertices[im2-2],vertices[im2-1]);
		GtTopProposed.push_back( Gt.getGtSz(vertices[im2-1],vertices[im2]) );
		for(int dim2=0; dim2<im2_new-im2-1; dim2++){
			GtTopProposed.push_back( GtBottomProposed.front() );
			GtBottomProposed.erase( GtBottomProposed.begin() );
		}
		GtBottomProposed.erase( GtBottomProposed.begin() );
		GtTopProposed.back() = Gt.getSzGtSz(vertices[im2_new-2],vertices[im2_new-1] );
	}

	extFlavorProposed = extFlavor;
	WProductProposed = WProduct;
	calculate_GtTopProductProposed();
	calculate_GtBottomProductProposed();

	list_v_to_W_proposed = list_v_to_W;
	list_W_to_vi_proposed = list_W_to_vi;
	list_W_to_vj_proposed = list_W_to_vj;
};
void configuration::updateShiftIm2(){
	im2 = im2_new;
	indexExtVertices[0] = im2-1;
	indexExtVertices[1] = im2;
	indexExtVertices[2] = vertices.size()-1;

	GtTop = GtTopProposed;
	GtTopProduct = GtTopProductProposed;

	GtBottom = GtBottomProposed;
	GtBottomProduct = GtBottomProductProposed;
};
void configuration::proposeShift(const int& iv){
	//cout<<"iv: "<<iv<<endl;
	iv_new = iv;
	double ti = vertices[iv-1].gettime();
	double tj;

	if(iv==vertices.size()-1) tj = parameters::beta;
	else tj = vertices[iv+1].gettime();
	coordinate v_next(tj);

	double dtau = tj - ti;
	//double tau_new = ti * dtau*unidist(mt);
	vi_new.genrand(vertices[iv-1], v_next);

	extFlavorProposed = extFlavor;

	if(iv<im2){
		//cout<<"top"<<endl;
		GtTopProposed = GtTop;
		if(iv==im2-1){
			GtTopProposed[iv-1] = Gt.getSzGtSz(vertices[iv-1],vi_new);
		}
		else{
			GtTopProposed[iv-1] = Gt.getGtSz(vertices[iv-1],vi_new);
			if(iv==im2-2) GtTopProposed[iv] = Gt.getSzGtSz(vi_new,vertices[iv+1]);
			else GtTopProposed[iv] = Gt.getGtSz(vi_new,vertices[iv+1]);
		}
	}
	else{
		//cout<<"bottom"<<endl;
		GtBottomProposed = GtBottom;
		//for(auto& it: GtBottomProposed) cout<<endl<<it<<endl;
		if(iv==vertices.size()-1){
			//cout<<"I think this is here"<<endl;
			//cout<<vertices[iv-1]<<vi_new<<endl;
			//cout<<endl<<Gt.getSzGtSz(vertices[iv-1], vi_new)<<endl;
			//cout<<iv-1-im2<<endl;

			GtBottomProposed[iv-im2-1] = Gt.getSzGtSz(vertices[iv-1], vi_new);
		}
		else if(iv==im2){
			if(iv==vertices.size()-2) GtBottomProposed[iv-im2] = Gt.getSzGtSz(vi_new,vertices[iv+1]);
			else GtBottomProposed[iv-im2] = Gt.getGtSz(vi_new,vertices[iv+1]);
		}
		else{
			GtBottomProposed[iv-1-im2] = Gt.getGtSz(vertices[iv-1],vi_new);
			if(iv==vertices.size()-2) GtBottomProposed[iv-im2] = Gt.getSzGtSz(vi_new,vertices[iv+1]);
			else GtBottomProposed[iv-im2] = Gt.getGtSz(vi_new,vertices[iv+1]);
		}
	}

	iW_new = list_v_to_W[iv];
	int jv = (list_W_to_vi[iW_new] == iv ? list_W_to_vj[iW_new]: list_W_to_vi[iW_new]);

	//cout<<vertices[iv]<<vertices[jv]<<endl;
	if(jv>iv) Wt_new = mWt.getWt(vi_new, vertices[jv]);
	else Wt_new = mWt.getWt(vertices[jv], vi_new);

	WProductProposed = WProduct*Wt_new/Ws[iW_new];
	if(iv<im2){
		calculate_GtTopProductProposed();
		GtBottomProductProposed = GtBottomProduct;
	}
	else{
		GtTopProductProposed = GtTopProduct;
		calculate_GtBottomProductProposed();
	}
};
void configuration::proposeShift(const int& iv, const double& t_n){
	iv_new = iv;
	vi_new.settime(t_n);

	extFlavorProposed = extFlavor;

	if(iv<im2){
		GtTopProposed = GtTop;
		if(iv==im2-1){
			GtTopProposed[iv-1] = Gt.getSzGtSz(vertices[iv-1],vi_new);
		}
		else{
			GtTopProposed[iv-1] = Gt.getGtSz(vertices[iv-1],vi_new);
			if(iv==im2-2) GtTopProposed[iv] = Gt.getSzGtSz(vi_new,vertices[iv+1]);
			else GtTopProposed[iv] = Gt.getGtSz(vi_new,vertices[iv+1]);
		}
	}
	else{
		GtBottomProposed = GtBottom;
		if(iv==vertices.size()-1){
			GtBottomProposed[iv-im2-1] = Gt.getSzGtSz(vertices[iv-1], vi_new);
		}
		else if(iv==im2){
			if(iv==vertices.size()-2) GtBottomProposed[iv-im2] = Gt.getSzGtSz(vi_new,vertices[iv+1]);
			else GtBottomProposed[iv-im2] = Gt.getGtSz(vi_new,vertices[iv+1]);
		}
		else{
			GtBottomProposed[iv-1-im2] = Gt.getGtSz(vertices[iv-1],vi_new);
			if(iv==vertices.size()-2) GtBottomProposed[iv-im2] = Gt.getSzGtSz(vi_new,vertices[iv+1]);
			else GtBottomProposed[iv-im2] = Gt.getGtSz(vi_new,vertices[iv+1]);
		}
	}

	iW_new = list_v_to_W[iv];
	int jv = (list_W_to_vi[iW_new] == iv ? list_W_to_vj[iW_new]: list_W_to_vi[iW_new]);

	if(jv>iv) Wt_new = mWt.getWt(vi_new, vertices[jv]);
	else Wt_new = mWt.getWt(vertices[jv], vi_new);

	WProductProposed = WProduct*Wt_new/Ws[iW_new];
	if(iv<im2){
		calculate_GtTopProductProposed();
		GtBottomProductProposed = GtBottomProduct;
	}
	else{
		GtTopProductProposed = GtTopProduct;
		calculate_GtBottomProductProposed();
	}
};

void configuration::updateShift(){
	//cout<<"* updateShift update *"<<endl;
	vertices[iv_new] = vi_new;
	Ws[iW_new] = Wt_new;
	WProduct = WProductProposed;

	if(iv_new<im2){
		if(iv_new==im2-1){
			GtTop[iv_new-1] = GtTopProposed[iv_new-1];
		}
		else{
			GtTop[iv_new-1] = GtTopProposed[iv_new-1];
			GtTop[iv_new] = GtTopProposed[iv_new];
		}
		GtTopProduct = GtTopProductProposed;
	}
	else{
		if(iv_new==im2){
			GtBottom.front() = GtBottomProposed.front();
		}
		else if(iv_new==vertices.size()-1){
			GtBottom.back() = GtBottomProposed.back();
		}
		else{
			GtBottom[iv_new-1-im2] = GtBottomProposed[iv_new-1-im2];
			GtBottom[iv_new-im2] = GtBottomProposed[iv_new-im2];
		}
		GtBottomProduct = GtBottomProductProposed;
	}
};
void configuration::proposeFlip(const int& iv){
	//cout<<"iv: "<<iv<<endl;
	iv_new = iv;

	extFlavorProposed = extFlavor;
	extFlavorProposed((iv/2), (iv%2)) ^= 1;

	WProductProposed = WProduct;
	GtTopProductProposed = GtTopProduct;
	GtBottomProductProposed = GtBottomProduct;
	
	/*double ti = vertices[iv-1].gettime();
	double tj;

	if(iv==vertices.size()-1) tj = parameters::beta;
	else tj = vertices[iv+1].gettime();
	coordinate v_next(tj);

	double dtau = tj - ti;
	double tau_new = ti * dtau*unidist(mt);
	vi_new.genrand(vertices[iv-1], v_next);*/

	/*if(iv<im2){
		cout<<"top"<<endl;
		GtTopProposed = GtTop;
		if(iv==im2-1){
			GtTopProposed[iv-1] = Gt.getSzGtSz(vertices[iv-1],vi_new);
		}
		else{
			GtTopProposed[iv-1] = Gt.getSzGt(vertices[iv-1],vi_new);
			GtTopProposed[iv] = Gt.getSzGt(vi_new,vertices[iv+1]);
		}
	}
	else{
		cout<<"bottom"<<endl;
		GtBottomProposed = GtBottom;
		for(auto& it: GtBottomProposed) cout<<endl<<it<<endl;
		if(iv==vertices.size()-1){
			//cout<<"I think this is here"<<endl;
			//cout<<vertices[iv-1]<<vi_new<<endl;
			//cout<<endl<<Gt.getSzGtSz(vertices[iv-1], vi_new)<<endl;
			//cout<<iv-1-im2<<endl;

			GtBottomProposed[iv-1-im2] = Gt.getSzGtSz(vertices[iv-1], vi_new);
		}
		else if(iv==im2){
			GtBottomProposed[iv-im2] = Gt.getSzGt(vi_new,vertices[iv+1]);
		}
		else{
			GtBottomProposed[iv-1-im2] = Gt.getSzGt(vertices[iv-1],vi_new);
			GtBottomProposed[iv-im2] = Gt.getSzGt(vi_new,vertices[iv+1]);
		}
	}

	iW = list_v_to_W[iv];
	int jv = (list_W_to_vi[iW] = iv ? list_W_to_vj[iW]: list_W_to_vi[iW]);

	Wt_new = mWt.getWt(vertices[iv], vertices[jv]);
	WProductProposed = WProduct*Wt_new/Ws[iW];
	if(iv<im2){
		calculate_GtTopProductProposed();
		GtBottomProductProposed = GtBottomProduct;
	}
	else{
		GtTopProductProposed = GtTopProduct;
		calculate_GtBottomProductProposed();
	}*/
};
void configuration::updateFlip(){
	//cout<<"* updateFlip update *"<<endl;
	extFlavor = extFlavorProposed;
};
bool configuration::checkProposedIrreducibility(){
	//cout<<"** check irreducibility module **"<<endl;
	int Nvertices = list_v_to_W_proposed.size();
	if(Nvertices<3) return true;
	else{
		graph topology;

		vector<vector<int> > adjLtmp(Nvertices);
		for(int i=0; i<Nvertices; i++) adjLtmp[i].reserve(3*Nvertices/2);

		for(int iv=0; iv<im2_new-1; iv++){
			adjLtmp[iv].push_back(iv+1);
			adjLtmp[iv+1].push_back(iv);
		}
		for(int iv=im2_new; iv<Nvertices-1; iv++){
			adjLtmp[iv].push_back(iv+1);
			adjLtmp[iv+1].push_back(iv);
		}
		for(int iW=0; iW<list_W_to_vi_proposed.size(); iW++){
			adjLtmp[ list_W_to_vi_proposed[iW] ].push_back( list_W_to_vj_proposed[iW] );
			adjLtmp[ list_W_to_vj_proposed[iW] ].push_back( list_W_to_vi_proposed[iW] );
		}

		//cout<<"adjLtmp"<<endl;
		//for(auto i: adjLtmp){
		//	for(auto j: i) cout<<j;
		//	cout<<endl;
		//}

		topology.set(adjLtmp.size(),adjLtmp);
		return topology.two_edge_connectivity();
	}
}
bool configuration::checkProposedCompactness(){
	//cout<<"** check compactness module **"<<endl;
	//cout<<endl<<"list_v_to_W_proposed: ";
	//for(auto& iterRef: list_v_to_W_proposed) cout<<iterRef;
	//cout<<endl
	//	<<"list_W_to_vi_proposed: ";
	//for(auto& iterRef: list_W_to_vi_proposed) cout<<iterRef;
	//cout<<endl
	//	<<"list_W_to_vj_proposed: ";
	//for(auto& iterRef: list_W_to_vj_proposed) cout<<iterRef;
	//cout<<endl<<endl;


	//printconfiguration(cout);
	//printProposed(cout);
	//cout<<"removing W-line : "<<*Wremove<<endl;
	//int Nvertices = verticesProposed.size();
	int Nvertices = list_v_to_W_proposed.size();
	//cout<<"Nvertices: "<<Nvertices<<endl<<flush;
	if(Nvertices<3) return true;
	else{
		vector<vector<int> > GadjLtmp, WadjLtmp;

		vector<int> supervertices(Nvertices);

		/*list<vertex>::iterator ref0 = vertices.begin();

		list<line>::iterator iter;

		list<vertex>::iterator vi_iter = Wremove->getconnected_vi();
		list<vertex>::iterator vj_iter = Wremove->getconnected_vj();

		list<line>::iterator Gi_in = vi_iter->getconnected_Gin();
		list<line>::iterator Gj_in = vj_iter->getconnected_Gin();
		list<line>::iterator Gi_out = vi_iter->getconnected_Gout();
		list<line>::iterator Gj_out = vj_iter->getconnected_Gout();

		int iv_remove = distance(ref0,vi_iter);
		int jv_remove = distance(ref0,vj_iter);

		int iv_prev = distance(ref0,Gi_in->getconnected_vi());
		int jv_prev = distance(ref0,Gj_in->getconnected_vi());
		int iv_next = distance(ref0,Gi_out->getconnected_vj());
		int jv_next = distance(ref0,Gj_out->getconnected_vj());*/


		//bool physical;
		//int iv, jv;
		int Nsupervertices;
		int isv_run;

		//cout<<"G topology construct"<<endl;
		three_edge Gtopology;

		isv_run = 0;
		for(int i=0; i<list_W_to_vi_proposed.size(); i++){
			supervertices[list_W_to_vi_proposed[i]] = isv_run;
			supervertices[list_W_to_vj_proposed[i]] = isv_run;
			isv_run++;
		}

		/*for(iter=Ws.begin(); iter!=Ws.end(); ++iter){
			if(iter!=Wremove){
				iv = distance(ref0,iter->getconnected_vi());
				jv = distance(ref0,iter->getconnected_vj());

				supervertices[iv] = isv_run;
				supervertices[jv] = isv_run;
				isv_run++;
			}
		}*/

		Nsupervertices = isv_run;
		//cout<<"Nsupervertices = "<<Nsupervertices<<endl;
		GadjLtmp.resize(Nsupervertices);
		//for(int i=0; i<Nsupervertices; i++) GadjLtmp[i].reserve(Nvertices);
		for(int i=0; i<Nvertices-1; i++){
			//cout<<"supervertices["<<i<<"] = "<<supervertices[i]<<endl;
			if(supervertices[i]!=supervertices[i+1]){
				GadjLtmp[supervertices[i]].push_back(supervertices[i+1]);
				GadjLtmp[supervertices[i+1]].push_back(supervertices[i]);
			}
		}

		/*for(iter=Gs.begin(); iter!=Gs.end(); ++iter){
			iv = distance(ref0,iter->getconnected_vi());

		//cout<<"G push_back"<<endl<<flush;
		iv = supervertices[iv_prev]; jv = supervertices[iv_next];
		if(iv!=jv){
			GadjLtmp[iv].push_back(jv);
			GadjLtmp[jv].push_back(iv);
		}
		iv = supervertices[jv_prev]; jv = supervertices[jv_next];
		if(iv!=jv){
			GadjLtmp[iv].push_back(jv);
			GadjLtmp[jv].push_back(iv);
		}
		for(iter=Gs.begin(); iter!=Gs.end(); ++iter){
			iv = distance(ref0,iter->getconnected_vi());
			jv = distance(ref0,iter->getconnected_vj());
			if((iv!=iv_remove) && (iv!=jv_remove) && (jv!=iv_remove) && (jv!=jv_remove) && (iv!=jv)){
				GadjLtmp[supervertices[iv]].push_back(supervertices[jv]);
				GadjLtmp[supervertices[jv]].push_back(supervertices[iv]);
			}
		}*/

		//cout<<endl<<"GadjLtmp"<<endl;;
		//for(auto i: GadjLtmp){
		//	for(auto j: i) cout<<j;
		//	cout<<endl;
		//}

		Gtopology.set(GadjLtmp.size(),GadjLtmp);
		bool Gcompactness = Gtopology.three_edge_connectivity();

		return Gcompactness;
	}
};

bool configuration::checkWProduct(){
	double WProductBruteForce = 1.0;
	for(auto& iterRef: Ws)
		WProductBruteForce *= iterRef;

	if( fabs( (WProductBruteForce - WProduct)/WProduct )<1.0e-10 ) return true;
	else{
		cout<<"WProductBruteForce = "<<WProductBruteForce<<endl;
		cout<<"WProduct = "<<WProduct<<endl;
		return false;
	}
};
bool configuration::checkGProduct(){
	Eigen::MatrixXd GtTopProductBruteForce = Identity2x2;
	Eigen::MatrixXd GtBottomProductBruteForce = Identity2x2;
	GtTopProductBruteForce *= Gt.getSzGtSz(vertices[im2-2],vertices[im2-1]);
	for(int i=im2-3; i>-1; i--)
		GtTopProductBruteForce *= Gt.getGtSz(vertices[i],vertices[i+1]);

	int Nv = vertices.size();
	GtBottomProductBruteForce *= Gt.getSzGtSz(vertices[Nv-2],vertices[Nv-1]);
	for(int i=Nv-3; i>im2-1; i--)
		GtBottomProductBruteForce *= Gt.getGtSz(vertices[i],vertices[i+1]);


	if( ( GtTopProductBruteForce - GtTopProduct ).norm()<1.0e-10 && ( GtBottomProductBruteForce - GtBottomProduct ).norm()<1.0e-10 ) return true;
	else{
		cout<<"GtTopProductBruteForce: "<<endl<<GtTopProductBruteForce<<endl;
		cout<<"GtTopProduct: "<<endl<<GtTopProduct<<endl<<endl;
		cout<<"GtBottomProductBruteForce: "<<endl<<GtBottomProductBruteForce<<endl;
		cout<<"GtBottomProduct: "<<endl<<GtBottomProduct<<endl;

		return false;
	}
};
bool configuration::checkGtSign(){
	bool running = true;

	for(int i=0; i<GtTop.size()-1; i++){
		//cout<<"Where?"<<endl;
		if( GtTop[i](0,0)<0.0 ) running = false;
		if( GtTop[i](0,1)<0.0 ) running = false;
		if( GtTop[i](1,0)>0.0 ) running = false;
		if( GtTop[i](1,1)>0.0 ) running = false;
	}

	if( GtTop.back()(0,0)<0.0 ) running = false;
	if( GtTop.back()(0,1)<0.0 ) running = false;
	if( GtTop.back()(1,0)<0.0 ) running = false;
	if( GtTop.back()(1,1)<0.0 ) running = false;

	for(int i=0; i<GtBottom.size()-1; i++){
		//cout<<"Here?"<<endl;
		if( GtBottom[i](0,0)<0.0 ) running = false;
		if( GtBottom[i](0,1)<0.0 ) running = false;
		if( GtBottom[i](1,0)>0.0 ) running = false;
		if( GtBottom[i](1,1)>0.0 ) running = false;
	}

	if( GtBottom.back()(0,0)<0.0 ) running = false;
	if( GtBottom.back()(0,1)<0.0 ) running = false;
	if( GtBottom.back()(1,0)<0.0 ) running = false;
	if( GtBottom.back()(1,1)<0.0 ) running = false;

	return running;
};
void configuration::print(ostream& ostr){
	ostr<<endl<<"* configuration *"<<endl;
	ostr<<"extFlavor: "<<endl
		<<extFlavor<<endl;

	ostr<<"vertices: ";
	for(int iv=0; iv<vertices.size(); iv++) ostr<<vertices[iv]<<"-W"<<list_v_to_W[iv]<<" ";
	ostr<<endl;

	ostr<<"im2: "<<im2<<endl;

	ostr<<endl<<"list_v_to_W: ";
	for(auto& iterRef: list_v_to_W) ostr<<iterRef;
	ostr<<endl
		<<"list_W_to_vi: ";
	for(auto& iterRef: list_W_to_vi) ostr<<iterRef;
	ostr<<endl
		<<"list_W_to_vj: ";
	for(auto& iterRef: list_W_to_vj) ostr<<iterRef;
	ostr<<endl<<endl;

	ostr<<"mWt list: ";
	for(int i=0; i<Ws.size(); i++) ostr<<"(v"<<list_W_to_vi[i]<<"-"<<Ws[i]<<"-v"<<list_W_to_vj[i]<<")";
	ostr<<endl;

	ostr<<"GtTop list: "<<endl;
	for(auto& it : GtTop) ostr<<it<<endl<<"---"<<endl;
	ostr<<endl;

	ostr<<"GtBottom list: "<<endl;
	for(auto& it : GtBottom) ostr<<it<<endl<<"---"<<endl;
	ostr<<endl;

	ostr<<"WProduct: "<<endl
		<<WProduct<<endl<<endl;
	ostr<<"GtTopProduct: "<<endl
		<<GtTopProduct<<endl<<endl;
	ostr<<"GtBottomProduct: "<<endl
		<<GtBottomProduct<<endl<<endl;

	ostr<<"MC weight: "<<getWeight()<<endl;

	ostr<<endl<<"*****************"<<endl;
};
void configuration::printProposed(ostream& ostr){
	ostr<<endl<<"* proposed configuration *"<<endl;
	ostr<<"selected vertex: ("<<iv_new<<","<<jv_new<<")"<<endl;
	ostr<<"new v: ("<<vi_new<<","<<vj_new<<")"<<endl;
	ostr<<"new im2: "<<im2_new<<endl;

	ostr<<"selected iW: "<<iW_new<<endl;
	ostr<<"new W: "<<Wt_new<<endl;

	ostr<<"extFlavorProposed: "<<endl
		<<extFlavorProposed<<endl;

	//ostr<<"vertices: ";
	//for(auto& it : vertices) ostr<<it<<" ";
	//ostr<<endl;

	//ostr<<"mWt list: ";
	//for(auto& it : Ws) ostr<<it<<"| ";
	//ostr<<endl;
	ostr<<endl<<"list_v_to_W_proposed: ";
	for(auto& iterRef: list_v_to_W_proposed) ostr<<iterRef;
	ostr<<endl
		<<"list_W_to_vi_proposed: ";
	for(auto& iterRef: list_W_to_vi_proposed) ostr<<iterRef;
	ostr<<endl
		<<"list_W_to_vj_proposed: ";
	for(auto& iterRef: list_W_to_vj_proposed) ostr<<iterRef;
	ostr<<endl<<endl;


	ostr<<"GtTopProposed list: "<<endl;
	for(auto& it : GtTopProposed) ostr<<it<<endl<<"---"<<endl;
	ostr<<endl;

	ostr<<"GtBottomProposed list: "<<endl;
	for(auto& it : GtBottomProposed) ostr<<it<<endl<<"---"<<endl;
	ostr<<endl;

	ostr<<"WProductProposed: "<<endl
		<<WProductProposed<<endl<<endl;
	ostr<<"GtTopProductProposed: "<<endl
		<<GtTopProductProposed<<endl<<endl;
	ostr<<"GtBottomProductProposed: "<<endl
		<<GtBottomProductProposed<<endl<<endl;

	ostr<<"proposed MC weight: "<<getProposedWeight()<<endl;
	ostr<<endl<<"**************************"<<endl;
};
