#include <iostream>
#include <cmath>
#include <bitset>
#include "parameters.hpp"
#include "element.hpp"
using namespace std;
using namespace parameters;

///////////////////////////////
// coordinate implementation //
///////////////////////////////

coordinate::coordinate(){ /*space.resize(Dspace);*/ };
coordinate::coordinate(const double& t_i){ time = t_i; };
void coordinate::settime(const double& t_i){ time = t_i; };
void coordinate::setref(){
	time = 0.0;
	//for(int i=0; i<Dspace; i++)
	//	space(i) = 0;
};
void coordinate::genrand(){
	time = parameters::beta*unidist(mt);
	//for(int i=0; i<Dspace; i++)
	//	space(i) = static_cast<int>(Lx*unidist(mt));
};
/*void coordinate::genrand(const coordinate& x_c){
	time = parameters::beta*unidist(mt);
	for(int i=0; i<Dspace; i++)
		space(i) = (Lx + x_c.space(i) + bidist(mt) - site_generate_radius/2)%Lx;
};*/
void coordinate::genrand(const coordinate& x_i, const coordinate& x_j){
	time = x_i.gettime() + fmod(x_j.gettime() - x_i.gettime() + parameters::beta,parameters::beta)*unidist(mt);
	//space = x_i.getspace();
};
/*double coordinate::getsitegenprobability(const coordinate& x_c){
	double prob = 1.;
	//cout<<"biprobability.size() = "<<biprobability.size()<<endl;
	for(int i=0; i<Dspace; i++){
		prob *= biprobability[min(abs(space(i) - x_c.space(i)),Lx-abs(space(i) - x_c.space(i))) + site_generate_radius/2];
		//cout<<i<<"th index : "<<min(abs(space(i) - x_c.space(i)),Lx-abs(space(i) - x_c.space(i))) + site_generate_radius/2<<endl;
	}
	return prob;
};*/
double coordinate::gettime() const { return time; };
//Eigen::VectorXi coordinate::getspace() const { return space; };
bool coordinate::operator == (const coordinate& x_i){
	const double eps = 1.e-10;
	bool tmp = true;
	tmp = tmp && (abs(time-x_i.time)<eps);
	//for(int i=0; i<space.size(); i++)
	//	tmp = tmp && (space(i)==x_i.space(i));
	return tmp;
};
bool coordinate::operator != (const coordinate& x_i){
	return !((*this)==x_i);
};
ostream& operator << (std::ostream& ostr, const coordinate& x){
	ostr<<"("<<x.time<<")";
	//for(int i=0; i<Dspace; i++) ostr<<x.space(i)<<",";
	//ostr<<")";
	return ostr;
};
bool coordinate::operator < (const coordinate& xi) const {
	return (time < xi.time);
};
coordinate mean(const coordinate& x_i, const coordinate& x_j){
	coordinate x_c;
	x_c.time = (x_i.time + x_j.time)/2.0;
	//for(int i=0; i<Dspace; i++) x_c.space(i) = (x_i.space(i) + x_j.space(i))/2;
	return x_c;
};
coordinate mean(const coordinate& x_1, const coordinate& x_2, const coordinate& x_3, const coordinate& x_4){
	coordinate x_c;
	x_c.time = (x_1.time + x_2.time + x_3.time + x_3.time)/4.0;
	//for(int i=0; i<Dspace; i++) x_c.space(i) = (x_1.space(i) + x_2.space(i) + x_3.space(i) + x_4.space(i))/4;
	return x_c;
};

///////////////////////////
// vertex implementation //
///////////////////////////

/*vertex::vertex(){};
vertex::vertex(unsigned flavor_i, coordinate x_i)
	: flavor(flavor_i), x(x_i){ connected_lines.resize(Nptvertex); };
void vertex::assign(const vertex& v_i){ 
	flavor = v_i.flavor; 
	x = v_i.x; 
};
void vertex::setflavor(const int& flavor_i){ flavor = flavor_i; };
void vertex::setx(const coordinate& x_i){ x = x_i; };
void vertex::setconnected_lines(std::list<line>::iterator W_iter, std::list<line>::iterator Gin_iter, std::list<line>::iterator Gout_iter){
	connected_lines[0] = W_iter;
	connected_lines[1] = Gin_iter;
	connected_lines[2] = Gout_iter;
};
void vertex::setconnected_Gin(std::list<line>::iterator Gin_iter){ connected_lines[1] = Gin_iter; };
void vertex::setconnected_W(std::list<line>::iterator W_iter){ connected_lines[0] = W_iter; };

void vertex::setconnected_loop(std::list<fermionicloop>::iterator loop_iter){ connected_loop = loop_iter; };

unsigned vertex::getflavor(){ return flavor; };
coordinate vertex::getx(){ return x; };
list<line>::iterator vertex::getconnected_W(){ return connected_lines[0]; };
list<line>::iterator vertex::getconnected_Gin(){ return connected_lines[1]; };
list<line>::iterator vertex::getconnected_Gout(){ return connected_lines[2]; };
list<fermionicloop>::iterator vertex::getconnected_loop(){ return connected_loop; };

ostream& operator << (ostream& ostr, const vertex& v){
	//ostr<<"{flavor:"<<v.flavor<<", coordinate:"<<v.x<<"}"<<endl;
	ostr<<"{flavor:"<<v.flavor<<", coordinate:"<<v.x<<"}"<<endl
		<<"connected W-line"<<(*v.connected_lines[0])<<endl
		<<"connected Gin-line"<<(*v.connected_lines[1])<<endl
		<<"connected Gout-line"<<(*v.connected_lines[2])<<endl;
	return ostr;
};*/


/////////////////////////
// line implementation //
/////////////////////////

/*
line::line(bool physical_i, int type_i, unsigned iflavor_i, unsigned jflavor_i, coordinate xi_i, coordinate xj_i)
	: physical(physical_i), type(type_i), iflavor(iflavor_i), jflavor(jflavor_i), xi(xi_i), xj(xj_i){
		if(physical){
			if(type==0){
				valueij = Gt[iflavor].getXvalue(xi,xj);
				valueji = Gt[iflavor].getXvalue(xj,xi);
			}
			else{
				valueij = -Wt[iflavor].getXvalue(xi,xj);
				valueji = -Wt[iflavor].getXvalue(xj,xi);
			}
		}
		else{
			valueij = 1.;
			valueji = 1.;
		}
};
line::line(const line& Line_i){
	physical = Line_i.physical;
	type = Line_i.type;
	//dress = Line_i.dress;
	iflavor = Line_i.iflavor;
	jflavor = Line_i.jflavor;
	xi = Line_i.xi;
	xj = Line_i.xj;
	valueij = Line_i.valueij;
	valueji = Line_i.valueji;
};
void line::assign(const line& line_i){
	physical = line_i.physical;
	type = line_i.type;
	//dress = line_i.dress;
	iflavor = line_i.iflavor;
	jflavor = line_i.jflavor;
	xi = line_i.xi;
	xj = line_i.xj;
	valueij = line_i.valueij;
	valueji = line_i.valueji;
};
void line::setphysical(const bool& physical_i){
       	physical = physical_i; 

	if(physical){
		if(type==0){
			valueij = Gt[jflavor].getXvalue(xi,xj);
			valueji = Gt[iflavor].getXvalue(xj,xi);
		}
		else{
			valueij = -Wt[iflavor].getXvalue(xi,xj);
			valueji = -Wt[iflavor].getXvalue(xj,xi);
		}
	}
	else{
		valueij = 1.;
		valueji = 1.;
	}
};*/
/*void line::setdress(const int& dress_i){
	dress = dress_i;
	if(physical){
		if(type==0){
			valueij = Gt[iflavor].getXvalue(xi,xj);
			valueji = Gt[iflavor].getXvalue(xj,xi);
		}
		else{
			valueij = -Wt[iflavor].getXvalue(xi,xj);
			valueji = -Wt[iflavor].getXvalue(xj,xi);
		}
	}
	else{
		valueij = 1.;
		valueji = 1.;
	}
};*/
/*void line::setflavor(const int& iflavor_i, const int& jflavor_i){ 
	iflavor =iflavor_i;
	jflavor =jflavor_i;

};
void line::setx(const coordinate& xi_i, const coordinate& xj_i){ 
	xi = xi_i; 
	xj = xj_i; 

	if(physical){
		if(type==0){
			valueij = Gt[iflavor].getXvalue(xi,xj);
			valueji = Gt[iflavor].getXvalue(xj,xi);
		}
		else{
			valueij = -Wt[iflavor].getXvalue(xi,xj);
			valueji = -Wt[iflavor].getXvalue(xj,xi);
		}
	}
	else{
		valueij = 1.;
		valueji = 1.;
	}
};
void line::setxi(const coordinate& xi_i){ 
	xi = xi_i; 

	if(physical){
		if(type==0){
			valueij = Gt[iflavor].getXvalue(xi,xj);
			valueji = Gt[iflavor].getXvalue(xj,xi);
		}
		else{
			valueij = -Wt[iflavor].getXvalue(xi,xj);
			valueji = -Wt[iflavor].getXvalue(xj,xi);
		}
	}
	else{
		valueij = 1.;
		valueji = 1.;
	}
};
void line::setxj(const coordinate& xj_i){ 
	xj = xj_i; 

	if(physical){
		if(type==0){
			valueij = Gt[iflavor].getXvalue(xi,xj);
			valueji = Gt[iflavor].getXvalue(xj,xi);
		}
		else{
			valueij = -Wt[iflavor].getXvalue(xi,xj);
			valueji = -Wt[iflavor].getXvalue(xj,xi);
		}
	}
	else{
		valueij = 1.;
		valueji = 1.;
	}
};
void line::setj(const unsigned& jflavor_i, const coordinate& xj_i){
	jflavor = jflavor_i;
	xj = xj_i;

	if(physical){
		if(type==0){
			valueij = Gt[iflavor].getXvalue(xi,xj);
			valueji = Gt[iflavor].getXvalue(xj,xi);
		}
		else{
			valueij = -Wt[iflavor].getXvalue(xi,xj);
			valueji = -Wt[iflavor].getXvalue(xj,xi);
		}
	}
	else{
		valueij = 1.;
		valueji = 1.;
	}
};
void line::setconnected_vertices(std::list<vertex>::iterator vi_iter, std::list<vertex>::iterator vj_iter){
	connected_vertices[0] = vi_iter;
	connected_vertices[1] = vj_iter;
};
void line::setconnected_vj(std::list<vertex>::iterator vj_iter){ connected_vertices[1] = vj_iter; };
//void line::setconnected_loops(std::list<fermionicloop>::iterator loopi_iter,std::list<fermionicloop>::iterator loopj_iter){
//	connected_loops[0] = loopi_iter;
//	connected_loops[1] = loopj_iter;
//};

bool line::getphysical() const { return physical; };
int line::gettype() const { return type; };
//bool line::getisG(){ return isG; };
//unsigned line::getdress() const { return dress; };
unsigned line::getiflavor() const { return iflavor; };
unsigned line::getjflavor() const { return jflavor; };
coordinate line::getxi() const { return xi; };
coordinate line::getxj() const { return xj; };
list<vertex>::iterator line::getconnected_vi(){ return connected_vertices[0]; };
list<vertex>::iterator line::getconnected_vj(){ return connected_vertices[1]; };
//list<fermionicloop>::iterator line::getconnected_loopi(){ return connected_loops[0]; };
//list<fermionicloop>::iterator line::getconnected_loopj(){ return connected_loops[1]; };
double line::getvalueij(){ return valueij; };
double line::getvalueji(){ return valueji; };

ostream& operator << (ostream& ostr, const line& l){
	ostr<<"(physical:"<<l.physical<<", valueij:"<<l.valueij<<", valueji:"<<l.valueji<<");  ";
	//ostr<<"(physical:"<<l.physical<<", dress:"<<l.dress<<", valueij:"<<l.valueij<<", valueji:"<<l.valueji<<");  ";
	if(l.type==0){
		ostr<<"{flavor:"<<l.iflavor<<", coordinate:"<<l.xi<<"} ----- {flavor:"<<l.jflavor<<", coordinate:"<<l.xj<<"}";
		//ostr<<"connected vertices"<<endl <<(*(l.connected_vertices[0])) <<(*(l.connected_vertices[1]))<<endl;
	}
	else{
		ostr<<"{flavor:"<<l.iflavor<<", coordinate:"<<l.xi<<"} = = = {flavor:"<<l.jflavor<<", coordinate:"<<l.xj<<"}";
		//ostr<<"connected vertices"<<endl <<(*l.connected_vertices[0]) <<(*l.connected_vertices[1])<<endl;
	}
	return ostr;
};
*/


///////////////////////////
// bubble implementation //
///////////////////////////

/*bubble::bubble(std::list<line>::iterator Gi_iter, std::list<line>::iterator Gj_iter){
	Gs[0] = Gi_iter;
	Gs[1] = Gj_iter;

	value = -1. * Gi_iter->getvalue() * Gj_iter->getvalue();
};
list<vertex>::iterator bubble::getconnected_vi(){ return Gs[0]->getconnected_vi(); };
list<vertex>::iterator bubble::getconnected_vj(){ return Gs[1]->getconnected_vi(); };
list<line>::iterator bubble::getGi(){ return Gs[0]; };
list<line>::iterator bubble::getGj(){ return Gs[1]; };
double bubble::getvalue(){ return value; };
ostream& operator << (ostream& ostr, const bubble& b){
	//ostr<<"{flavor:"<<v.flavor<<", coordinate:"<<v.x<<"}"<<endl;
	//ostr<<"{flavor:"<<v.flavor<<", coordinate:"<<v.x<<"}"<<endl
	ostr<<"connected Gi-line"<<(*b.Gs[0])<<endl
		<<"connected Gj-line"<<(*b.Gs[1])<<endl;
	return ostr;
};*/



//////////////////////////////////
// fermionicloop implementation //
//////////////////////////////////

/*
fermionicloop::fermionicloop(std::list<vertex>::iterator v_iter){ vs.push_back(v_iter); };
fermionicloop::fermionicloop(std::list<std::list<vertex>::iterator > vs_i)
	: vs(vs_i) {}
void fermionicloop::set(std::list<std::list<vertex>::iterator> vs_i){ vs = vs_i; };
void fermionicloop::setvalue(const double& valuerh_i, const double& valuelh_i){
	valuerh = valuerh_i;
	valuelh = valuelh_i;
}
void fermionicloop::calvalue(){
	double valuerh_tmp = 1.;
	double valuelh_tmp = 1.;

	list<list<vertex>::iterator>::iterator vs_iiter;
	list<line>::iterator G_iter;
	for(vs_iiter=vs.begin(); vs_iiter!=vs.end(); ++vs_iiter){
		G_iter = (*vs_iiter)->getconnected_Gout();
		valuerh_tmp *= G_iter->getvalueij();
		valuelh_tmp *= G_iter->getvalueji();
	}

	valuerh = valuerh_tmp;
	valuelh = valuelh_tmp;
};
void fermionicloop::push_back(const std::list<vertex>::iterator v_iter){ vs.push_back(v_iter); };
void fermionicloop::insert(const std::list<vertex>::iterator& v_prev, const std::list<vertex>::iterator& v_iter){
	list<list<vertex>::iterator>::iterator v_iiter = vs.begin();
	while(*v_iiter!=v_prev){ v_iiter++; };
	vs.insert(++v_iiter,v_iter);
};
void fermionicloop::erase(const std::list<vertex>::iterator& v_iter){
	list<list<vertex>::iterator>::iterator v_iiter = vs.begin();
	while(*v_iiter!=v_iter){ v_iiter++; };
	vs.erase(v_iiter);
};
int fermionicloop::getvssize(){ return vs.size(); };
std::list<std::list<vertex>::iterator> fermionicloop::getvs(){ return vs; };
std::list<std::list<vertex>::iterator> fermionicloop::getpath(const std::list<vertex>::iterator& vi, const std::list<vertex>::iterator& vj){
	// generate the path of range [vi,vj)
	std::list<std::list<vertex>::iterator> path;
	std::list<std::list<vertex>::iterator>::iterator v_iter;
	std::list<std::list<vertex>::iterator>::iterator vi_iter, vj_iter;

	for(v_iter=vs.begin(); v_iter!=vs.end(); ++v_iter){
		if(*v_iter==vi){
			vi_iter = v_iter;
			break;
		}
	}
	v_iter = vi_iter;
	do{
		path.push_back(*(v_iter++));
		if(v_iter==vs.end()) v_iter = vs.begin();
	} while( *v_iter!=vj );

	return path;
};
double fermionicloop::getvaluerh(){ return valuerh; };
double fermionicloop::getvaluelh(){ return valuelh; };
void fermionicloop::proposeloopspinflip(){
	//cout<<"* propose loop spin flip *"<<endl<<flush;
	list<vertex>::iterator vi_iter, vj_iter;
	list<line>::iterator G_iter, UW_iter;

	vs_modified.resize(vs.size());
	//Gs_modified.resize(vs.size());

	UWs_intra.clear(); UWs_intra.reserve(vs.size()/2);
	UWs_inter.clear(); UWs_inter.reserve(vs.size());
	UWs_modified_intra.clear(); UWs_modified_intra.reserve(vs.size()/2);
	UWs_modified_inter.clear(); UWs_modified_inter.reserve(vs.size());
	
	vector<bool> discovered(vs.size(),false);
	vector<bool> isxmodified(vs.size(),false);

	bool isinternal, itoi;
	int UW_type;
	int nflavor, iflavor, jflavor;
	coordinate x_tmp;

	loopflip_proposal = 1.0;
	loopflip_probability = 1.0;
	for(int iv=0; iv<vs.size(); iv++){
		//cout<<iv<<"th vertex"<<endl;
		vi_iter = vs[iv];
		nflavor = (vi_iter->getflavor()^1);
		vs_modified[iv] = *vi_iter;
		vs_modified[iv].setflavor(nflavor);
		if(!discovered[iv]){
			//cout<<"  not discoverd"<<endl;
			discovered[iv] = true;

			UW_iter = vi_iter->getconnected_W();
			UW_type = UW_iter->gettype();
			iflavor = UW_iter->getiflavor();
			jflavor = UW_iter->getjflavor();

			isinternal = false;
			itoi = UW_iter->getconnected_vi()==vi_iter;
			if(UW_type==2 && !(iflavor^jflavor)){
				vj_iter = (itoi ? (UW_iter->getconnected_vj()) : (UW_iter->getconnected_vi()));
				for(int jv=iv; jv<vs.size(); jv++){
					if(vj_iter==vs[jv]){
						discovered[jv] = true;
						isinternal = true;
						break;
					}
				}
				if(isinternal) UWs_intra.push_back(UW_iter);
				else UWs_inter.push_back(UW_iter);
			}
			else{ UWs_inter.push_back(UW_iter); }
			//cout<<"isinternal : "<<isinternal<<endl;
			
			if(UW_type==1){
				//cout<<"connect U-line"<<endl;
				coordinate x_tmp;
				x_tmp.genrand(vi_iter->getx());
				vs_modified[iv].setx(x_tmp);
				isxmodified[iv] = true;
				if(itoi){
					UWs_modified_inter.emplace_back(UW_iter->getphysical(),2,0,nflavor,jflavor,x_tmp,UW_iter->getxj());
				}
				else{
					UWs_modified_inter.emplace_back(UW_iter->getphysical(),2,0,iflavor,nflavor,UW_iter->getxi(),x_tmp);
				}
				loopflip_proposal *= parameters::beta/x_tmp.getsitegenprobability(vi_iter->getx())/2.;
				loopflip_probability *= UWs_modified_inter.back().getvalue()/UW_iter->getvalue();
			}
			else if(iflavor^jflavor){
				UWs_modified_inter.emplace_back(UW_iter->getphysical(),2,0,nflavor,nflavor,UW_iter->getxi(),UW_iter->getxj());
				loopflip_proposal *= 1./2.;
				loopflip_probability *= UWs_modified_inter.back().getvalue()/UW_iter->getvalue();
			}
			else if(!isinternal){
				int UW_type_n = 1+ static_cast<int>(2*unidist(mt));
				if(UW_type_n==1){
					if(itoi){
						x_tmp = UW_iter->getxj();
						vs_modified[iv].setx(x_tmp);
						UWs_modified_inter.emplace_back(UW_iter->getphysical(),UW_type_n,0,nflavor,jflavor,x_tmp,x_tmp);
					}
					else{
						x_tmp = UW_iter->getxi();
						vs_modified[iv].setx(x_tmp);
						UWs_modified_inter.emplace_back(UW_iter->getphysical(),UW_type_n,0,iflavor,nflavor,x_tmp,x_tmp);
					}
					isxmodified[iv] = true;
					loopflip_proposal *= 2./parameters::beta*x_tmp.getsitegenprobability(vi_iter->getx());
				}
				else{
					if(itoi){
						UWs_modified_inter.emplace_back(UW_iter->getphysical(),2,0,nflavor,jflavor,UW_iter->getxi(),UW_iter->getxj());
					}
					else{
						UWs_modified_inter.emplace_back(UW_iter->getphysical(),2,0,iflavor,nflavor,UW_iter->getxi(),UW_iter->getxj());
					}
					loopflip_proposal *= 2.;
				}
				loopflip_probability *= UWs_modified_inter.back().getvalue()/UW_iter->getvalue();
			}
			else{
				UWs_modified_intra.push_back(*UW_iter);
				UWs_modified_intra.back().setflavor(nflavor,nflavor);
				//UWs_modified_intra.emplace_back(UW_iter->getphysical(),2,0,nflavor,nflavor,UW_iter->getxi(),UW_iter->getxj());
			}
		}
	}

	Giters.clear(); Giters.reserve(vs.size());
	Gs_modified.clear(); Gs_modified.reserve(vs.size());
	for(int iv=0; iv<vs.size(); iv++){
		//cout<<"iv : "<<iv<<", (iv+1)%vs.size() : "<<(iv+1)%vs.size()<<endl;
		vi_iter = vs[iv];
		vj_iter = vs[((iv+1)%vs.size())];

		nflavor = (vi_iter->getflavor()^1);

		G_iter = vi_iter->getconnected_Gout();
		Giters.push_back(G_iter);
		Gs_modified.push_back(*G_iter);
		Gs_modified.back().setflavor(nflavor,nflavor);
		if(isxmodified[iv] || isxmodified[((iv+1)%vs.size())])
			Gs_modified.back().setx(vs_modified[iv].getx(),vs_modified[(iv+1)%vs.size()].getx());

		loopflip_probability *= Gs_modified[iv].getvalue()/G_iter->getvalue();
	}

	//cout<<"loopflip_proposal = "<<loopflip_proposal<<endl;
	loopflip_probability *= loopflip_proposal;
};*/
/*void fermionicloop::proposeloopspinflip(const vector<coordinate>& xs, const vector<int>& UW_types){
	//cout<<"* propose loop spin flip *"<<endl<<flush;
	list<vertex>::iterator vi_iter, vj_iter;
	list<line>::iterator G_iter, UW_iter;

	vs_modified.resize(vs.size());
	//Gs_modified.resize(vs.size());

	UWs_intra.clear();
	UWs_inter.clear();
	UWs_modified_intra.clear();
	UWs_modified_inter.clear();
	
	std::vector<bool> discovered(vs.size(),false);

	bool isinternal, itoi;
	int UW_type;
	int nflavor, iflavor, jflavor;
	coordinate x_tmp;

	loopflip_proposal = 1.0;
	loopflip_probability = 1.0;
	for(int iv=0; iv<vs.size(); iv++){
		//cout<<iv<<"th vertex"<<endl;
		vi_iter = vs[iv];
		nflavor = (vi_iter->getflavor()^1);
		vs_modified[iv] = *vi_iter;
		vs_modified[iv].setflavor(nflavor);
		if(!discovered[iv]){
			//cout<<"  not discoverd"<<endl;
			discovered[iv] = true;

			UW_iter = vi_iter->getconnected_W();
			UW_type = UW_iter->gettype();
			iflavor = UW_iter->getiflavor();
			jflavor = UW_iter->getjflavor();

			isinternal = false;
			itoi = UW_iter->getconnected_vi()==vi_iter;
			if(UW_type==2 && !(iflavor^jflavor)){
				vj_iter = (itoi ? (UW_iter->getconnected_vj()) : (UW_iter->getconnected_vi()));
				for(int jv=iv; jv<vs.size(); jv++){
					if(vj_iter==vs[jv]){
						discovered[jv] = true;
						isinternal = true;
						break;
					}
				}
				if(isinternal) UWs_intra.push_back(UW_iter);
				else UWs_inter.push_back(UW_iter);
			}
			else{ UWs_inter.push_back(UW_iter); }
			//cout<<"isinternal : "<<isinternal<<endl;
			
			if(UW_type==1){
				//cout<<"connect U-line"<<endl;
				coordinate x_tmp = xs[iv];
				//x_tmp.genrand(vi_iter->getx());
				vs_modified[iv].setx(x_tmp);
				if(itoi){
					UWs_modified_inter.emplace_back(UW_iter->getphysical(),2,0,nflavor,jflavor,x_tmp,UW_iter->getxj());
				}
				else{
					UWs_modified_inter.emplace_back(UW_iter->getphysical(),2,0,iflavor,nflavor,UW_iter->getxi(),x_tmp);
				}
				loopflip_proposal *= parameters::beta/x_tmp.getsitegenprobability(vi_iter->getx())/2.;
				loopflip_probability *= UWs_modified_inter.back().getvalue()/UW_iter->getvalue();
			}
			else if(iflavor^jflavor){
				UWs_modified_inter.emplace_back(UW_iter->getphysical(),2,0,nflavor,nflavor,UW_iter->getxi(),UW_iter->getxj());
				loopflip_proposal *= 1./2.;
				loopflip_probability *= UWs_modified_inter.back().getvalue()/UW_iter->getvalue();
			}
			else if(!isinternal){
				int UW_type_n = UW_types[iv];
				//int UW_type_n = 1+ static_cast<int>(2*unidist(mt));
				if(UW_type_n==1){
					if(itoi){
						x_tmp = xs[iv];
						//x_tmp = UW_iter->getxj();
						vs_modified[iv].setx(x_tmp);
						UWs_modified_inter.emplace_back(UW_iter->getphysical(),UW_type_n,0,nflavor,jflavor,x_tmp,x_tmp);
					}
					else{
						x_tmp = xs[iv];
						//x_tmp = UW_iter->getxi();
						vs_modified[iv].setx(x_tmp);
						UWs_modified_inter.emplace_back(UW_iter->getphysical(),UW_type_n,0,iflavor,nflavor,x_tmp,x_tmp);
					}
					loopflip_proposal *= 2./parameters::beta*x_tmp.getsitegenprobability(vi_iter->getx());
				}
				else{
					if(itoi){
						UWs_modified_inter.emplace_back(UW_iter->getphysical(),2,0,nflavor,jflavor,UW_iter->getxi(),UW_iter->getxj());
					}
					else{
						UWs_modified_inter.emplace_back(UW_iter->getphysical(),2,0,iflavor,nflavor,UW_iter->getxi(),UW_iter->getxj());
					}
					loopflip_proposal *= 2.;
				}
				loopflip_probability *= UWs_modified_inter.back().getvalue()/UW_iter->getvalue();
			}
			else{
				UWs_modified_intra.emplace_back(UW_iter->getphysical(),2,0,nflavor,nflavor,UW_iter->getxi(),UW_iter->getxj());
			}
		}
	}

	Giters.clear();
	Gs_modified.clear();
	for(int iv=0; iv<vs.size(); iv++){
		//cout<<"iv : "<<iv<<", (iv+1)%vs.size() : "<<(iv+1)%vs.size()<<endl;
		vi_iter = vs[iv];
		vj_iter = vs[((iv+1)%vs.size())];

		nflavor = (vi_iter->getflavor()^1);

		G_iter = vi_iter->getconnected_Gout();
		Giters.push_back(G_iter);
		Gs_modified.push_back(*G_iter);
		Gs_modified.back().setflavor(nflavor,nflavor);
		Gs_modified.back().setx(vs_modified[iv].getx(),vs_modified[(iv+1)%vs.size()].getx());

		loopflip_probability *= Gs_modified[iv].getvalue()/G_iter->getvalue();
	}

	//cout<<"loopflip_proposal = "<<loopflip_proposal<<endl;
	loopflip_probability *= loopflip_proposal;
};*/
/*vector<coordinate> fermionicloop::getxvector(){
	vector<coordinate> xs;
	for(auto vi: vs){
		xs.emplace_back(vi->getx());
	}
	return xs;
};*/
/*vector<int> fermionicloop::getUWtypes(){
	vector<int> UWtypes;
	list<line>::iterator UWiter;
	for(auto vi: vs){
		UWiter = vi->getconnected_W();
		UWtypes.emplace_back(UWiter->gettype());
	}
	return UWtypes;
};*/
//double fermionicloop::get_loopspinflip_probability(){ return loopflip_probability; };
/*void fermionicloop::loopspinflip(){
	for(int i=0; i<vs.size(); i++)
		vs[i]->assign(vs_modified[i]);
	for(int i=0; i<UWs_intra.size(); i++)
		UWs_intra[i]->assign(UWs_modified_intra[i]);
	for(int i=0; i<UWs_inter.size(); i++)
		UWs_inter[i]->assign(UWs_modified_inter[i]);
	for(int i=0; i<Giters.size(); i++)
		Giters[i]->assign(Gs_modified[i]);
};
ostream& operator << (ostream& ostr, const fermionicloop& fl){
	ostr<<"-- (valuerh: "<<fl.valuerh<<")"<<endl<<"|"<<endl;
	for(auto vi: fl.vs){ ostr<<*vi<<"|"<<endl; }
	ostr<<"-- (valuelh: "<<fl.valuelh<<")"<<endl;
	return ostr;
};*/



///////////////////////////////
// setUWinter implementation //
///////////////////////////////

/*setUWinter::setUWinter(){
	int Nreserve = 2*(Nordermax+1);

	isloopflavorfixed.reserve(Nreserve);
	loopflavor.reserve(Nreserve);
	UWinters.reserve(Nreserve);
	UWvalues.reserve(Nreserve);
	iloops.reserve(Nreserve);
	iloops.reserve(Nreserve);
};
setUWinter::setUWinter(const int& Nfermionicloops_i) 
	: Nfermionicloops(Nfermionicloops_i) {
	int Nreserve = 2*(Nordermax+1);

	isloopflavorfixed.reserve(Nreserve);
	loopflavor.reserve(Nreserve);
	UWinters.reserve(Nreserve);
	UWvalues.reserve(Nreserve);
	iloops.reserve(Nreserve);
	iloops.reserve(Nreserve);

	isloopflavorfixed.resize(Nfermionicloops);	
	loopflavor.resize(Nfermionicloops);	
};
setUWinter::setUWinter(const setUWinter& copy_i){
	int Nreserve = 2*(Nordermax+1);

	isloopflavorfixed.reserve(Nreserve);
	loopflavor.reserve(Nreserve);
	UWinters.reserve(Nreserve);
	UWvalues.reserve(Nreserve);
	iloops.reserve(Nreserve);
	iloops.reserve(Nreserve);

	this->Nfermionicloops = copy_i.Nfermionicloops;
	this->isloopflavorfixed = copy_i.isloopflavorfixed;
	this->loopflavor = copy_i.loopflavor;
	this->UWinters = copy_i.UWinters;
	this->iloops = copy_i.iloops;
	this->jloops = copy_i.jloops;
};

void setUWinter::setNfermionicloops(const int& Nfermionicloops_i){
	Nfermionicloops = Nfermionicloops_i;
	isloopflavorfixed.resize(Nfermionicloops);	
	loopflavor.resize(Nfermionicloops);	
};
void setUWinter::clear(){
	Nfermionicloops = 1;
	isloopflavorfixed.clear();
	loopflavor.clear();
	UWinters.clear();
	UWvalues.clear();
	iloops.clear();
	jloops.clear();
};
void setUWinter::push_back(const line& UWinter_i, const int& iloop, const int& jloop){
	UWinters.push_back(UWinter_i);
	iloops.push_back(iloop);
	jloops.push_back(jloop);

	Eigen::MatrixXd UWvalue(Nflavor,Nflavor);
	if(UWinter_i.getphysical()){
		if( UWinter_i.gettype()==1 ){
			UWvalue(0,0) = 0.;
			UWvalue(1,0) = -U;
			UWvalue(0,1) = -U;
			UWvalue(1,1) = 0.;
		}
		else{
			UWvalue(0,0) = -Wt[0][0][0].getXvalue(UWinter_i.getxi(),UWinter_i.getxj());
			UWvalue(1,0) = -Wt[0][0][1].getXvalue(UWinter_i.getxi(),UWinter_i.getxj());
			UWvalue(0,1) = UWvalue(1,0);
			UWvalue(1,1) = UWvalue(0,0);
		}
	}
	else{
		UWvalue(0,0) = 1.;
		UWvalue(1,0) = 1.;
		UWvalue(0,1) = 1.;
		UWvalue(1,1) = 1.;
	}
	UWvalues.push_back(UWvalue);
};
void setUWinter::releaseloopflavor(){
	for(int i=0; i<isloopflavorfixed.size(); i++){ isloopflavorfixed[i] = false; }
};
void setUWinter::fixloopflavor(const int& iloop, const int& flavor_i){
	isloopflavorfixed[iloop] = true;
	loopflavor[iloop] = flavor_i;
};
void setUWinter::setisloopflavorfixed(const std::vector<bool>& isloopflavorfixed_n){ isloopflavorfixed = isloopflavorfixed_n; };
void setUWinter::setloopflavor(const std::vector<int>& loopflavor_n){ loopflavor = loopflavor_n; };
bool setUWinter::getisloopflavorfixed(const int& iloop){
	return isloopflavorfixed[iloop];
};
int setUWinter::getloopflavor(const int& iloop){
	return loopflavor[iloop];
};
void setUWinter::setphysical(const line& line_i, const bool& physical_i){
	int iUWinter;
	vector<line>::iterator iter;
	for(iter=UWinters.begin(); iter!=UWinters.end(); ++iter){
		if(
			iter->gettype()==line_i.gettype()
			&& iter->getxi()==line_i.getxi()
			&& iter->getxj()==line_i.getxj()
		){
			iter->setphysical(physical_i); 
			iUWinter = distance(UWinters.begin(),iter);
			break;
		};
	};

	if(physical_i){
		if( UWinters[iUWinter].gettype()==1 ){
			UWvalues[iUWinter](0,0) = 0.;
			UWvalues[iUWinter](1,0) = -U;
			UWvalues[iUWinter](0,1) = -U;
			UWvalues[iUWinter](1,1) = 0.;
		}
		else{
			UWvalues[iUWinter](0,0) = -Wt[0][0][0].getXvalue(UWinters[iUWinter].getxi(),UWinters[iUWinter].getxj());
			UWvalues[iUWinter](1,0) = -Wt[0][0][1].getXvalue(UWinters[iUWinter].getxi(),UWinters[iUWinter].getxj());
			UWvalues[iUWinter](0,1) = UWvalues[iUWinter](1,0);
			UWvalues[iUWinter](1,1) = UWvalues[iUWinter](0,0);
		}
	}
	else{
		UWvalues[iUWinter](0,0) = 1.;
		UWvalues[iUWinter](1,0) = 1.;
		UWvalues[iUWinter](0,1) = 1.;
		UWvalues[iUWinter](1,1) = 1.;
	}
};
void setUWinter::erase(const line& line_i){
	int iUWinter;
	vector<line>::iterator iter;
	for(iter=UWinters.begin(); iter!=UWinters.end(); ++iter){
		if(
			iter->gettype()==line_i.gettype()
			&& iter->getxi()==line_i.getxi()
			&& iter->getxj()==line_i.getxj()
		){
			iUWinter = distance(UWinters.begin(),iter);
			break;
		};
	};
	UWinters.erase(next(UWinters.begin(),iUWinter));
	UWvalues.erase(next(UWvalues.begin(),iUWinter));
	iloops.erase(next(iloops.begin(),iUWinter));
	jloops.erase(next(jloops.begin(),iUWinter));
};
void setUWinter::transform(const line& line_o, const line& line_n){
	int iUWinter;
	vector<line>::iterator iter;
	for(iter=UWinters.begin(); iter!=UWinters.end(); ++iter){
		if(
			iter->gettype()==line_o.gettype()
			&& iter->getxi()==line_o.getxi()
			&& iter->getxj()==line_o.getxj()
		){
			*iter = line_n;
			iUWinter = distance(UWinters.begin(),iter);
			break;
		};
	};
	if(UWinters[iUWinter].getphysical()){
		if( UWinters[iUWinter].gettype()==1 ){
			UWvalues[iUWinter](0,0) = 0.;
			UWvalues[iUWinter](1,0) = -U;
			UWvalues[iUWinter](0,1) = -U;
			UWvalues[iUWinter](1,1) = 0.;
		}
		else{
			UWvalues[iUWinter](0,0) = -Wt[0][0][0].getXvalue(UWinters[iUWinter].getxi(),UWinters[iUWinter].getxj());
			UWvalues[iUWinter](1,0) = -Wt[0][0][1].getXvalue(UWinters[iUWinter].getxi(),UWinters[iUWinter].getxj());
			UWvalues[iUWinter](0,1) = UWvalues[iUWinter](1,0);
			UWvalues[iUWinter](1,1) = UWvalues[iUWinter](0,0);
		}
	}
	else{
		UWvalues[iUWinter](0,0) = 1.;
		UWvalues[iUWinter](1,0) = 1.;
		UWvalues[iUWinter](0,1) = 1.;
		UWvalues[iUWinter](1,1) = 1.;
	}
};
void setUWinter::calvalue(){
	//cout<<"** setUWinter::calvalue() **"<<endl;
	//cout<<"   terms: ";
	if(Nfermionicloops==1) (this->value) = 1;
	else{
		int power = 0;
		for(int iloop=0; iloop<isloopflavorfixed.size(); iloop++) if(!isloopflavorfixed[iloop]) power++;

		vector<bool>::iterator iter;
		int flavorset, irunning;
		double term;
		int iflavor, jflavor;

		(this->value) = 0.;
		for(int iconf=0; iconf<pow(2,power); iconf++){
			//cout<<"iconf = "<<iconf<<endl<<flush;
			flavorset = iconf;
			irunning = 0;
			for(int iloop=0; iloop<loopflavor.size(); iloop++){
				if(!isloopflavorfixed[iloop]){
					loopflavor[iloop] = static_cast<int>((flavorset&(1<<irunning))!=0);
					irunning++;
				}
			}

			term = 1.;
			for(int iUW=0; iUW<UWinters.size(); iUW++){
				iflavor = loopflavor[iloops[iUW]];
				jflavor = loopflavor[jloops[iUW]];

				if( abs(UWvalues[iUW](iflavor,jflavor))>1.0e-10 ) 
					term *= UWvalues[iUW](iflavor,jflavor);
				else{
					term = 0.;
					break;
				}
			}
			(this->value) += term;
			//cout<<term<<" + ";
		}
	}
	//cout<<endl;
	//cout<<"****************************"<<endl;
};
double setUWinter::calvalueafteraddvertices(const line& line_n, const int& iloop_n, const int& jloop_n){
	//cout<<"** setUWinter::calvalueafteraddvertices() **"<<endl<<flush;
	double valtmp;
	if(Nfermionicloops==1) valtmp = 1;
	else{
		Eigen::MatrixXd UWvalue_add(Nflavor,Nflavor);
		if(line_n.getphysical()){
			if( line_n.gettype()==1 ){
				UWvalue_add(0,0) = 0.;
				UWvalue_add(1,0) = -U;
				UWvalue_add(0,1) = -U;
				UWvalue_add(1,1) = 0.;
			}
			else{
				UWvalue_add(0,0) = -Wt[0][0][0].getXvalue(line_n.getxi(),line_n.getxj());
				UWvalue_add(0,1) = -Wt[0][0][1].getXvalue(line_n.getxi(),line_n.getxj());
				UWvalue_add(1,0) = UWvalue_add(0,1);
				UWvalue_add(1,1) = UWvalue_add(0,0);
			}
		}
		else{ 
			UWvalue_add(0,0) = 1.;
			UWvalue_add(1,0) = 1.;
			UWvalue_add(0,1) = 1.;
			UWvalue_add(1,1) = 1.;
		}

		int power = 0;
		for(int iloop=0; iloop<isloopflavorfixed.size(); iloop++) if(!isloopflavorfixed[iloop]) power++;

		vector<bool>::iterator iter;
		int flavorset, irunning;
		double term;
		int iflavor, jflavor;

		valtmp = 0.;
		for(int iconf=0; iconf<pow(2,power); iconf++){
			flavorset = iconf;
			irunning = 0;
			for(int iloop=0; iloop<loopflavor.size(); iloop++){
				if(!isloopflavorfixed[iloop]){
					loopflavor[iloop] = static_cast<int>((flavorset&(1<<irunning))!=0);
					irunning++;
				}
			}

			term = 1.;

			iflavor = loopflavor[iloop_n];
			jflavor = loopflavor[jloop_n];
			term *= UWvalue_add(iflavor,jflavor);

			for(int iUW=0; iUW<UWinters.size(); iUW++){
				iflavor = loopflavor[iloops[iUW]];
				jflavor = loopflavor[jloops[iUW]];

				if( abs(UWvalues[iUW](iflavor,jflavor))>1.0e-10 ) 
					term *= UWvalues[iUW](iflavor,jflavor);
				else{
					term = 0.;
					break;
				}
			}

			valtmp += term;
		}
	}

	//cout<<"   value = "<<valtmp<<endl;
	//cout<<"*******************************************"<<endl<<flush;

	return valtmp;
};
double setUWinter::calvalueafterremovevertices(const line& line_o){
	//cout<<"** setUWinter::calvalueafterremovevertices() **"<<endl<<flush;
	double valtmp;
	if(Nfermionicloops==1) valtmp = 1;
	else{
		int iremove;
		bool line_found = false;
		vector<line>::iterator UWiter;
		for(UWiter=UWinters.begin(); UWiter!=UWinters.end(); ++UWiter){
			if(
				UWiter->gettype()==line_o.gettype()
				&& UWiter->getxi()==line_o.getxi()
				&& UWiter->getxj()==line_o.getxj()
			){
				iremove = distance(UWinters.begin(),UWiter);
				line_found = true;
				break;
			}
		}
		if(!line_found){
			cout<<"   ALERT! TARGET LINE IS NOT FOUND!"<<endl;
			cout<<"   target line: "<<line_o<<endl;
			cout<<*this<<endl;
			exit(EXIT_FAILURE);
		}

		int power = 0;
		for(int iloop=0; iloop<isloopflavorfixed.size(); iloop++) if(!isloopflavorfixed[iloop]) power++;

		vector<bool>::iterator iter;
		int flavorset, irunning;
		double term;
		int iflavor, jflavor;

		valtmp = 0.;
		for(int iconf=0; iconf<pow(2,power); iconf++){
			flavorset = iconf;
			irunning = 0;
			for(int iloop=0; iloop<loopflavor.size(); iloop++){
				if(!isloopflavorfixed[iloop]){
					loopflavor[iloop] = static_cast<int>((flavorset&(1<<irunning))!=0);
					irunning++;
				}
			}

			term = 1.;
			for(int iUW=0; iUW<UWinters.size(); iUW++){
				iflavor = loopflavor[iloops[iUW]];
				jflavor = loopflavor[jloops[iUW]];

				if(iUW!=iremove){
					if( abs(UWvalues[iUW](iflavor,jflavor))>1.0e-10 ) 
						term *= UWvalues[iUW](iflavor,jflavor);
					else{
						term = 0.;
						break;
					}
				}
			}

			valtmp += term;
		}
	}

	//cout<<"   value = "<<valtmp<<endl;
	//cout<<"*******************************************"<<endl<<flush;

	return valtmp;
};
double setUWinter::calvalueafterinterswap(const line& line_o, const std::vector<bool>& isloopflavorfixed_n, std::vector<int> loopflavor_n){
	//cout<<"** setUWinter::calvalueafterremovevertices() **"<<endl<<flush;
	double valtmp;
	if(Nfermionicloops==1) valtmp = 1;
	else{
		int iinterswap;
		bool line_found = false;
		vector<line>::iterator UWiter;
		Eigen::MatrixXd UWvalue_interswap(Nflavor,Nflavor);
		for(UWiter=UWinters.begin(); UWiter!=UWinters.end(); ++UWiter){
			if(
				UWiter->gettype()==line_o.gettype()
				&& UWiter->getxi()==line_o.getxi()
				&& UWiter->getxj()==line_o.getxj()
			){
				iinterswap = distance(UWinters.begin(),UWiter);
				line_found = true;
				if(!UWiter->getphysical()){
					if( UWiter->gettype()==1 ){
						UWvalue_interswap(0,0) = 0.;
						UWvalue_interswap(1,0) = -U;
						UWvalue_interswap(0,1) = -U;
						UWvalue_interswap(1,1) = 0.;
					}
					else{
						UWvalue_interswap(0,0) = -Wt[0][0][0].getXvalue(UWiter->getxi(),UWiter->getxj());
						UWvalue_interswap(0,1) = -Wt[0][0][1].getXvalue(UWiter->getxi(),UWiter->getxj());
						UWvalue_interswap(1,0) = UWvalue_interswap(0,1);
						UWvalue_interswap(1,1) = UWvalue_interswap(0,0);
					}
				}
				else{ 
					UWvalue_interswap(0,0) = 1.;
					UWvalue_interswap(1,0) = 1.;
					UWvalue_interswap(0,1) = 1.;
					UWvalue_interswap(1,1) = 1.;
				}
				break;
			}
		}
		if(!line_found){
			cout<<"   ALERT! TARGET LINE IS NOT FOUND!"<<endl;
			cout<<"   target line: "<<line_o<<endl;
			cout<<*this<<endl;
			exit(EXIT_FAILURE);
		}

		int power = 0;
		for(int iloop=0; iloop<isloopflavorfixed_n.size(); iloop++) if(!isloopflavorfixed_n[iloop]) power++;

		vector<bool>::iterator iter;
		int flavorset, irunning;
		double term;
		int iflavor, jflavor;

		valtmp = 0.;
		for(int iconf=0; iconf<pow(2,power); iconf++){
			flavorset = iconf;
			irunning = 0;
			for(int iloop=0; iloop<loopflavor_n.size(); iloop++){
				if(!isloopflavorfixed_n[iloop]){
					loopflavor_n[iloop] = static_cast<int>((flavorset&(1<<irunning))!=0);
					irunning++;
				}
			}

			term = 1.;
			for(int iUW=0; iUW<UWinters.size(); iUW++){
				iflavor = loopflavor_n[iloops[iUW]];
				jflavor = loopflavor_n[jloops[iUW]];

				if(iUW==iinterswap){
					if( abs(UWvalue_interswap(iflavor,jflavor))>1.0e-10 )
						term *= UWvalue_interswap(iflavor,jflavor);
					else{
						term = 0.;
						break;
					}
				}
				else{
					if( abs(UWvalues[iUW](iflavor,jflavor))>1.0e-10 ) 
						term *= UWvalues[iUW](iflavor,jflavor);
					else{
						term = 0.;
						break;
					}
				}
			}

			valtmp += term;
		}
	}

	//cout<<"   value = "<<valtmp<<endl;
	//cout<<"*******************************************"<<endl<<flush;

	return valtmp;
};
double setUWinter::calvalueafterinterswap(const std::vector<bool>& isloopflavorfixed_n, std::vector<int> loopflavor_n){
	//cout<<"** setUWinter::calvalueafterremovevertices() **"<<endl<<flush;
	double valtmp;
	if(Nfermionicloops==1) valtmp = 1;
	else{
		int power = 0;
		for(int iloop=0; iloop<isloopflavorfixed_n.size(); iloop++) if(!isloopflavorfixed_n[iloop]) power++;

		vector<bool>::iterator iter;
		int flavorset, irunning;
		double term;
		int iflavor, jflavor;

		valtmp = 0.;
		for(int iconf=0; iconf<pow(2,power); iconf++){
			flavorset = iconf;
			irunning = 0;
			for(int iloop=0; iloop<loopflavor_n.size(); iloop++){
				if(!isloopflavorfixed_n[iloop]){
					loopflavor_n[iloop] = static_cast<int>((flavorset&(1<<irunning))!=0);
					irunning++;
				}
			}

			term = 1.;
			for(int iUW=0; iUW<UWinters.size(); iUW++){
				iflavor = loopflavor_n[iloops[iUW]];
				jflavor = loopflavor_n[jloops[iUW]];

				if( abs(UWvalues[iUW](iflavor,jflavor))>1.0e-10 ) 
					term *= UWvalues[iUW](iflavor,jflavor);
				else{
					term = 0.;
					break;
				}
			}

			valtmp += term;
		}
	}

	//cout<<"   value = "<<valtmp<<endl;
	//cout<<"*******************************************"<<endl<<flush;

	return valtmp;
};
double setUWinter::calvalueafterintraswap(const line& line_o1, const line& line_o2, const std::vector<bool>& isloopflavorfixed_n, std::vector<int> loopflavor_n){
	//cout<<"** setUWinter::calvalueafterremovevertices() **"<<endl<<flush;
	double valtmp;
	if(Nfermionicloops==1) valtmp = 1;
	else{
		int iinterswap;
		vector<line>::iterator UWiter;
		Eigen::MatrixXd UWvalue_interswap(Nflavor,Nflavor);
		for(UWiter=UWinters.begin(); UWiter!=UWinters.end(); ++UWiter){
			if(
				(UWiter->gettype()==line_o1.gettype()
				&& UWiter->getxi()==line_o1.getxi()
				&& UWiter->getxj()==line_o1.getxj())
				||
				(UWiter->gettype()==line_o2.gettype()
				&& UWiter->getxi()==line_o2.getxi()
				&& UWiter->getxj()==line_o2.getxj())
			){
				iinterswap = distance(UWinters.begin(),UWiter);
				if(!UWiter->getphysical()){
					if( UWiter->gettype()==1 ){
						UWvalue_interswap(0,0) = 0.;
						UWvalue_interswap(1,0) = -U;
						UWvalue_interswap(0,1) = -U;
						UWvalue_interswap(1,1) = 0.;
					}
					else{
						UWvalue_interswap(0,0) = -Wt[0][0][0].getXvalue(UWiter->getxi(),UWiter->getxj());
						UWvalue_interswap(0,1) = -Wt[0][0][1].getXvalue(UWiter->getxi(),UWiter->getxj());
						UWvalue_interswap(1,0) = UWvalue_interswap(0,1);
						UWvalue_interswap(1,1) = UWvalue_interswap(0,0);
					}
				}
				else{ 
					UWvalue_interswap(0,0) = 1.;
					UWvalue_interswap(1,0) = 1.;
					UWvalue_interswap(0,1) = 1.;
					UWvalue_interswap(1,1) = 1.;
				}
			}
		}

		int power = 0;
		for(int iloop=0; iloop<isloopflavorfixed_n.size(); iloop++) if(!isloopflavorfixed_n[iloop]) power++;

		vector<bool>::iterator iter;
		int flavorset, irunning;
		double term;
		int iflavor, jflavor;

		valtmp = 0.;
		for(int iconf=0; iconf<pow(2,power); iconf++){
			flavorset = iconf;
			irunning = 0;
			for(int iloop=0; iloop<loopflavor_n.size(); iloop++){
				if(!isloopflavorfixed_n[iloop]){
					loopflavor_n[iloop] = static_cast<int>((flavorset&(1<<irunning))!=0);
					irunning++;
				}
			}

			term = 1.;
			for(int iUW=0; iUW<UWinters.size(); iUW++){
				iflavor = loopflavor_n[iloops[iUW]];
				jflavor = loopflavor_n[jloops[iUW]];

				if(iUW==iinterswap){
					if( abs(UWvalue_interswap(iflavor,jflavor))>1.0e-10 )
						term *= UWvalue_interswap(iflavor,jflavor);
					else{
						term = 0.;
						break;
					}
				}
				else{
					if( abs(UWvalues[iUW](iflavor,jflavor))>1.0e-10 ) 
						term *= UWvalues[iUW](iflavor,jflavor);
					else{
						term = 0.;
						break;
					}
				}
			}

			valtmp += term;
		}
	}

	//cout<<"   value = "<<valtmp<<endl;
	//cout<<"*******************************************"<<endl<<flush;

	return valtmp;
};

double setUWinter::calvalueaftertransformUW(const line& line_o, const line& line_n){
	//cout<<"** setUWinter::calvalueaftertransformUW() **"<<endl<<flush;
	double valtmp;
	if(Nfermionicloops==1) valtmp = 1;
	else{
		int itransform;
		bool line_found = false;
		vector<line>::iterator UWiter;
		Eigen::MatrixXd UWvalue_trans(Nflavor,Nflavor);
		for(UWiter=UWinters.begin(); UWiter!=UWinters.end(); ++UWiter){
			int iUW = distance(UWinters.begin(),UWiter);
			if(
				UWiter->gettype()==line_o.gettype()
				&& UWiter->getxi()==line_o.getxi()
				&& UWiter->getxj()==line_o.getxj()
			){
				itransform = distance(UWinters.begin(),UWiter);
				line_found = true;

				if(line_n.getphysical()){
					if( line_n.gettype()==1 ){
						UWvalue_trans(0,0) = 0.;
						UWvalue_trans(1,0) = -U;
						UWvalue_trans(0,1) = -U;
						UWvalue_trans(1,1) = 0.;
					}
					else{
						UWvalue_trans(0,0) = -Wt[0][0][0].getXvalue(line_n.getxi(),line_n.getxj());
						UWvalue_trans(0,1) = -Wt[0][0][1].getXvalue(line_n.getxi(),line_n.getxj());
						UWvalue_trans(1,0) = UWvalue_trans(0,1);
						UWvalue_trans(1,1) = UWvalue_trans(0,0);
					}
				}
				else{ 
					UWvalue_trans(0,0) = 1.;
					UWvalue_trans(1,0) = 1.;
					UWvalue_trans(0,1) = 1.;
					UWvalue_trans(1,1) = 1.;
				}

				break;
			}
		}
		if(!line_found){
			cout<<"   ALERT! TARGET LINE IS NOT FOUND!"<<endl;
			cout<<"   target line: "<<line_o<<endl;
			cout<<*this<<endl;
			exit(EXIT_FAILURE);
		}
		//cout<<"   itransform = "<<itransform<<endl<<flush;
		//cout<<"   terms: ";

		int power = 0;
		for(int iloop=0; iloop<isloopflavorfixed.size(); iloop++) if(!isloopflavorfixed[iloop]) power++;

		//cout<<"power = "<<power<<endl<<flush;
		//cout<<"pow(2,power) = "<<pow(2,power)<<endl<<flush;
		//cout<<"UWinters.size() = "<<UWinters.size()<<endl<<flush;
		//cout<<"iloops.size() = "<<iloops.size()<<endl<<flush;
		//cout<<"jloops.size() = "<<jloops.size()<<endl<<flush;

		vector<bool>::iterator iter;
		int flavorset, irunning;
		double term;
		int iflavor, jflavor;

		valtmp = 0.;
		for(int iconf=0; iconf<pow(2,power); iconf++){
			//cout<<"iconf = "<<iconf<<endl<<flush;
			flavorset = iconf;
			irunning = 0;
			for(int iloop=0; iloop<loopflavor.size(); iloop++){
				if(!isloopflavorfixed[iloop]){
					loopflavor[iloop] = static_cast<int>((flavorset&(1<<irunning))!=0);
					irunning++;
				}
			}

			//for(int iloop=0; iloop<loopflavor.size(); iloop++) cout<<loopflavor[iloop]<<flush;
			//cout<<endl<<"point"<<endl<<flush;

			term = 1.;
			for(int iUW=0; iUW<UWinters.size(); iUW++){
				iflavor = loopflavor[iloops[iUW]];
				jflavor = loopflavor[jloops[iUW]];

				//cout<<iUW<<flush;
				//cout<<UWinters[iUW]<<endl;
				//cout<<"("<<iflavor<<","<<jflavor<<")"<<endl<<flush;

				if(iUW==itransform){
					if( abs(UWvalue_trans(iflavor,jflavor))>1.0e-10 ) 
						term *= UWvalue_trans(iflavor,jflavor);
					else{
						term = 0.;
						break;
					}
				}
				else{
					if( abs(UWvalues[iUW](iflavor,jflavor))>1.0e-10 ) 
						term *= UWvalues[iUW](iflavor,jflavor);
					else{
						term = 0.;
						break;
					}
				}


			}
			valtmp += term;
			//cout<<term<<" + ";
		}
	}

	//cout<<endl;
	//cout<<"   calculated value = "<<valtmp<<endl;
	//cout<<"****************************"<<endl<<flush;

	return valtmp;
};


void setUWinter::setvalue(const double& val){ value = val; };
double setUWinter::getvalue(){ return value; };
ostream& operator << (ostream& ostr, const setUWinter& SUI){
	ostr<<"** setUWinter info **"<<endl
		<<"Nfermionicloops: "<<SUI.Nfermionicloops<<endl;
	ostr<<"isloopflavorfixed:"<<endl;
	for(int i=0; i<SUI.isloopflavorfixed.size(); i++) ostr<<SUI.isloopflavorfixed[i];
	ostr<<endl<<"loopflavor:"<<endl;
	for(int i=0; i<SUI.loopflavor.size(); i++) ostr<<SUI.loopflavor[i];
	ostr<<endl<<"external lines"<<endl;
	for(int i=0; i<SUI.UWinters.size(); i++)
		ostr<<"("<<SUI.iloops[i]<<","<<SUI.jloops[i]<<"): "<<SUI.UWinters[i]<<endl;
	ostr<<endl<<"UW values"<<endl;
	for(int i=0; i<SUI.UWvalues.size(); i++)
		ostr<<SUI.UWvalues[i]<<endl;
	ostr<<endl<<endl<<"value : "<<SUI.value<<endl;
	ostr<<"*********************"<<endl;

	return ostr;
};*/
