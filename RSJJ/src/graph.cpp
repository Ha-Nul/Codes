#include <algorithm>
#include <iostream>
#include <iterator>
#include "graph.hpp"
using namespace std;

edge::edge(){};
edge::edge(const int& i_i,const int& j_i): type(0), i(i_i), j(j_i){};
edge::edge(const int& type_i, const int& i_i,const int& j_i): type(type_i), i(i_i), j(j_i){};

graph::graph(){};
void graph::set(const int& Nvertice_i, vector<vector<int> > adjL_i){
       	Nvertice = Nvertice_i; 
	adjL = adjL_i;
};
void graph::set(const int& Nvertice_i, vector<vector<int> > adjL_i, vector<vector<int> > adjLtype_i){
       	Nvertice = Nvertice_i; 
	adjL = adjL_i;
	adjLtype = adjLtype_i;
};
bool graph::edge_connectivity(){
	dfi.resize(Nvertice);
	discovered.resize(Nvertice);
	for(auto iter=discovered.begin(); iter!=discovered.end(); ++iter) (*iter)=false;
	int root=0, dfi_running=0;

	dfs(root, dfi_running);

	bool connectivity = true;
	for(auto iter=discovered.begin(); iter!=discovered.end(); ++iter){
		connectivity = (connectivity && (*iter));
		if(!connectivity) break;
	}

	return connectivity;
	//for(auto& it : bridge) cout<<"("<<it.i<<","<<it.j<<")"<<endl;
};
void graph::dfs(const int& v, int& dfi_running){
		dfi[v] = (++dfi_running);
		discovered[v] = true;

		for(int w: adjL[v]) if(!discovered[w]) dfs(w,dfi_running);
};
void graph::dfs_assign(vector<int>& supervertice, const int& assign_i, const int& v){
		discovered[v] = true;
		supervertice[v] = assign_i;
		for(int w: adjL[v]) if(!discovered[w]) dfs_assign(supervertice,assign_i,w);

};
vector<int> graph::assignsupervertice(int& Nsupervertice_m){
	vector<int> supervertice(Nvertice);
	discovered.resize(Nvertice);
	for(auto iter=discovered.begin(); iter!=discovered.end(); ++iter) (*iter)=false;
	int isv_run = 0;
	for(int i=0; i<Nvertice; i++)
		if(!discovered[i])
			dfs_assign(supervertice,isv_run++,i);

	Nsupervertice_m = isv_run;
	return supervertice;
};

int graph::one_two_edge_connectivity(){
	dfi.resize(Nvertice);
	discovered.resize(Nvertice);
	for(auto iter=discovered.begin(); iter!=discovered.end(); ++iter) (*iter)=false;
	bridge.resize(0);
	//int root=0, parent=-1, auxtype=-1,  dfi_running=0;
	int root=0, parent=-1,  dfi_running=0;

	//dfs_bridge(root, parent, auxtype, dfi_running);
	dfs_bridge(root, parent, dfi_running);

	if( !(dfi_running==Nvertice) ) return 0;
	else if( static_cast<bool>(bridge.size()) ) return 1;
	else return 2;

	//return !(static_cast<bool>(bridge.size()));
	//for(auto& it : bridge) cout<<"("<<it.i<<","<<it.j<<")"<<endl;
};

bool graph::two_edge_connectivity(){
	dfi.resize(Nvertice);
	discovered.resize(Nvertice);
	for(auto iter=discovered.begin(); iter!=discovered.end(); ++iter) (*iter)=false;
	bridge.resize(0);
	//int root=0, parent=-1, auxtype=-1,  dfi_running=0;
	int root=0, parent=-1,  dfi_running=0;

	//dfs_bridge(root, parent, auxtype, dfi_running);
	dfs_bridge(root, parent, dfi_running);

	//return ((dfi_running==Nvertice) && !(static_cast<bool>(bridge.size())));
	return !(static_cast<bool>(bridge.size()));
	//for(auto& it : bridge) cout<<"("<<it.i<<","<<it.j<<")"<<endl;
};
bool graph::two_edge_connectivity(const int& type_i){
	dfi.resize(Nvertice);
	discovered.resize(Nvertice);
	for(auto iter=discovered.begin(); iter!=discovered.end(); ++iter) (*iter)=false;
	bridge.resize(0);
	int root=0, parent=-1, auxtype=-1,  dfi_running=0;

	//cout<<"depth-first search"<<endl;
	dfs_bridge(root, parent, auxtype, dfi_running);

	int count = 0;
	for(auto& it : bridge){
		//cout<<"("<<it.type<<","<<it.i<<","<<it.j<<")"<<endl;
		if(it.type==type_i) count++;
	}

	//if(!(static_cast<bool>(count))) cout<<"two-edge-connected"<<endl;
	//else cout<<"NOT two-edge-connected"<<endl;;

	return !(static_cast<bool>(count));
};
int graph::getNloop(){
	dfi.resize(Nvertice);
	discovered.resize(Nvertice);
	for(auto iter=discovered.begin(); iter!=discovered.end(); ++iter) (*iter)=false;
	int ancestor = -1, auxtype = -1;
	int root, dfi_running;
	int Ncluster=0;
	for(auto iter=discovered.begin(); iter!=discovered.end(); ++iter){
		if(!(*iter)){
			bridge.clear();
			root=distance(discovered.begin(),iter);
			dfi_running=0;
			dfs_bridge(root, ancestor, auxtype, dfi_running);
			if(bridge.size()==0 && adjL[root].size()!=0) Ncluster++;
		}
	}

	return Ncluster;
};
int graph::dfs_bridge(const int& v, const int& parent, int& dfi_running){
	discovered[v] = true;
	dfi[v] = (++dfi_running);
	//cout<<"v:"<<v<<", dfi[v]:"<<dfi[v]<<" ("<<type<<","<<parent<<","<<v<<")"<<endl;

	int lowpt = dfi[v];
	bool parent_found = false;
	//for(int w: adjL[v]){
	int w;
	for(int i=0; i<adjL[v].size(); i++){
		w = adjL[v][i];
		if(!discovered[w]) lowpt = min(lowpt, dfs_bridge(w, v, dfi_running));
		else if(w==parent && (!parent_found)) parent_found = true;
		else if(w!=parent || parent_found) lowpt = min(lowpt, dfi[w]);
	}

	if(lowpt==dfi[v] && parent!=-1) bridge.emplace_back( edge(parent, v) );

	return lowpt;
};
int graph::dfs_bridge(const int& v, const int& parent, const int& type, int& dfi_running){
	discovered[v] = true;
	dfi[v] = (++dfi_running);
	//cout<<"v:"<<v<<", dfi[v]:"<<dfi[v]<<" ("<<type<<","<<parent<<","<<v<<")"<<endl;

	int lowpt = dfi[v];
	bool parent_found = false;
	//for(int w: adjL[v]){
	int w;
	for(int i=0; i<adjL[v].size(); i++){
		w = adjL[v][i];
		if(!discovered[w]) lowpt = min(lowpt, dfs_bridge(w, v, adjLtype[v][i], dfi_running));
		else if(w==parent && (!parent_found)) parent_found = true;
		else if(w!=parent || parent_found) lowpt = min(lowpt, dfi[w]);
	}

	if(lowpt==dfi[v] && parent!=-1) bridge.emplace_back( edge(type, parent, v) );

	return lowpt;
};

bool three_edge::three_edge_connectivity(){
	//NOTE THAT Nvertice SHOULD BE LARGER THAN 0
	//cout<<"* three_edge_connectivity module *"<<endl<<flush;
	//cout<<"adjL"<<endl;
	//for(auto i: adjL){
	//	for(auto j: i) cout<<j;
	//	cout<<endl;
	//}
	//cout<<endl;

	dfi.resize(Nvertice);
	discovered.resize(Nvertice);

	lowpt.resize(Nvertice);

	nd.resize(Nvertice);
	next.resize(Nvertice);

	sigma.resize(Nvertice);

	deg.resize(Nvertice);

	for(int v=0; v<Nvertice; v++){
		deg[v]=0; 
		sigma[v].clear();
		sigma[v].push_back(v);
		discovered[v] = false;
	}

	int root = 0, ancestor = -1, dfi_running = 1;
	three_edge_connect(root, ancestor, dfi_running);

	bool tmp = true;
	int cluster_size;
	for(int v=0; v<Nvertice; v++){
		cluster_size = sigma[v].size();
		if(cluster_size!=0 && cluster_size<Nvertice){
			tmp = false;
			break;
		}
	}
	//if(tmp) cout<<"it is three-edge-connected"<<endl<<flush;
	//else cout<<"it is NOT three-edge-connected"<<endl<<flush;

	return tmp;
};
/*bool three_edge::three_edge_connectivity(){
	//NOTE THAT Nvertice SHOULD BE LARGER THAN 0
	//cout<<"* three_edge_connectivity module *"<<endl<<flush;
	dfi.resize(Nvertice);
	discovered.resize(Nvertice);

	lowpt.resize(Nvertice);
	sigma.resize(Nvertice);
	path.resize(Nvertice);
	deg.resize(Nvertice);

	for(int v=0; v<Nvertice; v++){
		deg[v]=0; 
		sigma[v].clear();
		sigma[v].push_back(v);
		discovered[v] = false;
	}

	int root = 0, ancestor = -1, dfi_running = 1;
	three_edge_connect(root, ancestor, dfi_running);

	//cout<<"lowpt"<<endl;
	//for(auto i: lowpt) cout<<i<<" ";
	//cout<<endl;

	bool tmp = true;
	int cluster_size;
	for(int v=0; v<Nvertice; v++){
		cluster_size = sigma[v].size();
		if(cluster_size!=0 && cluster_size<Nvertice){
			tmp = false;
			break;
		}
	}
	//if(tmp) cout<<"it is three-edge-connected"<<endl<<flush;
	//else cout<<"it is NOT three-edge-connected"<<endl<<flush;

	return tmp;
};*/
void three_edge::three_edge_connect(const int& w, const int& v, int& dfi_running){
	discovered[w] = true;
	dfi[w] = (dfi_running++);
	lowpt[w] = dfi[w];
	next[w] = w;
	//cout<<"for w:"<<w<<","<<endl
	//       <<"    dfi[w]:"<<dfi[w]<<endl
	//       <<"    lowpt[w]:"<<lowpt[w]<<endl
	//       <<"    next[w]:"<<next[w]<<endl;

	int u, Pu_x0;
	bool parent_found = false;
	for(int u: adjL[w]){
		//cout<<"* adjL "<<u<<" from "<<w<<" *"<<flush;
		deg[w]++;
		if(!discovered[u]){
			//cout<<" : outgoing tree-edge to "<<u<<endl;
			three_edge_connect(u,w,dfi_running);

			//cout<<"now at "<<w<<", investigating "<<u<<endl<<flush;
			//cout<<"deg["<<u<<"] = "<<deg[u]<<endl;

			if(deg[u]==2){ Pu_x0 = next[u]; }
			else{ Pu_x0 = u; }

			if(!(lowpt[w]>lowpt[u])){
				//cout<<"!(lowpt[w]>lowpt[u]) case"<<endl;
				absorb_path(w,Pu_x0);
			}
			else{
				//cout<<"lowpt[w]>lowpt[u] case"<<endl;
				lowpt[w] = lowpt[u];
				absorb_path(w,next[w]);
				next[w] = Pu_x0;
			}
		}
		else if(u==v && (!parent_found)){
			//cout<<" : incoming tree-edge"<<endl;
		       	parent_found = true;
		}
		else{
			if(dfi[w]>dfi[u]){
				//cout<<" : outgoing backedge of "<<w<<endl;
				//adjL_outgoing_backedge[w].push_back(u);
				if(dfi[u]<lowpt[w]){
					//cout<<".. and update (dfi[u]:"<<dfi[u]<<", lowpt[w]:"<<lowpt[w]<<")"<<endl;
					absorb_path(w,next[w]);
					lowpt[w] = dfi[u];
					next[w] = u;
				}
				//else{
				//	cout<<"adjL_outgoing_backedge["<<w<<"] : ";
				//	for(auto i: adjL_outgoing_backedge[w]) cout<<i;
				//	cout<<endl;
				//}
				//cout<<"next["<<w<<"] = "<<next[w]<<endl;
			}
			else if(dfi[w]<dfi[u]){
				//cout<<" : incoming backedge of "<<w<<endl;
				deg[w] -= 2;
				absorb_path(w,next[w],u);
			}
			else{
				//cout<<" : self loop of "<<w<<endl;
				deg[w] -= 2;
			}
		}
	};
	nd[w] = dfi_running - dfi[w];
};
void three_edge::absorb_path(const int& x0, const int& x1){
	//cout<<"* full path absorbtion : "<<flush;
	//cout<<"("<<x0<<")";
	int xi = x1;
	while(dfi[xi]>dfi[x0]){
		//cout<<xi;
		deg[x0] += deg[xi] - 2;
		sigma[x0].insert(sigma[x0].end(),sigma[xi].begin(),sigma[xi].end());
		sigma[xi].clear();
		//adjL_outgoing_backedge[x0].insert(
		//		adjL_outgoing_backedge[x0].end(),
		//		adjL_outgoing_backedge[xi].begin(),
		//		adjL_outgoing_backedge[xi].end()
		//		);
		xi = next[xi];
	};
	//cout<<endl;
	//cout<<"sigma["<<x0<<"] : ";
	//for(auto i: sigma[x0]) cout<<i;
	//cout<<endl;
	//cout<<"adjL_outgoing_backedge["<<x0<<"] : ";
	//for(auto i: adjL_outgoing_backedge[x0]) cout<<i;
	//cout<<endl;
};
void three_edge::absorb_path(const int& x0, const int& x1, const int& xk){
	//cout<<"* finite path absorbtion : "<<flush;
	//cout<<"("<<x0<<")";
	if(dfi[x1]>dfi[x0]) {
		int xi = x1;
		while((dfi[xi]>dfi[x0]) && !(dfi[xi]>dfi[xk]) && !(dfi[xk]>dfi[xi]+nd[xi]-1)){
			//cout<<xi;
			deg[x0] += deg[xi] - 2;
			sigma[x0].insert(sigma[x0].end(),sigma[xi].begin(),sigma[xi].end());
			sigma[xi].clear();
			//adjL_outgoing_backedge[x0].insert(
			//		adjL_outgoing_backedge[x0].end(),
			//		adjL_outgoing_backedge[xi].begin(),
			//		adjL_outgoing_backedge[xi].end()
			//		);
			xi = next[xi];
		}
		next[x0] = xi;
	}

	//if(dfi[x1]>dfi[x0]) next[x0] = next[xk];
	//cout<<endl;
	//cout<<"sigma["<<x0<<"] : ";
	//for(auto i: sigma[x0]) cout<<i;
	//cout<<endl;
	//cout<<"adjL_outgoing_backedge["<<x0<<"] : ";
	//for(auto i: adjL_outgoing_backedge[x0]) cout<<i;
	//cout<<endl;
};
/*void three_edge::three_edge_connect(const int& w, const int& v, int& dfi_running){
	discovered[w] = true;
	dfi[w] = (dfi_running++);
	lowpt[w] = dfi[w];
	path[w].clear(); path[w].push_back(w);
	//cout<<"for w:"<<w<<", dfi[w]:"<<dfi[w]<<endl;

	int u;	
	bool parent_found = false;
	//for(int u: adjL[w]){
	for(int index=0; index<adjL[w].size(); index++){
		u = adjL[w][index];
		//cout<<"* adjL "<<u<<" from "<<w<<" *"<<flush;
		deg[w] += 1;
		if(!discovered[u]){
			//cout<<" : not discovered yet (out-going tree edge)"<<endl;
			three_edge_connect(u,w,dfi_running);
			//cout<<"now at "<<w<<", investigating "<<u<<endl<<flush;
			//cout<<"deg["<<u<<"] = "<<deg[u]<<endl;
			if(deg[u]==2){ 
				path[u].pop_back(); 
				edge_cut(w,u);
			}

			//cout<<"lowpt["<<w<<"] = "<<lowpt[w]<<", lowpt["<<u<<"] = "<<lowpt[u]<<endl;

			if(!(lowpt[w]>lowpt[u])){
				//cout<<"!(lowpt[w]>lowpt[u]) case"<<endl;
				path[u].push_back(w);
				absorb_path(path[u]);
				//path[w] = path[u];
				//path[w].push_back(w);
				//absorb_path(path[w]);
				//path[w] = path[u];
				//path[u].clear();
				//cout<<"after path update"<<endl;
				//for(auto i: path){
				//	for(auto j:i) cout<<j;
				//	cout<<endl;
				//}

			}
			else{
				//cout<<"lowpt[w]>lowpt[u] case"<<endl;
				lowpt[w] = lowpt[u];
				absorb_path(path[w]);
				path[w] = path[u];
				path[w].push_back(w);
				//cout<<"after path update"<<endl;
				//for(auto i: path){
				//	for(auto j:i) cout<<j;
				//	cout<<endl;
				//}
			}
		}
		else if(u==v && (!parent_found)){
			//cout<<" : incoming tree-edge"<<endl;
		       	parent_found = true;
		}
		else{
			if(dfi[w]>dfi[u]){
				//cout<<" : outgoing backedge of "<<w<<endl;
				if(dfi[u]<lowpt[w]){
					//cout<<".. and update"<<endl;
					absorb_path(path[w]);
					lowpt[w] = dfi[u];
					path[w].clear(); path[w].push_back(w);
				}
			}
			else if(dfi[w]<dfi[u]){
				//cout<<" : incoming backedge of "<<w<<endl;
				//cout<<"--during the incoming backedge update"<<endl;
				//cout<<"--deg["<<w<<"]="<<deg[w]<<endl;
				deg[w] -= 2;
				//cout<<"--deg["<<w<<"]="<<deg[w]<<endl;
				absorb_path(path[w],u);
				//cout<<"--deg["<<w<<"]="<<deg[w]<<endl;
			}
			else{
				//cout<<" : self-loop of "<<w<<endl;
				deg[w] -= 2;
				//selfloop_cut(w);
			}
		}
		//cout<<"deg["<<w<<"] = "<<deg[w]<<", after update "<<u<<endl<<flush;
		//cout<<"adjL: "<<endl;
		//for(auto i:adjL){
		//	for(auto j:i) cout<<j;
		//	cout<<endl;
		//}
	};

	//vector<int>::iterator itend_w = remove(adjL[w].begin(),adjL[w].end(),-1);
	//adjL[w].resize( distance(adjL[w].begin(),itend_w) );
	//itend_w = remove(adjL[w].begin(),adjL[w].end(),w);
	//adjL[w].resize( distance(adjL[w].begin(),itend_w) );

	//cout<<"after finishing update at "<<w<<endl;
	//cout<<"adjL: "<<endl;
	//for(auto i:adjL){
	//	for(auto j:i) cout<<j;
	//	cout<<endl;
	//}
	//cout<<"sigma: "<<endl;
	//for(auto i: sigma){
	//	for(auto j:i) cout<<j;
	//	cout<<endl;
	//}
	//cout<<"path"<<endl;
	//for(auto i: path){
	//	for(auto j:i) cout<<j;
	//	cout<<endl;
	//}
	//cout<<"deg["<<w<<"] = "<<deg[w]<<endl;
};*/
/*void three_edge::absorb_path(const int* w){
	//cout<<"* path absorbtion *"<<endl<<flush;
	x0 = w;
	xi = next[x0];
	while(dfi[xi]>dfi[x0]){
		deg[x0] += deg[xi] - 2;
		edge_cut(x0,xi);
		sigma[x0].insert(sigma[x0].end(),sigma[xi].begin(),sigma[xi].end());
		adjL_outgoing_backedge[x0].insert(
				adjL_outgoing_backedge[x0].end(),
				adjL_outgoing_backedge[xi].begin(),
				adjL_outgoing_backedge[xi].end()
				);
		xi = next[xi];
	};
	next[x0] = lowpt[x0];
};*/
/*void three_edge::absorb_path(vector<int>& path_i){
	//cout<<"* path absorbtion *"<<endl<<flush;
	//cout<<"before absorbtion: ";
	//for(int i: path_i) cout<<i;
	//cout<<endl;
	int x0 = path_i.back(); path_i.pop_back();
	int xi;
	//cout<<"("<<x0<<")";
	while(!(path_i.empty())){
		xi = path_i.back();
		//cout<<xi;
		deg[x0] += deg[xi] - 2;
		for(int xi_j: sigma[xi]){ sigma[x0].push_back(xi_j); }
		sigma[xi].clear();
		path_i.pop_back();
		edge_cut(x0,xi);
	};
	path_i.push_back(x0);
	//cout<<endl;
	//cout<<"after absorbtion: ";
	//for(int i: path_i) cout<<i;
	//cout<<endl;
};*/
/*void three_edge::absorb_path(const int& w, const int& u){
	//cout<<"* finite path absorbtion *"<<endl<<flush;
	if(dfi[next[w]]>dfi[w]){
		int dfi_u = dfi[u];

		int x0 = path_i.back(); 
		path_i.pop_back();
		int xi = next[w];
		while(dfi[xi]>dfi_u){
			deg[x0] += deg[xi] - 2;
			edge_cut(x0,xi);
			sigma[x0].insert(sigma[x0].end(),sigma[xi].begin(),sigma[xi].end());
			adjL_outgoing_backedge[x0].insert(
					adjL_outgoing_backedge[x0].end(),
					adjL_outgoing_backedge[xi].begin(),
					adjL_outgoing_backedge[xi].end()
					);
			xi = next[xi];
		};
		next[x0] = next[xi];
	}
};*/

/*void three_edge::absorb_path(vector<int>& path_i, const int& until){
	//cout<<"* finite path absorbtion *"<<endl<<flush;
	//cout<<"until = "<<until<<endl;
	//cout<<"sigma"<<endl;
	//for(auto i: sigma){
	//	for(auto j:i) cout<<j;
	//	cout<<endl;
	//}
	//cout<<"adjL"<<endl;
	//for(auto i: adjL){
	//	for(auto j:i) cout<<j;
	//	cout<<endl;
	//}

	//cout<<"before absorbtion"<<endl;
	//for(int i: path_i) cout<<i;
	//cout<<endl;

	bool flag, found, empty;
	int x0 = path_i.back(); 
	//cout<<"("<<x0<<",deg["<<x0<<"]="<<deg[x0]<<")";
	found = false;
	for(int x0_j: sigma[x0]){
		if(x0_j==until){
			found = true;
			break;
		}
	}
	if(found) return;

	path_i.pop_back();
	int xi = path_i.back();
	do{
		//cout<<xi<<"(deg["<<xi<<"]="<<deg[xi]<<")";

		found = false;
		for(int xi_j: sigma[xi]){
			if(xi_j==until){
				found = true;
				break;
			}
		}

		//cout<<"before deg update, deg["<<x0<<"] = "<<deg[x0]<<endl;
		deg[x0] += deg[xi] - 2; // ALERT! LINE OF DEBATE
		//cout<<"after deg update, deg["<<x0<<"] = "<<deg[x0]<<endl;
		for(int xi_j: sigma[xi]){ sigma[x0].push_back(xi_j); }
		sigma[xi].clear();
		path_i.pop_back();
		edge_cut(x0,xi);

		empty = path_i.empty();
		if(!empty) xi = path_i.back();
		flag = !(found||empty);
	} while(flag);

	//while((!path_i.empty())&&(until!=path_i.back())){
	//	xi = path_i.back();
	//	cout<<xi<<"(deg["<<xi<<"]="<<deg[xi]<<")";
	//	deg[x0] += deg[xi] - 2; // ALERT! LINE OF DEBATE
	//	for(int xi_j: sigma[xi]){ sigma[x0].push_back(xi_j); }
	//	sigma[xi].clear();
	//	path_i.pop_back();
	//};
	//if((!path_i.empty())&&(until==path_i.back())){
	//	xi = path_i.back();
	//	cout<<xi<<"(deg["<<xi<<"]="<<deg[xi]<<")";
	//	deg[x0] += deg[xi] - 2; // ALERT! LINE OF DEBATE
	//	for(int xi_j: sigma[xi]){ sigma[x0].push_back(xi_j); }
	//	sigma[xi].clear();
	//	path_i.pop_back();
	//}

	//bool until_found = false;
	//if((!path_i.empty())&&(until==path_i.back())){
	//	cout<<"HMMMMMMM"<<endl<<flush;
	//	until_found = true;
	//}
	////else until_found = false;
	//bool flag = false;
	//while(!(path_i.empty() || until_found)){
	//	xi = path_i.back();
	//	cout<<xi<<"(deg["<<xi<<"]="<<deg[xi]<<")";
	//	deg[x0] += deg[xi] - 2; // ALERT! LINE OF DEBATE
	//	for(int xi_j: sigma[xi]){ sigma[x0].push_back(xi_j); }
	//	sigma[xi].clear();
	//	path_i.pop_back();
	//	if((!path_i.empty())&&(until==path_i.back())) until_found = true;
	//};
	//if((!path_i.empty())&&(until_found)){
	//	xi = path_i.back();
	//	cout<<xi<<"(deg["<<xi<<"]="<<deg[xi]<<")";
	//	deg[x0] += deg[xi] - 2; // ALERT! LINE OF DEBATE
	//	for(int xi_j: sigma[xi]){ sigma[x0].push_back(xi_j); }
	//	sigma[xi].clear();
	//	path_i.pop_back();
	//}

	path_i.push_back(x0);
	//cout<<endl;

	//cout<<"after absorbtion"<<endl;
	//for(int i: path_i) cout<<i;
	//cout<<endl;
};*/
/*void three_edge::edge_cut(const int& w, const int& u){
	//cout<<"* edge_cut ("<<w<<","<<u<<")*"<<endl<<flush;
	//cout<<"before edge cut"<<endl;
	//for(auto i:adjL){
	//	for(auto j:i) cout<<j;
	//	cout<<endl;
	//}
	//for(vector<int>::iterator iter=adjL[w].begin(); iter!=adjL[w].end(); ++iter){
	//	if(*iter==u){
	//		//adjL[w].erase(iter);
	//		(*iter) = -1;
	//		break;
	//	}
	//}
	//cout<<"after put flag to tree edge cut"<<endl;
	//for(auto i:adjL){
	//	for(auto j:i) cout<<j;
	//	cout<<endl;
	//}

	//for(vector<int>::iterator iter=adjL[u].begin(); iter!=adjL[u].end(); ++iter){
	//	if(*iter==w){
	//		adjL[u].erase(iter);
	//	}
	//}

	//for(vector<int>::iterator iter=adjL[u].begin(); iter!=adjL[u].end(); ++iter){
	//	if(*iter==w) adjL[u].erase(iter);
	//}

	//auto itend_w = remove(adjL[w].begin(),adjL[w].end(),u);
	//adjL[w].resize( distance(adjL[w].begin(),itend_w) );
	//vector<int>::iterator itend_u = remove(adjL[u].begin(),adjL[u].end(),w);
	//adjL[u].resize( distance(adjL[u].begin(),itend_u) );
	//cout<<"after remove in u"<<endl;
	//for(auto i:adjL){
	//	for(auto j:i) cout<<j;
	//	cout<<endl;
	//}
	//adjL[w].insert(adjL[w].end(),adjL[u].begin(),adjL[u].end());
	adjL[u].clear();
	//cout<<"after attach embodiments to tree edge cut"<<endl;
	//for(auto i:adjL){
	//	for(auto j:i) cout<<j;
	//	cout<<endl;
	//}

	//for(vector<vector<int> >::iterator adjL_it=adjL.begin(); adjL_it!=adjL.end(); ++adjL_it)
	//	std::replace(adjL_it->begin(),adjL_it->end(),u,w);

	//for(auto i:adjL){
	//	for(auto j:i) cout<<j;
	//	cout<<endl;
	//}
	//for(auto it=adjL[w].begin(); it!=adjL[w].end(); ++it){
	//	if((*it)==u) adjL[w].erase(it);
	//}
	//cout<<"update adjL[w] completed"<<endl<<flush;
	//for(auto it=adjL[u].begin(); it!=adjL[u].end(); ++it){
	//	if(*it==w) adjL[u].erase(it);
	//}
	//cout<<"update adjL[u] completed"<<endl<<flush;

	//cout<<"after edge cut"<<endl;
	//for(auto i:adjL){
	//	for(auto j:i) cout<<j;
	//	cout<<endl;
	//}
};*/
/*void three_edge::edge_cut(const int& w, const int& u){
	//cout<<"* edge_cut ("<<w<<","<<u<<")*"<<endl<<flush;
	//cout<<"before edge cut"<<endl;
	//for(auto i:adjL){
	//	for(auto j:i) cout<<j;
	//	cout<<endl;
	//}
	for(vector<int>::iterator iter=adjL[w].begin(); iter!=adjL[w].end(); ++iter){
		if(*iter==u){
			//adjL[w].erase(iter);
			(*iter) = -1;
			break;
		}
	}
	cout<<"after put flag to tree edge cut"<<endl;
	for(auto i:adjL){
		for(auto j:i) cout<<j;
		cout<<endl;
	}

	for(vector<int>::iterator iter=adjL[u].begin(); iter!=adjL[u].end(); ++iter){
		if(*iter==w){
			adjL[u].erase(iter);
		}
	}

	//for(vector<int>::iterator iter=adjL[u].begin(); iter!=adjL[u].end(); ++iter){
	//	if(*iter==w) adjL[u].erase(iter);
	//}

	//auto itend_w = remove(adjL[w].begin(),adjL[w].end(),u);
	//adjL[w].resize( distance(adjL[w].begin(),itend_w) );
	//vector<int>::iterator itend_u = remove(adjL[u].begin(),adjL[u].end(),w);
	//adjL[u].resize( distance(adjL[u].begin(),itend_u) );
	//cout<<"after remove in u"<<endl;
	//for(auto i:adjL){
	//	for(auto j:i) cout<<j;
	//	cout<<endl;
	//}
	//adjL[w].insert(adjL[w].end(),adjL[u].begin(),adjL[u].end());
	adjL[u].clear();
	//cout<<"after attach embodiments to tree edge cut"<<endl;
	//for(auto i:adjL){
	//	for(auto j:i) cout<<j;
	//	cout<<endl;
	//}

	for(vector<vector<int> >::iterator adjL_it=adjL.begin(); adjL_it!=adjL.end(); ++adjL_it)
		std::replace(adjL_it->begin(),adjL_it->end(),u,w);

	//for(auto i:adjL){
	//	for(auto j:i) cout<<j;
	//	cout<<endl;
	//}
	//for(auto it=adjL[w].begin(); it!=adjL[w].end(); ++it){
	//	if((*it)==u) adjL[w].erase(it);
	//}
	//cout<<"update adjL[w] completed"<<endl<<flush;
	//for(auto it=adjL[u].begin(); it!=adjL[u].end(); ++it){
	//	if(*it==w) adjL[u].erase(it);
	//}
	//cout<<"update adjL[u] completed"<<endl<<flush;

	//cout<<"after edge cut"<<endl;
	//for(auto i:adjL){
	//	for(auto j:i) cout<<j;
	//	cout<<endl;
	//}
};*/
/*void three_edge::selfloop_cut(const int& w){
	//cout<<"* selfloop_cut ("<<w<<")"<<endl<<flush;
	//cout<<"before edge cut, adjL:"<<endl;
	//for(auto i:adjL){
	//	for(auto j:i) cout<<j;
	//	cout<<endl;
	//}
	for(vector<int>::iterator iter=adjL[w].begin(); iter!=adjL[w].end(); ++iter){
		if(*iter==w){
			(*iter) = -1;
			//adjL[w].erase(iter);
			break;
		}
	}
	//cout<<"after edge cut, adjL:"<<endl;
	//for(auto i:adjL){
	//	for(auto j:i) cout<<j;
	//	cout<<endl;
	//}
};*/
