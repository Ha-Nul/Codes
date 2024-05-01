#ifndef __graph_hpp__included__
#define __graph_hpp__included__

#include <vector>
#include <list>

class edge{
	public:
	int type, i, j;
	edge();
	edge(const int& i_i,const int& j_i);
	edge(const int& type_i, const int& i_i,const int& j_i);
};

class graph{
	protected:
	int Nvertice;

	std::vector<std::vector<int> > adjL;
	std::vector<std::vector<int> > adjLtype;

	std::vector<int> dfi;
	std::vector<bool> discovered;
	std::vector<edge> bridge;

	//std::vector<int> lowpt, deg;
	//std::stack<int> path;
	//std::vector<vector<int> > sigma
	//std::vector<edge> E;

	//bool connected, two_edge_connected, three_edge_connected;

	public:
	graph();
	void set(const int& Nvertice_i, std::vector<std::vector<int> > adjL_i);
	void set(const int& Nvertice_i, std::vector<std::vector<int> > adjL_i, std::vector<std::vector<int> > adjLtype_i);
	bool edge_connectivity();
	void dfs(const int& v, int& dfi_running);
	void dfs_assign(std::vector<int>& supervertice, const int& assign_i, const int& v);
	std::vector<int> assignsupervertice(int& Nsupervertice_m);
	int one_two_edge_connectivity();
	bool two_edge_connectivity();
	bool two_edge_connectivity(const int& type_i);
	int getNloop();
	int dfs_bridge(const int& v, const int& parent, int& dfi_running);
	int dfs_bridge(const int& v, const int& parent, const int& type, int& dfi_running);
	bool three_edge_connectivity();
	void three_edge_connect();
	//int dfs(int& v, int& dfi_running);
	//int dfs(int v, std::vector<bool> discovered_i, int parent, int& dfi_i);
};

class three_edge: public graph{
	std::vector<int> lowpt;

	std::vector<int> nd;
	//std::vector<int> parent;
	std::vector<int> next;

	std::vector<std::list<int> > sigma;
	//std::vector<std::list<int> > adjL_outgoing_backedge;

	std::vector<int> deg;

	//std::vector<std::vector<int> > path;

	public:
	//graph();
	//void set(const int& Nvertice_i, std::vector<std::vector<int> > adjL_i);
	//void set(const int& Nvertice_i, std::vector<std::vector<int> > adjL_i, std::vector<std::vector<int> > adjLtype_i);
	//bool edge_connectivity();
	//void dfs(const int& v, int& dfi_running);
	//bool two_edge_connectivity();
	//bool two_edge_connectivity(const int& type_i);
	//int dfs_bridge(const int& v, const int& parent, const int& type, int& dfi_running);
	bool three_edge_connectivity();
	void three_edge_connect(const int& w, const int& v, int& dfi_running);
	void absorb_path(const int& x0, const int& x1);
	void absorb_path(const int& x0, const int& x1, const int& xk);
	//void absorb_path(std::vector<int>& path_i);
	//void absorb_path(std::vector<int>& path_i, const int& until);
	//int dfs(int& v, int& dfi_running);
	//int dfs(int v, std::vector<bool> discovered_i, int parent, int& dfi_i);
	//void edge_cut(const int& w, const int& u);
	//void selfloop_cut(const int& w);
};

#endif
