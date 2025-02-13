#pragma once

#include <Eigen/Dense>
#include <list>
#include <iterator>

class coordinate{
	public:
	double time;
	//Eigen::VectorXi space;

	coordinate();
	coordinate(const double& t_i);
	void settime(const double& t_i);
	void setref();
	void genrand();
	//void genrand(const coordinate& x_c);
	void genrand(const coordinate& x_i, const coordinate& x_j);
	//double getsitegenprobability(const coordinate& x_c);
	double gettime() const;
	//Eigen::VectorXi getspace() const;
	bool operator == (const coordinate& x_i);
	bool operator != (const coordinate& x_i);
	bool operator < (const coordinate& xi) const;
	friend std::ostream& operator << (std::ostream& ostr, const coordinate& x);
};

/*
coordinate mean(const coordinate& x_i, const coordinate& x_j);
coordinate mean(const coordinate& x_1, const coordinate& x_2, const coordinate& x_3, const coordinate& x_4);

class line;
class fermionicloop;

class vertex{
	unsigned flavor;
	coordinate x;

	std::vector<std::list<line>::iterator> connected_lines;
	std::list<fermionicloop>::iterator connected_loop;

	public:
	vertex();
	vertex(unsigned flavor_i, coordinate x_i);
	void assign(const vertex& v_i);
	void setflavor(const int& flavor_i);
	void setx(const coordinate& x_i);
	void setconnected_lines(std::list<line>::iterator W_iter, std::list<line>::iterator Gin_iter, std::list<line>::iterator Gout_iter);
	void setconnected_Gin(std::list<line>::iterator Gin_iter);
	void setconnected_W(std::list<line>::iterator W_iter);
	void setconnected_loop(std::list<fermionicloop>::iterator loop_iter);
	unsigned getflavor();
	coordinate getx();
	std::list<line>::iterator getconnected_W();
	std::list<line>::iterator getconnected_Gin();
	std::list<line>::iterator getconnected_Gout();
	std::list<fermionicloop>::iterator getconnected_loop();
	friend std::ostream& operator << (std::ostream& ostr, const vertex& v);
};

class line{
	bool physical;
	int type; // 0:G, 1:W

	//unsigned dress;
	unsigned iflavor, jflavor;
	coordinate xi, xj;

	std::list<vertex>::iterator connected_vertices[2];
	//std::list<fermionicloop>::iterator connected_loops[2];

	double valueij, valueji;

	public:
	line(bool physical_i, int type_i, unsigned iflavor_i, unsigned jflavor_i, coordinate xi_i, coordinate xj_i);
	line(const line& Line_i);
	void assign(const line& line_i);
	void setphysical(const bool& physical_i);
	//void setdress(const int& dress_i);
	void setflavor(const int& iflavor_i, const int& jflavor_i);
	void setx(const coordinate& xi_i,const coordinate& xj_i);
	void setxi(const coordinate& xj_i);
	void setxj(const coordinate& xj_i);
	void setj(const unsigned& jflavor_i, const coordinate& xj_i);
	void setconnected_vertices(std::list<vertex>::iterator vi_iter,std::list<vertex>::iterator vj_iter);
	void setconnected_vj(std::list<vertex>::iterator vj_iter);
	//void setconnected_loops(std::list<fermionicloop>::iterator loopi_iter,std::list<fermionicloop>::iterator loopj_iter);
	bool getphysical() const;
	int gettype() const;
	//bool getisG();
	//unsigned getdress() const;
	unsigned getiflavor() const;
	unsigned getjflavor() const;
	coordinate getxi() const;
	coordinate getxj() const;
	std::list<vertex>::iterator getconnected_vi();
	std::list<vertex>::iterator getconnected_vj();
	std::list<fermionicloop>::iterator getconnected_loopi();
	//std::list<fermionicloop>::iterator getconnected_loopj();
	double getvalueij();
	double getvalueji();
	friend std::ostream& operator << (std::ostream& ostr, const line& l);
};*/

/*class bubble{
	std::list<line>::iterator Gs[2];
	double value;

	public:
	bubble(std::list<line>::iterator Gi_iter, std::list<line>::iterator Gj_iter);
	std::list<vertex>::iterator getconnected_vi();
	std::list<vertex>::iterator getconnected_vj();
	std::list<line>::iterator getGi();
	std::list<line>::iterator getGj();
	double getvalue();
	friend std::ostream& operator << (std::ostream& ostr, const bubble& v);
};*/

/*class fermionicloop{
	std::list<std::list<vertex>::iterator> vs;

	double valuerh, valuelh;
	//std::vector<std::list<line>::iterator> UWs_intra;
	//std::vector<std::list<line>::iterator> UWs_inter;

	//std::vector<vertex> vs_modified;
	//std::vector<line> UWs_modified_intra;
	//std::vector<line> UWs_modified_inter;

	//std::vector<std::list<line>::iterator> Giters;
	//std::vector<line> Gs_modified;

	//double loopflip_proposal;
	//double loopflip_probability;

	public:
	fermionicloop(std::list<vertex>::iterator v_iter);
	fermionicloop(std::list<std::list<vertex>::iterator> vs_i);
	void set(std::list<std::list<vertex>::iterator> vs_i);
	void setvalue(const double& valuerh_i, const double& valuelh_i);
	void calvalue();
	//std::vector<std::list<line>::iterator> getUWs_intra();
	//std::vector<std::list<line>::iterator> getUWs_inter();
	void push_back(const std::list<vertex>::iterator v_iter);
	void insert(const std::list<vertex>::iterator& v_pos, const std::list<vertex>::iterator& v_iter);
	void erase(const std::list<vertex>::iterator& v_iter);
	int getvssize();
	std::list<std::list<vertex>::iterator> getvs();
	std::list<std::list<vertex>::iterator> getpath(const std::list<vertex>::iterator& vi, const std::list<vertex>::iterator& vj);
	double getvaluerh();
	double getvaluelh();
	//void proposeloopspinflip();
	//void proposeloopspinflip(const std::vector<coordinate>& xs, const std::vector<int>& UW_types);
	//std::vector<coordinate> getxvector();
	//std::vector<int> getUWtypes();
	//double get_loopspinflip_probability();
	//void loopspinflip();
	friend std::ostream& operator << (std::ostream& ostr, const fermionicloop& fl);
};*/

/*class setUWinter{
	int Nfermionicloops;
	std::vector<bool> isloopflavorfixed;
	std::vector<int> loopflavor;

	std::vector<line> UWinters;
	std::vector<Eigen::MatrixXd> UWvalues;
	std::vector<int> iloops, jloops;

	double value;

	public:
	setUWinter();
	setUWinter(const int& Nfermionicloops_i);
	setUWinter(const setUWinter& copy_i);
	void setNfermionicloops(const int& Nfermionicloops_i);
	void clear();
	void push_back(const line& UWinter_i, const int& iloop, const int& jloop);
	void erase(const line& line_i);
	void setphysical(const line& line_i, const bool& physical_i);
	void transform(const line& line_o, const line& line_n);
	void releaseloopflavor();
	void fixloopflavor(const int& iloop, const int& flavor_i);
	void setisloopflavorfixed(const std::vector<bool>& isloopflavorfixed_n);
	void setloopflavor(const std::vector<int>& loopflavor_n);
	bool getisloopflavorfixed(const int& iloop);
	int getloopflavor(const int& iloop);
	void calvalue();
	double calvalueafteraddvertices(const line& line_n, const int& iloop_n, const int& jloop_n);
	double calvalueafterremovevertices(const line& line_o);
	double calvalueafterinterswap(const line& line_o, const std::vector<bool>& isloopflavorfixed_n, std::vector<int> loopflavor_n);
	double calvalueafterinterswap(const std::vector<bool>& isloopflavorfixed_n, std::vector<int> loopflavor_n);
	double calvalueafterintraswap(const line& line_o1, const line& line_o2, const std::vector<bool>& isloopflavorfixed_n, std::vector<int> loopflavor_n);
	double calvalueaftertransformUW(const line& line_o, const line& line_n);
	void setvalue(const double& val);
	double getvalue();
	friend std::ostream& operator << (std::ostream& ostr, const setUWinter& SUI);
};*/
