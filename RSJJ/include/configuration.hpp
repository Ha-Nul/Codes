#pragma once

#include <vector>
#include <iostream>
#include <Eigen/Dense>
#include "element.hpp"
#include "green.hpp"
//#include "graph.hpp"

class configuration{
	static const int Nflavor = 2;
	static const int NvExt = 4;

	static const Eigen::Matrix2d Identity2x2;
	static const Eigen::Matrix2d Zero2x2;

	interaction_t mWt;
	pseudogreen_t Gt;

	Eigen::Matrix2i extFlavor;
	std::vector<int> indexExtVertices;

	std::vector<coordinate> vertices;
	std::vector<double> Ws;
	std::vector<Eigen::Matrix2d> GtTop;
	std::vector<Eigen::Matrix2d> GtBottom;

	int im2;
	coordinate v_beta;

	std::vector<int> list_v_to_W;
	std::vector<int> list_W_to_vi;
	std::vector<int> list_W_to_vj;

	//int iv_measure, jv_reference;

	double WProduct;
	Eigen::Matrix2d GtTopProduct;
	Eigen::Matrix2d GtBottomProduct;

	// proposed quantity;
	bool is_on_top;
	bool is_vi_on_top, is_vj_on_top;

	int iW_new, iW_new2, iv_new, jv_new;
	coordinate vi_new, vj_new;
	double Wt_new, Wt_new2;
	int im2_new;


	double WProductProposed;
	Eigen::Matrix2i extFlavorProposed;
	Eigen::Matrix2d GtTopProductProposed;
	Eigen::Matrix2d GtBottomProductProposed;
	std::vector<coordinate> verticesProposed;
	std::vector<Eigen::Matrix2d> GtTopProposed;
	std::vector<Eigen::Matrix2d> GtBottomProposed;

	std::vector<int> list_v_to_W_proposed;
	std::vector<int> list_W_to_vi_proposed;
	std::vector<int> list_W_to_vj_proposed;

	//std::vector<double>::iterator Witer;
	//std::vector<Eigen::MatrixXd>::iterator SzGiter;

	//graph Wgraph, Ggraph;

	public:
	configuration();
	void init(const std::vector<Eigen::Matrix2d>& Gt_i);
	void init(const std::vector<double>& mWt_i,const std::vector<Eigen::Matrix2d>& Gt_i);
	int getVerticesSize();
	int getWsSize();
	int getGtTopSize();
	int getim2();
	int getiv(const int& iW);
	int getjv(const int& iW);
	int getiW(const int& iv);
	int getiv_partner(const int& iv);
	double get_tinterval(const int& iv_i);
	double getWProduct();
	coordinate getVertex(const int& iv_i);
	coordinate getProposedVi();
	coordinate getProposedVj();
	bool getIsOnTop();
	//int getWsSize();
	//int getGsSize();
	//void fillWs();
	//void fillPGtop();
	//void fillPGbottom();
	void calculate_WProduct();
	void calculate_GtTopProduct();
	void calculate_GtTopProductProposed();
	void calculate_GtBottomProduct();
	void calculate_GtBottomProductProposed();
	double getWeight();
	double getProposedWeight();
	void proposeInsertVertices(const int& iv, const int& jv);
	void proposeInsertVertices(const int& iv, const coordinate& vi_new_i, const int& jv, const coordinate& vj_new_i, const bool& is_on_top_i);
	void updateInsertVertices();
	void proposeRemoveVertices(const int& iW_i);
	void updateRemoveVertices();
	//void proposeSwapWline(const int& iv_i, const int& jv_i);
	//void updateSwapWline();
	void proposeShiftIm2(const int& im2_new_i);
	void updateShiftIm2();
	void proposeShift(const int& iv);
	void proposeShift(const int& iv,const double& t_n);
	void updateShift();
	void proposeFlip(const int& iv);
	void updateFlip();
	bool checkProposedIrreducibility();
	bool checkProposedCompactness();
	std::vector<coordinate> getExtVertices();
	std::vector<int> getFlavorExt();
	//void initProposal();
	////double getInsertionRatio();
	////double getRemovalRatio();
	////void copy( const configuration& Conf_i );
	bool checkWProduct();
	bool checkGProduct();
	bool checkGtSign();
	void print(std::ostream& ostr);
	void printProposed(std::ostream& ostr);
};
