#include "export2d_tecplot.hpp"
#include <fstream>

using namespace HM2D;

namespace{
std::string default_bfun(int i){
	return std::string("boundary") + std::to_string(i);
}
}

void Export::GridTecplot(const GridData& g, std::string fn, vector<int> bndindex, BNamesFun bnames){
	//prepare data
	//edges
	auto& eds = g.vedges;
	//bzones
	std::map<int, vector<Edge*>> bzones;
	for (int i=0; i<eds.size(); ++i) if (bndindex[i]>=0){
		int b = bndindex[i];
		auto emp = bzones.emplace(b, vector<Edge*>());
		emp.first->second.push_back(eds[i].get());
	}

	//write to file
	std::ofstream of(fn);
	of.precision(10);
	//main header
	of<<"TITLE=\"Tecplot Export\""<<std::endl;
	of<<"VARIABLES=\"X\" \"Y\""<<std::endl;
	of<<"ZONE T=\"Grid\""<<std::endl;
	of<<"Nodes="<<g.vvert.size()<<std::endl;
	of<<"Faces="<<eds.size()<<std::endl;
	of<<"Elements="<<g.vcells.size()<<std::endl;
	of<<"ZONETYPE=FEPOLYGON"<<std::endl;
	of<<"DATAPACKING=BLOCK"<<std::endl;
	of<<"NumConnectedBoundaryFaces=0, TotalNumBoundaryConnections=0"<<std::endl;
	//main points
	for (int i=0; i<g.vvert.size(); ++i) of<<g.vvert[i]->x<<std::endl;
	for (int i=0; i<g.vvert.size(); ++i) of<<g.vvert[i]->y<<std::endl;
	//main edges
	aa::enumerate_ids_pvec(g.vvert);
	aa::enumerate_ids_pvec(g.vcells);
	for (int i=0; i<eds.size(); ++i) of<<eds[i]->first()->id+1<<" "<<eds[i]->last()->id+1<<std::endl;
	//main adjacents
	for (int i=0; i<eds.size(); ++i) {
		int lid = (eds[i]->has_left_cell()) ? eds[i]->left.lock()->id + 1 : 0;
		of<<lid<<std::endl;
	}
	for (int i=0; i<eds.size(); ++i){
		int rid = (eds[i]->has_right_cell()) ? eds[i]->right.lock()->id + 1 : 0;
		of<<rid<<std::endl;
	}

	//boundaries
	for (auto& v: bzones){
		std::string name = bnames(v.first);
		of<<"ZONE T=\""<<name<<"\""<<std::endl;
		of<<"D=(1, 2)"<<std::endl;
		of<<"E="<<v.second.size()<<std::endl;
		of<<"ZONETYPE=FELINESEG"<<std::endl;
		for (auto e: v.second) of<<e->first()->id + 1<<" "<<e->last()->id + 1<<std::endl;
	}

	of.close();
}


void Export::GridTecplot(const GridData& g, std::string fn, vector<int> bndindex){
	return GridTecplot(g, fn, bndindex, default_bfun);
}
