#include <fstream>
#include "tecplot_export_grid2d.hpp"

namespace{
std::string default_bfun(int i){
	return std::string("boundary") + std::to_string(i);
}
}

void GGeom::Export::GridTecplot(const GridGeom& g, std::string fn, vector<int> bndindex, BFun bnames){
	//prepare data
	//edges
	std::set<Edge> eds_set = g.get_edges();
	vector<Edge> eds(eds_set.begin(), eds_set.end());
	//bzones
	std::map<int, vector<Edge*>> bzones;
	for (int i=0; i<eds.size(); ++i) if (bndindex[i]>=0){
		int b = bndindex[i];
		auto emp = bzones.emplace(b, vector<Edge*>());
		emp.first->second.push_back(&eds[i]);
	}

	//write to file
	std::ofstream of(fn);
	of.precision(10);
	//main header
	of<<"TITLE=\"Tecplot Export\""<<std::endl;
	of<<"VARIABLES=\"X\" \"Y\""<<std::endl;
	of<<"ZONE T=\"Grid\""<<std::endl;
	of<<"Nodes="<<g.n_points()<<std::endl;
	of<<"Faces="<<eds.size()<<std::endl;
	of<<"Elements="<<g.n_cells()<<std::endl;
	of<<"ZONETYPE=FEPOLYGON"<<std::endl;
	of<<"DATAPACKING=BLOCK"<<std::endl;
	of<<"NumConnectedBoundaryFaces=0, TotalNumBoundaryConnections=0"<<std::endl;
	//main points
	for (int i=0; i<g.n_points(); ++i) of<<g.get_point(i)->x<<std::endl;
	for (int i=0; i<g.n_points(); ++i) of<<g.get_point(i)->y<<std::endl;
	//main edges
	for (int i=0; i<eds.size(); ++i) of<<eds[i].p1+1<<" "<<eds[i].p2+1<<std::endl;
	//main adjacents
	for (int i=0; i<eds.size(); ++i) of<<eds[i].cell_left+1<<std::endl;
	for (int i=0; i<eds.size(); ++i) of<<eds[i].cell_right+1<<std::endl;

	//boundaries
	for (auto& v: bzones){
		std::string name = bnames(v.first);
		of<<"ZONE T=\""<<name<<"\""<<std::endl;
		of<<"D=(1, 2)"<<std::endl;
		of<<"E="<<v.second.size()<<std::endl;
		of<<"ZONETYPE=FELINESEG"<<std::endl;
		for (auto e: v.second) of<<e->p1 + 1<<" "<<e->p2 + 1<<std::endl;
	}

	of.close();
}


void GGeom::Export::GridTecplot(const GridGeom& g, std::string fn, vector<int> bndindex){
	return GridTecplot(g, fn, bndindex, default_bfun);
}
