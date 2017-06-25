#include "export2d_vtk.hpp"
#include <fstream>
#include "contour.hpp"
#include "contour_tree.hpp"

using namespace HM2D;

namespace{

template<class C>
void add_vtk_data(const vector<C>& data, std::string name, std::string fn, bool is_cell, bool is_first){
	std::ofstream f(fn, std::ios_base::app);
	if (is_first && is_cell) f<<"CELL_DATA "<<data.size()<<std::endl;
	if (is_first && !is_cell) f<<"POINT_DATA "<<data.size()<<std::endl;
	std::string tp = (std::is_same<C, int>::value) ? " int 1" : " float 1";
	f<<"SCALARS "<<name<<tp<<std::endl;
	f<<"LOOKUP_TABLE default"<<std::endl;
	for (auto v: data) f<<v<<std::endl;
	f.close();
}

}


void Export::ContourVTK(const EdgeData& ec, std::string fn){
	std::ofstream fs(fn);
	fs<<"# vtk DataFile Version 3.0"<<std::endl;
	fs<<"EdgesData"<<std::endl;
	fs<<"ASCII"<<std::endl;
	VertexData pts = AllVertices(ec);
	aa::enumerate_ids_pvec(pts);
	//Points
	fs<<"DATASET UNSTRUCTURED_GRID"<<std::endl;
	fs<<"POINTS "<<pts.size()<< " float"<<std::endl;
	for (int i=0; i<pts.size(); ++i){
		fs<<float(pts[i]->x)<<" "<<float(pts[i]->y)<<" 0"<<std::endl;
	}
	//Cells
	fs<<"CELLS  "<<ec.size()<<"   "<<3*ec.size()<<std::endl;
	for (int i=0; i<ec.size(); ++i){
		fs<<2<<" "<<ec[i]->first()->id<<" "<<ec[i]->last()->id<<std::endl;
	}
	fs<<"CELL_TYPES  "<<ec.size()<<std::endl;
	for (int i=0;i<ec.size();++i) fs<<3<<std::endl;

	//boundary conditions if any
	vector<int> bcond(ec.size());
	for (int i=0; i<ec.size(); ++i) bcond[i] = ec[i]->boundary_type;
	if (std::all_of(bcond.begin(), bcond.end(), [](int a){ return a==0; })){
		return;
	}
	fs<<"CELL_DATA "<<bcond.size()<<std::endl;
	fs<<"SCALARS boundary_type int 1"<<std::endl;
	fs<<"LOOKUP_TABLE default"<<std::endl;
	for (auto v: bcond) fs<<v<<std::endl;

	fs.close();
}

void Export::VerticesVTK(const VertexData& data, std::string fn){
	std::ofstream fs(fn);
	fs<<"# vtk DataFile Version 3.0"<<std::endl;
	fs<<"VertexData"<<std::endl;
	fs<<"ASCII"<<std::endl;
	//Points
	fs<<"DATASET UNSTRUCTURED_GRID"<<std::endl;
	fs<<"POINTS "<<data.size()<< " float"<<std::endl;
	for (int i=0; i<data.size(); ++i){
		fs<<float(data[i]->x)<<" "<<float(data[i]->y)<<" 0"<<std::endl;
	}
	//Cells
	fs<<"CELLS  "<<data.size()<<"   "<<2*data.size()<<std::endl;
	for (int i=0; i<data.size(); ++i){
		fs<<"1 "<<i<<std::endl;
	}
	fs<<"CELL_TYPES  "<<data.size()<<std::endl;
	for (int i=0;i<data.size();++i) fs<<1<<std::endl;
	fs.close();
}

void Export::GridVTK(const GridData& g, std::string fn){
	std::ofstream fs(fn);
	fs<<"# vtk DataFile Version 3.0"<<std::endl;
	fs<<"HybMesh Grid 2D"<<std::endl;
	fs<<"ASCII"<<std::endl;
	//Points
	fs<<"DATASET UNSTRUCTURED_GRID"<<std::endl;
	fs<<"POINTS "<<g.vvert.size()<< " float"<<std::endl;
	for (int i=0;i<g.vvert.size();++i){
		auto p = g.vvert[i];
		fs<<(float)p->x<<" "<<(float)p->y<<" 0"<<std::endl;
	}
	//Cells
	int totc = g.vcells.size();
	for (auto c: g.vcells) totc+=c->edges.size();
	fs<<"CELLS  "<<g.vcells.size()<<"   "<<totc<<std::endl;
	aa::enumerate_ids_pvec(g.vvert);
	for (int i=0;i<g.vcells.size();++i){
		auto c = g.vcells[i];
		fs<<c->edges.size()<<"  ";
		auto op = HM2D::Contour::OrderedPoints(c->edges);
		op.resize(op.size()-1);
		for (auto p: op){ fs<<p->id<<" "; }
		fs<<std::endl;
	}
	fs<<"CELL_TYPES  "<<g.vcells.size()<<std::endl;
	for (int i=0;i<g.vcells.size();++i) fs<<7<<std::endl;
	fs.close();
}

void Export::GridCDataVTK(const GridData& g, const vector<double>& dt, std::string fn){
	assert(dt.size() == g.vcells.size());
	GridVTK(g, fn);
	add_vtk_data(dt, "cell_data", fn, true, true);
}
void Export::GridVDataVTK(const GridData& g, const vector<double>& dt, std::string fn){
	assert(dt.size() == g.vvert.size());
	GridVTK(g, fn);
	add_vtk_data(dt, "vertex_data", fn, false, true);
}

void Export::BoundaryVTK(const GridData& g, std::string fn){
	//save contour
	auto ct = HM2D::Contour::Tree::GridBoundary(g).alledges();
	HM2D::Export::ContourVTK(ct, fn.c_str());
	/*
	vector<int> bcond(ct.size());
	for (int i=0; i<ct.size(); ++i) bcond[i] = ct[i]->boundary_type;
	if (std::all_of(bcond.begin(), bcond.end(), [](int a){ return a==0; })){
		return;
	}
	aa::constant_ids_pvec(g.vedges, 0);
	for (int i=0; i<bcond.size(); ++i){
		if (i==g.vedges.size()) break;
		if (!g.vedges[i]->is_boundary()) continue;
		g.vedges[i]->id = bcond[i];
	}
	vector<int> convertedbc(ct.size(), 0);
	for (int i=0; i<ct.size(); ++i){
		convertedbc[i]=ct[i]->id;
	}

	//add to output file
	add_vtk_data(convertedbc, "boundary_type", fn, true);
	*/
}
