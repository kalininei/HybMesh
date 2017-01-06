#include <fstream>
#include "vtk_export_grid2d.hpp"
#include "procgrid.h"
#include "vtk_export2d.hpp"

namespace hme = GGeom::Export;
namespace{

template<class C>
void add_vtk_data(const vector<C>& data, std::string name, std::string fn, bool is_first = false){
	std::ofstream f(fn, std::ios_base::app);
	if (is_first) f<<"CELL_DATA "<<data.size()<<std::endl;
	std::string tp = (std::is_same<C, int>::value) ? " int 1" : " float 1";
	f<<"SCALARS "<<name<<tp<<std::endl;
	f<<"LOOKUP_TABLE default"<<std::endl;
	for (auto v: data) f<<v<<std::endl;
	f.close();
}

}

void hme::GridVTK(const GridGeom& g, std::string fn){
	std::ofstream fs(fn);
	fs<<"# vtk DataFile Version 3.0"<<std::endl;
	fs<<"HybMesh Grid 2D"<<std::endl;
	fs<<"ASCII"<<std::endl;
	//Points
	fs<<"DATASET UNSTRUCTURED_GRID"<<std::endl;
	fs<<"POINTS "<<g.n_points()<< " float"<<std::endl;
	for (int i=0;i<g.n_points();++i){
		auto p = g.get_point(i);
		fs<<p->x<<" "<<p->y<<" 0"<<std::endl;
	}
	//Cells
	fs<<"CELLS  "<<g.n_cells()<<"   "<<g.n_cells()+g.n_cellsdim()<<std::endl;
	for (int i=0;i<g.n_cells();++i){
		auto c = g.get_cell(i);
		fs<<c->dim()<<"  ";
		for (int j=0;j<c->dim();++j)
			fs<<c->get_point(j)->get_ind()<<" ";
		fs<<std::endl;
	}
	fs<<"CELL_TYPES  "<<g.n_cells()<<std::endl;
	for (int i=0;i<g.n_cells();++i) fs<<7<<std::endl;
	fs.close();
}

void hme::BoundaryVTK(const GridGeom& g, std::string fn, const vector<int>& bcond){
	//save contour
	auto ct = GGeom::Info::Contour(g).alledges();
	HM2D::Export::ContourVTK(ct, fn.c_str());
	if (bcond.size() == 0) return;

	//add boundary condition
	//assemble
	std::map<Edge, int> eind; 
	int i=0;
	for (auto& e: g.get_edges()) eind.emplace(e, i++);
	vector<int> convertedbc(eind.size(), 0);

	std::set<GridPoint> bp;
	for (auto it: g.get_bnd_points()) bp.insert(*it);

	i = 0;
	for (auto e: ct){
		int ind1 = bp.find(*e->first())->get_ind();
		int ind2 = bp.find(*e->last())->get_ind();
		//int ind1 = static_cast<GridPoint*>(e->pstart)->get_ind();
		//int ind2 = static_cast<GridPoint*>(e->pend)->get_ind();
		Edge fedge(ind1, ind2);
		auto fnd = eind.find(fedge);
		assert(fnd!=eind.end());
		if (fnd != eind.end()) convertedbc[i] = bcond[fnd->second];
		++i;
	}
	convertedbc.resize(i);

	//add to output file
	add_vtk_data(convertedbc, "boundary_type", fn, true);
}



