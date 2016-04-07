#ifndef NDEBUG

#include "debug_grid2d.h"
#include <fstream>
#include "vtk_export_grid2d.hpp"

using namespace GGeom;

double Debug::hash(const GridGeom& grid){
	double sum;
	for (int i=0; i<grid.n_cells(); ++i){
		auto c = grid.get_cell(i);
		for (int j=0; j<c->dim(); ++j){
			auto vert = c->get_point(j);
			sum += 0.333*sin(i*vert->x+3);
			sum -= 0.333*cos(i*vert->y+4);
		}
	}
	return sum;
}

namespace{

void add_vtk_point_data(const vector<double>& data, const char* name, const char* fn, bool is_first = false){
	std::ofstream f(fn, std::ios_base::app);
	if (is_first) f<<"POINT_DATA "<<data.size()<<std::endl;
	f<<"SCALARS "<<name<<" float 1"<<std::endl;
	f<<"LOOKUP_TABLE default"<<std::endl;
	for (auto v: data) f<<v<<std::endl;
	f.close();
}
void add_vtk_cell_data(const vector<double>& data, const char* name, const char* fn, bool is_first = false){
	std::ofstream f(fn, std::ios_base::app);
	if (is_first) f<<"CELL_DATA "<<data.size()<<std::endl;
	f<<"SCALARS "<<name<<" float 1"<<std::endl;
	f<<"LOOKUP_TABLE default"<<std::endl;
	for (auto v: data) f<<v<<std::endl;
	f.close();
}

}

void Debug::save_vtk(const GridGeom* g, const char* fn){
	GGeom::Export::GridVTK(*g, fn);
}

void Debug::save_vtk(const GridGeom* g, const vector<double>& data, const char* fn){
	Debug::save_vtk(g, fn);
	if (data.size() == g->n_cells()) add_vtk_cell_data(data, "data", fn, true);
	else add_vtk_point_data(data, "data", fn, true);
}

void Debug::save_vtk(const PContour* c, const char* fn){
	std::ofstream fs(fn);
	fs<<"# vtk DataFile Version 3.0"<<std::endl;
	fs<<"Contour 2D"<<std::endl;
	fs<<"ASCII"<<std::endl;
	//Points
	fs<<"DATASET UNSTRUCTURED_GRID"<<std::endl;
	fs<<"POINTS "<<c->n_points()<< " float"<<std::endl;
	for (int i=0;i<c->n_points();++i){
		auto p = c->get_point(i);
		fs<<p->x<<" "<<p->y<<" 0"<<std::endl;
	}
	//Cells
	fs<<"CELLS  "<<c->n_points()<<"   "<<3*c->n_points()<<std::endl;
	for (int i=0;i<c->n_points();++i){
		int inext = (i+1) % c->n_points();
		fs<<2<<" "<<i<<" "<<inext<<std::endl;
	}
	fs<<"CELL_TYPES  "<<c->n_points()<<std::endl;
	for (int i=0;i<c->n_points();++i) fs<<3<<std::endl;
	fs.close();
}

void Debug::save_vtk(const vector<PContour>& c, const char* fn){
	std::ofstream fs(fn);
	fs<<"# vtk DataFile Version 3.0"<<std::endl;
	fs<<"Contour 2D"<<std::endl;
	fs<<"ASCII"<<std::endl;
	//Points
	int numpts = 0;
	for (auto cont: c) numpts+=cont.n_points();
	fs<<"DATASET UNSTRUCTURED_GRID"<<std::endl;
	fs<<"POINTS "<<numpts<< " float"<<std::endl;
	for (auto cont: c){
		for (int i=0;i<cont.n_points();++i){
			auto p = cont.get_point(i);
			fs<<p->x<<" "<<p->y<<" 0"<<std::endl;
		}
	}
	//Cells
	fs<<"CELLS  "<<numpts<<"   "<<3*numpts<<std::endl;
	int globN = 0;
	for (auto cont: c){
		for (int i=0;i<cont.n_points();++i){
			int inext = (i+1) % cont.n_points();
			fs<<2<<" "<<i+globN<<" "<<inext+globN<<std::endl;
		}
		globN+=cont.n_points();
	}
	fs<<"CELL_TYPES  "<<numpts<<std::endl;
	for (int i=0;i<numpts;++i) fs<<3<<std::endl;
	fs.close();
	std::vector<double> cont_id;
	int id = 0;
	for (auto cont: c){
		for (int i=0;i<cont.n_points();++i) cont_id.push_back(id);
		++id;
	}
	add_vtk_point_data(cont_id, "Contours_ID", fn, true);
}

void Debug::save_vtk(const vector<PContour>& c, const vector<double>& data, const char* fn){
	save_vtk(c, fn);
	add_vtk_point_data(data, "custom_data", fn);
}

void Debug::save_vtk(const PtsGraph* g, const char* fn){
	std::ofstream fs(fn);
	fs<<"# vtk DataFile Version 3.0"<<std::endl;
	fs<<"PtsGraph"<<std::endl;
	fs<<"ASCII"<<std::endl;
	//Points
	fs<<"DATASET UNSTRUCTURED_GRID"<<std::endl;
	fs<<"POINTS "<<g->Nnodes()<< " float"<<std::endl;
	for (int i=0;i<g->Nnodes();++i){
		auto p = g->get_point(i);
		fs<<p->x<<" "<<p->y<<" 0"<<std::endl;
	}
	//Cells
	fs<<"CELLS  "<<g->Nlines()<<"   "<<3*g->Nlines()<<std::endl;
	for (int i=0;i<g->Nlines();++i){
		auto ln = g->get_line(i);
		fs<<2<<" "<<ln.first<<" "<<ln.second<<std::endl;
	}
	fs<<"CELL_TYPES  "<<g->Nlines()<<std::endl;
	for (int i=0;i<g->Nlines();++i) fs<<3<<std::endl;
	fs.close();
}

void Debug::save_vtk(const ContoursCollection& c, const char* fn){ save_vtk(c.contours_list(), fn); }
void Debug::save_vtk(const ContoursCollection* c, const char* fn){ save_vtk(c->contours_list(), fn); }
void Debug::save_vtk(const ContoursCollection& c, const vector<double>& pdata, const char* fn){
	save_vtk(c.contours_list(), pdata, fn);
}
void Debug::save_vtk(const GridGeom& g, const char* fn){ Debug::save_vtk(&g, fn);}
void Debug::save_vtk(const GridGeom& g, const vector<double>& data, const char* fn){ save_vtk(&g, data, fn);}
void Debug::save_vtk(const PtsGraph& g, const char* fn){ save_vtk(&g, fn);}


#endif
