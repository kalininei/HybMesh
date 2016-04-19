#ifndef NDEBUG

#include "debug_grid3d.hpp"
#include "stdarg.h"
#include "hmgrid3d.hpp"

using namespace HMGrid3D;

void Debug::info_gridedge(const Edge& e){
	Print("+++ Edge at %p. Length = %10.6f. %i vertices. id=%i\n", &e, e.length(), e.vertices.size(), e.id);
	tabs+=1;
	for (auto p: e.vertices){
		Print("Vertex at %p: (%10.6f, %10.6f, %10.6f). id=%i\n", p.get(), p->x, p->y, p->z, p->id);
	}
	tabs-=1;
	Cout()<<"+++++++++++++++++++++++++++++++++++++"<<std::endl;
}

void Debug::info_gridedge(const SGrid& g, int edge_index){
	info_gridedge(*g.vedges[edge_index]);
}

void Debug::info_gridface(const Face& f){
	Print("+++ Face at %p. %i Edges. id=%i\n", &f, f.edges.size(), f.id);
	Print("+++ Left cell: %p; Right cell: %p\n", f.left.lock().get(), f.right.lock().get());
	tabs+=1;
	for (auto e: f.edges) info_gridedge(*e);
	tabs-=1;
	Cout()<<"+++++++++++++++++++++++++++++++++++++"<<std::endl;
}

void Debug::info_gridface(const SGrid& grid, int face_index){
	info_gridface(*grid.vfaces[face_index]);
}

void Debug::info_gridcell(const Cell& c){
	Print("+++ Cell at %p. %i Faces. id=%i\n", &c, c.faces.size(), c.id);
	tabs+=1;
	for (auto e: c.faces) info_gridface(*e);
	tabs-=1;
	Cout()<<"+++++++++++++++++++++++++++++++++++++"<<std::endl;
}
void Debug::info_gridcell(const SGrid& grid, int cell_index){
	info_gridcell(*grid.vcells[cell_index]);
}

void Debug::save_bnd_vtk(const SGrid& grid){
	HMGrid3D::Export::BoundaryVTK(grid, "_dbgout.vtk");
}


#endif
