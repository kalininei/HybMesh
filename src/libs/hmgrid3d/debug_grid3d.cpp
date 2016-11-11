#ifndef NDEBUG

#include "debug_grid3d.hpp"
#include "stdarg.h"
#include "hmgrid3d.hpp"
#include <fstream>

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
	Export::BoundaryVTK(grid, "_dbgout.vtk");
}
void Debug::save_bnd_vtk(const GridData& grid){
	Export::BoundaryVTK(grid, "_dbgout.vtk");
}

void Debug::save_vtk(const SGrid& grid, const char* fn){
	Export::GridVTK(grid, fn);
	std::ofstream fs(fn, std::ios_base::app);
	//add x, y, z fields to make cuts easier
	fs<<"POINT_DATA "<<grid.n_vert<<std::endl;
	fs<<"SCALARS xcoords float "<<std::endl;
	fs<<"LOOKUP_TABLE default"<<std::endl;
	for (int i=0; i<grid.n_vert; ++i) fs<<grid.vert[3*i]<<std::endl;
	fs<<"SCALARS ycoords float "<<std::endl;
	fs<<"LOOKUP_TABLE default"<<std::endl;
	for (int i=0; i<grid.n_vert; ++i) fs<<grid.vert[3*i+1]<<std::endl;
	fs<<"SCALARS zcoords float "<<std::endl;
	fs<<"LOOKUP_TABLE default"<<std::endl;
	for (int i=0; i<grid.n_vert; ++i) fs<<grid.vert[3*i+2]<<std::endl;
}
void Debug::save_vtk(const GridData& grid, const char* fn){
	SGrid sg(grid);
	save_vtk(sg, fn);
}

void Debug::save_vtk(const SGrid& grid){ save_vtk(grid, "_dbgout.vtk"); }
void Debug::save_vtk(const GridData& grid){ save_vtk(grid, "_dbgout.vtk"); }


void Debug::save_vtk(const Surface& srf, const char* fn){
	Surface cubic = srf;
	ReallocateAll(cubic.faces);

	SGrid g;
	g.vcells.emplace_back(new HMGrid3D::Cell());
	shared_ptr<HMGrid3D::Cell> c0 = g.vcells[0];
	for (auto f: cubic.faces) {f->left = c0; f->right.reset();}
	g.vcells[0]->faces = cubic.faces;
	g.fill_from_cells();
	g.actualize_serial_data();
	HMGrid3D::Export::BoundaryVTK(g, fn);
}

void Debug::save_vtk(const Surface& srf){ save_vtk(srf, "_dbgout.vtk"); }

#endif
