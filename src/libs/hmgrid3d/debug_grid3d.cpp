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
	Print("+++ Left cell: %p (%i); Right cell: %p (%i)\n",
			f.left.lock().get(), (f.has_left_cell() ? f.left.lock()->id : -1),
			f.right.lock().get(), (f.has_right_cell() ? f.right.lock()->id : -1));
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

void Debug::save_grid_vtk(const GridData& grid){
	//save current enumeration
	RestoreIds<VertexData> s1(grid.vvert);
	RestoreIds<EdgeData> s2(grid.vedges);
	RestoreIds<FaceData> s3(grid.vfaces);
	RestoreIds<CellData> s4(grid.vcells);

	HMGrid3D::Export::GridVTK(grid, "_dbgout.vtk");
}
void Debug::save_cells_vtk(const CellData& grid){
	GridData g;
	g.vcells = grid;
	g.vfaces = all_faces(g.vcells);
	g.vedges = all_edges(g.vfaces);
	g.vvert = all_vertices(g.vedges);

	//temporary remove nonexisting cells from face-cell connectivity
	ShpVector<Cell> oldfc;
	for (auto f: g.vfaces){
		oldfc.push_back(f->left.lock());
		oldfc.push_back(f->right.lock());
		if (std::find(grid.begin(), grid.end(), f->left.lock()) == grid.end()) f->left.reset();
		if (std::find(grid.begin(), grid.end(), f->right.lock()) == grid.end()) f->right.reset();
	}
	//push single faced cells
	for (auto c: g.vcells) if (c->faces.size() == 1){
		auto fc = c->faces[0];
		if (!fc->has_left_cell()) fc->left = c;
		else fc->right = c;
	}

	save_grid_vtk(g);

	//place connectivity back
	auto it = oldfc.begin();
	for (auto f: g.vfaces){
		f->left = *it++;
		f->right = *it++;
	}
}

void Debug::save_cells_faces_vtk(const CellData& grid){
	FaceData fc = all_faces(grid);
	save_faces_vtk(fc);
}

void Debug::save_cell_vtk(shared_ptr<Cell> c){
	return save_cells_vtk(CellData{c});
}

void Debug::save_edges_vtk(const EdgeData& edges){
	VertexData vert = all_vertices(edges);
	RestoreIds<VertexData> s2(vert);
	enumerate_ids_pvec(vert);
	std::ofstream fs("_dbgout.vtk");
	fs<<"# vtk DataFile Version 3.0"<<std::endl;
	fs<<"3D Grid Edges"<<std::endl;
	fs<<"ASCII"<<std::endl;
	fs<<"DATASET UNSTRUCTURED_GRID"<<std::endl;
	fs<<"POINTS "<<vert.size()<< " float"<<std::endl;
	for (int i=0; i<vert.size(); ++i)
		fs<<vert[i]->x<<" "<<vert[i]->y<<" "<<vert[i]->z<<std::endl;

	int nffull = 0;
	for (auto e: edges) nffull += (1 + e->vertices.size());
	fs<<"CELLS  "<<edges.size()<<"   "<<nffull<<std::endl;
	for (auto e: edges){
		fs<<e->vertices.size();
		for (auto v: e->vertices) fs<<" "<<v->id;
		fs<<std::endl;
	}
	fs<<"CELL_TYPES  "<<edges.size()<<std::endl;
	for (auto& e: edges) fs<<'4'<<std::endl;
	fs.close();
}
void Debug::save_faces_vtk(const FaceData& faces){
	VertexData vert = all_vertices(faces);
	RestoreIds<VertexData> s2(vert);
	enumerate_ids_pvec(vert);
	std::ofstream fs("_dbgout.vtk");
	fs<<"# vtk DataFile Version 3.0"<<std::endl;
	fs<<"3D Grid Faces"<<std::endl;
	fs<<"ASCII"<<std::endl;
	fs<<"DATASET UNSTRUCTURED_GRID"<<std::endl;
	fs<<"POINTS "<<vert.size()<< " float"<<std::endl;
	for (int i=0; i<vert.size(); ++i)
		fs<<vert[i]->x<<" "<<vert[i]->y<<" "<<vert[i]->z<<std::endl;

	int nffull = 0;
	for (auto f: faces) nffull += (1 + f->edges.size());
	fs<<"CELLS  "<<faces.size()<<"   "<<nffull<<std::endl;
	for (auto f: faces){
		fs<<f->edges.size();
		auto av = f->sorted_vertices();
		for (auto v: av) fs<<" "<<v->id;
		fs<<std::endl;
	}
	fs<<"CELL_TYPES  "<<faces.size()<<std::endl;
	for (auto& f: faces) fs<<'7'<<std::endl;
	fs.close();
}

void Debug::save_surf_vtk(const Surface& srf){
	auto ae = srf.alledges();
	auto av = srf.allvertices();

	RestoreIds<VertexData> s1(av);
	RestoreIds<EdgeData> s2(ae);
	RestoreIds<FaceData> s3(srf.faces);

	HMGrid3D::Export::SurfaceVTK(srf, "_dbgout.vtk");
}

VertexData Debug::all_vertices(const CellData& a){
	return all_vertices(all_faces(a));
}
VertexData Debug::all_vertices(const FaceData& a){
	return all_vertices(all_edges(a));
}
VertexData Debug::all_vertices(const EdgeData& a){
	VertexData ret;
	for (auto edge: a)
	for (auto v: edge->vertices) ret.push_back(v);
	ret = aa::no_dublicates(ret);
	return ret;
}
EdgeData Debug::all_edges(const CellData& a){
	return all_edges(all_faces(a));
}
EdgeData Debug::all_edges(const FaceData& a){
	EdgeData ret;
	for (auto face: a)
	for (auto e: face->edges) ret.push_back(e);
	ret = aa::no_dublicates(ret);
	return ret;
}
FaceData Debug::all_faces(const CellData& a){
	FaceData ret;
	for (auto cell: a)
	for (auto f: cell->faces) ret.push_back(f);
	ret = aa::no_dublicates(ret);
	return ret;
}

#endif
