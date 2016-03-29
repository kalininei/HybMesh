#ifndef NDEBUG

#include "debug_grid3d.hpp"
#include "stdarg.h"

using namespace HMGrid3D;

int Debug::tabs = 0;

void Debug::Print(const char* fmt, ...){
	std::string s;
	for (int i=0; i<tabs; ++i) s+="        ";
	s+=fmt;
	va_list va;
	const char* format = s.c_str();
	va_start(va, fmt);
	vprintf(format, va);
	va_end(va);
}

std::ostream& Debug::Cout(){
	for (int i=0; i<tabs; ++i) std::cout<<"        ";
	return std::cout;
}

void Debug::info_gridedge(const Edge& e){
	Print("+++ Edge at %p. Length = %10.6f. %i vertices.\n", &e, e.length(), e.vertices.size());
	tabs+=1;
	for (auto p: e.vertices){
		Print("Vertex at %p: (%10.6f, %10.6f, %10.6f)\n", p.get(), p->x, p->y, p->z);
	}
	tabs-=1;
	Cout()<<"+++++++++++++++++++++++++++++++++++++"<<std::endl;
}

void Debug::info_gridedge(const Grid& g, int edge_index){
	auto alled = g.alledges();
	info_gridedge(*alled[edge_index]);
}

void Debug::info_gridface(const Face& f){
	Print("+++ Face at %p. %i Edges\n", &f, f.edges.size());
	Print("+++ Left cell: %p; Right cell: %p\n", f.left.get(), f.right.get());
	tabs+=1;
	for (auto e: f.edges) info_gridedge(*e);
	tabs-=1;
	Cout()<<"+++++++++++++++++++++++++++++++++++++"<<std::endl;
}

void Debug::info_gridface(const Grid& grid, int face_index){
	auto allfc = grid.allfaces();
	info_gridface(*allfc[face_index]);
}

void Debug::info_gridcell(const Cell& c){
	Print("+++ Cell at %p. %i Faces\n", &c, c.faces.size());
	tabs+=1;
	for (auto e: c.faces) info_gridface(*e);
	tabs-=1;
	Cout()<<"+++++++++++++++++++++++++++++++++++++"<<std::endl;
}
void Debug::info_gridcell(const Grid& grid, int cell_index){
	info_gridcell(*grid.cells[cell_index]);
}

#endif
