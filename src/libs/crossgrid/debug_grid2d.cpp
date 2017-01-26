#ifndef NDEBUG

#include "debug_grid2d.hpp"
#include "debug2d.hpp"
#include "buildgrid.hpp"
#include <fstream>

using namespace HM2D;

void HM2D::Grid::Debug::save_wf_vtk(const Grid::Impl::PtsGraph& c){
	auto ecol = c.toedges();
	HM2D::Debug::save_edges_vtk(ecol);

}

void HM2D::Grid::Debug::save_bbfinder_vtk(const BoundingBoxFinder& bf){
	Point p1 = bf.pmin();
	Point p2 = bf.pmax();
	auto rg = Grid::Constructor::RectGrid(p1, p2, bf.nx(), bf.ny());
	HM2D::Debug::save_grid_vtk(rg);
	std::ofstream fs("_dbgout.vtk", std::ios::app);
	fs<<"CELL_DATA "<<bf.nsqr()<<std::endl;
	fs<<"SCALARS entries_num int 1"<<std::endl;
	fs<<"LOOKUP_TABLE default"<<std::endl;
	for (int i=0; i<bf.nsqr(); ++i) fs<<bf.sqr_entries(i).size()<<std::endl;
}
void HM2D::Grid::Debug::save_bbfinder_vtk(const BoundingBoxFinder& bf, const vector<int>& dt){
	save_bbfinder_vtk(bf);
	std::ofstream fs("_dbgout.vtk", std::ios::app);
	fs<<"SCALARS data int 1"<<std::endl;
	fs<<"LOOKUP_TABLE default"<<std::endl;
	for (auto a: dt) fs<<a<<std::endl;
}


#endif
