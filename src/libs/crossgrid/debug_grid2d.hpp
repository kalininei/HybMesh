#ifndef HYBMESH_GRID2D_DEBUG_HPP
#define HYBMESH_GRID2D_DEBUG_HPP

#ifndef NDEBUG

#include "primitives2d.hpp"
#include "debug2d.hpp"
#include "wireframegrid.hpp"
#include "finder2d.hpp"

namespace HM2D{ namespace Grid{

struct Debug: public HMDebug{
	static void report_grid_problems(const GridData& g);
	static void report_grid_problems(const GridData& g, double maxskew);

	static void save_wf_vtk(const Grid::Impl::PtsGraph& c);

	static void save_bbfinder_vtk(const BoundingBoxFinder& bf);
	static void save_bbfinder_vtk(const BoundingBoxFinder& bf, const vector<int>& dt);
	static void save_raster_vtk(const Finder::RasterizeEdges& rs);
	static void save_raster_vtk(shared_ptr<Finder::RasterizeEdges> rs);

};

}}

#endif
#endif
