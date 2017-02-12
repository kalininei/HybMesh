#ifndef HYBMESH_GRID2D_DEBUG_HPP
#define HYBMESH_GRID2D_DEBUG_HPP

#ifndef NDEBUG

#include "primitives2d.hpp"
#include "debug2d.hpp"
#include "wireframegrid.hpp"

namespace HM2D{ namespace Grid{

struct Debug: public HMDebug{
	static void save_wf_vtk(const Grid::Impl::PtsGraph& c);

	static void save_bbfinder_vtk(const BoundingBoxFinder& bf);
	static void save_bbfinder_vtk(const BoundingBoxFinder& bf, const vector<int>& dt);
};

}}

#endif
#endif
