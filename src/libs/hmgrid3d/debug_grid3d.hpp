#ifndef HYBMESH_GRID3D_DEBUG_HPP
#define HYBMESH_GRID3D_DEBUG_HPP

#ifndef NDEBUG

#include "primitives_grid3d.hpp"

namespace HMGrid3D{

struct Debug{
	static void info_gridedge(const Grid& grid, int edge_index);
	static void info_gridface(const Grid& grid, int face_index);
	static void info_gridcell(const Grid& grid, int cell_index);
	static void info_gridedge(const Edge& edge);
	static void info_gridface(const Face& face);
	static void info_gridcell(const Cell& cell);
	

	//printing with tabulation
	static int tabs;
	static void Print(const char* fmt, ...);
	static std::ostream& Cout();

};

}

#endif
#endif

