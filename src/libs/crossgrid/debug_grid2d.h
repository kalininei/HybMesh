#ifndef HYBMESH_GRID2D_DEBUG_H
#define HYBMESH_GRID2D_DEBUG_H

#ifndef NDEBUG

#include "hmdebug.hpp"
#include "grid.h"
#include "bgeom2d.h"
#include "wireframegrid.h"

namespace GGeom{

struct Debug: public HMDebug{
	static double hash(const GridGeom& grid);

	template<class G>
	static void outhash(G&& grid, std::string prefix=""){
		if (prefix.size() > 0) prefix = " ("+prefix+")";
		std::cout<<"Hash for grid2d"<<prefix<<": "<<hash(grid)<<std::endl;
	}

	static void save_vtk(const GridGeom* g, const char* fn);
	static void save_vtk(const GridGeom* g, const vector<double>& pdata, const char* fn);
	static void save_vtk(const PContour* c, const char* fn);
	static void save_vtk(const vector<PContour>& c, const char* fn);
	static void save_vtk(const vector<PContour>& c, const vector<double>& data, const char* fn);
	static void save_vtk(const PtsGraph* g, const char* fn);

	static void save_vtk(const ContoursCollection& c, const char* fn);
	static void save_vtk(const ContoursCollection* c, const char* fn);
	static void save_vtk(const ContoursCollection& c, const vector<double>& pdata, const char* fn);
	static void save_vtk(const GridGeom& g, const char* fn);
	static void save_vtk(const GridGeom& g, const vector<double>& pdata, const char* fn);
	static void save_vtk(const PtsGraph& g, const char* fn);
};



}

#endif
#endif
