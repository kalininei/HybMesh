#ifndef HYBMESH_GRID3D_DEBUG_HPP
#define HYBMESH_GRID3D_DEBUG_HPP

#ifndef NDEBUG

#include "primitives_grid3d.hpp"
#include "hmdebug.hpp"

namespace HMGrid3D{

struct Debug: public HMDebug{
	static void info_gridedge(const Grid& grid, int edge_index);
	static void info_gridface(const Grid& grid, int face_index);
	static void info_gridcell(const Grid& grid, int cell_index);
	static void info_gridedge(const Edge& edge);
	static void info_gridface(const Face& face);
	static void info_gridcell(const Cell& cell);
	
	template<class G>
	static double hash(const G& grid){
		double sum;
		int i=0;
		for (auto vert: grid.allvertices()){
			sum += 0.333*sin(i*vert->x+3);
			sum -= 0.333*cos(i*vert->y+4);
			sum += 0.333*sin(i*vert->z+5);
			++i;
		}
		return sum;
	}
	template<class G>
	static void outhash(const G& grid, std::string prefix=""){
		if (prefix.size() > 0) prefix = " ("+prefix+")";
		std::cout<<"Hash for grid3d"<<prefix<<": "<<hash(grid)<<std::endl;
	}
};

}

#endif
#endif

