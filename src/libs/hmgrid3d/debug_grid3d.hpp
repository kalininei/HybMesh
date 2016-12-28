#ifndef HYBMESH_GRID3D_DEBUG_HPP
#define HYBMESH_GRID3D_DEBUG_HPP

#ifndef NDEBUG

#include "serialize_grid3d.hpp"
#include "hmdebug.hpp"

namespace HMGrid3D{

struct Debug: public HMDebug{
	static void info_gridedge(const SGrid& grid, int edge_index);
	static void info_gridface(const SGrid& grid, int face_index);
	static void info_gridcell(const SGrid& grid, int cell_index);
	static void info_gridedge(const Edge& edge);
	static void info_gridface(const Face& face);
	static void info_gridcell(const Cell& cell);
	
	template<class G>
	static double hash(const G& grid){
		double sum;
		int i=0;
		for (auto vert: grid.vvert){
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

	static void save_bnd_vtk(const SGrid& grid);
	static void save_grid_vtk(const GridData& grid);
	static void save_surf_vtk(const Surface& grid);
	static void save_cell_vtk(shared_ptr<Cell> c);
	static void save_cells_vtk(const CellData& grid);
	static void save_cells_faces_vtk(const CellData& grid);
	static void save_edges_vtk(const EdgeData& edges);
	static void save_faces_vtk(const FaceData& faces);

	//extracts primitives without ids change
	static VertexData all_vertices(const CellData&);
	static VertexData all_vertices(const FaceData&);
	static VertexData all_vertices(const EdgeData&);
	static EdgeData all_edges(const CellData&);
	static EdgeData all_edges(const FaceData&);
	static FaceData all_faces(const CellData&);
};

}

#endif
#endif

