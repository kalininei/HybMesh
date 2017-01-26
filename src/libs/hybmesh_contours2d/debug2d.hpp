#ifndef HYBMESH_CONTOURS2D_DEBUG_HPP
#define HYBMESH_CONTOURS2D_DEBUG_HPP

#ifndef NDEBUG

#include "primitives2d.hpp"
#include "contour.hpp"
#include "tree.hpp"
#include "hmdebug.hpp"

namespace HM2D{

struct Debug: public HMDebug{
	static void info_contour(const EdgeData& c);
	static void info_tree(const Contour::Tree& c);
	static void info_tree_node(const Contour::Tree::TNode& c);

	//contour in geogebra script format
	static void geogebra_contour(const EdgeData& c);
	static void geogebra_ecollection(const EdgeData& ecol);
	static void geogebra_tree(const Contour::Tree& c);
	static void geogebra_box(const BoundingBox& c);

	//to vtk format
	static void save_grid_vtk(const GridData& c);
	static void save_cells_vtk(const CellData& c);
	static void save_cells_vtk_id(const CellData& c);

	static void save_edges_vtk(const EdgeData& c);
	static void save_edges_vtk(const EdgeData& c, const std::map<Point*, double>& data);
	static void save_edges_vtk_id(const EdgeData& c);

	static void save_vertices_vtk(const VertexData& c);
	static void save_vertices_vtk_id(const VertexData& c);
	static void save_vertices_vtk(const EdgeData& c);
	static void save_vertices_vtk_id(const EdgeData& c);

	//extract primitives without id change
	static VertexData allvertices(const EdgeData& c);
	static EdgeData alledges(const CellData& c);
};

}

#endif
#endif
