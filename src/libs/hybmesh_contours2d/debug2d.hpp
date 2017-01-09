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

	//contour in vtk format
	static void save_edges_vtk(const EdgeData& c);
	static void save_edges_vtk(const EdgeData& c, const std::map<Point*, double>& data);
	static void save_vertices_vtk(const VertexData& c);

	static double hash(const EdgeData& ecol);

	template<class G>
	static void outhash(const G& grid, std::string prefix=""){
		if (prefix.size() > 0) prefix = " ("+prefix+")";
		std::cout<<"Hash for contour 2d"<<prefix<<": "<<hash(grid)<<std::endl;
	}

	//extract primitives without id change
	static VertexData allvertices(const EdgeData& c);
};

}

#endif
#endif
