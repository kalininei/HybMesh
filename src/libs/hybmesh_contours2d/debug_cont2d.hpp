#ifndef HYBMESH_CONTOURS2D_DEBUG_HPP
#define HYBMESH_CONTOURS2D_DEBUG_HPP

#ifndef NDEBUG

#include "collections.hpp"
#include "contour.hpp"
#include "tree.hpp"

namespace HMCont2D{

struct Debug{
	static void info_contour(const Contour& c);
	static void info_ecollection(const ECollection& c);
	static void info_tree(const ContourTree& c);
	static void info_tree_node(const ContourTree::TreeNode& c);

	//contour in geogebra script format
	static void geogebra_contour(const Contour& c);
	static void geogebra_tree(const ContourTree& c);
	static void geogebra_etree(const ExtendedTree& c);
	static void geogebra_box(const BoundingBox& c);

	//contour in vtk format
	static void vtk_contours(const ShpVector<HMCont2D::Contour>& c, const char* fn="_debug.vtk");

	//printing with tabulation
	static int tabs;
	static void Print(const char* fmt, ...);
	static std::ostream& Cout();

	//aux data
	Debug(): ptr(0), i(0) { pre(); }
	void* ptr;
	int i;

	void pre();
};
extern Debug _dbg;

}

#endif
#endif
