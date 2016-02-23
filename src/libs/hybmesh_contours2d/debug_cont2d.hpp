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
	static void geogebra_ecollection(const ECollection& ecol);

	//contour in vtk format
	static void vtk_contours(const ShpVector<HMCont2D::Contour>& c, const char* fn="_debug.vtk");
	static void vtk_contour(const ECollection& col, std::map<Point*, double>& pdata,
			double defval=0, const char* fn="_debug.vtk");
	static void vtk_contour(const ECollection& col, std::map<Edge*, double>& edata,
			double defval=0, const char* fn="_debug.vtk");
	template<class C1, class C2>
	static void vtk_contour(const ECollection& col, std::map<C1, C2>& edata,
			double defval=0, const char* fn="_debug.vtk"){
		std::map<Edge*, double> edata1;
		for (auto e: edata) edata1[const_cast<Edge*>(e.first)] = (double)e.second;
		vtk_contour(col, edata1, defval, fn);
	}

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
