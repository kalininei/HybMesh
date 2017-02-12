#ifndef HYBMESH_CONTOURS2D_GPC_CORE_HPP
#define HYBMESH_CONTOURS2D_GPC_CORE_HPP
#include "primitives2d.hpp"
#include "contour_tree.hpp"

extern "C"{
#include "gpc.h"
}

namespace HM2D{ namespace Impl{

class GpcTree{
	gpc_polygon poly;
	GpcTree(): poly{0, 0, 0}{}
public:
	GpcTree(const EdgeData& inp);
	GpcTree(const Contour::Tree& inp);
	GpcTree(const GpcTree& other);
	GpcTree(GpcTree&& other) noexcept;
	GpcTree& operator=(const GpcTree& other);
	~GpcTree();

	Contour::Tree ToContourTree() const;
	
	static GpcTree Union(const GpcTree& c1, const GpcTree& c2);
	static GpcTree Intersect(const GpcTree& c1, const GpcTree& c2);
	static GpcTree Substract(const GpcTree& c1, const GpcTree& c2);
	static GpcTree Xor(const GpcTree& c1, const GpcTree& c2);
};

}}
#endif
