#ifndef HYBMESH_CONTOURS2D_GPC_CORE_HPP
#define HYBMESH_CONTOURS2D_GPC_CORE_HPP

#include "hybmesh_contours2d.hpp"

extern "C"{
#include "gpc.h"
}

namespace HMCont2D{ namespace Impl{

class GpcTree{
	gpc_polygon poly;
	GpcTree(): poly{0, 0, 0}{}
public:
	GpcTree(const HMCont2D::Contour& inp);
	GpcTree(const HMCont2D::ContourTree& inp);
	GpcTree(const GpcTree& other);
	GpcTree(GpcTree&& other) noexcept;
	GpcTree& operator=(const GpcTree& other);
	~GpcTree();

	HMCont2D::Container<HMCont2D::ContourTree> ToContourTree() const;
	
	static GpcTree Union(const GpcTree& c1, const GpcTree& c2);
	static GpcTree Intersect(const GpcTree& c1, const GpcTree& c2);
	static GpcTree Substract(const GpcTree& c1, const GpcTree& c2);
};

}}
#endif
