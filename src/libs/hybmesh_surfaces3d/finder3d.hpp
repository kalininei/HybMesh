#ifndef HYBMESH_FINDER2D_HPP
#define HYBMESH_FINDER2D_HPP
#include "primitives3d.hpp"

namespace HM3D{ namespace Finder{

std::tuple<
	int,         //index of closest vertex within vec
	double       //measure to closest vertex
> ClosestPoint(const VertexData& vec, const Point3& v);


}}
#endif
