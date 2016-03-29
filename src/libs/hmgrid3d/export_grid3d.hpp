#ifndef EXPORT_GRID_3D_HPP
#define EXPORT_GRID_3D_HPP
#include "hmgrid3d.hpp"

namespace HMGrid3D{namespace Export{


std::tuple<
	std::vector<double>,   //points
	std::vector<int>       //
>
Serialize(const HMGrid3D::Grid& g);


void GridVTK(const HMGrid3D::Grid& g, const char* fn);


void BoundaryVTK(const HMGrid3D::Grid& g, const char* fn);



}}
#endif
