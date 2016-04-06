#ifndef VTK_EXPORT_GRID2D_HPP
#define VTK_EXPORT_GRID2D_HPP
#include "grid.h"

namespace GGeom{ namespace Export{

void GridVTK(const GridGeom& g, std::string fn);

void BoundaryVTK(const GridGeom& g, std::string fn, const vector<int>& bcond={});

}}

#endif

