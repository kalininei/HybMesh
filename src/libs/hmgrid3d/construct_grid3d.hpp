#ifndef  CONSTRUCT_GRID3D_HPP
#define  CONSTRUCT_GRID3D_HPP

#include "hmgrid3d.hpp"
#include "grid.h"

namespace HMGrid3D{namespace Constructor{

//build a cuboid in [0, 0, 0]x[lx, ly, lz] and translate it to leftp.
//boundary types of resulting grid are (for unit cube):
//    x = 0 -> bt = 1
//    x = 1 -> bt = 2
//    y = 0 -> bt = 3
//    y = 1 -> bt = 4
//    z = 0 -> bt = 5
//    z = 1 -> bt = 6
HMGrid3D::Grid Cuboid(HMGrid3D::Vertex leftp, double lx, double ly, double lz, int nx, int ny, int nz);

HMGrid3D::Grid SweepGrid2D(const GridGeom& g2d, const vector<double>& zcoords);


}}

#endif
