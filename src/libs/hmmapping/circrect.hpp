#ifndef HMMAPPING_CIRCRECT_HPP
#define HMMAPPING_CIRCRECT_HPP

#include "grid.h"

namespace HMGMap{

//fills circle with rectangular cells: linear algorithm
//8 should be a divisor of n
//a - side of an internal square normalized by radius
//hcoef - radial step near the outer side of contour will be hcoef*2*pi*rad/n;
//hcoef = 1 gives square cells
//hcoef < 1 gives refinement towards outer boundary
GridGeom Circ4Prototype(Point center, double rad, int n, double a=1.0, double hcoef=1.0);

GridGeom Circ4Prototype2(Point center, double rad, int n, double a=1.0, double hcoef=1.0);

}

#endif
