#ifndef HMMAPPING_CIRCRECT_HPP
#define HMMAPPING_CIRCRECT_HPP

#include "grid.h"

namespace HMMap{

//fills circle with rectangular cells: linear algorithm
//8 should be a divisor of n
//a - side of an internal square normalized by radius
//hcoef - radial step near the outer side of contour will be hcoef*2*pi*rad/n;
//hcoef = 1 gives square cells
//hcoef < 1 gives refinement towards outer boundary
//algos : "linear", "laplace", "orthogonal-circ", "orthogonal-rect"
GridGeom Circ4Prototype(Point center, double rad, int n, std::string algo,
	double a=1.0, double hcoef=1.0);
}

#endif
