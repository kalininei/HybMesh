#ifndef HYBMESH_LAPLACE_BEM2D_HPP
#define HYBMESH_LAPLACE_BEM2D_HPP

#include "contour_tree.hpp"
#include "nodes_compare.h"

namespace HMBem{

//2D solution using constant elements
class LaplaceCE2D{
	vector<double> dirvals, neuvals;
	vector<int> iprev;
	vector<bool> isdir;
	vector<double> len, len2, loglen, nx, ny, xs, ys, xm, ym;
	mutable CoordinateMap2D<int> vfinder;
	mutable bool no_bmap;
	int N;
	int where_is(double x, double y) const;
	//-> integral_Ck[F(x, y)], integral_Ck[dF/dn(x, y)]
	std::pair<double, double> Ffun(int k, double x, double y) const;
	double B(int k, double x, double y) const;
	double E(int k, double x, double y) const;
public:
	LaplaceCE2D(const HM2D::Contour::Tree& area);
	void Solve();
	//set boundary values
	//default is df/dn=0
	void dirichlet_value(int i, double value){ isdir[i]=true; dirvals[i]=value; }
	void neumann_value(int i, double value){ isdir[i]=false; neuvals[i]=value; }
	std::pair<double, double> boundary_data(int i) const { return std::make_pair(dirvals[i], neuvals[i]); }
	double value_at(double x, double y) const;
	double internal_value_at(double x, double y) const;
	double boundary_value_at(double x, double y) const;
};

}

#endif
