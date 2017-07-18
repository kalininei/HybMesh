#ifndef HYBMESH_SIZEFUN_HPP
#define HYBMESH_SIZEFUN_HPP
#include "contour_tree.hpp"
#include "femgrid43.hpp"
#include "piecewise.hpp"

namespace HM2D{ namespace Grid{

class SizeFun{
public:
	virtual double sz(const Point& p) const=0;
	virtual double sz_proj(const Point& p) const=0;
	virtual vector<double> sz(const vector<Point>& p) const;
	virtual HMMath::LinearPiecewise sz_boundary_segment(const HM2D::EdgeData& pos) const=0;
	virtual double minstep() const = 0;
	virtual double maxstep() const = 0;
};

//build size function using source edges with id=1 and given points as sources
shared_ptr<SizeFun> BuildSizeFunction(const Contour::Tree& source,
		const vector<std::pair<Point, double>>& psrc={});

// If source has any edges marked with id!=1 then a size function will be built
// and those edges will be resegmented.
//
// While this resegmentation vertices with id=1 will not be touched. All primitives
// which were not changed will not be reallocated.
//
// If force_sizefun=true then a size function will always be built and returned,
// otherwise it will be constructed only if resegmentation occured.
shared_ptr<SizeFun> ApplySizeFunction(Contour::Tree& source,
		const vector<std::pair<Point, double>>& psrc={},
		bool force_sizefun=false);



}}

#endif

