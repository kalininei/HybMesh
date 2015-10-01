#ifndef HYBMESH_HMFEM_FEMGRID43
#define HYBMESH_HMFEM_FEMGRID43
#include "grid.h"

namespace HMFem{ 

//grid shares points with another gridgeom object
//but has connectivity with only 3/4 nodes elements
class Grid43: public GridGeom{
public:
	//build grid by simplification of
	//arbitraty cell grid
	static shared_ptr<Grid43>
	Build(GridGeom* orig);

	//simple triangulatrion of singly-connected area
	static shared_ptr<Grid43>
	Build3(const vector<Point>& p, double h);

	//build a triangle grid which shares points whith
	//common grid
	static shared_ptr<Grid43>
	Build3(GridGeom* orig);

	class Approximator;
	shared_ptr<Approximator> GetApprox() const;
};

class Grid43::Approximator{
	const Grid43* grid;
	shared_ptr<GGeom::Info::CellFinder> cfinder;
public:
	Approximator(const Grid43* g, int nx=20, int ny=20):
		grid(g),
		cfinder(new GGeom::Info::CellFinder(g, nx, ny)){};

	//function value calculator.
	//If p is not within grid throws OutOfArea
	double Val(Point p, const vector<double>& fun) const;
	vector<double> Vals(Point, const vector<const vector<double>*>& funs) const;
};



}
#endif
