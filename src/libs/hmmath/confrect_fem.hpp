#ifndef HYBMESH_HMMATH_CONFRECT_FEM_HPP
#define HYBMESH_HMMATH_CONFRECT_FEM_HPP
#include "hmconformal.hpp"
#include "femgrid43.hpp"

namespace HMMath{ namespace Conformal{ namespace Impl{
namespace ConfFem{

class ToRect: public HMMath::Conformal::Rect{
	// ==== data filled by constructor
	//fem grid in path area
	shared_ptr<HMFem::Grid43> grid;
	shared_ptr<HMFem::Grid43::Approximator> approx;
	//indexes of original points in grid
	//order is the same with given path
	vector<int> origs;
	//Contours made of original points. shares points with grid.
	//direction coinsides with canonic axes: u[0->1], v[0->1]
	HMCont2D::Contour left, bottom, right, top;
	//discrete mapping function to canonic domain
	vector<double> u, v;
	//conform module or -1 if construction has failed
	double _module;
	
	//fem grid in canonic area
	shared_ptr<HMFem::Grid43> inv_grid;
	shared_ptr<HMFem::Grid43::Approximator> inv_approx;
	//discrete mapping function to physical domain
	vector<double> inv_u, inv_v;

	//h - linear size of fem grid element.
	//If fails -> _module = -1;
	ToRect(const vector<Point>& path, int i1, int i2, int i3, double h);

	//==== Constructor subroutines
	//build fem grid: fills grid, approx, origs, left/bottom/right/top
	void BuildGrid(const vector<Point>& path, int i1, int i2, int i3, double h);
	//main costruction procedure: fills u, v, module
	void DoMapping();
	//fills inv_* data
	void BuildInverse();
	//get boundary points set as:
	//[0]->left, [1]->bottom, [2]->right, [3]->top
	//all points resulting points are sorted by contour weights
	std::array<vector<const GridPoint*>, 4> BndPoints() const;


	//estimates linear size of triangle grid
	//Nmax - maximum number of nodes
	//hrec - recommended linear size
	static double HEstimate(const vector<Point>& path, int segn, int nmax);
public:
	static shared_ptr<ToRect>
	Build(const vector<Point>& path, int i1, int i2, int i3, const Options& opt=Options());

	//conformal module
	double module() const override { return _module; }
	//Points from rectangle to polygon
	vector<Point> MapToPolygon(const vector<Point>& input) const override;
	//Points from polygon to rectangle
	vector<Point> MapToRectangle(const vector<Point>& input) const override;
	//Mapped Rectangle
	vector<Point> RectPoints() const override;
};


}}}}
#endif

