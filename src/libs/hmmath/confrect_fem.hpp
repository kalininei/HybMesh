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
	//HMCont2D::Contour left, bottom, right, top;
	vector<int> ileft, ibottom, iright, itop;
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

class ToAnnulus: public HMMath::Conformal::Annulus{
	//data filled by constructor
	double _module;
	shared_ptr<HMFem::Grid43> grid;
	shared_ptr<HMFem::Grid43::Approximator> approx;
	vector<double> u, v;   //radius, phi in canonical domain
	Point pzero;  //inner point in (x,y) with phi = 0;

	//indicies of grid points of original points
	vector<int> outer_origs, inner_origs;
	vector<int> outer, inner;
	vector<int> raz_io, raz_oi;

	//fem grid in canonic area
	shared_ptr<HMFem::Grid43> inv_grid;
	shared_ptr<HMFem::Grid43::Approximator> inv_approx;
	//discrete mapping function to physical domain
	vector<double> inv_u, inv_v;

	//constructor
	ToAnnulus(const vector<Point>& outer_path,
			const vector<Point>& inner_path, double h);

	//Constructor subroutines
	//builds doubly connected grid. Temporary fills grid, approx, outer, inner
	void BuildGrid1(const vector<Point>& outer_path,
			const vector<Point>& inner_path, double h);

	//fills u (temporary), _module
	void DoMappingU();

	//returns curve of steepest descent of u from inner to outer contours
	HMCont2D::Container<HMCont2D::Contour> SteepestDescentUCurve();

	//fills vgrid, vapprox, v_outer/inner_origs, v_raz_io/oi
	void BuildGrid2(const vector<Point>& outer_path,
			const vector<Point>& inner_path,
			double h,
			const HMCont2D::Contour& razor);

	//fills u, v
	void DoMapping();

	//fills inv_* data
	void BuildInverse();

	//estimates linear size of triangle grid
	//Nmax - maximum number of nodes
	//hrec - recommended linear size
	static double HEstimate(const vector<Point>& outer_path,
			const vector<Point>& inner_path, int segn, int nmax);

	// ===================== mapping subroutins and data
	mutable shared_ptr<HMCont2D::Contour> _inv_cont;
	const HMCont2D::Contour* InvGridContour() const;

public:
	static shared_ptr<ToAnnulus>
	Build(const vector<Point>& outer_path, const vector<Point>& inner_path,
			const Options& opt=Options());

	//===== overriden
	//= RadInner < 1.0
	double module() const override { return _module; }
	//get original points mapped to inner circle
	vector<Point> InnerCirclePoints() const override;
	//get original points mapped to outer circle
	vector<Point> OuterCirclePoints() const override;
	//mapping functions
	vector<Point> MapToOriginal(const vector<Point>& input) const override;
	vector<Point> MapToAnnulus(const vector<Point>& input) const override;

	//calculate angle at which i-th point is mapped to annulus
	double PhiInner(int i) const override;
	double PhiOuter(int i) const override;
};

}}}}
#endif

