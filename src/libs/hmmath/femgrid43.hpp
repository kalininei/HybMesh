#ifndef HYBMESH_HMFEM_FEMGRID43
#define HYBMESH_HMFEM_FEMGRID43
#include "grid.h"
#include "procgrid.h"

namespace HMFem{ 

//grid shares points with another gridgeom object
//but has connectivity with only 3/4 nodes elements
class Grid43: public GridGeom{
public:
	//build grid by simplification of
	//arbitraty cell grid
	static shared_ptr<Grid43>
	Build(GridGeom* orig);

	//triangulatrion of singly-connected area
	static shared_ptr<Grid43>
	Build3(const vector<Point>& p, double h);

	//triangulation of multiply-connected area
	static shared_ptr<Grid43>
	Build3(const vector<vector<Point>>& cnts, double h);

	//triangulation of multiply-connected area
	//with constraint lines
	static shared_ptr<Grid43>
	Build3(const HMCont2D::ContourTree& conts,
		const ShpVector<HMCont2D::Contour>& constraints,
		double h);

	//build a triangle grid which shares points whith
	//common grid
	static shared_ptr<Grid43>
	Build3(GridGeom* orig);

	class Approximator;
	shared_ptr<Approximator> GetApprox() const;
	class IsolineBuilder;

	//first and last points  in pts should be equal to existing grid points
	//connected by boundary edge. All middle points will be added between them,
	//cell will be triangulated.
	void AddSegments(const vector<vector<Point>>& pts);
};

//build auxiliary triangle grid for elliptic problems fem solution
//all points of tree will present in resulting grid
//nrec, nmax - recommended and maximum allowed number of resulting grid vertices
struct TAuxGrid3: public HMCallback::ExecutorBase{
	HMCB_SET_PROCNAME("Auxiliary Triangulation");
	HMCB_SET_DEFAULT_DURATION(100);

	//minumum allowed length of section normalized by recommended length
	static constexpr double ZEROANGLE = M_PI/4;
	static constexpr double CORRECTION_FACTOR = 1.13;

	Grid43 _run(const HMCont2D::ContourTree& _tree,
			const vector<HMCont2D::Contour>& _constraints,
			int nrec, int nmax);
	Grid43 _run(const HMCont2D::ContourTree& tree, int nrec, int nmax);
	Grid43 _run(const HMCont2D::Contour& cont, int nrec, int nmax);
private:
	HMCont2D::ContourTree tree;
	ShpVector<HMCont2D::Contour> constraints;
	HMCont2D::PCollection pcol;

	double step_estimate(const HMCont2D::ContourTree& tree, int nrec);
	void adopt_boundary(HMCont2D::ContourTree& tree, double h,
			vector<vector<Point>>& lost);
	void adopt_contour(HMCont2D::Contour& cont, double h,
			vector<vector<Point>>& lost);
	void input(const HMCont2D::ContourTree& _tree, const vector<HMCont2D::Contour>& _constraints);
	vector<Point*> gather_section(const vector<Point*>& ordered, int start, double h);

	std::unordered_set<Point*> mandatory_points;
	void mandatory_corner_points(const HMCont2D::Contour& cont);
	void mandatory_intersections(HMCont2D::Contour& c1, HMCont2D::Contour& c2);
	void clear();
};
extern HMCallback::FunctionWithCallback<TAuxGrid3> AuxGrid3;

class Grid43::Approximator{
	const Grid43* grid;
	shared_ptr<GGeom::Info::CellFinder> cfinder;

	static Point LocalCoordinates(const Cell* c, Point p);
	static Point LocalCoordinates3(const Cell* c, Point p);
	static Point LocalCoordinates4(const Cell* c, Point p);
	static double Interpolate(const Cell* c, Point ksieta, const vector<double>& fun);
	static double Interpolate3(const Cell* c, Point ksieta, const vector<double>& fun);
	static double Interpolate4(const Cell* c, Point ksieta, const vector<double>& fun);
public:
	Approximator(const Grid43* g, int nx=20, int ny=20):
		grid(g),
		cfinder(new GGeom::Info::CellFinder(g, nx, ny)){};

	//function value calculator.
	//If p is not within grid throws OutOfArea
	double Val(Point p, const vector<double>& fun) const;
	vector<double> Vals(Point, const vector<const vector<double>*>& funs) const;
};

class Grid43::IsolineBuilder{
	shared_ptr<Grid43> grid;
	shared_ptr<Grid43::Approximator> approx;
	
	//if isoline fun=value passes cell c then 
	//adds isoline section to Ecollection (returns true)
	//else do nothing (returns false)
	bool AddLine(const Cell* c, HMCont2D::ECollection&, HMCont2D::PCollection&,
			double value, const vector<double>& fun) const;
	//builds a segment using weighted grid points.
	//adds it to ecol.
	void BuildLine(int i0, int i1, double wi,
			int j0, int j1, double wj,
			HMCont2D::ECollection& ecol,
			HMCont2D::PCollection& pcol) const;
public:
	IsolineBuilder(shared_ptr<Grid43>& grid43, shared_ptr<Grid43::Approximator> approx43=0);

	//draw isoline from point in both directions till
	//closing or boundary.
	HMCont2D::Container<HMCont2D::Contour>
	FromPoint(Point p, const vector<double>& fun) const;
};



}
#endif
