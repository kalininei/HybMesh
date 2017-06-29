#ifndef HYBMESH_HMFEM_FEMGRID43
#define HYBMESH_HMFEM_FEMGRID43
#include "primitives2d.hpp"
#include "contour_tree.hpp"
#include "hmcallback.hpp"

namespace HMFem{ 

namespace Grid43{

struct Approximator{
	const HM2D::GridData* grid;
	shared_ptr<BoundingBoxFinder> cfinder;
	vector<std::array<HM2D::Vertex*, 4>> cellvert;
	vector<std::array<int, 4>> icellvert;
	vector<bool> is3; //cell array
	vector<bool> isbnd; //vertex array
	//try to find point amoung cells with positive ordering.
	//if fails->searches amoung others
	//If given point is out of area but close to the boundary,
	//returns nearest cell with correct ksieta (out of [0, 1] triangle)
	//if fails->throws EOutOfArea
	int FindPositive(const Point& p, Point& ksieta) const;
	void FillJ3(std::array<double, 5>& J, int ic) const;
	void FillJ4(std::array<double, 5>& J, const Point& p, int ic) const;

	Point LocalCoordinates(int ic, Point p) const;
	Point LocalCoordinates3(int ic, Point p) const;
	Point LocalCoordinates4(int ic, Point p) const;
	double Interpolate(int ic, Point ksieta, const vector<double>& fun) const;
	double Interpolate3(int ic, Point ksieta, const vector<double>& fun) const;
	double Interpolate4(int ic, Point ksieta, const vector<double>& fun) const;

	struct NodePos{
		static const int Out = 0;
		static const int BndVertex = 1;
		static const int BndEdge = 2;
		static const int InternalVertex = 3;
		static const int InternalEdge = 4;
		static const int Internal = 5;

		int pos;
		int nvert, nve1, nve2, ncell;
		Point ksieta;
	};
	NodePos FindNodePos(Point pos) const;
	
	// -> node1, node2, ksi
	std::tuple<int, int, double> BndCoordinates(Point p) const;
//public:
	struct EOutOfArea: public std::runtime_error{
		EOutOfArea(): std::runtime_error("out of area"){}
	};
	Approximator(const HM2D::GridData* g, int n=40);

	//function value calculator.
	//If p is not within grid throws OutOfArea
	double Val(Point p, const vector<double>& fun) const;
	vector<double> Vals(Point, const vector<const vector<double>*>& funs) const;

	double BndVal(const Point& p, const vector<double>& fun) const;
	vector<double> BndVals(const Point& p, const vector<const vector<double>*>& funs) const;

	friend class IsolineBuilder;
};

class IsolineBuilder{
	shared_ptr<HM2D::GridData> grid;
	shared_ptr<Grid43::Approximator> approx;
	
	//if isoline fun=value passes cell c then 
	//adds isoline section to Ecollection (returns true)
	//else do nothing (returns false)
	bool AddLine(int c, HM2D::EdgeData&,
	             double value, const vector<double>& fun) const;
	//builds a segment using weighted grid points.
	//adds it to ecol.
	void BuildLine(int i0, int i1, double wi,
	               int j0, int j1, double wj,
	               HM2D::EdgeData& ecol) const;
public:
	IsolineBuilder(shared_ptr<HM2D::GridData>& grid43, shared_ptr<Grid43::Approximator> approx43=0);

	//draw isoline from point in both directions till
	//closing or boundary.
	HM2D::EdgeData FromPoint(Point p, const vector<double>& fun) const;
};



//first and last points  in pts should be equal to existing grid points
//connected by boundary edge. All middle points will be added between them,
//cell will be triangulated.
void AddSegments(HM2D::GridData& grid, const vector<vector<Point>>& pts);

//places a pts to grid if pts is within grid area.
//adds a new vertex with respective cell split if pts!=existing grid vertex
//additional pos output is a position of pts before insertion.
//return null if pts is outside grid otherwise pointer to existing (old or new) grid vertex.
shared_ptr<HM2D::Vertex> GuaranteePoint(HM2D::GridData& grid, Point pts,
		HMFem::Grid43::Approximator& approx);
shared_ptr<HM2D::Vertex> GuaranteePoint(HM2D::GridData& grid, Point pts,
		HMFem::Grid43::Approximator& approx,
		HMFem::Grid43::Approximator::NodePos& pos);
};

//build auxiliary triangle grid for elliptic problems fem solution
//all points of tree will present in resulting grid
//detached contours of input tree will be ignored. Use constraints argument.
//nrec - recommended number of resulting grid vertices
//hrec - recommended step size. Resulting step is a maximum of both condition.
//nrec or hrec could be -1. If so the corresponding condition is ignored.
//if resulting N>nmax exception will be raised
struct TAuxGrid3: public HMCallback::ExecutorBase{
	HMCB_SET_PROCNAME("Auxiliary triangulation");
	HMCB_SET_DEFAULT_DURATION(100);

	//minumum allowed length of section normalized by recommended length
	static constexpr double ZEROANGLE = M_PI/4;
	static constexpr double CORRECTION_FACTOR = 1.13;

	HM2D::GridData _run(const HM2D::Contour::Tree& tree,
			const vector<HM2D::EdgeData>& constraints,
			int nrec, int nmax, double hrec=-1);
	HM2D::GridData _run(const HM2D::Contour::Tree& tree, int nrec, int nmax, double hrec=-1);
	HM2D::GridData _run(const HM2D::EdgeData& cont, int nrec, int nmax, double hrec=-1);
private:
	HM2D::Contour::Tree tree;
	ShpVector<HM2D::EdgeData> constraints;
	HM2D::Contour::Tree ttree;

	double step_estimate(const HM2D::Contour::Tree& tree, int nrec, double hrec);
	void adopt_boundary(HM2D::Contour::Tree& tree, double h,
			vector<vector<Point>>& lost);
	void adopt_contour(HM2D::EdgeData& cont, double h,
			vector<vector<Point>>& lost);
	void adopt_complicated_connections(HM2D::Contour::Tree& cont,
			vector<vector<Point>>& lost);
	void input(const HM2D::Contour::Tree& _tree, const vector<HM2D::EdgeData>& _constraints);
	vector<Point*> gather_section(const vector<Point*>& ordered, int start, double h);

	std::set<shared_ptr<HM2D::Vertex>> mandatory_points;
	void mandatory_intersections(HM2D::EdgeData& c1, HM2D::EdgeData& c2);
	static bool angle_check(const vector<Point*>& line);
	void clear();
};
extern HMCallback::FunctionWithCallback<TAuxGrid3> AuxGrid3;


}
#endif
