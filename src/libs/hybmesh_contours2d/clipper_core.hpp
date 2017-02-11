#ifndef HYBMESH_CONTOURS2D_CLIPPER_CORE_HPP
#define HYBMESH_CONTOURS2D_CLIPPER_CORE_HPP

#include "clipper.hpp"
#include "contour.hpp"
#include "tree.hpp"
#include "algos.hpp"

namespace HM2D{ namespace Impl{
class ClipperPath;
class ClipperTree;

//area of integer arithmetic computing is within [-Resolution, Resolution] square
//const ClipperLib::cInt CLIPPER_RESOLUTION = 1e8;
const ClipperLib::cInt CLIPPER_RESOLUTION = 100.0/geps;
//normalized to [0, 1] distance between true and approximated arcs
const double ClipperArcTolerance = 0.0001;

// ================================ ClipperPath
struct ClipperObject{
	//BoundingBox in real geometry
	BoundingBox bbox;
	//scaling:
	//   real_geometry = p0 + integer_geometry / factor
	//   integer_geometry = factor * (real_geomtery - p0)
	Point p0;
	long double factor;

	ClipperObject(): bbox(0.0, 0.0, 0.0, 0.0), p0(0,0), factor(1.0){}
	ClipperObject(const BoundingBox& bb) { ApplyBoundingBox(bb);}

	//scaling methods
	ClipperLib::IntPoint ToIntGeom(const Point& p) const;
	Point ToRealGeom(const ClipperLib::IntPoint& p) const;
	virtual void ApplyBoundingBox(const BoundingBox& newbbox);
};


struct ClipperPath: public ClipperObject{
	//Clipper::Path object
	ClipperLib::Path data;

	//constructor
	ClipperPath(): ClipperObject(), data(){};
	ClipperPath(const EdgeData&);

	//scaling methods
	void ApplyBoundingBox(const BoundingBox& newbbox) override;

	// === Modify geometry
	void SetPoints(const ShpVector<Point>& p);
	void AddPointToEnd(const Point& p);

	// === Methods
	Contour::Tree Offset(double delta, HM2D::Contour::Algos::OffsetTp tp) const;
	//-1 - point is on polygon
	// 0 - point is outside polygon
	// 1 - point is in polygon
	int WhereIs(Point p) const;

	EdgeData ToHMContainer();
	static EdgeData HMContainer(ClipperLib::Path&, const BoundingBox& bbox, bool is_closed=true);

	//intersection of closed contours.
	//if embedded1/2 is true then contours in сorresponding vector is treated as a tree.
	static Contour::Tree Intersect(
			vector<ClipperPath>& pths1,
			vector<ClipperPath>& pths2,
			bool embedded1,
			bool embedded2);

	//union of closed contours. 
	//if embedded1/2 is true then contours in сorresponding vector are treated as a tree.
	static Contour::Tree Union(
			vector<ClipperPath>& pths1,
			vector<ClipperPath>& pths2,
			bool embedded1,
			bool embedded2);

	//substruction of closed contours. 
	//if embedded1/2 is true then contours in сorresponding vector are treated as a tree.
	static Contour::Tree Substruct(
			vector<ClipperPath>& pths1,
			vector<ClipperPath>& pths2,
			bool embedded1,
			bool embedded2);

	//cuts lines by area intersection
	static Contour::Tree CutLines(
			vector<ClipperPath>& area,
			vector<ClipperPath>& lines,
			bool embedded_area);

};


// ================================ ClipperTree =========================
struct ClipperTree: public ClipperObject{
	//overriden
	void ApplyBoundingBox(const BoundingBox& newbbox) override;
	//from tree
	static ClipperTree Build(const Contour::Tree& tree);
	//to container
	static Contour::Tree HMContainer(const ClipperLib::PolyTree&, const BoundingBox&);

	//sorting physical points: INSIDE/OUTSIDE/BOUND for each point
	//BOUND result is not reliable.
	//!!! Use this only if it is known that pts do not lie on boundaries.
	vector<int> SortOutPoints(const vector<Point>& pts) const;

private:
	vector<ClipperPath> data;
	vector<bool> isopen;
};


}}


#endif

