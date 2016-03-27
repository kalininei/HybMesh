#ifndef HYBMESH_CONTOURS2D_CLIPPER_CORE_HPP
#define HYBMESH_CONTOURS2D_CLIPPER_CORE_HPP

#include "hybmesh_contours2d.hpp"
#include "clipper.hpp"

namespace HMCont2D{ namespace Impl{
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
	ClipperPath(const Contour&);

	//scaling methods
	void ApplyBoundingBox(const BoundingBox& newbbox) override;

	// === Modify geometry
	void SetPoints(const ShpVector<Point>& p);
	void AddPointToEnd(const Point& p);

	// === Methods
	Container<ContourTree> Offset(double delta, HMCont2D::OffsetTp tp) const;
	//-1 - point is on polygon
	// 0 - point is outside polygon
	// 1 - point is in polygon
	int WhereIs(Point p) const;

	Container<Contour> ToHMContainer();
	static Container<Contour> HMContainer(ClipperLib::Path&, const BoundingBox& bbox, bool is_closed=true);

	//intersection of closed contours.
	//if embedded1/2 is true then contours in сorresponding vector is treated as a tree.
	static Container<ContourTree> Intersect(
			vector<ClipperPath>& pths1,
			vector<ClipperPath>& pths2,
			bool embedded1,
			bool embedded2);

	//union of closed contours. 
	//if embedded1/2 is true then contours in сorresponding vector are treated as a tree.
	static Container<ContourTree> Union(
			vector<ClipperPath>& pths1,
			vector<ClipperPath>& pths2,
			bool embedded1,
			bool embedded2);

	//substruction of closed contours. 
	//if embedded1/2 is true then contours in сorresponding vector are treated as a tree.
	static Container<ContourTree> Substruct(
			vector<ClipperPath>& pths1,
			vector<ClipperPath>& pths2,
			bool embedded1,
			bool embedded2);

	//cuts lines by area intersection
	static Container<ExtendedTree> CutLines(
			vector<ClipperPath>& area,
			vector<ClipperPath>& lines,
			bool embedded_area);

};


// ================================ ClipperTree =========================
struct ClipperTree: public ClipperObject{
	//overriden
	void ApplyBoundingBox(const BoundingBox& newbbox) override;
	//from tree
	static ClipperTree Build(const ContourTree& tree);
	//to container
	static Container<ContourTree> HMContainer(const ClipperLib::PolyTree&, const BoundingBox&);
	//sorting physical points: INSIDE/OUTSIDE/BOUND for each point
	vector<int> SortOutPoints(const vector<Point>& pts) const;

private:
	vector<ClipperPath> data;
	vector<bool> isopen;
};


}}



#endif

