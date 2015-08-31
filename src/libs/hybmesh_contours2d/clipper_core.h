#ifndef HYBMESH_CONTOURS2D_CLIPPER_CORE_H
#define HYBMESH_CONTOURS2D_CLIPPER_CORE_H
#include "bgeom2d.h"
#include "clipper.hpp"


namespace HMCont2DCore{

//10^16
//area of computing is within [-Resolution, Resolution] square
const ClipperLib::cInt CLIPPER_RESOLUTION = 10000000000000000;

class Path{
	//scaling:
	//   real_geometry = p0 + integer_geometry / factor
	//   integer_geometry = factor * (real_geomtery - p0)
	Point p0;
	long double factor;
	ClipperLib::IntPoint ToIntGeom(const Point& p);
	Point ToRealGeom(const ClipperLib::IntPoint& p);

	//Clipper::Path object
	mutable ClipperLib::Path data;

	//BoundingBox in real geometry
	BoundingBox bbox;
	void ApplyBoundingBox(const BoundingBox& newbbox);
public:
	Path();
	
	// === Modify geometry
	void SetPoints(const ShpVector<Point>& p);
	void AddPointToEnd(const Point& p);
	void Reverse();

	// === Geometry Info
	int ClosedDirection() const; //INSIDE or OUTSIDE
	//returns true if point is in polygon
	//Does orientation matter ?????
	bool IsWithin(const Point& p); 
	//Orientation independent area
	double Area() const;
	//Orientation dependent area
	double SignedArea() const;
};



}


#endif
