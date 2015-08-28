#ifndef HYBMESH_CONTOURS2D_HPP
#define HYBMESH_CONTOURS2D_HPP
#include "hmproject.h"
#include "bgeom2d.h"

namespace HMCont2D{

//2D contour base class
class Contour{
protected:
	//path
	ShpVector<Point> pts;

	//reallocates all pts
	//This is used for making contours with unique dataset
	void PtsReallocate();
public:
	//=== Constructors
	Contour(){};

	//=== Set Geometry
	void Clear();
	void SetPointList(const vector<double>& xy);
	void SetPointList(const vector<Point>& xy);
	void AddPointToEnd(double x, double y);
	void AddPointToEnd(Point xy);

	//=== Modify Geometry
	void Simplify();

	//=== GeometryInfo
	int NumPoints() const { return pts.size(); }
	virtual bool IsClosed() const = 0;
	//true if no self crosses
	virtual bool IsValid() const = 0;
	//bounding box
	BoundingBox BuildBoundingBox() const;

	//=== Operations
	//if two contours intersects each other
	static bool HasCrosses(const Contour& c1, const Contour& c2);

};

//Closed contour
//   path in anticlockwise direction
//   self crossing is not allowed
//   endpoints are not doubled
class ClosedContour: public Contour{
public:
	//=== Constructors
	ClosedContour(): Contour(){}
	ClosedContour DeepCopy() const;

	//=== GeometryInfo
	bool IsClosed() const override {return false;}
	bool IsValid() const override;

};

//Path as sequence of points
//   self crossing is not allowed
class OpenContour: public Contour{
public:
	//=== Constructors
	OpenContour(): Contour(){}
	OpenContour DeepCopy() const;

	//=== GeometryInfo
	bool IsClosed() const override {return true;}
	bool IsValid() const override;
};

//Contours structure
class ContourTree{
	ShpVector<ClosedContour> conts;
	void RebuildStructure();
public:
	//=== Constructors
	ContourTree(){}

	//=== AddContours
	//Contours are copied to ContourTree structure
	//throws std::runtime_error if something went wrong
	void AddContour(const ClosedContour& cont);

	//=== Geometry Info
	BoundingBox BuildBoundingBox() const;
};

}//HMCont2D

#endif
