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
	//boundary type of point
	std::map<const Point*, int> btype;

	//reallocates all pts
	//This is used for making contours with unique dataset
	void PtsReallocate();
public:
	//=== Constructors
	Contour(){};

	//=== Set Geometry
	void Clear();
	void AddPointToEnd(double x, double y, int b=0);
	void AddPointToEnd(Point xy, int b=0);

	//=== GeometryInfo
	int NumPoints() const { return pts.size(); }
	virtual int NumEdges() const { return pts.size() - 1; }
	virtual const Point* Pnt(int i) const { return pts[i].get(); }
	std::tuple<const Point*, const Point*, int> Edge(int i) const;

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

	//build diffent contours from connected points
	//reads only closed contours
	static vector<ClosedContour> FromConnectedPoints(const vector<Point>& xy,
		const vector<int>& i0, const vector<int>& i1,
		const vector<int>& bnd = vector<int>());


	//=== GeometryInfo
	bool IsClosed() const override {return false;}
	bool IsValid() const override;
	int NumEdges() const override { return pts.size(); }
	const Point* Pnt(int i) const override { return  (i>=pts.size()) ? Pnt(i-NumPoints()) : pts[i].get(); }

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
	int NumPoints() const;
	int NumEdges() const;
	int NumContours() const { return conts.size(); }
	const ClosedContour* Cont(int i) const { return conts[i].get(); }
	BoundingBox BuildBoundingBox() const;

};

void SaveVtk(const ContourTree& Tree, const char* fname);


}//HMCont2D

#endif
