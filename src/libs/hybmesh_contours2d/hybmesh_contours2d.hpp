#ifndef HYBMESH_CONTOURS2D_HPP
#define HYBMESH_CONTOURS2D_HPP
#include "hmproject.h"
#include "bgeom2d.h"

namespace HMCont2D{

//2D contour base class
class Contour{
	std::unique_ptr<int> __not_copiable_object;
protected:
	// =========== Data
	//path
	ShpVector<Point> pts;
	//boundary type of point
	std::map<const Point*, int> btype;

	//object which implements procedures from clipper lib
	//HMCont2DCore::Path
	mutable void* _core;

	// ============ Methods
	//get data from another contour.
	//object c is not longer valid after this
	void GetData(Contour& c);

	//make a deep copy of contour c data into this
	void CopyData(const Contour& c);
public:
	//=== Constructors
	Contour();
	Contour(Contour&& other) { GetData(other); }
	Contour(const Contour& other) { CopyData(other); }
	Contour& operator=(const Contour& other){ if (this != &other) CopyData(other); return *this; }
	virtual ~Contour();

	//=== Set Geometry
	void Clear();
	void AddPointToEnd(double x, double y, int b=0);
	void AddPointToEnd(Point xy, int b=0);

	//=== Modify Geometry
	virtual void Reverse();

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
//   path in counter clockwise direction
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


	//=== Modify Geometry
	void ForceDirection(bool is_inner);

	//=== GeometryInfo
	bool IsClosed() const override {return false;}
	bool IsValid() const override;
	int NumEdges() const override { return pts.size(); }
	const Point* Pnt(int i) const override { return  (i>=pts.size()) ? Pnt(i-NumPoints()) : pts[i].get(); }

	//does point lie strictly within contour not regarding to its orientation
	bool IsWithinGeom(const Point& p) const;

	//INSIDE or OUTSIDE
	int Direction() const;

	//SignedArea < 0 for Direction() == Outside. Area always > 0.
	double Area() const;
	double SignedArea() const;
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

	struct Entry{
		Entry(ClosedContour* c){
			self = c;
			parent = 0;
		}
		ClosedContour* self;
		Entry* parent;
		vector<Entry*> children;
	};
	void ForceMultiplicity(Entry* e, bool is_inner);
	Entry* FindEntry(ClosedContour*);
	ShpVector<Entry> entries;
	vector<Entry*> top_level;

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
	double Area() const;
};

void SaveVtk(const ContourTree& Tree, const char* fname);


}//HMCont2D

#endif
