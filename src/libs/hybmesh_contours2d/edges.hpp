#ifndef HYBMESH_CONTOUR2D_EDGES_H
#define HYBMESH_CONTOUR2D_EDGES_H

#include "hmproject.h"
#include "bgeom2d.h"

namespace HMCont2D{

//Basic Edge Type
struct Edge{
	Point* pstart;
	Point* pend;

	Edge(Point* p1=0, Point* p2=0): pstart(p1), pend(p2){}

	bool contains(const Point* p) const { return (pstart == p || pend == p); }
	double meas() const { return Point::meas(*pstart, *pend); }
	double length() const { return sqrt(meas()); }

	Point* sibling(const Point* p) const {
		if (p == pstart) return pend;
		else if (p == pend) return pstart;
		else throw std::runtime_error("Invalid edge point");
	}

	static bool AreConnected(const Edge& e1, const Edge& e2){
		return (e1.pend == e2.pstart ||
			e1.pend == e2.pend   ||
			e1.pstart == e2.pstart ||
			e1.pstart == e2.pend);
	}

	static std::array<Point*, 3> PointOrder(const Edge& e1, const Edge& e2){
		if (e1.pend == e2.pstart) return { e1.pstart, e1.pend, e2.pend };
		else if (e1.pend == e2.pend) return { e1.pstart, e1.pend, e2.pstart };
		else if (e1.pstart == e2.pstart) return { e1.pend, e1.pstart, e2.pend };
		else if (e1.pstart == e2.pend) return { e1.pend, e1.pstart, e2.pstart };
		else throw std::runtime_error("Failed to build edges connection");
	}
};

struct BEdge: public Edge{
	int tp;
};


}


#endif
