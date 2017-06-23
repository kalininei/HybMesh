#ifndef HMCONT2D_CONSTRUCTOR_HPP
#define HMCONT2D_CONSTRUCTOR_HPP

#include "bgeom2d.h"
#include "primitives2d.hpp"

//Construction routines (unlike assembler routines) always return deep copies of input data
namespace HM2D{ namespace Contour{ namespace Constructor{

EdgeData Circle(int N, double rad, Point cnt);
EdgeData Circle(int N, Point cnt, Point point_on_curve);

EdgeData FromPoints(const vector<double>& pnt, bool force_closed=false);
EdgeData FromPoints(const vector<Point>& pnt, bool force_closed=false);
EdgeData Rectangle(Point botleft, Point topright);

//returns deep copies
EdgeData FromContours(const vector<EdgeData>& conts,
	bool last_close=true, bool fullshift=true, std::set<int> fixed={});

//Contour cut from another contours with deep copy
EdgeData CutContour(const EdgeData& cont,
		const Point& pstart, const Point& pend);
EdgeData CutContour(const EdgeData& cont,
		const Point& pstart, int direction, double len);

//builds parametric spline which passes given points
EdgeData Spline(const vector<Point>& pnt, int nedges, bool force_closed=false);

//resulting contours have no crosses and complicated connections.
//equal resulting contours vertices are shared
vector<EdgeData> ExtendedSeparate(const EdgeData& ecol);

}}

namespace ECol{ namespace Constructor{

EdgeData FromRaw(int npnt, int neds, double* pnt, int* eds);


}}

}

#endif
