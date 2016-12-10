#ifndef HMCONT2D_ALGOS_HPP
#define HMCONT2D_ALGOS_HPP
#include "contour.hpp"
#include "tree.hpp"
#include "cont_partition.hpp"
#include "cont_repart.hpp"

namespace HMCont2D{

//Offset contour
enum class OffsetTp{
	//rounding at corners
	RC_CLOSED_POLY,
	RC_OPEN_ROUND,
	RC_OPEN_BUTT,
	//square corners
	SC_CLOSED_POLY,
};

namespace Algos{
// ================================== Offset
//takes into account direction of source and sign of delta:
//all positives -> offsets to the left etc.
Container<ContourTree> Offset(const Contour& source, double delta, OffsetTp tp);
//forces singly connected output contour. tp = CLOSED_POLY or OPEN_ROUND
Container<Contour> Offset1(const Contour& source, double delta);


// ============================ Crosses and intersections
//finds first cross (with respect to length of c1) of contours c1, c2.
//returns <0>: if cross was found
//        <1>: cross point
//        <2,3>: normalized length coordinate of intersection
std::tuple<bool, Point, double, double> 
Cross(const Contour& c1, const Contour& c2);

vector<std::tuple<bool, Point, double, double>>
CrossAll(const Contour& c1, const Contour& c2);

//returns true if c1 and c2 have common area (not a point, but may be an edge)
bool DoIntersect(const Contour& c1, const Contour& c2);
bool DoIntersect(const ContourTree& t1, const Contour& c2);

//returns true if c1 and c2 have common area (not a point, not an edge)
//is not realiable if  Area(c1) >> Area(c2) 
bool DoReallyIntersect(const Contour& c1, const Contour& c2);

//Calculate points position with respect to contour.
//Contour direction is taken into account.
//Builds vector with INSIDE/OUTSIDE/BOUND for each point.
vector<int> SortOutPoints(const Contour& t1, const vector<Point>& pnt);
vector<int> SortOutPoints(const ContourTree& t1, const vector<Point>& pnt);

// =========================== Simplifications
//remove points which lie on the same edge
//removes zero length edges
Contour Simplified(const Contour& cont);
ContourTree Simplified(const ContourTree& t1);
ExtendedTree Simplified(const ExtendedTree& t1);
ECollection Simplified(const ECollection& ecol, double degree_angle, bool no_break_at_id_change=false);

Container<ECollection> NoCrosses(const ECollection& ecol);
void MergePoints(ECollection& ecol);
void DeleteUnusedPoints(Container<ECollection>& econt);

// =========================== Smoothing
//calculates vector representing direction of contour at point p smoother by lengh len
//if p doesn't lie on c -> project it to c and calculate
Vect SmoothedDirection(const Contour& c, Point* p, int direction, double len);

//Gives a vector defining a director of contour in given point p which should be amoung contour points.
//direction = 1 -> consider contour as it is
//direction = -1 -> revert contour before procedure
//len_forward, len_backward - length of smoothing.
//   *_forward and *_backward are steps along given contour direction
//   before possible reversion due to `direction=-1` option.
//result is a unit vector
Vect SmoothedDirection2(const Contour& c, const Point* p, int direction, double len_forwad, double len_backward);


}
};


#endif
