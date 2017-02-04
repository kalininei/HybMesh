#ifndef HMCONT2D_ALGOS_HPP
#define HMCONT2D_ALGOS_HPP

#include "primitives2d.hpp"
#include "tree.hpp"

namespace HM2D{

//algorithms for consequitive contours and
//contour trees
namespace Contour{ namespace Algos{

//Offset contour
enum class OffsetTp{
	//rounding at corners
	RC_CLOSED_POLY,
	RC_OPEN_ROUND,
	RC_OPEN_BUTT,
	//square corners
	SC_CLOSED_POLY,
};

//takes into account direction of source and sign of delta:
//offsets to the right side if delta>0.
Tree Offset(const EdgeData& source, double delta, OffsetTp tp);
//forces singly connected output contour. tp = CLOSED_POLY or OPEN_ROUND
EdgeData Offset1(const EdgeData& source, double delta);

//Calculate points position with respect to contour.
//Contour direction is taken into account.
//Builds vector with INSIDE/OUTSIDE/BOUND for each point.
vector<int> SortOutPoints(const EdgeData& t1, const vector<Point>& pnt);
vector<int> SortOutPoints(const Tree& t1, const vector<Point>& pnt);

// ========= Simplifications
//remove points which lie on the same edge
//removes zero length edges
EdgeData Simplified(const EdgeData& cont);
Contour::Tree Simplified(const Tree& t1);

//Gives a vector defining a director of contour in given point p which should be amoung contour points.
//direction = 1 -> consider contour as it is
//direction = -1 -> revert contour before procedure
//len_forward, len_backward - length of smoothing.
//   *_forward and *_backward are steps along given contour direction
//   before possible reversion due to `direction=-1` option.
//result is a unit vector
Vect SmoothedDirection2(const EdgeData& c, const Point* p, int direction, double len_forwad, double len_backward);

//removes points by indicies
//if inner point is removed then a new edge will be created.
//It will keep all features of the edge previous to deleted point
void RemovePoints(EdgeData& data, vector<int> ipnt);

//get boundary types from contour coordinate
vector<int> BTypesFromWeights(const EdgeData&, const vector<double>&);
}}

// ============================ algorithms for shattered edges collections
namespace ECol{ namespace Algos{
//Does simplifications of simple contours assembled from ecol shattered edges.
//makes deep copies of input edges (not points).
//if angle<0 returns everything back making copies of edges but not points.
EdgeData Simplified(const EdgeData& ecol, double degree_angle, bool no_break_at_bt_change=false,
		const VertexData& keep={});

//deep copies everything
EdgeData NoCrosses(const EdgeData& ecol);

void MergePoints(EdgeData& ecol);

//set 'to' boundary types from closest 'from' edges.
void AssignBTypes(const EdgeData& from, EdgeData& to);
}}


}
#endif
