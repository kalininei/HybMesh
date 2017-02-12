#ifndef HMCONT2D_CONT_REPART_HPP
#define HMCONT2D_CONT_REPART_HPP

#include "contour.hpp"

namespace HM2D{namespace Contour{ namespace Algos{

//takes contours and makes partition by creating vertex points so that:
//1) all sections were less than h
//2) all subsegments points lie in a diamond built
//   on first and last points with angle=cirt_angle at it
//3) mandatory points were among vertex points
//4) any segment were no more than 3 times longer than adjacent ones
//5) initial contour points were snapped to vertex points if they lie
//   closer than snap_dist * its subcontour length

//returns vector of vector were each first point is the vertex point
//and all others are contour points which lie within a segment
vector<VertexData> Coarsening(const EdgeData& cont, 
		const std::set<shared_ptr<Vertex>>& mandatory,
		double h, double crit_angle, double snap_dist);


//Gives a vector defining a director of contour in given point p which should be amoung contour points.
//direction = 1 -> consider contour as it is
//direction = -1 -> revert contour before procedure
//len_forward, len_backward - length of smoothing.
//   *_forward and *_backward are steps along given contour direction
//   before possible reversion due to `direction=-1` option.
//result is a unit vector
Vect SmoothedDirection2(const EdgeData& c, const Point* p, int direction, double len_forwad, double len_backward);

}}}

#endif
