#ifndef HMCONT2D_CONT_REPART_HPP
#define HMCONT2D_CONT_REPART_HPP

#include "contour.hpp"

namespace HMCont2D{namespace Algos{

//takes contours and makes partition by creating vertex points of so that:
//1) all sections were less than h
//2) all subsegments points lie in a diamond built
//   on first and last points with angle=cirt_angle at it
//3) mandatory points were among vertex points
//4) any segment were no more than 3 times longer than adjacent ones
//5) initial contour points were snapped to vertex points if they lie
//   closer than snap_dist * its subcontour length

//returns vector of vector were each first point is the vertex point
//and all others are contour points which lie within a segment
vector<vector<Point*>> Coarsening(const HMCont2D::Contour& cont, 
		const std::set<Point*>& mandatory,
		HMCont2D::PCollection& pcol,
		double h, double crit_angle, double snap_dist);

}}

#endif
