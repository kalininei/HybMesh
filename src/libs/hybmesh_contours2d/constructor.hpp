#ifndef HMCONT2D_CONSTRUCTOR_HPP
#define HMCONT2D_CONSTRUCTOR_HPP

#include "hybmesh_contours2d.hpp"

namespace HMCont2D{ namespace Constructor{

//Shaped contours
HMCont2D::Container<Contour> Circle(int N, double rad, Point cnt);

//Contour from points
HMCont2D::Contour ContourFromPoints(const vector<Point*>& pnt, bool force_closed=false);
HMCont2D::Contour ContourFromPoints(const HMCont2D::PCollection& dt, bool force_closed=false);
HMCont2D::Container<Contour> ContourFromPoints(vector<double> pnt, bool force_closed=false);
HMCont2D::Container<Contour> ContourFromPoints(vector<Point> pnt, bool force_closed=false);

//Contour from another contours with deep copy
HMCont2D::Container<Contour> CutContour(const HMCont2D::Contour& cont, const Point& pstart, int direction, double len);


}};

#endif
