#ifndef HMCONT2D_CONSTRUCTOR_HPP
#define HMCONT2D_CONSTRUCTOR_HPP

#include "hybmesh_contours2d.hpp"

namespace HMCont2D{ namespace Constructor{

HMCont2D::Container<Contour> Circle(int N, double rad, Point cnt);

HMCont2D::Contour ContourFromPoints(const vector<Point*>& pnt, bool forse_closed=false);

HMCont2D::Container<Contour> ContourFromPoints(vector<double> pnt, bool forse_closed=false);
HMCont2D::Container<Contour> ContourFromPoints(vector<Point> pnt, bool force_closed=false);


}};

#endif
