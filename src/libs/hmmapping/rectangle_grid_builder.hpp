#ifndef RECTANGLE_GRID_BUILDER_HPP
#define RECTANGLE_GRID_BUILDER_HPP
#include "grid.h"
#include "hybmesh_contours2d.hpp"

namespace HMGMap{

//Contours points will be moved to form a rectangle. No reverses.
GridGeom LinearRectGrid(HMCont2D::Contour& left, HMCont2D::Contour& bot,
	HMCont2D::Contour& right, HMCont2D::Contour& top);

GridGeom ConformalRectGrid(HMCont2D::Contour& left, HMCont2D::Contour& bot,
	HMCont2D::Contour& right, HMCont2D::Contour& top);

GridGeom LaplasRectGrid(HMCont2D::Contour& left, HMCont2D::Contour& bot,
	HMCont2D::Contour& right, HMCont2D::Contour& top);

GridGeom FDMLaplasRectGrid(HMCont2D::Contour& left, HMCont2D::Contour& bot,
	HMCont2D::Contour& right, HMCont2D::Contour& top);

}
#endif

