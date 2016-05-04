#ifndef RECTANGLE_GRID_BUILDER_HPP
#define RECTANGLE_GRID_BUILDER_HPP
#include "grid.h"
#include "hybmesh_contours2d.hpp"
#include "gridmap.hpp"

namespace HMGMap{

//Contours points will be moved to form a rectangle. No reverses.
GridGeom LinearRectGrid(HMCont2D::Contour& left, HMCont2D::Contour& bot,
	HMCont2D::Contour& right, HMCont2D::Contour& top);

struct TOrthogonalRectGrid: public HMCallback::ExecutorBase{
	HMCB_SET_PROCNAME("Orthogonal custom rectangular grid");
	HMCB_SET_DEFAULT_DURATION(100);

	GridGeom _run(HMCont2D::Contour& left, HMCont2D::Contour& bot,
		HMCont2D::Contour& right, HMCont2D::Contour& top);
};
extern HMCallback::FunctionWithCallback<TOrthogonalRectGrid> OrthogonalRectGrid;


//algo is similar to HMGMap::Options::algo = {inverse-laplace, direct-laplace}
struct TLaplaceRectGrid: public HMCallback::ExecutorBase{
	HMCB_SET_PROCNAME("Custom rectangular grid");
	HMCB_SET_DEFAULT_DURATION(100);

	GridGeom _run(HMCont2D::Contour& left, HMCont2D::Contour& bot,
		HMCont2D::Contour& right, HMCont2D::Contour& top, std::string algo);
};
extern HMCallback::FunctionWithCallback<TLaplaceRectGrid> LaplaceRectGrid;

GridGeom FDMLaplasRectGrid(HMCont2D::Contour& left, HMCont2D::Contour& bot,
	HMCont2D::Contour& right, HMCont2D::Contour& top);

}
#endif

