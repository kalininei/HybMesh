#ifndef RECTANGLE_GRID_BUILDER_HPP
#define RECTANGLE_GRID_BUILDER_HPP
#include "grid.h"
#include "primitives2d.hpp"
#include "gridmap.hpp"

namespace HMMap{

//Contours points will be moved to form a rectangle. No reverses.
GridGeom LinearRectGrid(HM2D::EdgeData& left, HM2D::EdgeData& bot,
	HM2D::EdgeData& right, HM2D::EdgeData& top);

struct TOrthogonalRectGrid: public HMCallback::ExecutorBase{
	HMCB_SET_PROCNAME("Orthogonal custom rectangular grid");
	HMCB_SET_DEFAULT_DURATION(100);

	GridGeom _run(HM2D::EdgeData& left, HM2D::EdgeData& bot,
		HM2D::EdgeData& right, HM2D::EdgeData& top);
};
extern HMCallback::FunctionWithCallback<TOrthogonalRectGrid> OrthogonalRectGrid;


//algo is similar to HMGMap::Options::algo = {inverse-laplace, direct-laplace}
struct TLaplaceRectGrid: public HMCallback::ExecutorBase{
	HMCB_SET_PROCNAME("Custom rectangular grid");
	HMCB_SET_DEFAULT_DURATION(100);

	GridGeom _run(HM2D::EdgeData& left, HM2D::EdgeData& bot,
		HM2D::EdgeData& right, HM2D::EdgeData& top, std::string algo);
};
extern HMCallback::FunctionWithCallback<TLaplaceRectGrid> LaplaceRectGrid;

GridGeom FDMLaplasRectGrid(HM2D::EdgeData& left, HM2D::EdgeData& bot,
	HM2D::EdgeData& right, HM2D::EdgeData& top);

GridGeom LinearTFIRectGrid(HM2D::EdgeData& left, HM2D::EdgeData& bot,
		HM2D::EdgeData& right, HM2D::EdgeData& top);

GridGeom LinearTFIRectGrid2(HM2D::EdgeData& left, HM2D::EdgeData& bot,
		HM2D::EdgeData& right, HM2D::EdgeData& top);

//c = [cleft, cbot, cright, ctop]
GridGeom CubicTFIRectGrid(HM2D::EdgeData& left, HM2D::EdgeData& bot,
		HM2D::EdgeData& right, HM2D::EdgeData& top, std::array<double, 4> c);

}
#endif

