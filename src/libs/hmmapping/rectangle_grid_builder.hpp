#ifndef RECTANGLE_GRID_BUILDER_HPP
#define RECTANGLE_GRID_BUILDER_HPP
#include "primitives2d.hpp"
#include "gridmap.hpp"

namespace HMMap{

//Contours points will be moved to form a rectangle. No reverses.
HM2D::GridData LinearRectGrid(HM2D::EdgeData& left, HM2D::EdgeData& bot,
	HM2D::EdgeData& right, HM2D::EdgeData& top);

struct TOrthogonalRectGrid: public HMCallback::ExecutorBase{
	HMCB_SET_PROCNAME("Orthogonal custom rectangular grid");
	HMCB_SET_DEFAULT_DURATION(100);

	HM2D::GridData _run(HM2D::EdgeData& left, HM2D::EdgeData& bot,
		HM2D::EdgeData& right, HM2D::EdgeData& top);
};
extern HMCallback::FunctionWithCallback<TOrthogonalRectGrid> OrthogonalRectGrid;


//algo is similar to HMGMap::Options::algo = {inverse-laplace, direct-laplace}
struct TLaplaceRectGrid: public HMCallback::ExecutorBase{
	HMCB_SET_PROCNAME("Custom rectangular grid");
	HMCB_SET_DEFAULT_DURATION(100);

	HM2D::GridData _run(HM2D::EdgeData& left, HM2D::EdgeData& bot,
		HM2D::EdgeData& right, HM2D::EdgeData& top, std::string algo);
};
extern HMCallback::FunctionWithCallback<TLaplaceRectGrid> LaplaceRectGrid;

HM2D::GridData FDMLaplaceRectGrid(HM2D::EdgeData& left, HM2D::EdgeData& bot,
	HM2D::EdgeData& right, HM2D::EdgeData& top);

HM2D::GridData LinearTFIRectGrid(HM2D::EdgeData& left, HM2D::EdgeData& bot,
		HM2D::EdgeData& right, HM2D::EdgeData& top);

HM2D::GridData LinearTFIRectGrid2(HM2D::EdgeData& left, HM2D::EdgeData& bot,
		HM2D::EdgeData& right, HM2D::EdgeData& top);

//c = [cleft, cbot, cright, ctop]
HM2D::GridData CubicTFIRectGrid(HM2D::EdgeData& left, HM2D::EdgeData& bot,
		HM2D::EdgeData& right, HM2D::EdgeData& top, std::array<double, 4> c);

}
#endif

