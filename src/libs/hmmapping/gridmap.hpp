#ifndef HYBMESH_GRIDMAP_HPP
#define HYBMESH_GRIDMAP_HPP

#include "hmproject.h"
#include "primitives2d.hpp"
#include "hmcallback.hpp"

namespace HMMap{

class MapException: public std::runtime_error{
public:
	MapException(std::string m) noexcept: std::runtime_error(
			std::string("Grid mapping exception: ") + m){};
};

class EInvalidGrid: public MapException{
public:
	HM2D::GridData invalid_grid;
	EInvalidGrid(HM2D::GridData&& g) noexcept:
		MapException("Resulting grid is not valid"),
		invalid_grid(std::move(g)){};
};

struct Options{
	int fem_nmax;
	int fem_nmin;
	int fem_nrec;

	std::string snap;   //NO, ADD_VERTICES, SHIFT_VERTICES
	std::string algo;   //direct-laplace, inverse-laplace
	bool btypes_from_contour; //or from grid

	Options(std::string _algo="inverse-laplace", std::string _snap="NO",
			bool btypes_from_contour=true):
		fem_nmax(100000), fem_nmin(100), fem_nrec(1000),
		snap(_snap), algo(_algo),
		btypes_from_contour(btypes_from_contour){}
};

struct TMapGrid: public HMCallback::ExecutorBase{
	HMCB_SET_PROCNAME("Grid mapping");
	HMCB_SET_DEFAULT_DURATION(110);

	HM2D::GridData _run(const HM2D::GridData& base,
			const HM2D::EdgeData& area,
			vector<Point> base_points,
			vector<Point> mapped_points,
			bool reversed,
			Options opt);
};
extern HMCallback::FunctionWithCallback<TMapGrid> MapGrid;

}
#endif
