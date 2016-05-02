#ifndef HYBMESH_GRIDMAP_HPP
#define HYBMESH_GRIDMAP_HPP

#include "hmproject.h"
#include "hybmesh_contours2d.hpp"
#include "grid.h"

namespace HMGMap{

class MapException: public std::runtime_error{
public:
	MapException(std::string m) noexcept: std::runtime_error(
			std::string("Grid mapping exception: ") + m){};
};

struct Options{
	int fem_nmax;
	int fem_nmin;
	int fem_nrec;
	int fem_nedge;

	std::string snap;   //NO, ADD_VERTICES, SHIFT_VERTICES
	std::string algo;   //direct-laplace, inverse-laplace

	Options(std::string _algo="inverse-laplace", std::string _snap="NO"):
		fem_nmax(100000), fem_nmin(100), fem_nrec(1000), fem_nedge(3),
		snap(_snap), algo(_algo){}
};

struct TMapGrid: public HMCallback::ExecutorBase{
	HMCB_SET_PROCNAME("Domain mapping");
	HMCB_SET_DEFAULT_DURATION(110);

	GridGeom _run(const GridGeom& base,
			const HMCont2D::ECollection& area,
			vector<Point> base_points,
			vector<Point> mapped_points,
			Options opt);
};
extern HMCallback::FunctionWithCallback<TMapGrid> MapGrid;

}
#endif
