#ifndef HYBMESH_UNITE_GRIDS_HPP
#define HYBMESH_UNITE_GRIDS_HPP

#include "primitives2d.hpp"
#include "hmcallback.hpp"
#include "contour_tree.hpp"

namespace HM2D{ namespace Grid{ namespace Algos{

struct OptUnite{
	OptUnite(double buffer_size=0,
		bool preserve_bp=false,
		bool empty_holes=false,
		double angle0=0,
		int filler=0): buffer_size(buffer_size), preserve_bp(preserve_bp),
			       empty_holes(empty_holes), angle0(angle0), filler(filler){}
	double buffer_size;
	bool preserve_bp;
	bool empty_holes;
	double angle0;
	int filler;  //0 - triangles, 1 - recombined.
};


struct TUniteGrids: public HMCallback::ExecutorBase{
	HMCB_SET_PROCNAME("Unite grids");
	HMCB_SET_DEFAULT_DURATION(100);

	GridData _run(const GridData& base, const GridData& sec, const OptUnite& opt);
};
extern HMCallback::FunctionWithCallback<TUniteGrids> UniteGrids;


struct TCombineGrids: public HMCallback::ExecutorBase{
	HMCB_SET_PROCNAME("Combine grids");
	HMCB_SET_DEFAULT_DURATION(100);

	GridData _run(const GridData& g1, const GridData& g2, bool keep_g2_holes=false);
};
extern HMCallback::FunctionWithCallback<TCombineGrids> CombineGrids;
}}}

#endif
