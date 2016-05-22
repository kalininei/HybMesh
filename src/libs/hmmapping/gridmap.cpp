#include "gridmap.hpp"
#include "domapping.hpp"

HMCallback::FunctionWithCallback<HMMap::TMapGrid> HMMap::MapGrid;

GridGeom HMMap::TMapGrid::_run(const GridGeom& base,
		const HMCont2D::ECollection& area,
		vector<Point> base_points,
		vector<Point> mapped_points,
		Options opt){
	//1) scale
	callback->step_after(5, "Initializing");
	//2) set input data
	shared_ptr<HMMap::Impl::DoMapping> dm;
	if (opt.algo == "inverse-laplace"){
		dm.reset(new HMMap::Impl::InverseMapping(opt));
	} else if (opt.algo == "direct-laplace"){
		dm.reset(new HMMap::Impl::DirectMapping(opt));
	} else {
		assert(false);
	}
	dm->set_grid(base);
	dm->set_contour(area);
	dm->set_points(base_points, mapped_points);
	//3) build
	GridGeom ret = dm->run(*callback);
	//4) unscale and return
	callback->step_after(5, "Finilizing");
	return ret;
}
