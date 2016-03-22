#include "hmmapping.hpp"
#include "gridmap.hpp"
using namespace HMGMap;

namespace{
GridGeom scale_base_grid(const GridGeom& b, ScaleBase& sc){
	GridGeom ret = GGeom::Constructor::DeepCopy(b);
	sc = ret.do_scale();
	return ret;
}

HMCont2D::Container<HMCont2D::ECollection>
scale_mapped_domain(const HMCont2D::ECollection& a, ScaleBase& sc){
	HMCont2D::Container<HMCont2D::ECollection> ret =
		HMCont2D::Container<HMCont2D::ECollection>::DeepCopy(a);
	sc = HMCont2D::ECollection::Scale01(ret);
	return ret;
}
}

GridGeom HMGMap::MapGrid(const GridGeom& base, const HMCont2D::ECollection& area,
		vector<Point> base_points,
		vector<Point> mapped_points,
		Options opt){
	//1) scale
	ScaleBase bscale, cscale;
	GridGeom g = scale_base_grid(base, bscale);
	HMCont2D::Container<HMCont2D::ECollection> c = scale_mapped_domain(area, cscale);
	bscale.scale(base_points.begin(), base_points.end());
	cscale.scale(mapped_points.begin(), mapped_points.end());
	//2) set input data
	HMGMap::Impl::DoMapping dm(opt);
	dm.set_grid(g);
	dm.set_contour(c);
	dm.set_points(base_points, mapped_points);
	//3) build
	GridGeom ret = dm.run();
	//4) unscale and return
	ret.undo_scale(cscale);
	return ret;
}
