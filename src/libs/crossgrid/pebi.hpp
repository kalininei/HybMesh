#ifndef HYBMESH_PEBI_HPP
#define HYBMESH_PEBI_HPP
#include "hmcallback.hpp"
#include "primitives2d.hpp"
#include "tree.hpp"

namespace HM2D{

namespace Mesher{

struct TUnstructuredPebi: public HMCallback::ExecutorBase{
	HMCB_SET_PROCNAME("Pebi grid building");
	HMCB_SET_DEFAULT_DURATION(110);

	GridData _run(const Contour::Tree& source);
	GridData _run(const Contour::Tree& source, const std::map<Point, double>& embedded);
};
extern HMCallback::FunctionWithCallback<TUnstructuredPebi> UnstructuredPebi;

}

namespace Grid{ namespace Constructor{

GridData TriToPebi(const GridData& g);
GridData RegularHexagonal(Point cnt, double area_rad, double cell_rad, bool strict_area=false);
GridData RegularHexagonal(Point cnt1, Point cnt2, double cell_rad, bool strict_area=false);

}}

}
#endif

