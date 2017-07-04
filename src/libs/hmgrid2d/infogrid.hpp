#ifndef HYBMESH_INFOGRID_HPP
#define HYBMESH_INFOGRID_HPP
#include "primitives2d.hpp"
#include "contour_tree.hpp"
#include "hmcallback.hpp"

namespace HM2D{ namespace Grid{

//area of grid bounding tree
double Area(const GridData& grid);

//areas cell by cell
vector<double> CellAreas(const GridData& grid);

//extracts cells which are fully inside (what = INSIDE) or
//fully outside (what = OUTSIDE) of given domain
//if (what == BOUND) then cells which cross domain will be extracted.
//!!! Contact here is not a cross.
struct TExtractCells: public HMCallback::ExecutorBase{
	HMCB_SET_PROCNAME("Extract cells");
	HMCB_SET_DEFAULT_DURATION(100);

	CellData _run(const GridData& grid, const Contour::Tree& domain, int what);
	CellData _run(const GridData& grid, const EdgeData& domain, int what);
};
extern HMCallback::FunctionWithCallback<TExtractCells> ExtractCells;

//calculate skewness
vector<double> Skewness(const GridData& grid);

}}
#endif
