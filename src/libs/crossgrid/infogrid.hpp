#ifndef HYBMESH_INFOGRID_HPP
#define HYBMESH_INFOGRID_HPP
#include "primitives2d.hpp"
#include "tree.hpp"

namespace HM2D{ namespace Grid{

//area of grid bounding tree
double Area(const GridData& grid);

//areas cell by cell
vector<double> CellAreas(const GridData& grid);

//extracts cells which are fully inside (what = INSIDE) or outside (what = OUTSIDE) of given domain
CellData ExtractCells(const GridData& grid, const Contour::Tree& domain, int what);
CellData ExtractCells(const GridData& grid, const EdgeData& domain, int what);

//calculate skewness
vector<double> Skewness(const GridData& grid);

}}
#endif
