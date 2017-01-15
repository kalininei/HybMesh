#ifndef HYBMESH_INFOGRID_HPP
#define HYBMESH_INFOGRID_HPP
#include "primitives2d.hpp"

namespace HM2D{ namespace Grid{

//area of grid bounding tree
double Area(const GridData& grid);

vector<double> CellAreas(const GridData& grid);

//calculate skewness
vector<double> Skewness(const GridData& grid);

}}
#endif
