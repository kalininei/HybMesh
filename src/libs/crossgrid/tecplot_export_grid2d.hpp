#ifndef TECPLOT_EXPORT_GRID2D_HPP
#define TECPLOT_EXPORT_GRID2D_HPP

#include "fluent_export_grid2d.hpp"

namespace GGeom{ namespace Export{ 

void GridTecplot(const GridGeom& g, std::string fn, vector<int> bndindex);

void GridTecplot(const GridGeom& g, std::string fn, vector<int> bndindex, BFun bnames);

}}

#endif
