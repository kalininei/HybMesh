#ifndef TECPLOT_EXPORT_GRID2D_HPP
#define TECPLOT_EXPORT_GRID2D_HPP

#include "export2d_fluent.hpp"

namespace HM2D{ namespace Export{ 

void GridTecplot(const GridData& g, std::string fn, vector<int> bndindex);

void GridTecplot(const GridData& g, std::string fn, vector<int> bndindex, BNamesFun bnames);

}}

#endif

