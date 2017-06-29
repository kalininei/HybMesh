#ifndef HYBMESH_HEALGRID_HPP
#define HYBMESH_HEALGRID_HPP

#include "primitives2d.hpp"

namespace HM2D{ namespace Grid{ namespace Algos{

//checks for negative cell areas, self intersections
bool Check(const GridData& g);

//merges coincident vertices, checks rotation
void Heal(GridData& from);

//deletes all nullptr vcells entries and
//recalculates vvert and vedges with respect to it;
void RestoreFromCells(GridData&);

//angle0 - minimum not-allowed angle
void NoConcaveCells(GridData& grid, double angle0=180, bool ignore_bnd=false);


void RemoveShortEdges(GridData& grid, double ref_len);

}}}




#endif
