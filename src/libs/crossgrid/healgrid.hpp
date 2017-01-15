#ifndef HYBMESH_HEALGRID_HPP
#define HYBMESH_HEALGRID_HPP

#include "primitives2d.hpp"

namespace HM2D{ namespace Grid{ namespace Algos{

//checks for negative cell areas, self intersections
bool Check(const GridData& g);

//if boundary edge points (p1, p2) lie on contour then
//all significant contour points between (p1, p2) will present in grid
//grid should be located to the left of the contour
//snap_nodes is a list of points which will be snapped to contour before procedure starts
void SnapToContour(GridData& grid, const EdgeData& cont,
		const VertexData& snap_nodes);

//shifts boundary grid node to significant contour vertex
//if it is non-significant by itself. Otherwise does nothing
void ShiftToContour(GridData& grid, const EdgeData& cont,
		const VertexData& snap_nodes);

//merges coincident vertices, checks rotation
void Heal(GridData& from);

//deletes all nullptr vcells entries and
//recalculates vvert and vedges with respect to it;
void RestoreFromCells(GridData&);

//angle 0 - minimum not-allowed angle
void NoConcaveCells(GridData& grid, double angle0=180);


void RemoveShortEdges(GridData& grid, double ref_len);

}}}




#endif
