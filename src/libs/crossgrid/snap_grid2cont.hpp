#ifndef HYBMESH_SNAPGRID_HPP
#define HYBMESH_SNAPGRID_HPP

#include "primitives2d.hpp"

namespace HM2D{ namespace Grid{ namespace Algos{
// If boundary edge points (p1, p2) lie on contour then
// (all/significant) contour points between (p1, p2) will present in grid
// snap_nodes is a list of points which will be snapped to contour before procedure starts
// !!! @grid should be located to the left of the @cont
void SnapToContour(GridData& grid, const EdgeData& cont,
		const VertexData& snap_nodes, bool only_significant=true);

//shifts boundary grid node to significant contour vertex
//if it is non-significant by itself. Otherwise does nothing
void ShiftToContour(GridData& grid, const EdgeData& cont,
		const VertexData& snap_nodes);


// takes grid boundary section from gp1 to gp2 and snaps it to contour.
// gp1 point will be shifted to contour first point;
// gp2 point will be shifted to contour last point;
// intermediate points will be snapped according to their weight coordinates.
// if snap_strategy == 'shift', intermediate points will be then shifted to closest
//   contour vertex.
// gp1, gp2 are @grid boundary point pointers.
// if gp1 == gp2 then cont should be a closed contour else open.
//
// Grid contour will be treated counterclockwise regardless its nesting.
GridData SnapBndSection(const GridData& grid, const EdgeData& cont,
	Vertex* gp1, Vertex* gp2, std::string snap_strategy);
}}}

#endif
