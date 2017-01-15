#ifndef HYBMESH_MODGRID_HPP
#define HYBMESH_MODGRID_HPP
#include "primitives2d.hpp"
#include "tree.hpp"

namespace HM2D{ namespace Grid{ namespace Algos{

//adds 'from' data to 'to'.
//merges coincident boundary primitives.
//duplicate primitives will be taken from `from` grid
void MergeTo(const GridData& from, GridData& to);

//adds non-repeating (by pointer) vertices, edges and cells
void ShallowAdd(const GridData& from, GridData& to);

//splits i-th cell into two cells connecting lnode1-lnode2 where
//lnode1, lnode2 local vertex indicies
//if lnode1 = -1 then calculates it by largest angle
//if lnode2 = -1 then finds best fitted second node
//adds new cell and new edge to the end of grids vcells and vedges table.
//if impossible -> returns false
bool SplitCell(GridData& grid, int icell, int lnode1=-1, int lnode2=-1);

//splits ith edge pasting new points. i-th edge object will not be deleted from the grid but shrinked.
//writes result to grid only if no cell self-intersections were found until force=true
//returns: if edge was splitted
bool SplitEdge(GridData& grid, int iedge, const vector<Point>& apoints, bool force=false);

//no complicated boundary cell edges
//angle is between [0, 180]
//0 -- deletes only non-significant edge points
//180 -- deletes all intermediate points
void SimplifyBoundary(GridData& grid, double angle);

//all cells with dimensin higher than maxdim will be splitted
void CutCellDims(GridData& grid, int maxdim);

//extracts cells which are fully inside (what = INSIDE) or outside (what = OUTSIDE) of given domain
CellData ExtractCells(const GridData& grid, const Contour::Tree& domain, int what);
//removes cells by index
void RemoveCells(GridData& grid, const vector<int>& icells);
//removes cells which has points lying inside (what = INSIDE) or outside (what = OUTSIDE) of given domain
void RemoveCells(GridData& grid, const Contour::Tree& domain, int what);

}}}
#endif
