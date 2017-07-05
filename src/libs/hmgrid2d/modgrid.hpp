#ifndef HYBMESH_MODGRID_HPP
#define HYBMESH_MODGRID_HPP
#include "primitives2d.hpp"
#include "contour_tree.hpp"

namespace HM2D{ namespace Grid{ namespace Algos{

//sort all vertices by coordinates
//edges -> by lowest vertex index
//cells -> by lowest edge index
//make all edges start with lowest indexed vertex
//make all cell->edges start with lowest indexed edge
void UniqueRearrange(GridData& from);

//adds 'from' data to 'to'.
//merges coincident boundary primitives.
//duplicate primitives will be taken from `from` grid
void MergeTo(const GridData& from, GridData& to);

//does the same as MergeTo but also places all 'from' boundary nodes
//to 'to' boundary if they lie on the boundary but do not equal 'to' vertex
//Unlike MergeTo this doesn't copy 'from' primitives but creates deep copied ones.
void MergeBoundaries(const GridData& from, GridData& to);

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

//all cells with dimension higher than maxdim will be splitted
void CutCellDims(GridData& grid, int maxdim);

//removes edges by index and restore grid topology
//for inner edges left cell will be preserved and concatenated with the right one.
//for boundary edges removes whole adjacent cell
void RemoveEdges(GridData& grid, vector<int> iedges);

//removes cells by index and restore grid topology
void RemoveCells(GridData& grid, vector<int> icells);
void RemoveCellsById(GridData& grid, int id);
//removes cells which has points lying inside (what = INSIDE) or outside (what = OUTSIDE) of given domain
void RemoveCells(GridData& grid, const Contour::Tree& domain, int what);

}}}
#endif
