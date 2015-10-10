#ifndef HYBMESH_BGRID_IMPOSE_HPP
#define HYBMESH_BGRID_IMPOSE_HPP

#include "bgrid.hpp"

namespace HMBlay{namespace Impl{

//Main procedure.
//Does self imposition of a grid using priority function.
//Weeps out all cells which lie to the right of the source contour.
shared_ptr<BGrid> BGridImpose(shared_ptr<BGrid> grid, std::function<double(const Cell*)> prifun,
		const HMCont2D::Contour& source);

//Builds a contour from cell
HMCont2D::Contour Cell2Cont(const Cell* c);

class NoHangingNodes{
	BGrid* grid;
	GGeom::Info::CellFinder cfinder;

	//Finds cell which should contain (but doesn't) handing node p.
	//Returns c - found cell, ind - local index of cell edge
	//which should be splitted by this node
	void FindCellForNode(GridPoint* p, Cell*& c, int& ind);

	//places point p to a cell at edge ind.
	//new node local index will be ind+1
	void PlaceNode(GridPoint* p, Cell* c, int ind);

	//splits cell connecting points p1 and p2
	void CellSplit(Cell* cell, const GridPoint* p1, const GridPoint* p2);

	//divides cell by connecting ind-th node
	//with closest not collinear
	void SimpleSplit(Cell* cell, int ind);

	//divides cell into to using ind-th point
	//as a start point
	void SimplifyCell(Cell* cell, int ind);

	//processes hanging node p
	void AddHangingNode(GridPoint* p);

	//for c - pseudo contour made of hanging nodes
	//adds all such nodes to list
	void AnalyzeSingularContour(const HMCont2D::Contour& c, std::list<GridPoint*>& lst);

	//finds list of all hanging nodes
	std::list<GridPoint*> Find();
public:
	//Modifies _grid removing hanging nodes:
	//nodes which don't belong to all adjacent cells
	NoHangingNodes(BGrid& _grid);
};


}}



#endif
