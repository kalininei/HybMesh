#ifndef HYBMESH_BGRID_IMPOSE_HPP
#define HYBMESH_BGRID_IMPOSE_HPP

#include "bgrid.hpp"

namespace HMBlay{namespace Impl{

//Main procedure.
//Does self imposition of a grid using priority function.
//Weeps out all cells which lie to the right of the source contour.
void BGridImpose(BGrid& grid, std::function<double(const HM2D::Cell*)> prifun,
		const HM2D::EdgeData& source);

/*
class NoHangingNodes{
	BGrid* grid;
	std::unique_ptr<BoundingBoxFinder> cfinder;

	//Finds cell which should contain (but doesn't) handing node p.
	//Returns c - found cell, ind - local index of cell edge
	//which should be splitted by this node
	void FindCellForNode(HM2D::Vertex* p, HM2D::Cell*& c, int& ind);

	//places point p to a cell at edge ind.
	//new node local index will be ind+1
	//void PlaceNode(HM2D::Vertex* p, HM2D::Cell* c, int ind);

	//splits cell connecting points p1 and p2
	void CellSplit(shared_ptr<HM2D::Cell> cell, shared_ptr<HM2D::Vertex> p1, shared_ptr<HM2D::Vertex> p2);

	//divides cell by connecting ind-th node
	//with closest not collinear
	void SimpleSplit(shared_ptr<HM2D::Cell> cell, int ind);

	//divides cell into to using ind-th point
	//as a start point
	void SimplifyCell(shared_ptr<HM2D::Cell> cell, int ind);

	//processes hanging node p
	void AddHangingNode(HM2D::Vertex* p);

	//for c - pseudo contour made of hanging nodes
	//adds all such nodes to list
	void AnalyzeSingularContour(const HM2D::EdgeData& c, std::list<HM2D::Vertex*>& lst);

	//finds list of all hanging nodes
	std::list<HM2D::Vertex*> Find();
public:
	//Modifies _grid removing hanging nodes:
	//nodes which don't belong to all adjacent cells
	NoHangingNodes(BGrid& _grid);
};

*/

}}



#endif
