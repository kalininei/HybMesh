#ifndef HYBMESH_TREVERTER2D_HPP
#define HYBMESH_TREVERTER2D_HPP
#include "contour.hpp"
#include "tree.hpp"

namespace HM2D{
	
namespace Contour{ namespace R{

//this reverses edges and
//forces correct direction of points within edges
class ReallyRevert{
	EdgeData* obj;
	std::vector<bool> reverted_edges;
	bool permanent;
public:
	//delete default constructors to avoid premature reversion.
	ReallyRevert(const ReallyRevert&) = delete;

	ReallyRevert(const EdgeData& ed);
	~ReallyRevert();

	void make_permanent(){permanent = true; }

	static void Permanent(EdgeData& ed){
		ReallyRevert a(ed);
		a.make_permanent();
	}
};

class ReallyDirect{
	EdgeData* obj;
	std::vector<bool> reverted_edges;
	bool permanent;
public:
	//delete default constructors to avoid premature reversion.
	ReallyDirect(const ReallyDirect&) = delete;

	ReallyDirect(const EdgeData& ed);
	~ReallyDirect();

	void make_permanent(){permanent = true; }

	static void Permanent(EdgeData& ed){
		ReallyDirect a(ed);
		a.make_permanent();
	}
};

//sets first point of contour closest to given one.
//for open contours chooses only between first and last points.
//makes reversions of all edges.
class ForceFirst{
	EdgeData* obj;
	int oldstart;
	std::unique_ptr<ReallyDirect> really_direct;
	std::unique_ptr<ReallyRevert> really_revert;
public:
	//delete default constructors to avoid premature reversion.
	ForceFirst(const ForceFirst&) = delete;

	ForceFirst(const EdgeData& ed, Point p0);
	~ForceFirst();

	void make_permanent(){
		if (really_direct) really_direct->make_permanent();
		if (really_revert) really_revert->make_permanent();
		oldstart=0;
	}

	static void Permanent(EdgeData& ed, Point p0){
		ForceFirst a(ed, p0);
		a.make_permanent();
	}
};

//forces clockwise or counterclockwise direction of closed contour
//for open contours makes direct reverse
class Clockwise{
	std::unique_ptr<ReallyDirect> really_direct;
	std::unique_ptr<ReallyRevert> really_revert;
public:
	Clockwise(const Clockwise&) = delete;

	//direct = true -> clockwise;
	//direct = false -> counterclockwise
	Clockwise(const EdgeData& ed, bool direct);
	~Clockwise(){}

	void make_permanent(){
		if (really_direct) really_direct->make_permanent();
		if (really_revert) really_revert->make_permanent();
	}

	static void Permanent(EdgeData& ed, bool direct){
		Clockwise c(ed, direct);
		c.make_permanent();
	}
};

//forces counterclockwise for even level tree nodes
//and clockwise for odd level tree nodes
class RevertTree{
	std::list<std::unique_ptr<Clockwise>> really_clockwise;
	std::list<std::unique_ptr<ReallyDirect>> really_direct;
public:
	RevertTree(const RevertTree&) = delete;

	RevertTree(const Tree& tree);
	~RevertTree(){}

	void make_permanent(){
		for (auto& it:really_clockwise) it->make_permanent();
		for (auto& it:really_direct) it->make_permanent();
	}

	static void Permanent(Tree& tree){
		RevertTree r(tree);
		r.make_permanent();
	}

};

//all edges which has single cell connection will be
//reverted so that it is a left cell
class LeftCells{
	vector<int> reverted_edges_inds;
	EdgeData* obj;
public:
	LeftCells(const LeftCells&) = delete;

	LeftCells(const EdgeData& ed);
	~LeftCells();

	void make_permanent(){
		reverted_edges_inds.clear();
	}

	static void Permanent(EdgeData& ed){
		LeftCells r(ed);
		r.make_permanent();
	}
};

}}}
#endif
