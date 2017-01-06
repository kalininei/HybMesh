#ifndef HYBMESH_TREVERTER2D_HPP
#define HYBMESH_TREVERTER2D_HPP
#include "contour.hpp"

namespace HM2D{ namespace Contour{

//this reverses edges and
//forces correct direction of points within edges
class ReallyRevert{
	EdgeData* obj;
	std::vector<bool> reverted_edges;
	bool permanent;
public:
	//delete default constructors to avoid premature reversion.
	ReallyRevert(const ReallyRevert&) = delete;
	ReallyRevert(ReallyRevert&&) = delete;
	ReallyRevert& operator=(const ReallyRevert&) = delete;

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
	ReallyDirect(ReallyDirect&&) = delete;
	ReallyDirect& operator=(const ReallyDirect&) = delete;

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
	ForceFirst(ForceFirst&&) = delete;
	ForceFirst& operator=(const ForceFirst&) = delete;

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

}}
#endif
