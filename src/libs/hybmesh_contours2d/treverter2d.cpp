#include "treverter2d.hpp"
using namespace HM2D;
using namespace HM2D::Contour;
using namespace HM2D::Contour::R;

ReallyRevert::ReallyRevert(const EdgeData& ed){
	permanent = false;
	obj = const_cast<EdgeData*>(&ed);
	auto av = OrderedPoints(ed);
	//reverse edges
	reverted_edges.resize(ed.size(), false);
	for (int i=0; i<ed.size(); ++i){
		if (ed[i]->first() == av[i]){
			reverted_edges[i] = true;
			ed[i]->reverse();
		}
	}
	//reverse contour
	Reverse(*obj);
}

ReallyRevert::~ReallyRevert(){
	if (!permanent){
		Reverse(*obj);
		for (int i=0; i<obj->size(); ++i) if (reverted_edges[i]){
			(*obj)[i]->reverse();
		}
	}
}

ReallyDirect::ReallyDirect(const EdgeData& ed){
	permanent = false;
	obj = const_cast<EdgeData*>(&ed);
	auto av = OrderedPoints(ed);
	//reverse edges
	reverted_edges.resize(ed.size(), false);
	for (int i=0; i<ed.size(); ++i){
		if (ed[i]->first() != av[i]){
			reverted_edges[i] = true;
			ed[i]->reverse();
		}
	}
}

ReallyDirect::~ReallyDirect(){
	if (!permanent){
		for (int i=0; i<obj->size(); ++i) if (reverted_edges[i]){
			(*obj)[i]->reverse();
		}
	}
}

ForceFirst::ForceFirst(const EdgeData& ed, Point p0){
	obj = const_cast<EdgeData*>(&ed);
	if (IsOpen(ed)){
		double d0 = Point::meas(p0, *First(ed));
		double d1 = Point::meas(p0, *Last(ed));
		if (d0 < d1){
			really_direct.reset(new ReallyDirect(ed));
		} else {
			really_revert.reset(new ReallyRevert(ed));
		}
		oldstart = 0;
	} else {
		really_direct.reset(new ReallyDirect(ed));
		auto av = OrderedPoints(ed);
		int fc = std::get<0>(FindClosestNode(av, p0));
		std::rotate(obj->begin(), obj->begin() + fc, obj->end());
		oldstart = ed.size() - fc;
	}
}

ForceFirst::~ForceFirst(){
	std::rotate(obj->begin(), obj->begin()+oldstart, obj->end());
}

Clockwise::Clockwise(const EdgeData& ed, bool direct){
	assert(IsClosed(ed));
	bool dirnow = Contour::Area(ed) < 0;
	if (dirnow == direct){
		really_direct.reset(new ReallyDirect(ed));
	} else {
		really_revert.reset(new ReallyRevert(ed));
	}
}

RevertTree::RevertTree(const Tree& tree){
	for (auto& n: tree.nodes){
		if (n->isdetached()){
			really_direct.emplace_back(
				new ReallyDirect(n->contour)
			);
		} else {
			really_clockwise.emplace_back(
				new Clockwise(n->contour, n->level % 2 == 1)
			);
		}
	}
}
