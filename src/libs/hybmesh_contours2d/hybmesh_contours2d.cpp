#include "hybmesh_contours2d.h"
#include "hybmesh_contours2d.hpp"

int hybmesh_contours2d_ping(int a){
	return 2*a;
}

using namespace HMCont2D;

void Contour::AddPointToEnd(double x, double y){
	shared_ptr<Point> newp(new Point(x, y));
	pts.push_back(newp);
}

void Contour::AddPointToEnd(Point xy){
	std::cout<<"DUMMY AddPointToEnd"<<std::endl;
}

bool Contour::HasCrosses(const Contour& c1, const Contour& c2){
	std::cout<<"DUMMY HasCrosses"<<std::endl;
	return false;
}
void Contour::PtsReallocate(){
	std::cout<<"DUMMY PtsReallocate"<<std::endl;
}

BoundingBox Contour::BuildBoundingBox() const{
	if (NumPoints() == 0) { return BoundingBox(0, 0, 0, 0); }
	auto ret = BoundingBox(pts[0]->x, pts[0]->y, pts[0]->x, pts[0]->y);
	for (auto p: pts) ret.WidenWithPoint(*p);
	return ret;
}

ClosedContour ClosedContour::DeepCopy() const{
	ClosedContour c2(*this);
	c2.PtsReallocate();
	return c2;
}

bool ClosedContour::IsValid() const{
	std::cout<<"DUMMY IsValid()"<<std::endl;
	return true;
}

void ContourTree::AddContour(const ClosedContour& cont){
	//check for contour validity
	if (!cont.IsValid()) throw std::runtime_error("Invalid contour was added to ContourTree");
	//check if contour intersects other contours
	for (auto& c: conts){
		if (Contour::HasCrosses(cont, *c))
			throw std::runtime_error("Impossible to add contour to ContourTree due to intersections");
	}
	//add to contours list
	aa::add_shared(conts, cont.DeepCopy());
	//rebuild structure
	RebuildStructure();
}

void ContourTree::RebuildStructure(){
	std::cout<<"DUMMY RebuildStructure"<<std::endl;
}

BoundingBox ContourTree::BuildBoundingBox() const{
	vector<BoundingBox> bd;
	for (auto c: conts) bd.push_back(c->BuildBoundingBox());
	return BoundingBox(bd);
}

