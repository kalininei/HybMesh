#include <iostream>
#include "hybmesh_contours2d.h"
#include "hybmesh_contours2d.hpp"

using namespace HMCont2D;

void add_check(bool ex, const char* info = 0){
	if (info==0){
		std::cout<<"\tunknown check: ";
	} else{
		std::cout<<"\t"<<info;
	}
	if (ex){
		std::cout<<": True"<<std::endl;
	} else {
		std::cout<<": False <<<<<<<<<<<<<<<<<<<"<<std::endl;
	}
};

void test1(){
	std::cout<<"ContourTree building"<<std::endl;
	Point p1(0,0), p2(1,0), p3(1,1), p4(0,1);
	Point p5(0.3,0.3), p6(0.5,0.3), p7(0.5,0.5), p8(0.3,0.5);
	auto a = ClosedContour::FromConnectedPoints(
			{p1, p2, p3, p4, p5, p6, p7, p8},
			{0,1,2,3,4,5,6,7},
			{1,2,3,0,5,6,7,4});
	add_check(a.size() == 2, "number of contours");
	ContourTree tree;
	for (auto& c: a) tree.AddContour(c);
	add_check(tree.NumContours() == 2, "number of contours in the tree");
	add_check(tree.NumPoints() == 8, "number of points in the tree");
	add_check(fabs(tree.Area() - 0.096) < 1e-12, "Tree area");

	std::cout<<a[0].SignedArea()<<std::endl;
	std::cout<<a[1].SignedArea()<<std::endl;
	if (a[1].IsWithinGeom(Point(0.4, 0.4))){
		std::cout<<"within"<<std::endl;
	} else {
		std::cout<<"without"<<std::endl;
	}
	a[1].Reverse();
	std::cout<<a[1].SignedArea()<<std::endl;
	if (a[1].IsWithinGeom(Point(0.9, 0.9))){
		std::cout<<"within"<<std::endl;
	} else {
		std::cout<<"without"<<std::endl;
	}

	SaveVtk(tree, "cont2d_test1.vtk");
}

int main(){
	std::cout<<"hybmesh_contours2d testing"<<std::endl;
	if (hybmesh_contours2d_ping(1) == 2) 
		std::cout<<"Ping OK"<<std::endl;
	test1();
}
