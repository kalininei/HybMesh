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
	Point p1(0,0);
	Point p2(1,0);
	Point p3(1,1);
	Point p4(0,1);
	auto a = ClosedContour::FromConnectedPoints(
			{p1, p2, p3, p4},
			{0,1,2,3},
			{1,2,3,0});
	add_check(a.size() == 1, "number of contours");
	ContourTree tree;
	for (auto& c: a){
		tree.AddContour(c);
	}
	add_check(tree.NumPoints() == 4, "number of points in a tree");
	SaveVtk(tree, "cont2d_test1.vtk");
}

int main(){
	std::cout<<"hybmesh_contours2d testing"<<std::endl;
	if (hybmesh_contours2d_ping(1) == 2) 
		std::cout<<"Ping OK"<<std::endl;
	test1();
}
