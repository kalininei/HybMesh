#include "hmgrid3d.hpp"
#include "procgrid.h"

int FAILED_CHECKS = 0;

void add_check(bool ex, std::string info){
	if (info.size()==0){
		std::cout<<"\tunknown check: ";
	} else{
		std::cout<<"\t"<<info;
	}
	if (ex){
		std::cout<<": True"<<std::endl;
	} else {
		++FAILED_CHECKS;
		std::cout<<": False <<<<<<<<<<<<<<<<<<<"<<std::endl;
	}
};

void test01(){
	std::cout<<"export import cuboid"<<std::endl;
	auto g1 = HMGrid3D::Constructor::Cuboid({0, 0, 0}, 1, 2, 5, 3, 3, 3);
	HMGrid3D::Export::BoundaryVTK(g1, "c1.vtk");
	HMGrid3D::Export::GridVTK(g1, "g1.vtk");
	add_check(g1.n_vertices() == 64 && g1.n_cells() == 27 &&
			g1.n_edges() == 144 && g1.n_faces() == 108,
			"cuboid primitives number"); 
}

void test02(){
	std::cout<<"rectangular grid parallel sweep"<<std::endl;
	auto g2d = GGeom::Constructor::RectGrid01(10, 13);
	auto g3d = HMGrid3D::Constructor::SweepGrid2D(g2d, {0.3, 0.4, 0.8});
	HMGrid3D::Export::BoundaryVTK(g3d, "c1.vtk");
	HMGrid3D::Export::GridVTK(g3d, "g1.vtk");
}

int main(){
	test01();
	test02();
	
	if (FAILED_CHECKS ==1){
		std::cout<<FAILED_CHECKS<<" test failed <<<<<<<<<<<<<<<<<<<"<<std::endl;
	} else if (FAILED_CHECKS > 1) {
		std::cout<<FAILED_CHECKS<<" tests failed <<<<<<<<<<<<<<<<<<<"<<std::endl;
	} else {
		std::cout<<"All tests passed"<<std::endl;
	}
	std::cout<<"DONE"<<std::endl;
}
