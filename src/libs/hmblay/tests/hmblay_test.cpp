#include "hmblay.hpp"
#include "fileproc.h"

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
	std::cout<<"01. Boundary Layers grid from circle"<<std::endl;
	std::string cn;
	// 1. Build Edges Collection 
	auto col = HMCont2D::Constructor::Circle(8, 2.0, Point(0, 0));
	// 2. set inp
	HMBlay::Input inp;
	inp.edges = &col;
	inp.direction = HMBlay::DirectionFromString("INNER");
	inp.bnd_step_method = HMBlay::MethFromString("NO");
	inp.partition = {0.0, 0.2, 0.4, 1.5};
	inp.bnd_step = 0.1;
	inp.round_off = true;
	inp.start=inp.end=Point(0,0);
	inp.sharp_angle = inp.corner_angle = 0.0;
	inp.regular_angle = 300;
	inp.start = inp.end = Point(0,0);

	cn = "Inner full circle";
	GridGeom Ans1 = HMBlay::BuildBLayerGrid({inp});
	add_check(fabs(Ans1.area() - 10.606602)<0.00001, cn);
	save_vtk(Ans1, "_dbgout.vtk");
	HMCont2D::SaveVtk(col, "orig.vtk");

	cn = "Outer full circle";
	inp.direction = HMBlay::DirectionFromString("OUTER");
	GridGeom Ans2 = HMBlay::BuildBLayerGrid({inp});
	add_check(fabs(Ans2.area() - 23.3345237) <1e-6, cn);
}

void test02(){
	//cn = "Smoothing normals";
	//inp.bnd_step_method = HMBlay::MethFromString("KEEP_ORIGIN");
	//GridGeom Ans4 = HMBlay::BuildBLayerGrid({inp});
	//add_check(fabs(Ans4.area() - 23.3345237) <1e-6, cn);
	
	//GridGeom Ans5 = HMBlay::BuildBLayerGrid({inp});
	//add_check(fabs(Ans2.area() - 23.3345237) <1e-6, cn);
	//cn = "Circle with different partitions";
	//HMBlay::Input inp2(inp);
	//inp2.partition = {0.0, 0.2, 0.3};
	//inp2.start = inp.end = Point(1,0);
	//inp.bnd_step_method = inp2.bnd_step_method = HMBlay::MethFromString("KEEP_ORIGIN");
	//GridGeom Ans6 = HMBlay::BuildBLayerGrid({inp, inp2});
	//save_vtk(Ans6, "test21.vtk");
};


int main(){
	test01();

	if (FAILED_CHECKS ==1){
		std::cout<<FAILED_CHECKS<<" test failed <<<<<<<<<<<<<<<<<<<"<<std::endl;
	} else if (FAILED_CHECKS > 1) {
		std::cout<<FAILED_CHECKS<<" tests failed <<<<<<<<<<<<<<<<<<<"<<std::endl;
	} else {
		std::cout<<"All tests passed"<<std::endl;
	}
	std::cout<<"DONE"<<std::endl;
}
