#include <iostream>
#include "crossgrid.h"
#include <vector>

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

//build a rectangular structured grid 
Grid* rectangular_grid(double x0, double y0,
		double x1, double y1, int Nx, int Ny){
	double hx = (x1 - x0)/Nx;
	double hy = (y1 - y0)/Ny;
	//points
	std::vector<double> pts;
	for (int j=0; j<Ny+1; ++j){
		for (int i=0;i<Nx+1;++i){
			pts.push_back(i*hx+x0);
			pts.push_back(j*hy+y0);
		}
	}
	//cells
	auto pts_ind = [Nx, Ny](int i, int j){
		return j*(Nx+1)+i;
	};
	std::vector<int> cls;
	for (int j=0; j<Ny; ++j){
		for (int i=0; i<Nx; ++i){
			cls.push_back(4);
			cls.push_back(pts_ind(i,j));
			cls.push_back(pts_ind(i+1,j));
			cls.push_back(pts_ind(i+1,j+1));
			cls.push_back(pts_ind(i,j+1));
		}
	}
	return grid_construct((Nx+1)*(Ny+1), Nx*Ny, &pts[0], &cls[0]);
}

void test1(){
	std::cout<<"1. grid creation and communication"<<std::endl;
	double points[] = {0,0, 1,0, 1,1, 0,1, 1,0.5};
	int cells[] = {3,0,1,4, 4,0,4,2,3};
	Grid* g = grid_construct(5,2,points,cells);
	add_check(grid_npoints(g)==5, "number of points");
	add_check(grid_ncells(g)==2, "number of cells");
	grid_save_vtk(g, "out1.vtk");
	grid_free(g);
}

void test2(){
	std::cout<<"2. merging two grids. Secondary grid lies within the main"<<std::endl;
	Grid* gmain = rectangular_grid(0,0, 1,1, 10, 10);
	Grid* gsec  = rectangular_grid(0.3,0.3, 0.6, 0.6, 30, 30);
	Grid* res = cross_grids(gmain, gsec, 0.05);
	grid_save_vtk(gmain,"out_main2.vtk");
	grid_save_vtk(gsec,"out_sec2.vtk");
	grid_save_vtk(res,"out_res2.vtk");
	grid_free(gmain);
	grid_free(gsec);
	grid_free(res);
}

void test3(){
	std::cout<<"3. merging two grids. Secondary grid crosses area of the main"<<std::endl;
	Grid* gmain = rectangular_grid(0,0, 1,1, 10, 10);
	Grid* gsec  = rectangular_grid(0.3,-0.1, 0.6, 0.2, 30, 30);
	Grid* res = cross_grids(gmain, gsec, 0.15);
	grid_save_vtk(gmain,"out_main3.vtk");
	grid_save_vtk(gsec,"out_sec3.vtk");
	grid_save_vtk(res,"out_res3.vtk");
	grid_free(gmain);
	grid_free(gsec);
	grid_free(res);
}

void test4(){
	std::cout<<"4. Secondary grid covers the corner of main"<<std::endl;
	Grid* gmain = rectangular_grid(0,0, 1,1, 10, 10);
	Grid* gsec  = rectangular_grid(-0.3,-0.3, 0.5, 0.5, 30, 30);
	Grid* res = cross_grids(gmain, gsec, 0.2);
	grid_save_vtk(gmain,"out_main4.vtk");
	grid_save_vtk(gsec,"out_sec4.vtk");
	grid_save_vtk(res,"out_res4.vtk");

	grid_free(gmain);
	grid_free(gsec);
	grid_free(res);
}

void test5(){
	std::cout<<"5. Diamond within a square grid"<<std::endl;
	Grid* gmain = rectangular_grid(-5,-5, 5,5, 20, 20);
	double pnt2[] = {
		0,0,0,-0.5,0.5,0,0,0.5,-0.5,0,0,-1,1,0,0,1,-1,0
	};
	int cls2[] ={
		3,1,0,4,  3,4,0,3,  3,0,2,3,  3,1,2,0,
		4,1,4,8,5,  4,8,4,3,7,  4,3,2,6,7,  4,1,5,6,2
	};
	Grid* gsec  = grid_construct(9, 8, pnt2, cls2);

	Grid* res = cross_grids(gmain, gsec, 2.0);
	grid_save_vtk(gmain,"out_main5.vtk");
	grid_save_vtk(gsec,"out_sec5.vtk");
	grid_save_vtk(res,"out_res5.vtk");

	grid_free(gmain);
	grid_free(gsec);
	grid_free(res);
}

int main(){
	crossgrid_internal_tests();
	test1();
	test2();
	test3();
	test4();
	test5();
	std::cout<<"DONE"<<std::endl;
}
