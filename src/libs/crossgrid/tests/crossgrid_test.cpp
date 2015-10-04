#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include "crossgrid.h"
#include "fileproc.h"
#include "hmcport.h"
#include "trigrid.h"

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

Cont* uniform_polygon(double xc, double yc, int N, double rad){
	//points
	std::vector<double> pts;
	for (int i=0; i<N; ++i){
		double ang = i*2*M_PI/N;
		pts.push_back(rad*cos(ang) + xc);
		pts.push_back(rad*sin(ang) + yc);
	}
	//edges
	std::vector<int> ed;
	for (int i=0; i<N; ++i){
		ed.push_back(i);
		ed.push_back(i+1);
	}
	ed.back() = 0;
	return contour_construct(N, N, &pts[0], &ed[0]);
}

Cont* unite_contours(Cont* c1, Cont* c2){
	std::vector<double> pfin;
	std::vector<int> efin;

	double* pts; int* eds;
	int Npnt; int Neds;
	//first
	contour_get_info(c1, &Npnt, &Neds, &pts, &eds);
	std::copy(eds, eds+2*(Neds), std::back_inserter(efin));
	std::copy(pts, pts+2*(Npnt), std::back_inserter(pfin));
	contour_free_info(&pts, &eds);
	//second
	int Npnt_old = Npnt;
	contour_get_info(c2, &Npnt, &Neds, &pts, &eds);
	std::for_each(eds, eds+2*Neds, [&Npnt_old](int& e){ e+=Npnt_old;});
	std::copy(eds, eds+2*Neds, std::back_inserter(efin));
	std::copy(pts, pts+2*Npnt, std::back_inserter(pfin));
	contour_free_info(&pts, &eds);
	//build a contour
	int Np = pfin.size()/2; int Ne = efin.size()/2;
	return contour_construct(Np, Ne, &pfin[0], &efin[0]);
}

std::pair<int, int> number_of_nonconv_cells(Grid* g){
	if (g==NULL) return std::make_pair(-1, -1);
	//-> non-convex cells, cells with hanging nodes
	auto res = std::make_pair(0, 0);
	auto& nonconv = res.first;
	auto& hangn = res.second;

	auto trisgn = [](double* p0, double* p1, double* p2)->double{
		double ax = p0[0]-p2[0], ay = p0[1]-p2[1];
		double bx = p1[0]-p2[0], by = p1[1]-p2[1];
		return ax*by-ay*bx;
	};

	double* pts; int *cls;
	int nc = grid_ncells(g);
	grid_get_points_cells2(g, &pts, &cls);
	int* c1 = cls;
	for (int i=0; i<nc; ++i){
		int n = *c1++;
		if (n<3) { throw std::runtime_error("cell with <3 nodes"); }
		if (n==3){ c1+=3; continue; }
		bool is_nonconv = false;
		bool is_hang = false;
		for (int j=0; j<n; ++j){
			auto pp = pts + 2*(c1[ j==0 ? n-1 : j-1 ]);
			auto pc = pts + 2*(c1[j]);
			auto pn = pts + 2*(c1[ j==n-1 ? 0 : j+1 ]);
			double a = trisgn(pp, pc, pn);
			if (fabs(a) < 1e-8) is_hang = true;
			else if (a<0) is_nonconv = true;
		}
		if (is_nonconv) ++nonconv;
		else if (is_hang) ++hangn;
		c1+=n;
	}
	
	delete[] pts;
	delete[] cls;
	return res;
}

bool check_convexity(Grid* g, int nconv_cells, int hang_cells){
	auto r = number_of_nonconv_cells(g);
	return (r.first == nconv_cells && r.second == hang_cells);
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
	Grid* res = cross_grids(gmain, gsec, 0.05, 1, 0);
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
	Grid* res = cross_grids(gmain, gsec, 0.15, 1, 0);
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
	Grid* res = cross_grids(gmain, gsec, 0.2, 1, 0);
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

	Grid* res = cross_grids(gmain, gsec, 2.0, 1, 0);
	grid_save_vtk(gmain,"out_main5.vtk");
	grid_save_vtk(gsec,"out_sec5.vtk");
	grid_save_vtk(res,"out_res5.vtk");

	grid_free(gmain);
	grid_free(gsec);
	grid_free(res);
}

void test6(){
	std::cout<<"6. Different density"<<std::endl;
	Grid* gmain = rectangular_grid(0,0, 1,1, 10, 10);
	Grid* gsec  = rectangular_grid(0.5,0.5, 0.6, 0.6, 30, 30);
	grid_save_vtk(gmain,"out_main6.vtk");
	grid_save_vtk(gsec,"out_sec6.vtk");

	Grid* res = cross_grids(gmain, gsec, 0.3, 1, 0);
	grid_save_vtk(res,"out_res6_1.vtk");
	grid_free(res);

	res = cross_grids(gmain, gsec, 0.3, 1, 0);
	grid_save_vtk(res,"out_res6_5.vtk");
	grid_free(res);

	res = cross_grids(gmain, gsec, 0.3, 1, 0);
	grid_save_vtk(res,"out_res6_9.vtk");
	grid_free(res);

	grid_free(gmain);
	grid_free(gsec);
}

void test7(){
	std::cout<<"7. Merging non crossing areas"<<std::endl;
	Grid* gmain = rectangular_grid(0,0, 1,1, 10, 10);
	Grid* gsec  = rectangular_grid(2,0, 3,1, 10, 10);
	Grid* gsec2  = rectangular_grid(1,1, 2,2, 10, 10);
	Grid* gsec3  = rectangular_grid(2,1.05, 3, 2.05, 10, 10);
	Grid* res = cross_grids(gmain, gsec, 0.2, 1, 0);
	Grid* res2 = cross_grids(res, gsec2, 0.2, 1, 0);
	Grid* res3 = cross_grids(res2, gsec3, 0.2, 1, 0);
	grid_save_vtk(res2,"out_res7.vtk");
	add_check(grid_npoints(res)==242 && grid_ncells(res)==200, "merge non crossing");
	add_check(grid_npoints(res2)==361 && grid_ncells(res2)==300, "merge grids with congruent point");
	add_check(grid_npoints(res3)>482 && grid_ncells(res3)>400, "merge grids with tangent edges");
	grid_free(gmain);
	grid_free(gsec);
	grid_free(gsec2);
	grid_free(gsec3);
	grid_free(res);
	grid_free(res2);
	grid_free(res3);
}

void test8(){
	std::cout<<"8. Merging areas with complicated intersections"<<std::endl;
	Grid* gmain = rectangular_grid(0,0, 1,1, 10, 10);
	Grid* gsec  = rectangular_grid(4,0, 5,1, 10, 10);
	Grid* gsec2  = rectangular_grid(-0.5,0.3, 5.5,0.6, 100, 10);
	Grid* res = cross_grids(gmain, gsec, 0.2, 1, 0);
	Grid* res2 = cross_grids(res, gsec2, 0.2, 1, 0);
	add_check(check_convexity(res2, 0, 4), "hanging nodes number");
	grid_save_vtk(res2,"out_res8.vtk");
	grid_free(gmain);
	grid_free(gsec);
	grid_free(gsec2);
	grid_free(res);
	grid_free(res2);
}

void test9(){
	std::cout<<"9. Boundary points control"<<std::endl;
	Grid* gmain = rectangular_grid(0,0, 7,7, 7, 7);
	Grid* gsec  = rectangular_grid(2.5,-1, 4.99, 1, 30, 30);
	Grid* res = cross_grids(gmain, gsec, 0.2, 1, 0);
	Grid* res2 = cross_grids(gmain, gsec, 0.2, 0, 0);
	std::vector<double> pts(grid_npoints(res)*2), pts2(grid_npoints(res2)*2);
	grid_get_points_cells(res, &pts[0], 0);
	grid_get_points_cells(res2, &pts2[0], 0);
	auto find_pt = [](double x, double y, const std::vector<double>& p)->bool{
		auto it = p.begin();
		while (it!=p.end()){
			bool t1 = (fabs(*it++ - x)<1e-6);
			bool t2 = (fabs(*it++ - y)<1e-6);
			if (t1 && t2) return true;
		}
		return false;
	};
	add_check(find_pt(5, 0, pts), "boundary point was set");
	add_check(!find_pt(5, 0, pts2), "boundary point was ignored");
	grid_save_vtk(res,"out_res9_1.vtk");
	grid_save_vtk(res2,"out_res9_2.vtk");
	grid_free(gmain);
	grid_free(gsec);
	grid_free(res);
	grid_free(res2);
}

void test10(){
	std::cout<<"10. Buffer zone is bigger then outer grid"<<std::endl;
	Grid* gmain = rectangular_grid(0,0, 1,1, 10, 10);
	Grid* gsec  = rectangular_grid(0.31,-0.3, 0.695, 0.5, 10, 10);

	Grid* res = cross_grids(gmain, gsec, 1.0, 1, 0);
	Grid* res2 = cross_grids(gmain, gsec, 1.0, 0, 0);

	grid_save_vtk(res,"out_res10_1.vtk");
	grid_save_vtk(res2,"out_res10_2.vtk");
	grid_free(gmain);
	grid_free(gsec);
	grid_free(res);
	grid_free(res2);
}

void test11(){
	std::cout<<"11. Grid combine with not single connected result"<<std::endl;
	Grid* gmain = rectangular_grid(0,0, 1,1, 10, 10);

	std::vector<double> pts_sec = {0.2, -0.2, 0.8, -0.2, 
		0.8, 0.2, 0.2, 0.2, 0.3, -0.1, 0.7, -0.1, 0.7, 0.1, 0.3, 0.1};
	std::vector<int> cls_sec = {4,0,1,5,4, 4,1,2,6,5, 4,6,2,3,7, 4,3,0,4,7};
	Grid* gsec = grid_construct(8, 4, &pts_sec[0], &cls_sec[0]);
	Grid* res = cross_grids(gmain, gsec, 0, 0, 0);
	grid_save_vtk(res, "out_res11.vtk");
	add_check(grid_npoints(res) == 125 && grid_ncells(res) == 96, "resulting topology");
	grid_free(gmain);
	grid_free(gsec);
	grid_free(res);
}

void test12(){
	//!!! Different results depending on i don't know what
	std::cout<<"12. Big buffer for the polygon with hole imposition"<<std::endl;
	Grid* gmain = rectangular_grid(0,0, 1,1, 10, 10);
	std::vector<double> pts_sec = {0.2, -0.2, 0.8, -0.2, 
		0.8, 0.2, 0.2, 0.2, 0.3, -0.1, 0.7, -0.1, 0.7, 0.1, 0.3, 0.1};
	for (size_t i=1; i<pts_sec.size(); i+=2) pts_sec[i]+=0.5;
	std::vector<int> cls_sec = {4,0,1,5,4, 4,1,2,6,5, 4,6,2,3,7, 4,3,0,4,7};
	Grid* gsec = grid_construct(8, 4, &pts_sec[0], &cls_sec[0]);
	Grid* res = cross_grids(gmain, gsec, 0.2, 1, 0);
	grid_save_vtk(res, "out_res12.vtk");
	//frontal
	//add_check(grid_ncells(res) == 132 && grid_npoints(res) == 99, "resulting topology");
	//Delaunay
	//add_check(grid_ncells(res) == 122 && grid_npoints(res) == 94, "resulting topology");

	grid_free(gmain);
	grid_free(gsec);
	grid_free(res);
}

void test13(){
	std::cout<<"13. Big data processing"<<std::endl;
	Grid* gmain = rectangular_grid(0,0, 1,1, 100, 100);
	Grid* gsec = rectangular_grid(-2,0, -1,1, 10, 10);
	Grid* res = cross_grids(gmain, gsec, 0.2, 1, 0);
	grid_save_vtk(res, "out_res13.vtk");

	grid_free(gmain);
	grid_free(gsec);
	grid_free(res);
}

void test14(){
	std::cout<<"14. Grid minus contour: basic"<<std::endl;
	auto g = rectangular_grid(0, 0, 10, 10, 12, 12);
	double pts[] = {
		2,2, 5,2, 6,7, 2,6, 3,3, 4,3, 4,4, 3,5
	};
	int edges[] = {
		0,1, 3,0, 2,3, 2,1, 7,6, 6,5, 4,5, 7,4
	};
	auto c = contour_construct(8, 8, pts, edges);
	auto res1 = grid_exclude_cont(g, c, true);
	auto res2 = grid_exclude_cont(g, c, false);
	add_check(grid_ncells(res1) == 138 && grid_npoints(res1) == 184, "outer resulting topology");
	add_check(grid_ncells(res2) == 35 && grid_npoints(res2) == 55, "inner resulting topology");
	
	grid_save_vtk(res2, "out_res14.vtk");
	grid_free(res1);
	grid_free(res2);
	grid_free(g);
	cont_free(c);
}

void test15(){
	std::cout<<"15. Grid minus contour: non-trivial topology"<<std::endl;
	auto g = rectangular_grid(0, 0, 1, 1, 10, 10);
	double points[] = {
		3,3, 4,3, 4,4, 3,4
	};
	int edges[] = {
		0,1, 1,2, 2,3, 3,0
	};
	auto c = contour_construct(4, 4, points, edges);

	auto res = grid_exclude_cont(g, c, true);
	add_check(grid_ncells(res) == 100 && grid_npoints(res) == 121, "fully inside");
	delete res;

	res = grid_exclude_cont(g, c, false);
	add_check(grid_ncells(res) == 0 && grid_npoints(res) == 0, "fully outside");
	delete res; 

	double points2[] = {
		0.03,0.03, 0.04,0.03, 0.04,0.04, 0.03,0.04
	};
	delete c; c = contour_construct(4, 4, points2, edges);
	res = grid_exclude_cont(g, c, false);
	add_check(grid_ncells(res) == 1 && grid_npoints(res) == 4, "inner contour within cell");
	delete res;

	res = grid_exclude_cont(g, c, true);
	add_check(grid_ncells(res) == 101 && grid_npoints(res) == 125, "outer contour within cell (intrusion)");
	grid_save_vtk(res, "out_res15.vtk");
	delete res;

	delete c;
	delete g;
}

void test16(){
	std::cout<<"16. Contours edges correlation"<<std::endl;
	double p1[] = {
		0,0, 4,2, 4,4, 0,6,
		1,2, 1,4, 3,4, 2,1
	};
	double p2[] = {
		1,2, 2,3, 3,4, 1,3.5, 1,4,
		0,6, 6,3, 0,0
	};
	int e1[] = {
		0,7, 2,1, 0,3, 4,5, 2,3, 6,5, 6,4, 1,7
	};
	int e2[] = {
		7,6, 5,6, 2,1, 2,4, 0,1, 0,3, 4,3, 5,7
	};
	auto c1 = contour_construct(8, 8, p1, e1);
	auto c2 = contour_construct(8, 8, p2, e2);
	int v1[] = {0,1,2,3,4,5,6,7};
	int v2[] = {0,0,0,0,0,0,0,0};
	int ans[] = {7,4,6,5,6,3,3,2};
	
	add_contour_bc(c1, c2, v1, v2, -10);
	
	bool good=true;
	for (int i=0; i<8; ++i) if (v2[i]!=ans[i]) good=false;
	add_check(good, "check vector values");

	cont_free(c1);
	cont_free(c2);
}

void test17(){
	std::cout<<"17. Large scale differences"<<std::endl;
	auto bigg=rectangular_grid(-20, -10, 100, 10, 12, 2);
	auto smallg=rectangular_grid(0, 0, 0.1, 0.1, 10, 10);
	Grid* unig = cross_grids(bigg, smallg, 0.2, 1, 0);

	grid_save_vtk(unig, "out_res17.vtk");
	grid_free(bigg);
	grid_free(smallg);
	grid_free(unig);

}

void test18(){
	std::cout<<"18. Excluding of multiple contours"<<std::endl;
	auto c1 = uniform_polygon(0, 0, 5, 1);
	auto c2 = uniform_polygon(3, 2.7, 12, 0.2);
	auto g = rectangular_grid(-2,-2,6,4,30,30);

	//1) excluding one by one
	auto g2 = grid_exclude_cont(g, c1, 1);
	auto g3 = grid_exclude_cont(g2, c2, 1);
	add_check(grid_ncells(g3) == 869 && grid_npoints(g3)==962, "one by one exclusion");

	//2) simultaneous exclusion
	auto c3=unite_contours(c1, c2);
	auto g4 = grid_exclude_cont(g, c3, 1);
	add_check(grid_ncells(g4) == 869 && grid_npoints(g4)==962, "simulataneous exclusion");
	grid_save_vtk(g4, "out_res18.vtk");

	grid_free(g);
	grid_free(g2);
	grid_free(g3);
	grid_free(g4);
	cont_free(c1);
	cont_free(c2);
	cont_free(c3);
}

void test19(){
	std::cout<<"19. No hanging nodes"<<std::endl;
	auto g1 = rectangular_grid(0,0,1,1,10,10);
	auto g2 = rectangular_grid(0.5,0.3,1.5,0.65,13,10);
	auto g3 = cross_grids(g1, g2, 0.1, 0, 0); 
	add_check(check_convexity(g3, 0, 0), "no hanging nodes check");
	grid_save_vtk(g3, "out_res19.vtk");
	grid_free(g1);
	grid_free(g2);
	grid_free(g3);
};

void test20(){
	std::cout<<"20. Empty holes"<<std::endl;
	auto g1 = rectangular_grid(0,0,1,1,6,6);
	auto c1 = uniform_polygon(0.5, 0.5, 10, 0.3);
	auto g2 = grid_exclude_cont(g1, c1, 1);
	auto g3 = rectangular_grid(0.45, 0.45, 0.55, 0.55, 3,3);
	auto g4 = cross_grids(g2, g3, 0, 0, 0);
	auto g5 = rectangular_grid(-5, -5, 5, 5, 19, 21);
	auto g6 = cross_grids(g5, g4, 0.1, 0, 1);

	double a = grid_area(g5)-contour_area(c1)+grid_area(g3);
	add_check( fabs(grid_area(g6) - a)<1e-6, "resulting grid with a hole");

	grid_save_vtk(g6, "out_res20.vtk");
	grid_free(g1);
	grid_free(g2);
	grid_free(g4);
	grid_free(g3);
	grid_free(g5);
	grid_free(g6);
	cont_free(c1);
};

void test21(){
	std::cout<<"21. Constrained triangulation"<<std::endl;
	auto no_edge_intersections = [](Point p0, Point p1, const GridGeom& g){
		auto edges = g.get_edges();
		double ksieta[2];
		for (auto& e: edges){
			auto gp1 = g.get_point(e.p1);
			auto gp2 = g.get_point(e.p2);
			if (SectCross(p0, p1, *gp1, *gp2, ksieta)){
				if (ksieta[1] > geps && ksieta[1] < 1-geps){
					return false;
				}
			}
		}
		return true;
	};
	// === 
	vector<Point> outer1 { Point(0,0), Point(1,0), Point(1,1) , Point(0.6, 1.0), Point (0, 1)};
	vector<Point> cons1 { Point(0.1, 0.1), Point(0.8, 0.8)};
	auto g1 = TriGrid::TriangulateAreaConstrained({outer1}, {cons1}, 0.05);
	add_check(no_edge_intersections(cons1[0], cons1[1], *g1),
			"One line constraint");
	// === 
	vector<Point> cons2 { Point(0.1, 0.1), Point(0.8, 0.8), Point (0.6, 1.0)};
	auto g2 = TriGrid::TriangulateAreaConstrained({outer1}, {cons2}, 0.05);
	add_check(no_edge_intersections(cons2[0], cons2[1], *g2) &&
		  no_edge_intersections(cons2[1], cons2[2], *g2),
			"Two lines constraint, common point");
	// === 
	vector<Point> outer3 { Point(0,0), Point(1,0), Point(1,1), Point (0, 1)};
	vector<Point> cons3 { Point(0.1, 0.1), Point(0.8, 0.8), Point(0.6, 1.0)};
	auto g3 = TriGrid::TriangulateAreaConstrained({outer3}, {cons3}, 0.05);
	add_check(no_edge_intersections(cons3[0], cons3[1], *g3) &&
		  no_edge_intersections(cons3[1], cons3[2], *g3),
			"Two lines constraint, no common point");
	// ===
	vector<Point> outer4 { Point(0.3, 0.3), Point(0.6, 0.3), Point(0.6, 0.6), Point(0.3, 0.6)};
	vector<Point> cons4 { Point(0.3, 0.6), Point(0.85, 1) };
	auto g4 = TriGrid::TriangulateAreaConstrained({outer3, outer4}, {cons4}, 0.05);
	add_check(no_edge_intersections(cons4[0], cons4[1], *g4),
			"Doubly connected contour");
	
};

int main(){
	crossgrid_silent_callback();
	crossgrid_internal_tests();
	test1();
	test2();
	test3();
	test4();
	test5();
	test6();
	test7();
	test8();
	test9();
	test10();
	test11();
	test12();
	test13();
	test14();
	test15();
	test16();
	test17();
	test18();
	test19();
	test20();
	test21();


	if (FAILED_CHECKS ==1){
		std::cout<<FAILED_CHECKS<<" test failed <<<<<<<<<<<<<<<<<<<"<<std::endl;
	} else if (FAILED_CHECKS > 1) {
		std::cout<<FAILED_CHECKS<<" tests failed <<<<<<<<<<<<<<<<<<<"<<std::endl;
	} else {
		std::cout<<"All tests passed"<<std::endl;
	}
	std::cout<<"DONE"<<std::endl;
}
