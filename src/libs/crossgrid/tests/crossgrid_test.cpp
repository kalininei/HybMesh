#include "hmtesting.hpp"
#include "cport_grid2d.h"
#include "cport_cont2d.h"
#include "trigrid.hpp"
#include "pebi.hpp"
#include "buildcont.hpp"
#include "buildgrid.hpp"
#include "healgrid.hpp"
#include "infogrid.hpp"
#include "unite_grids.hpp"
#include "modgrid.hpp"
#include "wireframegrid.hpp"
#include "assemble2d.hpp"
#include "debug2d.hpp"
#include "debug_grid2d.hpp"
#include "export2d_vtk.hpp"
#include "export2d_fluent.hpp"
#include "export2d_hm.hpp"
#include "export2d_tecplot.hpp"
#include "import2d_hm.hpp"
#include "snap_grid2cont.hpp"
#include "inscribe_grid.hpp"
#include "partcont.hpp"

using HMTesting::add_check;
using HMTesting::add_file_check;
using HMTesting::add_filesize_check;

int nocb(const char*, const char*, double, double){ return 0; }

double maxskew(const HM2D::GridData& g){
	vector<double> s=HM2D::Grid::Skewness(g);
	return *max_element(s.begin(), s.end());
}

std::pair<int, int> number_of_nonconv_cells(void* g){
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

	int dims[3];
	g2_dims(g, dims);

	int *cz=new int[dims[2]], *cls=0, clssize=0;
	double* pts=new double[dims[0]*2];
	g2_tab_cellsizes(g, cz);
	g2_tab_cellvert(g, &clssize, &cls);
	g2_tab_vertices(g, pts);
	int* c1=cls;
	for (int i=0; i<dims[2]; ++i){
		int n = cz[i];
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
	delete[] cz;
	delete[] cls;
	return res;
}

bool check_convexity(void* g, int nconv_cells, int hang_cells){
	auto r = number_of_nonconv_cells(g);
	return (r.first == nconv_cells && r.second == hang_cells);
}

bool has_edges(const HM2D::EdgeData& gd, const HM2D::VertexData& pts, const std::vector<int>& edges){
	vector<bool> ret(edges.size()/2, false);

	for (auto& e: gd){
		for (int i=0; i<ret.size(); ++i) if (!ret[i]){
			if ((*e->pfirst()==*pts[edges[2*i]] && *e->plast()==*pts[edges[2*i+1]]) ||
			    (*e->pfirst()==*pts[edges[2*i+1]] && *e->plast()==*pts[edges[2*i]])){
				ret[i] = true;
				break;
			}

		}
	}
	return std::all_of(ret.begin(), ret.end(), [](bool a){return a;});
}

bool edge_len_within(const HM2D::EdgeData& ed, double minlen, double maxlen){
	double l0 = ed[0]->length(), l1 = ed[0]->length();
	for (auto e: ed){
		double len = e->length();
		if (len < l0) l0 = len;
		if (len > l1) l1 = len;
	}
	//std::cout<<l0<<" "<<l1<<std::endl;
	if (l0<minlen || l1>maxlen) return false;
	else return true;
}

void* grid_construct(int npoints, int ncells, double* points, int* cells){
	int* cellsizes = new int[ncells];
	int* it=cells;
	int sum=0;
	for (int i=0; i<ncells; ++i){
		cellsizes[i] = *it;
		sum += cellsizes[i];
		it += cellsizes[i]+1;
	}
	int* cellvert = new int[sum];
	it = cells;
	int* it2 = cellvert;
	for (int i=0; i<ncells; ++i){
		int n = *it++;
		for (int j=0; j<n; ++j){
			*it2++ = *it++;
		}
	}
	void* ret;
	g2_from_points_cells(
		npoints, points, ncells, cellsizes, cellvert,
		0, 0, &ret);
	delete[] cellsizes;
	delete[] cellvert;
	return ret;
}

void test0(){
	using namespace HM2D;
	using namespace HM2D::Grid;
	using namespace HM2D::Grid::Impl;
	std::cout<<"0. wireframe grid tests"<<std::endl;
	{
		double pts[] = {0,0,9,2,8,6,4,1};
		int cls[] = {3,0,3,2, 3,3,1,2};
		GridData G = Constructor::FromRaw(4, 2, pts, cls, -1);
		PtsGraph gr(G);

		vector<double> vv {2, -1, 9, -1, 4, 5, 2, 5};
		EdgeData ed = Contour::Constructor::FromPoints(vv, true);
		Contour::Tree tree = Contour::Tree::Assemble(ed);
		PtsGraph newgraph = PtsGraph::cut(gr, tree, OUTSIDE);
		GridData newgeom =  newgraph.togrid();
		add_check(newgeom.vcells.size()==4 && newgeom.vvert.size()==10, "PtsGraph cut");
	}
	{
		double pts1[] = {0,0, 2,2, 0,4, -2,2, 0,2};
		double pts2[] = {-1,1, 1,1, 1, 2.5, -1,3};
		int cls1[] = {3,0,1,4, 3,4,1,2, 3,3,4,2, 3,3,0,4};
		int cls2[] = {4,0,1,2,3};
		GridData gmain = Constructor::FromRaw(5,4, pts1, cls1, -1);
		GridData gsec = Constructor::FromRaw(4,1, pts2, cls2, -1);
		GridData comb = Algos::CombineGrids(gmain, gsec);
		add_check(comb.vcells.size()==8 && comb.vvert.size()==12, "combine grids");
	}
	{
		GridData grid1 = Constructor::RectGrid01(10, 10);
		GridData grid2 = Constructor::RectGrid(Point(1, 1), Point(2, 2), 10, 10);
		GridData grid3 = Constructor::RectGrid(Point(2, 1.85), Point(3, 2.85), 10, 10);
		Algos::OptUnite opt(0.0, true);
		auto cross1 = Algos::UniteGrids(grid1, grid2, opt);
		auto cross2 = Algos::UniteGrids(cross1, grid3, opt);
		add_check(cross2.vvert.size()==362 && cross2.vcells.size()==300, "grid subdivide1");
		Export::GridVTK(cross2, "g2.vtk");
		auto div = SplitData(cross2.vcells);
		add_check(div.size()==2, "grid subdive2");
		auto g1 = Constructor::InvokeGrid::Permanent(div[0]);
		auto g2 = Constructor::InvokeGrid::Permanent(div[1]);
		add_check(g1.vvert.size()==121 && g1.vcells.size()==100, "grid subdivide3");
		add_check(g2.vvert.size()==242 && g2.vcells.size()==200, "grid subdivide4");
	}
	//PtsGraph::togrid with intrusion
	{ //case1
		auto c1 = Contour::Constructor::Circle(6, 4, Point(0, 0));
		auto c2 = Contour::Constructor::Circle(4, 2, Point(0, 0));
		PtsGraph g1(c1);
		g1.add_edges(c2);
		auto grid = g1.togrid();
		add_check(fabs(Contour::Area(c1)-Area(grid)) < 1e-6, "grid from two polys");
		HM2D::Export::GridVTK(grid, "g1.vtk");
	}
	{ //case 2
		auto c1 = Contour::Constructor::Circle(7, 5, Point(3, 2));
		auto c2 = Contour::Constructor::Circle(4, 2, Point(3, 3));
		auto c3 = Contour::Constructor::Circle(3, 0.5, Point(3.2, 3.1));
		PtsGraph g1(c1); g1.add_edges(c2); g1.add_edges(c3);
		auto grid = g1.togrid();
		add_check(fabs(Contour::Area(c1)-Area(grid)) < 1e-6, "grid from 3 nested polys");
	}
	{ //case 3
		auto c0 = Contour::Constructor::Circle(5, 0.5, Point(3, 3));
		auto c1 = Contour::Constructor::Circle(5, 1, Point(6, 6));
		auto c2 = Contour::Constructor::Circle(8, 2, Point(2, 3));
		auto c3 = Contour::Constructor::Circle(4, 0.1, Point(2.1, 3));
		auto c4 = Contour::Constructor::Circle(5, 0.3, Point(2, 3));
		auto c5 = Contour::Constructor::Circle(3, 0.1, Point(1, 4));
		PtsGraph g1(c0);
		g1.add_edges(c1); g1.add_edges(c2);
		g1.add_edges(c3); g1.add_edges(c4); g1.add_edges(c5);
		auto grid = g1.togrid();
		double a = Contour::Area(c1) + Contour::Area(c2);
		add_check(fabs(a-Area(grid)) < 1e-6 &&
			ECol::Assembler::GridBoundary(grid).size() == 13,
			"grid from complicated nested structure");
		Export::GridVTK(grid, "g1.vtk");
	}
}
void test1(){
	std::cout<<"1. grid creation and communication"<<std::endl;
	double points[] = {0,0, 1,0, 1,1, 0,1, 1,0.5};
	int cells[] = {3,0,1,4, 4,0,4,2,3};
	void* g = grid_construct(5,2,points,cells);
	int dims[3];
	g2_dims(g, dims);
	add_check(dims[0]==5, "number of points");
	add_check(dims[2]==2, "number of cells");
	g2_free(g);
}

void test2(){
	std::cout<<"2. merging two grids. Secondary grid lies within the main"<<std::endl;
	HM2D::GridData gmain = HM2D::Grid::Constructor::RectGrid01(10, 10);
	HM2D::GridData gsec = HM2D::Grid::Constructor::RectGrid(Point(0.3, 0.3), Point(0.6, 0.6), 30, 30);
	void* res;
	g2_unite_grids(&gmain, &gsec, 0.05, 1, 0, 0, "3", &res, nocb);
	g2_free(res);
}

void test3(){
	std::cout<<"3. merging two grids. Secondary grid crosses area of the main"<<std::endl;
	HM2D::GridData gmain = HM2D::Grid::Constructor::RectGrid(Point(0, 0), Point(1, 1), 10, 10);
	HM2D::GridData gsec = HM2D::Grid::Constructor::RectGrid(Point(0.3, -0.1), Point(0.6, 0.2), 30, 30);
	void* res;
	g2_unite_grids(&gmain, &gsec, 0.15, 1, 0, 0, "3", &res, nocb);
	g2_free(res);
}

void test4(){
	std::cout<<"4. Secondary grid covers the corner of main"<<std::endl;
	HM2D::GridData gmain = HM2D::Grid::Constructor::RectGrid(Point(0, 0), Point(1, 1), 10, 10);
	HM2D::GridData gsec = HM2D::Grid::Constructor::RectGrid(Point(-0.3, -0.3), Point(0.5, 0.5), 30, 30);
	void* res;
	g2_unite_grids(&gmain, &gsec, 0.2, 1, 0, 0, "3", &res, nocb);
	g2_free(res);
}

void test5(){
	std::cout<<"5. Diamond within a square grid"<<std::endl;
	HM2D::GridData gmain = HM2D::Grid::Constructor::RectGrid(Point(-5, -5), Point(5, 5), 20, 20);
	double pnt2[] = {
		0,0,0,-0.5,0.5,0,0,0.5,-0.5,0,0,-1,1,0,0,1,-1,0
	};
	int cls2[] ={
		3,1,0,4,  3,4,0,3,  3,0,2,3,  3,1,2,0,
		4,1,4,8,5,  4,8,4,3,7,  4,3,2,6,7,  4,1,5,6,2
	};
	void* gsec  = grid_construct(9, 8, pnt2, cls2);

	void* res;
	g2_unite_grids(&gmain, gsec, 2.0, 1, 0, 0, "3", &res, nocb);
	add_check(ISEQ(HM2D::Grid::Area(*static_cast<HM2D::GridData*>(res)), 100), "area check");

	g2_free(gsec);
	g2_free(res);
}

void test6(){
	std::cout<<"6. Different density"<<std::endl;
	HM2D::GridData gmain = HM2D::Grid::Constructor::RectGrid(Point(0, 0), Point(1, 1), 10, 10);
	HM2D::GridData gsec = HM2D::Grid::Constructor::RectGrid(Point(0.5, 0.5), Point(0.6, 0.6), 30, 30);

	void* res;
	g2_unite_grids(&gmain, &gsec, 0.3, 1, 0, 0, "3", &res, nocb);
	g2_free(res);

	g2_unite_grids(&gmain, &gsec, 0.3, 1, 0, 0, "3", &res, nocb);
	g2_free(res);

	g2_unite_grids(&gmain, &gsec, 0.3, 1, 0, 0, "3", &res, nocb);
	g2_free(res);

}

void test7(){
	std::cout<<"7. Merging non crossing areas"<<std::endl;
	HM2D::GridData gmain = HM2D::Grid::Constructor::RectGrid(Point(0, 0), Point(1, 1), 10, 10);
	HM2D::GridData gsec = HM2D::Grid::Constructor::RectGrid(Point(2, 0), Point(3, 1), 10, 10);
	HM2D::GridData gsec2 = HM2D::Grid::Constructor::RectGrid(Point(1, 1), Point(2, 2), 10, 10);
	HM2D::GridData gsec3 = HM2D::Grid::Constructor::RectGrid(Point(2, 1.05), Point(3, 2.05), 10, 10);

	void *res, *res2, *res3;
	g2_unite_grids(&gmain, &gsec, 0.2, 1, 0, 0, "3", &res, nocb);
	g2_unite_grids(res, &gsec2, 0.2, 1, 0, 0, "3", &res2, nocb);
	g2_unite_grids(res2, &gsec3, 0.2, 1, 0, 0, "3", &res3, nocb);

	int dim[3], dim2[3], dim3[3];
	g2_dims(res, dim);
	g2_dims(res2, dim2);
	g2_dims(res3, dim3);
	add_check(dim[0]==242 && dim[2]==200, "merge non crossing");
	add_check(dim2[0]==361 && dim2[2]==300, "merge grids with congruent point");
	add_check(dim3[0]>482 && dim3[2]>400, "merge grids with tangent edges");
	HM2D::Export::GridVTK(*static_cast<HM2D::GridData*>(res3), "g3.vtk");

	g2_free(res);
	g2_free(res2);
	g2_free(res3);
}

void test8(){
	std::cout<<"8. Merging areas with complicated intersections"<<std::endl;
	HM2D::GridData gmain = HM2D::Grid::Constructor::RectGrid(Point(0, 0), Point(1, 1), 10, 10);
	HM2D::GridData gsec = HM2D::Grid::Constructor::RectGrid(Point(4, 0), Point(5, 1), 10, 10);
	HM2D::GridData gsec2 = HM2D::Grid::Constructor::RectGrid(Point(-0.5, 0.3), Point(5.5, 0.6), 100, 10);

	void *res, *res2;
	g2_unite_grids(&gmain, &gsec, 0.2, 1, 0, 0, "3", &res, nocb);
	g2_unite_grids(res, &gsec2, 0.2, 1, 0, 0, "3", &res2, nocb);
	add_check(check_convexity(res2, 0, 4), "hanging nodes number");
	g2_free(res);
	g2_free(res2);
}

void test9(){
	std::cout<<"9. Boundary points control"<<std::endl;
	HM2D::GridData gmain = HM2D::Grid::Constructor::RectGrid( Point(0, 0), Point(7, 7), 7, 7);
	HM2D::GridData gsec = HM2D::Grid::Constructor::RectGrid( Point(2.5, -1), Point(4.99, 1), 30, 30);
	void* res, *res2;
	g2_unite_grids(&gmain, &gsec, 0.2, 1, 0, 0, "3", &res, nocb);
	g2_unite_grids(&gmain, &gsec, 0.2, 0, 0, 0, "3", &res2, nocb);
	int dim[3], dim2[3];
	g2_dims(res, dim);
	g2_dims(res2, dim2);
	std::vector<double> pts(dim[0]*2), pts2(dim2[0]*2);

	g2_tab_vertices(res, &pts[0]);
	g2_tab_vertices(res2, &pts2[0]);

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
	g2_free(res);
	g2_free(res2);
}

void test10(){
	std::cout<<"10. Buffer zone is bigger then outer grid"<<std::endl;
	HM2D::GridData gmain = HM2D::Grid::Constructor::RectGrid(Point(0,0), Point(1,1), 10, 10);
	HM2D::GridData gsec = HM2D::Grid::Constructor::RectGrid(Point(0.31,-0.3), Point(0.695,0.5), 10, 10);

	void* res, *res2;
	g2_unite_grids(&gmain, &gsec, 1.0, 1, 0, 0, "3", &res, nocb);
	g2_unite_grids(&gmain, &gsec, 1.0, 0, 0, 0, "3", &res2, nocb);

	g2_free(res);
	g2_free(res2);
}

void test11(){
	std::cout<<"11. Grid combine with not single connected result"<<std::endl;
	HM2D::GridData gmain = HM2D::Grid::Constructor::RectGrid(Point(0,0), Point(1,1), 10, 10);
	
	std::vector<double> pts_sec = {0.2, -0.2, 0.8, -0.2, 
		0.8, 0.2, 0.2, 0.2, 0.3, -0.1, 0.7, -0.1, 0.7, 0.1, 0.3, 0.1};
	std::vector<int> cls_sec = {4,0,1,5,4, 4,1,2,6,5, 4,6,2,3,7, 4,3,0,4,7};
	void *gsec, *res;
	gsec = grid_construct(8, 4, &pts_sec[0], &cls_sec[0]);
	g2_unite_grids(&gmain, gsec, 0, 0, 0, 0, "3", &res, nocb);
	int dim[3];
	g2_dims(res, dim);
	add_check(dim[0] == 125 && dim[2] == 96, "resulting topology");
	g2_free(gsec);
	g2_free(res);
}

void test12(){
	//!!! Different results depending on i don't know what
	std::cout<<"12. Big buffer for the polygon with hole imposition"<<std::endl;
	HM2D::GridData gmain = HM2D::Grid::Constructor::RectGrid(Point(0,0), Point(1,1), 10, 10);
	std::vector<double> pts_sec = {0.2, -0.2, 0.8, -0.2, 
		0.8, 0.2, 0.2, 0.2, 0.3, -0.1, 0.7, -0.1, 0.7, 0.1, 0.3, 0.1};
	for (size_t i=1; i<pts_sec.size(); i+=2) pts_sec[i]+=0.5;
	std::vector<int> cls_sec = {4,0,1,5,4, 4,1,2,6,5, 4,6,2,3,7, 4,3,0,4,7};
	void *gsec, *res;
	gsec = grid_construct(8, 4, &pts_sec[0], &cls_sec[0]);
	g2_unite_grids(&gmain, gsec, 0.2, 1, 0, 0, "3", &res, nocb);
	//frontal
	//add_check(grid_ncells(res) == 132 && grid_npoints(res) == 99, "resulting topology");
	//Delaunay
	//add_check(grid_ncells(res) == 122 && grid_npoints(res) == 94, "resulting topology");

	add_check(maxskew(*static_cast<HM2D::GridData*>(res))<0.8, "skewness check");
	g2_free(gsec);
	g2_free(res);
}

void test13(){
	std::cout<<"13. Big data processing"<<std::endl;
	HM2D::GridData gmain = HM2D::Grid::Constructor::RectGrid(Point(0,0), Point(1,1), 100, 100);
	HM2D::GridData gsec = HM2D::Grid::Constructor::RectGrid(Point(-2, 0), Point(-1, 1), 10, 10);
	void *res;
	g2_unite_grids(&gmain, &gsec, 0.2, 1, 0, 0, "3", &res, nocb);
	g2_free(res);
}

void test14(){
	std::cout<<"14. Grid minus contour: basic"<<std::endl;
	{
		HM2D::GridData g = HM2D::Grid::Constructor::RectGrid(Point(0,0), Point(10,10), 12, 12);
		double pts[] = {
			2,2, 5,2, 6,7, 2,6, 3,3, 4,3, 4,4, 3,5
		};
		int edges[] = {
			0,1, 3,0, 2,3, 2,1, 7,6, 6,5, 4,5, 7,4
		};
		auto c = HM2D::ECol::Constructor::FromRaw(8, 8, pts, edges);

		void *res1, *res2;
		g2_exclude_cont(&g, &c, true, &res1, nocb);
		g2_exclude_cont(&g, &c, false, &res2, nocb);
		int dim1[3], dim2[3];
		g2_dims(res1, dim1);
		g2_dims(res2, dim2);
		add_check(dim1[2] == 138 && dim1[0] == 184, "outer resulting topology");
		add_check(dim2[2] == 35 && dim2[0] == 55, "inner resulting topology");
		
		g2_free(res1);
		g2_free(res2);
	}
	{
		HM2D::GridData g = HM2D::Grid::Constructor::RectGrid01(10, 12);
		auto c = HM2D::Contour::Constructor::FromPoints({0, 0.5, 0.5, 0, 1.1, 0.5, 0.5, 1.1}, true);
		void *res3;
		g2_exclude_cont(&g, &c, true, &res3, nocb);
		HM2D::Export::GridVTK(*static_cast<HM2D::GridData*>(res3), "res3.vtk");
		int dim[3];
		g2_dims(res3, dim);
		add_check(dim[2] == 67 && dim[0] == 102, "cut contour touches grid contour");
		g2_free(res3);
	}
}

void test15(){
	std::cout<<"15. Grid minus contour: non-trivial topology"<<std::endl;
	HM2D::GridData g = HM2D::Grid::Constructor::RectGrid(Point(0, 0), Point(1, 1), 10, 10);
	double points[] = { 3,3, 4,3, 4,4, 3,4 };
	int edges[] = { 0,1, 1,2, 2,3, 3,0 };
	auto c = HM2D::ECol::Constructor::FromRaw(4, 4, points, edges);

	void* res;
	g2_exclude_cont(&g, &c, true, &res, nocb);
	int dim[3];
	g2_dims(res, dim);
	add_check(dim[2] == 100 && dim[0] == 121, "fully inside");
	g2_free(res);

	g2_exclude_cont(&g, &c, false, &res, nocb);
	g2_dims(res, dim);
	add_check(dim[2] == 0 && dim[0] == 0, "fully outside");
	g2_free(res); 

	double points2[] = {
		0.03,0.03, 0.04,0.03, 0.04,0.04, 0.03,0.04
	};
	c = HM2D::ECol::Constructor::FromRaw(4, 4, points2, edges);
	g2_exclude_cont(&g, &c, false, &res, nocb);
	g2_dims(res, dim);
	add_check(dim[2] == 1 && dim[0] == 4, "inner contour within cell");
	g2_free(res);

	g2_exclude_cont(&g, &c, true, &res, nocb);
	double arres, arg, arecol;
	g2_area(res, &arres);
	g2_area(&g, &arg);
	c2_area(&c, &arecol);
	add_check(ISZERO(arres - arg + arecol), "outer contour within cell (intrusion)");
	HM2D::Export::GridVTK(*static_cast<HM2D::GridData*>(res), "g1.vtk");
	g2_free(res);
}

void test16(){
	//std::cout<<"16. Contours edges correlation"<<std::endl;
	//double p1[] = {
	//        0,0, 4,2, 4,4, 0,6,
	//        1,2, 1,4, 3,4, 2,1
	//};
	//double p2[] = {
	//        1,2, 2,3, 3,4, 1,3.5, 1,4,
	//        0,6, 6,3, 0,0
	//};
	//int e1[] = {
	//        0,7, 2,1, 0,3, 4,5, 2,3, 6,5, 6,4, 1,7
	//};
	//int e2[] = {
	//        7,6, 5,6, 2,1, 2,4, 0,1, 0,3, 4,3, 5,7
	//};
	//auto c1 = HM2D::ECol::Constructor::FromRaw(8, 8, p1, e1);
	//auto c2 = HM2D::ECol::Constructor::FromRaw(8, 8, p2, e2);
	//int v1[] = {0,1,2,3,4,5,6,7};
	//int v2[] = {0,0,0,0,0,0,0,0};
	//int ans[] = {7,4,6,5,6,3,3,2};
	
	//set_ecollection_bc_force(&c1, &c2, v1, v2, 3);
	
	//bool good=true;
	//for (int i=0; i<8; ++i) if (v2[i]!=ans[i]) good=false;
	//add_check(good, "check vector values");
}

void test17(){
	std::cout<<"17. Large scale differences"<<std::endl;
	HM2D::GridData bigg = HM2D::Grid::Constructor::RectGrid(Point(-20, -10), Point(100, 10), 12, 2);
	HM2D::GridData smallg = HM2D::Grid::Constructor::RectGrid(Point(0, 0), Point(0.1, 0.1), 10, 10);
	void* unig;
	g2_unite_grids(&bigg, &smallg, 0.2, 1, 0, 0, "3", &unig, nocb);

	g2_free(unig);
}

void test18(){
	std::cout<<"18. Excluding of multiple contours"<<std::endl;
	auto c1 = HM2D::Contour::Constructor::Circle(5, 1.0, Point(0, 0)); 
	auto c2 = HM2D::Contour::Constructor::Circle(12, 0.2, Point(3, 2.7)); 
	auto g = HM2D::Grid::Constructor::RectGrid(Point(-2, -2), Point(6, 4), 30, 30);

	//1) excluding one by one
	void *g2, *g3;
	g2_exclude_cont(&g, &c1, 1, &g2, nocb);
	g2_exclude_cont(g2, &c2, 1, &g3, nocb);
	int dim[3];
	g2_dims(g3, dim);
	add_check(dim[2] == 869 && dim[0]==962, "one by one exclusion");

	//2) simultaneous exclusion
	c2.insert(c2.end(), c1.begin(), c1.end());
	void *g4;
	g2_exclude_cont(&g, &c2, 1, &g4, nocb);
	g2_dims(g4, dim);
	add_check(dim[2] == 869 && dim[0]==962, "simulataneous exclusion");

	g2_free(g2);
	g2_free(g3);
	g2_free(g4);
}

void test19(){
	std::cout<<"19. No hanging nodes"<<std::endl;
	auto g1 = HM2D::Grid::Constructor::RectGrid01(10, 10);
	auto g2 = HM2D::Grid::Constructor::RectGrid(Point(0.5, 0.3), Point(1.5, 0.65), 13, 10);
	void *g3;
	g2_unite_grids(&g1, &g2, 0.1, 0, 0, 0, "3", &g3, nocb);
	add_check(check_convexity(g3, 0, 0), "no hanging nodes check");
	g2_free(g3);
};

void test20(){
	std::cout<<"20. Empty holes"<<std::endl;
	auto g1 = HM2D::Grid::Constructor::RectGrid01(6, 6);
	auto c1 = HM2D::Contour::Constructor::Circle(10, 0.3, Point(0.5, 0.5));
	void *g2;
	g2_exclude_cont(&g1, &c1, 1, &g2, nocb);

	auto g3 = HM2D::Grid::Constructor::RectGrid(Point(0.45, 0.45), Point(0.55, 0.55), 3, 3);
	void* g4;
	g2_unite_grids(g2, &g3, 0, 0, 0, 0, "3", &g4, nocb);

	auto g5 = HM2D::Grid::Constructor::RectGrid(Point(-5, -5), Point(5, 5), 19, 21);
	void* g6;
	g2_unite_grids(&g5, g4, 0.1, 0, 1, 0, "3", &g6, nocb);

	double a5, a1, a3, a6;
	g2_area(&g5, &a5);
	g2_area(&g3, &a3);
	c2_area(&c1, &a1);
	g2_area(g6, &a6);
	double a = a5-a1+a3;
	add_check( fabs(a6 - a) < 1e-6, "resulting grid with a hole");

	g2_free(g2);
	g2_free(g4);
	g2_free(g6);
};

namespace{
HM2D::Contour::Tree tree_w_constraints(const vector<vector<Point>>& outer,
		const vector<vector<Point>>& constr, double defsize){
	//outer
	HM2D::Contour::Tree ret;
	for (auto& c: outer){
		ret.add_contour(HM2D::Contour::Constructor::FromPoints(c, true));
	}
	for (auto& c: constr){
		ret.add_detached_contour(HM2D::Contour::Constructor::FromPoints(c, false));
	}
	ret = HM2D::Mesher::PrepareSource(ret, defsize);
	return ret;
}
};

void test21(){
	std::cout<<"21. Constrained triangulation"<<std::endl;
	auto no_edge_intersections = [](Point p0, Point p1, const HM2D::GridData& g){
		double ksieta[2];
		for (auto& e: g.vedges){
			auto gp1 = e->first();
			auto gp2 = e->last();
			if (SectCross(p0, p1, *gp1, *gp2, ksieta)){
				if (ksieta[1] > geps && ksieta[1] < 1-geps){
					return false;
				}
			}
		}
		return true;
	};
	auto no_edges_intersections = [&](vector<Point> c, const HM2D::GridData& g){
		for (int i=0; i<c.size()-1; ++i){
			Point p1 = c[i];
			Point p2 = c[i+1];
			if (!no_edge_intersections(p1, p2, g)) return false;
		}
		return true;
	};
	auto has_point = [](Point p, const HM2D::GridData& g){
		for(auto& v: g.vvert){ if (*v == p) return true; }
		return false;
	};

	// === 
	vector<Point> outer1 { Point(0,0), Point(1,0), Point(1,1) , Point(0.6, 1.0), Point (0, 1)};
	vector<Point> cons1 { Point(0.1, 0.1), Point(0.8, 0.8)};
	auto t1 = tree_w_constraints({outer1}, {cons1}, 0.05);
	auto g1 = HM2D::Mesher::UnstructuredTriangle(t1);
	HM2D::Export::GridVTK(g1, "g1.vtk");
	add_check(no_edge_intersections(cons1[0], cons1[1], g1),
			"One line constraint");
	// === 
	vector<Point> cons2 { Point(0.1, 0.1), Point(0.8, 0.8), Point (0.6, 1.0)};
	auto t2 = tree_w_constraints({outer1}, {cons2}, 0.05);
	auto g2 = HM2D::Mesher::UnstructuredTriangle(t2);
	add_check(no_edge_intersections(cons2[0], cons2[1], g2) &&
		  no_edge_intersections(cons2[1], cons2[2], g2),
			"Two lines constraint, common point");
	// === 
	vector<Point> outer3 { Point(0,0), Point(1,0), Point(1,1), Point (0, 1)};
	vector<Point> cons3 { Point(0.1, 0.1), Point(0.8, 0.8), Point(0.6, 1.0)};
	auto t3 = tree_w_constraints({outer3}, {cons3}, 0.05);
	auto g3 = HM2D::Mesher::UnstructuredTriangle(t3);
	add_check(no_edge_intersections(cons3[0], cons3[1], g3) &&
		  no_edge_intersections(cons3[1], cons3[2], g3),
		  "Two lines constraint, no common point");
	// ===
	vector<Point> outer4 { Point(0.3, 0.3), Point(0.6, 0.3), Point(0.6, 0.6), Point(0.3, 0.6)};
	vector<Point> cons4 { Point(0.3, 0.6), Point(0.85, 1) };
	auto t4 = tree_w_constraints({outer3, outer4}, {cons4}, 0.05);
	auto g4 = HM2D::Mesher::UnstructuredTriangle(t4);
	add_check(no_edge_intersections(cons4[0], cons4[1], g4),
			"Doubly connected contour");
	// ===
	vector<Point> outer5 {Point(0, 0), Point(1,0), Point(1,1), Point(0,1) };
	vector<Point> outer6 {Point(0.2, 0.2), Point(0.6,0.2), Point(0.6,0.6), Point(0.2,0.6) };
	vector<Point> outer7 {Point(1.1, 0.1), Point(2,0.1), Point(2,1.1), Point(1.1,1.1) };
	vector<Point> cons8;
	for (int i=0; i<17; ++i){
		Point c(0.63*sin(double(i)/16*2*M_PI), 0.63*cos(double(i)/16*2*M_PI));
		cons8.push_back(c+Point(1, 0.5));
	}
	auto t5 = tree_w_constraints({outer5, outer6, outer7}, {cons8}, 0.05);
	auto g5 = HM2D::Mesher::UnstructuredTriangle(t5);
	add_check(no_edges_intersections(cons8, g5),
			"Constrant is not fully inside meshing zone");
	// ===
	auto t6 = tree_w_constraints({outer5}, {cons8}, 0.03);
	auto g6 = HM2D::Mesher::UnstructuredTriangleRecomb(t6);
	add_check(no_edges_intersections(cons8, g6),
			"Recombined mesh with polyline constraint 1");

	// === 
	auto t7 = tree_w_constraints({outer5, outer6, outer7}, {cons8}, 0.03);
	auto g7 = HM2D::Mesher::UnstructuredTriangleRecomb(t7);
	add_check(no_edges_intersections(cons8, g7),
			"Recombined mesh with polyline constraint 2");
	// ===
	auto t8 = tree_w_constraints({outer5}, {}, 0.1);
	CoordinateMap2D<double> emb;
	emb.add(0.3, 0.3, 0.01);
	emb.add(0.7, 0.7, 0.2);
	emb.add(1.3, 1.3, 0.01);
	auto g8 = HM2D::Mesher::UnstructuredTriangle(t8, emb);
	add_check(has_point(Point(0.3, 0.3), g8) && has_point(Point(0.7, 0.7), g8),
			"Embedded points");
};

void test22(){
	std::cout<<"22. Snap and shift boundaries"<<std::endl;
	auto cont = HM2D::Contour::Constructor::FromPoints(
		{0,0, 6,0, 9,1, 9,5, 4,5, 2,4, 0,2}, true);
	double pts[] = {0,0, 9,2.5, 8,5, 3.8,4.9, 2.2,4.1, 1.3,3.3};
	int cls[] = {0,1,2,3,4,5};

	auto ans1 = HM2D::Grid::Constructor::FromRaw(6, 1, pts, cls, 6);
	HM2D::Grid::Algos::SnapToContour(ans1, cont, {});
	add_check(fabs(HM2D::Contour::Area(cont) - HM2D::Grid::Area(ans1))<1e-12,
		"snapping of single cell");

	auto ans2 = HM2D::Grid::Constructor::FromRaw(6, 1, pts, cls, 6);
	HM2D::Export::GridVTK(ans2, "a22.vtk");
	HM2D::Grid::Algos::ShiftToContour(ans2, cont, {});
	add_check(fabs(35.5 - HM2D::Grid::Area(ans2))<1e-12,
		"shifting vertices of single cell");

	HM2D::Export::GridVTK(ans2, "t22.vtk");
	HM2D::Export::ContourVTK(cont, "c22.vtk");
}

void test23(){
	std::cout<<"23. Fluent exporting"<<std::endl;
	auto bfun1 = [](double x, double y)->int{
		if (x>1 || y>1) return 1;
		else if (ISEQ(x, 1) || ISEQ(y, 1) || ISZERO(x) || ISZERO(y)) return 2;
		else return 3;
	};
	auto build_bcond = [](HM2D::GridData& g, std::function<int(double, double)> fun){
		for (const auto& e: g.vedges){
			int val = 0;
			if (e->is_boundary()){
				auto pc = e->center();
				val = fun(pc.x, pc.y);
			}
			e->boundary_type = val;
		}
	};
	auto bfun2 = [](double x, double y)->int{
		if (ISZERO(x) && y<=0.5) return 1;
		if (ISZERO(x-1) && y>=0.5) return 2;
		if (ISZERO(y-1)) return 3;
		if (ISZERO(y)) return 4;
		return 0;
	};
	auto g1 = HM2D::Grid::Constructor::RectGrid01(20, 20);
	auto g2 = HM2D::Grid::Constructor::Ring(Point{0.8,0.8}, 0.3, 0.1, 20, 10);
	{
		HM2D::Grid::Algos::OptUnite uopt;
		uopt.empty_holes = true; uopt.angle0 = 10;
		HM2D::GridData g3 = HM2D::Grid::Algos::UniteGrids(g1, g2, uopt);
		build_bcond(g3, bfun1);

		HM2D::Export::GridVTK(g3, "g1.vtk");
		HM2D::Export::BoundaryVTK(g3, "c1.vtk");
		add_file_check(12480458951093668814U, "g1.vtk", "grid to vtk");
		add_file_check(9640670790211601721U, "c1.vtk", "grid contour to vtk");

		bool hasfailed = false;
		try{
			HM2D::Export::GridMSH(g3, "g1.msh");
		} catch (...){
			hasfailed = true;
		}
		add_check(hasfailed, "improper grid for fluent export");
	}
	{
		HM2D::Grid::Algos::OptUnite uopt(0.1);
		uopt.empty_holes = true; 
		uopt.preserve_bp=false; uopt.angle0=30;
		auto g3 = HM2D::Grid::Algos::UniteGrids(g1, g2, uopt);
		build_bcond(g3, bfun1);
		HM2D::Export::GridVTK(g3, "g1.vtk");
		HM2D::Export::GridMSH(g3, "g1.msh");
		add_file_check(17424596082284710940U, "g1.msh", "grid to fluent");

		HM2D::Export::GridMSH(g3, "g2.msh", [](int i){ return (i==2)?"sqr":"circ";});
		add_file_check(2210112131490853668U, "g2.msh", "grid to fluent with bnd names");
	}
	{
		auto g3 = HM2D::Grid::Constructor::RectGrid01(4, 4);
		build_bcond(g3, bfun2);
		HM2D::Export::PeriodicData dt;

		dt.clear(); dt.add_data(1, 2, false);
		HM2D::Export::GridMSH(g3, "g1.msh", dt);
		add_file_check(10436800471085896365U, "g1.msh", "fluent simple direct periodic");

		dt.clear(); dt.add_data(1, 2, true);
		HM2D::Export::GridMSH(g3, "g2.msh", dt);
		add_file_check(5529949316670468060U, "g2.msh", "fluent simple reversed periodic");

		dt.add_data(3, 4, true);
		HM2D::Export::GridMSH(g3, "g3.msh", dt);
		add_file_check(6149640790734922004U, "g3.msh", "fluent with 2 periodic boundaries");
	}
	{
		HM2D::GridData g4 = HM2D::Grid::Constructor::Circle(Point(0.5, 0.5), 0.1, 10, 3, true);
		HM2D::Grid::Algos::OptUnite uopt(0.1);
		uopt.empty_holes=true; uopt.preserve_bp=true;
		auto g3 = HM2D::Grid::Algos::UniteGrids(g1, g4, uopt);
		build_bcond(g3, bfun2);
		HM2D::Export::PeriodicData dt;
		dt.add_data(1, 2, true);
		dt.add_data(3, 4, true);
		HM2D::Export::GridMSH(g3, "g1.msh", [](int i)->std::string{
					if (i==1 || i==2) return "periodic-short";
					if (i==3 || i==4) return "periodic-long";
					return "no-periodic";
				}, dt);
		add_file_check(16736560883794732383U, "g1.msh", "periodic with complicated mesh");
	}
}

void test24(){
	std::cout<<"24. Pebi grid building"<<std::endl;
	{
		auto tree = tree_w_constraints({{Point(0, 0), Point(1, 0), Point(1, 1), Point(0, 1)}}, {}, 0.1);
		auto trig = HM2D::Mesher::UnstructuredTriangle(tree);
		auto g1 = HM2D::Grid::Constructor::TriToPebi(trig);
		auto els = HM2D::Contour::ELengths(g1.vedges);
		add_check(maxskew(g1)<0.5, "skewness check");
		add_check(*min_element(els.begin(), els.end())>0.005, "edges lengths");
	}
	{
		vector<Point> tt {Point(0, 0), Point(1, 0), Point(1, 1), Point(0, 1), Point(0, 0.02)};
		auto tree = tree_w_constraints({tt}, {}, 0.1);
		auto trig = HM2D::Mesher::UnstructuredTriangle(tree);
		auto g1 = HM2D::Grid::Constructor::TriToPebi(trig);
		HM2D::Export::GridVTK(trig, "g1.vtk");
		HM2D::Export::GridVTK(g1, "g2.vtk");
		add_check(maxskew(g1)<0.7, "pebi points out of area");
	}
}
void test25(){
	std::cout<<"25. Import/Export"<<std::endl;

	auto g1 = HM2D::Grid::Constructor::RectGrid01(3, 3);
	auto g2 = HM2D::Grid::Constructor::Circle(Point(0, 0), 10, 10, 4, false);

	//mixed
	{
		auto writer = HM2D::Export::GridHMG(g1, "g1.hmg", "grid1", "ascii");
		writer.AddCellVertexConnectivity();
		writer.AddCellEdgeConnectivity();
		std::vector<int> somedata1(g1.vcells.size(), 2);
		std::vector<int> somedata2(g1.vcells.size(), 1);
		writer.AddCellData("data1", somedata1, false);
		writer.AddCellData("data1", somedata2, false);
		writer.AddEdgeData("data3", vector<double>(24, 0.123), true);
		writer.AddVertexData("data2", vector<char>(g1.vvert.size(), 12), true);
		writer.AddEdgeData("vecdata", vector<vector<float>>(24, vector<float>(5, -1./6.)), false);
		vector<vector<char>> cd(g1.vcells.size(), vector<char>(2, 23));
		cd[0].resize(3, 0);
		writer.AddCellData("vecdata1", cd, false);
		writer.AddCellData("vecdata2", cd, true);
		writer.Flush();
		add_filesize_check(3655, "g1.hmg", "mixed grid output");
	}

	//pure binary
	{
		auto writer2 = HM2D::Export::GridHMG(g1, "g2.hmg", "grid1", "bin");
		writer2.AddCellEdgeConnectivity();
		writer2.Flush();
		add_filesize_check(1444, "g2.hmg", "binary grid output");
	}

	//multiple grids output
	//ascii
	{
		auto writer3 = HM2D::Export::MultipleGridsHMG({&g1, &g2}, {"grid1", "grid2"}, "g3.hmg", "ascii");
		writer3.sub(0)->AddCellEdgeConnectivity();
		writer3.sub(1)->AddCellVertexConnectivity();
		writer3.sub(1)->AddCellVertexConnectivity();
		writer3.Flush();
		add_file_check(11365718584822709592U, "g3.hmg", "multiple grid output ascii");
	}

	//binary
	{
		auto writer4 = HM2D::Export::MultipleGridsHMG({&g1, &g2}, {"grid1", "grid2"}, "g4.hmg", "bin");
		writer4.sub(0)->AddCellEdgeConnectivity();
		writer4.sub(1)->AddCellVertexConnectivity();
		writer4.sub(0)->AddCellData("somedata", vector<vector<float>>(9, vector<float>(2, 1.)), false);
		writer4.Flush();
		add_filesize_check(4568, "g4.hmg", "multiple grid output binary");
	}

	//bin floats
	{
		auto writer5 = HM2D::Export::MultipleGridsHMG({&g1, &g2}, {"grid1", "grid2"}, "g5.hmg", "binfloat");
		writer5.sub(0)->AddCellEdgeConnectivity();
		writer5.sub(1)->AddCellVertexConnectivity();
		writer5.sub(0)->AddCellData("somedata", vector<vector<float>>(9, vector<float>(2, 1.)), false);
		writer5.Flush();
		add_filesize_check(3579, "g5.hmg", "multiple grid output binary floats");
	}

	// ================================ READER
	{
		auto rd1 = HM2D::Import::GridHMG("g1.hmg");
		auto& r1 = rd1.result;
		HM2D::Grid::Algos::UniqueRearrange(g1);
		HM2D::Grid::Algos::UniqueRearrange(*r1);
		HM2D::Export::GridVTK(g1, "g1.vtk");
		HM2D::Export::GridVTK(*r1, "g2.vtk");
		add_file_check("g1.vtk", "g2.vtk", "read ascii grid");

		auto r2 = HM2D::Import::GridHMG("g2.hmg").result;
		HM2D::Grid::Algos::UniqueRearrange(*r2);
		HM2D::Export::GridVTK(*r2, "g3.vtk");
		add_file_check("g1.vtk", "g3.vtk", "read binary grid");

		auto vf = rd1.vertices_fields();
		auto ef = rd1.edges_fields();
		auto cf = rd1.cells_fields();
		add_check(vf.size()==1 && ef.size()==2 && cf.size()==5, "read field names");

		add_check(
			[&](){
				vector<int> v = rd1.read_vertices_field<int>("data2"); 
				for (auto it: v) if (it != 12) return false;
				return true; }(),
			"reading char binary field with converting to int");

		add_check(
			[&](){
				vector<vector<long>> v = rd1.read_cells_vecfield<long>("__cell_vertices__");
				for (auto it: v) if (it.size()!=4) return false;
				if (v[3][1]!=5) return false;
				if (v[1][3]!=5) return false;
				return true; }(),
			"reading int vector[4] field with converting to long");
		add_check(
			[&](){
				vector<vector<double>> v = rd1.read_cells_vecfield<double>("vecdata2");
				if (v.size() != 9) return false;
				if (v[0].size() != 3) return false;
				if (v[8].size() != 2) return false;
				if (v[0][2] != 0) return false;
				for (int i=0; i<v.size(); ++i)
					if (v[i][0] != 23 || v[i][1] != 23) return false;
				return true; }(),
			"reading variable length vector");
	}

}

void test26(){
	std::cout<<"26. Assemble grid boundary"<<std::endl;
	{
		auto g1 = HM2D::Grid::Constructor::RectGrid01(10, 10);
		HM2D::Grid::Algos::RemoveCells(g1, {43, 44, 53, 54, 65, 66, 75, 76});
		vector<HM2D::EdgeData> conts = HM2D::Contour::Assembler::GridBoundary(g1);
		add_check(conts.size() == 3 &&
		          ISEQ(HM2D::Contour::Length(conts[0]), 4) &&
		          ISEQ(HM2D::Contour::Length(conts[1]), 0.8) &&
		          ISEQ(HM2D::Contour::Length(conts[2]), 0.8), "two inner areas");
	}
	{
		auto g1 = HM2D::Grid::Constructor::RectGrid01(10, 10);
		HM2D::Grid::Algos::RemoveCells(g1, {35,  43, 44, 53, 54, 65, 66, 75, 76});
		vector<HM2D::EdgeData> conts = HM2D::Contour::Assembler::GridBoundary(g1);
		add_check(conts.size() == 4 &&
		          ISEQ(HM2D::Contour::Length(conts[0]), 4) &&
		          ISEQ(HM2D::Contour::Length(conts[1]), 0.4) &&
		          ISEQ(HM2D::Contour::Length(conts[2]), 0.8) &&
	                  ISEQ(HM2D::Contour::Length(conts[3]), 0.8), "three inner areas");
	}
	{
		auto g1 = HM2D::Grid::Constructor::RectGrid01(10, 10);
		HM2D::Grid::Algos::RemoveCells(g1, {35, 46, 57, 43, 44, 53, 54, 65, 66, 75, 76});
		auto tree = HM2D::Contour::Tree::GridBoundary(g1);
		add_check(tree.nodes.size() == 3 &&
		          ISEQ(tree.area(), 0.89), "multiple connected inner areas");
	}
}

void test27(){
	std::cout<<"27. Cells substraction"<<std::endl;
	{
		auto g1 = HM2D::Grid::Constructor::RectGrid01(100, 100);
		auto c1 = HM2D::Contour::Constructor::FromPoints(
				{0.333,0.4354, 0.8,0.6, 1.1,1.8}, true);
		HM2D::Contour::Tree t1;
		t1.add_contour(c1);

		auto g2 = HM2D::Grid::Algos::SubstractCells(
			g1, t1, HM2D::Grid::Algos::SubstractCellsAlgo::PARTLY_INSIDE);
		add_check(g2.vcells.size() == 8347 && g2.vvert.size() == 8636, "partly inside");
		auto g3 = HM2D::Grid::Algos::SubstractCells(
			g1, t1, HM2D::Grid::Algos::SubstractCellsAlgo::PARTLY_OUTSIDE);
		add_check(g3.vcells.size() == 1463 && g3.vvert.size() == 1573, "partly outside");

		auto g4 = HM2D::Grid::Algos::SubstractCells(
			g1, t1, HM2D::Grid::Algos::SubstractCellsAlgo::FULLY_INSIDE);
		add_check(g4.vcells.size() == 8537 && g4.vvert.size() == 8824, "fully inside");

		auto g5 = HM2D::Grid::Algos::SubstractCells(
			g1, t1, HM2D::Grid::Algos::SubstractCellsAlgo::CROSS);
		add_check(g5.vcells.size() == 9789 && g5.vvert.size() == 10188, "cross");

		auto g6 = HM2D::Grid::Algos::SubstractCells(
			g1, t1, HM2D::Grid::Algos::SubstractCellsAlgo::NO_CROSS);
		add_check(g6.vcells.size() == 211 && g6.vvert.size() == 413, "no cross");
	}
	{
		auto g1 = HM2D::Grid::Constructor::RectGrid01(124, 102);
		auto c1 = HM2D::Contour::Constructor::FromPoints(
				{0.2,15./102., 0.5,15./102., 0.9,0.2, 0.8,0.8, 0.1,0.84}, true);
		auto c2 = HM2D::Contour::Constructor::Circle(16, 0.1, Point(0.5,0.5));
		auto c3 = HM2D::Contour::Constructor::Rectangle(Point(-0.3,-0.3), Point(0.15, 0.15));
		auto c4 = HM2D::Contour::Constructor::FromPoints({0.46,0.5, 0.5,0.5, 0.53,0.47}, false);
		HM2D::Contour::Tree t1;
		t1.add_contour(c1); t1.add_contour(c2);
		t1.add_contour(c3); t1.add_contour(c4);

		auto g2 = HM2D::Grid::Algos::SubstractCells(
			g1, t1, HM2D::Grid::Algos::SubstractCellsAlgo::NO_CROSS);
		add_check(g2.vcells.size() == 521 && g2.vvert.size() == 1003, "multilevel, no_cross");

		auto g3 = HM2D::Grid::Algos::SubstractCells(
			g1, t1, HM2D::Grid::Algos::SubstractCellsAlgo::FULLY_OUTSIDE);
		add_check(g3.vcells.size() == 5947 && g3.vvert.size() == 6199, "multilevel, fully outside");

		for (auto& pnt: g1.vvert) pnt->set(vecRotate(*pnt, 37/180.0*M_PI));
		for (auto& pnt: HM2D::AllVertices(t1.alledges())) pnt->set(vecRotate(*pnt, 37/180.0*M_PI));

		g2 = HM2D::Grid::Algos::SubstractCells(
			g1, t1, HM2D::Grid::Algos::SubstractCellsAlgo::NO_CROSS);
		add_check(g2.vcells.size() == 521 && g2.vvert.size() == 1003, "multilevel, no_cross, rotated");

		g3 = HM2D::Grid::Algos::SubstractCells(
			g1, t1, HM2D::Grid::Algos::SubstractCellsAlgo::FULLY_OUTSIDE);
		add_check(g3.vcells.size() == 5947 && g3.vvert.size() == 6199, "multilevel, fully outside, rotated");

	}
}

void test28(){
	std::cout<<"28. Inscribe grids"<<std::endl;
	//{
	//        auto g1 = HM2D::Grid::Constructor::RectGrid01(10, 10);
	//        auto c1 = HM2D::Contour::Constructor::Circle(32, 0.44, Point(0.72, 0.61));
	//        HM2D::Contour::Tree t1;
	//        t1.add_contour(c1);
	//        auto g2 = HM2D::Grid::Algos::InscribeGrid(g1, t1,
	//                HM2D::Grid::Algos::OptInscribe(0.2, true, 0, true, 0.0));
	//        add_check(fabs(HM2D::Contour::Area(c1) - HM2D::Grid::Area(g2))<1e-12 &&
	//                has_edges(g2.vedges, g1.vvert, {83,84, 84,73, 73,74}),
	//                        "square grid into intersecting circle");

	//        auto g3 = HM2D::Grid::Algos::InscribeGrid(g1, t1,
	//                HM2D::Grid::Algos::OptInscribe(0.2, false, 0, true, 0.0));
	//        add_check(fabs(0.48329 - HM2D::Grid::Area(g3))<1e-2 &&
	//                has_edges(g3.vedges, g1.vvert, {0,1, 99,100}) &&
	//                has_edges(g3.vedges, HM2D::AllVertices(c1), {11,12, 19,20}),
	//                        "square grid outside intersecting circle");
	//        HM2D::Export::GridVTK(g3, "g2.vtk");
	//}
	//{
	//        auto g1 = HM2D::Grid::Constructor::RegularHexagonal(Point(-1.1, -1.0), Point(2.04, 2.11), 0.05, true);
	//        auto gc1 = HM2D::Grid::Constructor::RectGrid01(11, 11);
	//        auto c1 = HM2D::Contour::Tree::GridBoundary(gc1);

	//        auto g2 = HM2D::Grid::Algos::InscribeGrid(g1, c1,
	//                HM2D::Grid::Algos::OptInscribe(0.1, true, 0, true, 0.0));
	//        add_check(fabs(1.0 - HM2D::Grid::Area(g2))<1e-12 &&
	//                has_edges(g2.vedges, g1.vvert, {1779,1692, 1922,1963}) &&
	//                has_edges(g2.vedges, HM2D::AllVertices(c1.alledges()), {1,2, 33,34}),
	//                        "inside, keep source");
	//        auto g3 = HM2D::Grid::Algos::InscribeGrid(g1, c1,
	//                HM2D::Grid::Algos::OptInscribe(0.1, true, 0, false, 0.0));
	//        add_check(fabs(1.0 - HM2D::Grid::Area(g3))<1e-12 &&
	//                has_edges(g3.vedges, g1.vvert, {1779,1692, 1922,1963}) &&
	//                edge_len_within(g3.vedges, 0.03, 0.075),
	//                        "inside, don't keep source");
	//        HM2D::Contour::Tree t2; t2.add_contour(HM2D::Contour::Constructor::Rectangle(Point(0,0), Point(1,1)));
	//        auto g4 = HM2D::Grid::Algos::InscribeGrid(g1, t2,
	//                HM2D::Grid::Algos::OptInscribe(0.1, false, 1, false, 0.0));
	//        add_check(fabs(HM2D::Grid::Area(g1) - 1. - HM2D::Grid::Area(g4))<1e-12 &&
	//                has_edges(g4.vedges, g1.vvert, {1770,1769}) &&
	//                edge_len_within(g4.vedges, 0.015, 0.15),
	//                        "outside, don't keep source");

	//}
	{
		auto g1 = HM2D::Grid::Constructor::Ring(Point(0, 0), 3, 0.5, 64, 20);
		HM2D::Contour::Tree t1;
		t1.add_contour(HM2D::Contour::Constructor::Rectangle(Point(-1.7, -1.7), Point(1.81, 1.83)));
		t1.add_contour(HM2D::Contour::Constructor::Circle(32, 0.3, Point(0.8, 0.8)));
		t1.add_contour(HM2D::Contour::Constructor::Rectangle(Point(-1.2, -1.2), Point(-0.3, -0.3)));
		for (auto e: t1.nodes[0]->contour) e->boundary_type = 1;
		for (auto e: t1.nodes[1]->contour) e->boundary_type = 2;
		for (auto e: t1.nodes[2]->contour) e->boundary_type = 3;
		for (auto e: HM2D::ECol::Assembler::GridBoundary(g1)) e->boundary_type = 4;

		//auto g2 = HM2D::Grid::Algos::InscribeGrid(g1, t1,
		//        HM2D::Grid::Algos::OptInscribe(0.3, true, 0, false, 30.0));
		//auto cg2 = HM2D::Contour::Assembler::GridBoundary(g2);
		//std::sort(cg2.begin(), cg2.end(), [](const HM2D::EdgeData& a, const HM2D::EdgeData& b){ return HM2D::Contour::Length(a)<HM2D::Contour::Length(b); });
		//add_check(fabs(t1.area() - HM2D::Grid::Area(g2))<1e-1 &&
		//        has_edges(g2.vedges, g1.vvert, {1056,1120, 648,584}) &&
		//        cg2.size() == 3 &&
		//        std::all_of(cg2[0].begin(), cg2[0].end(), [](const shared_ptr<HM2D::Edge>& e){ return e->boundary_type == 2; }) &&
		//        std::all_of(cg2[1].begin(), cg2[1].end(), [](const shared_ptr<HM2D::Edge>& e){ return e->boundary_type == 3; }) &&
		//        std::all_of(cg2[2].begin(), cg2[2].end(), [](const shared_ptr<HM2D::Edge>& e){ return e->boundary_type == 1; }),
		//                "multilevel, inside");

		auto g3 = HM2D::Grid::Algos::InscribeGrid.WTimer(g1, t1,
			HM2D::Grid::Algos::OptInscribe(0.1, false, 0, false, 30.0));
		add_check(fabs(16.9189 - HM2D::Grid::Area(g3))<1e-2, "multilevel, outside");
	}
}

void test29(){
	//{
	//        auto g1 = HM2D::Grid::Constructor::RectGrid01(14, 14);
	//        auto c1 = HM2D::Contour::Constructor::FromPoints({0.4,0.2, 0.8,0.8}, false);
	//        std::map<double, double> mp {{0, 0.05}, {1, 0.1}};
	//        auto c2 = HM2D::Contour::Algos::WeightedPartition(mp, c1);
	//        HM2D::Contour::Tree t2; t2.add_detached_contour(c2);

	//        auto g2 = HM2D::Grid::Algos::InscribeGrid(g1, t2,
	//                HM2D::Grid::Algos::OptInscribe(0.2, true, 0, true, 0.0));

	//        HM2D::Export::GridVTK(g1, "g1.vtk");
	//        HM2D::Export::ContourVTK(t2.alledges(), "c1.vtk");
	//        HM2D::Export::GridVTK(g2, "g2.vtk");
	//}
}

int main(){
	//test0();
	//test1();
	//test2();
	//test3();
	//test4();
	//test5();
	//test6();
	//test7();
	//test8();
	//test9();
	//test10();
	//test11();
	//test12();
	//test13();
	//test14();
	//test15();
	//test16();
	//test17();
	//test18();
	//test19();
	//test20();
	//test21();
	//test22();
	//test23();
	//test24();
	//test25();
	//test26();
	//test27();
	test28();
	//test29();

	HMTesting::check_final_report();
	std::cout<<"DONE"<<std::endl;
}
