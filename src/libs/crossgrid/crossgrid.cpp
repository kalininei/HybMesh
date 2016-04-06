#include "crossgrid.h"
#include "grid.h"
#include "wireframegrid.h"
#include "Gmsh.h"

namespace{

//gmesh library initialization
struct _gmsh_initializer{
	_gmsh_initializer(){GmshInitialize();}
	~_gmsh_initializer(){GmshFinalize();}
};
_gmsh_initializer gi;

} //namespace

// ========================== testing
namespace{

void add_check(bool ex, const char* info = 0){
	if (info==0) std::cout<<"\tunknown check: ";
	else std::cout<<"\t"<<info;
	if (ex) std::cout<<": True"<<std::endl;
	else std::cout<<": False <<<<<<<<<<<<<<<<<<<"<<std::endl;
};

//build a rectangular structured grid
GridGeom rectangular_grid(double x0, double y0,
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
	return GridGeom((Nx+1)*(Ny+1), Nx*Ny, &pts[0], &cls[0]);
}

PointsContoursCollection uniform_polygon(double xc, double yc, int N, double rad){
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
	return PointsContoursCollection(pts, ed);
}

void test1(){
	std::cout<<"contours collection tests"<<std::endl;
	Contour c1, c2, c3, c4, c5;
	c1.add_point(0,0); c1.add_point(1,0);
	c1.add_point(1,1); c1.add_point(0,1);
	c2.add_point(0.2, 0.3); c2.add_point(0.1, 0.1);
	c2.add_point(0.1, 0.9); c2.add_point(0.9, 0.9);
	c3.add_point(0.3, 0.8); c3.add_point(0.2, 0.8);
	c3.add_point(0.2, 0.7); c3.add_point(0.3, 0.7);
	c4.add_point(0.5, 0.1); c4.add_point(0.8, 0.1);
	c4.add_point(0.8, 0.4);

	ContoursCollection cc;
	cc.add_contour(c2);
	cc.add_contour(c1);
	cc.add_contour(c4);
	cc.add_contour(c3);

	add_check(cc.is_inside(Point(0,0))==0, "point within");
	add_check(cc.is_inside(Point(-1,0))==-1, "point within");
	add_check(cc.is_inside(Point(0.3,0.1))==1, "point within");
	add_check(cc.is_inside(Point(0.6,0.1))==0, "point within");
	add_check(cc.is_inside(Point(0.3,0.6))==-1, "point within");
	add_check(cc.is_inside(Point(0.25,0.75))==1, "point within");

	c5.add_point(0.9, 0.05); c5.add_point(0.9, 0.6);
	c5.add_point(0.3, 0.05);
	add_check(cc.is_inside(Point(0.7, 0.2))==-1, "point within");
	cc.add_contour(c5);
	add_check(cc.is_inside(Point(0.7, 0.2))== 1, "point within");
}

void test2(){
	std::cout<<"PtsGraph cut test"<<std::endl;
	double pts[] = {0,0,9,2,8,6,4,1};
	int cls[] = {3,0,3,2, 3,3,1,2};
	GridGeom G(4, 2, pts, cls);
	PtsGraph gr(G);

	Contour c;
	c.add_point(2,-1); c.add_point(9,-1);
	c.add_point(4, 5); c.add_point(2, 5);
	ContoursCollection cc({c});
	PtsGraph newgraph = PtsGraph::cut(gr, cc, OUTSIDE);
	GridGeom newgeom =  newgraph.togrid();
	//save_vtk(&gr, "graph.vtk");
	//save_vtk(&newgraph, "cutgraph.vtk");
	//save_vtk(&newgeom, "cutgrid.vtk");
	add_check(newgeom.n_cells()==4 && newgeom.n_points()==10, "grid geometry");
}

void test3(){
	std::cout<<"Grid combine test"<<std::endl;
	double pts1[] = {0,0, 2,2, 0,4, -2,2, 0,2};
	double pts2[] = {-1,1, 1,1, 1, 2.5, -1,3};
	int cls1[] = {3,0,1,4, 3,4,1,2, 3,3,4,2, 3,3,0,4};
	int cls2[] = {4,0,1,2,3};
	GridGeom gmain(5,4, pts1, cls1);
	GridGeom gsec(4,1, pts2, cls2);
	GridGeom* comb = GridGeom::combine(&gmain, &gsec);
	add_check(comb->n_cells()==8 && comb->n_points()==12, "grid geometry");
	//save_vtk(comb, "combined_grid.vtk");
	delete comb;
}

void test4(){
	std::cout<<"Grid subdivide"<<std::endl;
	auto grid1 = rectangular_grid(0, 0, 1, 1, 10, 10);
	auto grid2 = rectangular_grid(1, 1, 2, 2, 10, 10);
	auto grid3 = rectangular_grid(2,1.85, 3, 2.85, 10, 10);
	auto cross1 = GridGeom::cross_grids(&grid1, &grid2, 0.0, 0.5, true, false, 0);
	auto cross2 = GridGeom::cross_grids(cross1, &grid3, 0.0, 0.5, true, false, 0);
	add_check(cross2->n_points()==362 && cross2->n_cells()==300, "combined grid topology");
	auto div = cross2->subdivide();
	add_check(div.size()==2, "number of single connected grids");
	add_check(div[0]->n_points()==121 && div[0]->n_cells()==100, "grid1 topology");
	add_check(div[1]->n_points()==242 && div[1]->n_cells()==200, "grid2 topology");
	delete cross1;
	delete cross2;
}

void test5(){
	std::cout<<"Node finder"<<std::endl;
	NodeFinder finder(Point(3,3), 1.0, 1.0, 5, 10, 0.01);
	Point p[] = {
		Point(-2, 23),
		Point(3.25, 3.42),
		Point(3.2503, 3.4203),
		Point(3.261, 3.239),
		Point(2.999, 3.0),
		Point(3.2, 3.2),
		Point(3.203, 3.203),
		Point(3.197, 3.203),
		Point(3.197, 3.197),
		Point(3.203, 3.197),
		Point(4, 4)
	};
	add_check(finder.add(p) == 0, "out of rectangle");
	add_check(finder.add(p+1) == p+1, "first addition");
	add_check(finder.add(p+2) == p+1, "equal point");
	add_check(finder.add(p+3) == p+3, "not equal point");
	add_check(finder.add(p+4) == p+4, "point on the bottom left");
	add_check(finder.add(p+10) == p+10, "point on the top right");
	finder.add(p+5);
	add_check(finder.add(p+6) == p+5, "first square equal");
	add_check(finder.add(p+7) == p+5, "second square equal");
	add_check(finder.add(p+8) == p+5, "third square equal");
	add_check(finder.add(p+9) == p+5, "fourth square equal");
}

void test6(){
	std::cout<<"Contours Collection constructing"<<std::endl;
	vector<Point> p = {
		Point(0, 0),
		Point(1, 0),
		Point(1, 1),
		Point(0, 1),
		Point(0.3, 0.3),
		Point(0.6, 0.3),
		Point(0.6, 0.6),
		Point(0.3, 0.6)
	};
	vector<int> e = {
		7, 4,
		0, 1,
		2, 1,
		3, 0,
		2, 3,
		6, 5,
		7, 6,
		5, 4
	};

	auto c = PointsContoursCollection(p, e);
	add_check(fabs(c.area() - 0.91) < 1e-6, "collection area");
}

void test7(){
	std::cout<<"PtsGraph::togrid() with intrusion subprocedure"<<std::endl;
	auto getar = [](const GridGeom& g){
		double ar = 0;
		for (int i=0; i<g.n_cells(); ++i) ar += g.get_cell(i)->area();
		return ar;
	};
{ //case1
	auto c1 = uniform_polygon(0, 0, 6, 4);
	auto c2 = uniform_polygon(0,0,4,2);
	PtsGraph g1(c1);
	g1.add_edges(c2.contour(0));
	auto grid = g1.togrid();
	add_check(grid.n_cells() == 3 && fabs(c1.area()-getar(grid)) < 1e-6, "grid from two polys");
}
{ //case 2
	auto c1 = uniform_polygon(3,2,7, 5);
	auto c2 = uniform_polygon(3,3,4, 2);
	auto c3 = uniform_polygon(3.2,3.1,3, 0.5);
	PtsGraph g1(c1); g1.add_edges(c2.contour(0)); g1.add_edges(c3.contour(0));
	auto grid = g1.togrid();
	add_check(grid.n_cells() == 5 && fabs(c1.area()-getar(grid)) < 1e-6, "grid from 3 nested polys");
}

{ //case 3
	auto c0 = uniform_polygon(3,3,5,0.5);
	auto c1 = uniform_polygon(6,6,5,1);
	auto c2 = uniform_polygon(2,3,8,2);
	auto c3 = uniform_polygon(2.1,3,4,0.1);
	auto c4 = uniform_polygon(2,3,5,0.3);
	auto c5 = uniform_polygon(1,4,3,0.1);
	PtsGraph g1(c0); g1.add_edges(c1.contour(0)); g1.add_edges(c2.contour(0));
	g1.add_edges(c3.contour(0));g1.add_edges(c4.contour(0));g1.add_edges(c5.contour(0));
	auto grid = g1.togrid();
	double a = c1.area() + c2.area();
	add_check(grid.n_cells() == 10 && fabs(a-getar(grid)) < 1e-6, "grid from complicated nested structure");
}


}//test7

}//namespace

void crossgrid_internal_tests(){
	std::cout<<"crossgrid shared library internal tests ================"<<std::endl;
	test1();
	//TODO: check what's wrong with test2.
	//      It stopped working after I added boundary segments to delete list
	//      in PtsGraph::cut() procedure
	//test2();
	test3();
	test4();
	test5();
	test6();
	test7();
	std::cout<<"crossgrid shared library internal tests: DONE =========="<<std::endl;
}

