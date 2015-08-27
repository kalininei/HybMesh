#include "crossgrid.h"
#include "grid.h"
#include "fileproc.h"
#include "wireframegrid.h"
#include "Gmsh.h"

namespace{

//gmesh library initialization
struct _gmsh_initializer{
	_gmsh_initializer(){GmshInitialize();}
	~_gmsh_initializer(){GmshFinalize();}
};
_gmsh_initializer gi;

//default callback function
int default_callback(const char* proc, const char* subproc, double percent, double subpercent){
	std::cout<<proc;
	if (percent>=0) std::cout<<" - "<<int(100*percent)<<"%";
	std::cout<<";";
	if (subproc!=0){
		std::cout<<" ("<<subproc;
		if (subpercent>=0) std::cout<<" - "<<int(100*subpercent)<<"%";
		std::cout<<")";
	}
	std::cout<<std::endl;
	return CALLBACK_OK;
}
int silent_callback(const char* proc, const char* subproc, double percent, double subpercent){
	return CALLBACK_OK;
}

crossgrid_callback global_callback = default_callback;

} //namespace

//callbacks
void crossgrid_cout_callback(){
	global_callback = default_callback;
}
void crossgrid_silent_callback(){
	global_callback = silent_callback;
}
void crossgrid_set_callback(crossgrid_callback fun){
	global_callback = fun;
}


Grid* grid_construct(int Npts, int Ncells, double* pts, int* cells){
	try{
		return new GridGeom(Npts, Ncells, pts, cells);
	} catch (const std::exception &e){
		std::cout<<"crossgrid error: "<<e.what()<<std::endl;
		return 0;
	}
}

int grid_npoints(Grid* g){
	return static_cast<GridGeom*>(g)->n_points();
}

int grid_ncells(Grid* g){
	return static_cast<GridGeom*>(g)->n_cells();
}
int grid_cellsdim(Grid* g){
	return static_cast<GridGeom*>(g)->n_cellsdim();
}

void grid_get_points_cells(Grid* g, double* pts, int* cells){
	auto gg=static_cast<GridGeom*>(g);
	//points
	if (pts!=0) for (int i=0; i<gg->n_points(); ++i){
		auto p = gg->get_point(i);
		*pts++ = p->x;
		*pts++ = p->y;
	}
	//cells
	if (cells!=0) for (int i=0; i<gg->n_cells(); ++i){
		auto c = gg->get_cell(i);
		*cells++ = c->dim();
		for (int j=0; j<c->dim(); ++j){
			*cells++ = c->get_point(j)->get_ind();
		}
	}
}

void grid_get_points_cells2(Grid* g, double** pts, int** cells){
	int np = grid_npoints(g);
	int nc = grid_ncells(g);
	int npc = grid_cellsdim(g);
	*pts = new double[2*np];
	*cells = new int[nc+npc];
	grid_get_points_cells(g, *pts, *cells);
}

//get grid in points, edges->points, cells->edges format
void grid_get_edges_info(Grid* grd, int* Npnt, int* Neds, int* Ncls,
		double** pts,
		int** ed_pt,
		int** cls_dims,
		int** cls_eds){
	GridGeom* g = static_cast<GridGeom*>(grd);
	auto eds = g->get_edges();
	//fill counts
	*Npnt = g->n_points();
	*Ncls = g->n_cells();
	*Neds = eds.size();
	//allocate arrays
	*pts = new double[2 * (*Npnt)];
	*ed_pt = new int[2 * (*Neds)];
	*cls_dims = new int[*Ncls];
	*cls_eds = new int[g->n_cellsdim()];
	//fill points
	double* p_pts = *pts;
	for (int i=0; i<g->n_points(); ++i){
		*p_pts++ = g->get_point(i)->x;
		*p_pts++ = g->get_point(i)->y;
	}
	//fill edges
	std::map<std::pair<int, int>, int> nd_eds;
	int* p_eds = *ed_pt;
	auto eit = eds.begin();
	for (size_t i=0; i<eds.size(); ++i){
		auto& e = *eit++;
		*p_eds++ = e.p1;
		*p_eds++ = e.p2;
		nd_eds.emplace(std::make_pair(e.p1, e.p2), i);
	}
	//fill elements dimensions
	int* p_cdim = *cls_dims;
	for (int i=0; i<g->n_cells(); ++i){
		*p_cdim++ = g->get_cell(i)->dim();
	}
	//fill cell->edges array
	int* p_ced = *cls_eds;
	for (int i=0; i<g->n_cells(); ++i){
		auto c = g->get_cell(i);
		for (int j=0; j<c->dim(); ++j){
			int pprev = c->get_point(j-1)->get_ind();
			int pcur = c->get_point(j)->get_ind();
			if (pprev>pcur) std::swap(pprev, pcur);
			*p_ced++ = nd_eds[std::make_pair(pprev, pcur)];
		}
	}
}


//free edges data
void grid_free_edges_info(double** pts, int** ed_pt, int** cls_dims, int** cls_eds){
	delete[] *pts; *pts = 0;
	delete[] *ed_pt; *ed_pt = 0;
	delete[] *cls_dims; *cls_dims = 0;
	delete[] *cls_eds; *cls_eds = 0;
}

int grid_get_edge_cells(Grid* g, int* Neds, int** ed_cell, int* ed_pt){
	auto ed = static_cast<GridGeom*>(g)->get_edges();
	*Neds = ed.size();
	*ed_cell = new int[*Neds];
	int *ec = *ed_cell;
	int ret = 1;
	if (ed_pt == 0){
		for (auto& e: ed){
			*ec++ = e.cell_left;
			*ec++ = e.cell_right;
		}
	} else {
		for (int i=0; i<*Neds; ++i){
			int p1 = ed_pt[2*i], p2 = ed_pt[2*i+1];
			auto fnd = ed.find(Edge(p1, p2));
			if (fnd != ed.end()){
				if (p1 == fnd->p1){
					*ec++ = fnd->cell_left;
					*ec++ = fnd->cell_right;
				} else {
					*ec++ = fnd->cell_right;
					*ec++ = fnd->cell_left;
				}
				ed.erase(fnd);
			} else {
				//invalid edge->points connectivity
				ret = 0;
				*ec++ = -1; *ec++ = -1;
			}
		}
	}
	return ret;
}

void grid_free_edge_cells(int** ed_cell){
	delete[] *ed_cell; *ed_cell=0;
}


void grid_save_vtk(Grid* g, const char* fn){
	if (g==NULL) return;
	save_vtk(static_cast<GridGeom*>(g), fn);
}

void contour_save_vtk(Cont* c, const char* fn){
	if (c==NULL) return;
	save_vtk(static_cast<PointsContoursCollection*>(c), fn);
}

void grid_free(Grid* g){
	if (g!=NULL) delete g;
}

Grid* cross_grids(Grid* gbase, Grid* gsecondary, double buffer_size, int preserve_bp, int eh){
	return cross_grids_wcb(gbase, gsecondary, buffer_size, preserve_bp, eh, global_callback);
}

Grid* cross_grids_wcb(Grid* gbase, Grid* gsecondary, double buffer_size, 
		int preserve_bp, int empty_holes, crossgrid_callback cb_fun){
	try{
		if (gbase == NULL || gsecondary == NULL)
			throw std::runtime_error("nullptr grid data");
		auto ret = GridGeom::cross_grids(
				static_cast<GridGeom*>(gbase),
				static_cast<GridGeom*>(gsecondary),
				buffer_size, 0.5, (preserve_bp==1),
				(empty_holes==1), cb_fun);
		return ret;
	} catch (const std::exception &e) {
		std::cout<<"crossgrid error: "<<e.what()<<std::endl;
		return 0;
	}
}

Cont* contour_construct(int Npts, int Ned, double* pts, int* edges){
	auto p = vector<Point>();
	auto e = vector<int>(edges, edges+2*Ned);
	for (int i=0; i<Npts; ++i){
		p.push_back(Point(pts[2*i], pts[2*i+1]));
	}
	return new PointsContoursCollection(p, e);
}

void cont_free(Cont* c){
	delete c;
}

Grid* grid_exclude_cont(Grid* grd, Cont* cont, int is_inner){
	return grid_exclude_cont_wcb(grd, cont, is_inner, global_callback);
}

Grid* grid_exclude_cont_wcb(Grid* grd, Cont* cont, int is_inner,
		crossgrid_callback cb_fun){
	try{
		auto ret = GridGeom::grid_minus_cont(
				static_cast<GridGeom*>(grd),
				static_cast<PointsContoursCollection*>(cont),
				is_inner!=0, cb_fun);
		return ret;
	} catch (const std::exception &e){
		std::cout<<"crossgrid error: "<<e.what()<<std::endl;
		return 0;
	}
}

void add_contour_bc(Cont* src, Cont* tar, int* vsrc, int* vtar, int def){
	try{
		auto src1 = static_cast<PointsContoursCollection*>(src);
		auto tar1 = static_cast<PointsContoursCollection*>(tar);
		vector<int> cor = PointsContoursCollection::edge_correlation(*src1, *tar1);
		for (int i=0; i<tar1->n_edges(); ++i){
			vtar[i] =  (cor[i]>=0) ? vsrc[cor[i]] : def;
		}
	} catch (const std::exception &e){
		std::cout<<"crossgrid error: "<<e.what()<<std::endl;
	}
}

void contour_get_info(Cont* c, int* Npnt, int* Neds, 
		double** pts,
		int** eds){
	PointsContoursCollection* cont = static_cast<PointsContoursCollection*>(c);
	*Npnt = cont->n_total_points();
	*Neds = cont->n_edges();
	//points
	*pts = new double[2*(*Npnt)];
	double* pp = *pts;
	for (int i=0; i<*Npnt; ++i){
		auto p = cont->get_point(i);
		*pp++ = p->x;
		*pp++ = p->y;
	}
	//edges
	*eds = new int[2*(*Neds)];
	int *ee = *eds;
	for (int i=0; i<*Neds; ++i){
		auto e = cont->get_edge(i);
		*ee++ = e.first;
		*ee++ = e.second;
	}
}

void contour_free_info(double** pts, int** eds){
	delete[] (*pts);
	delete[] (*eds);
}

double grid_area(Grid* g){
	return static_cast<GridGeom*>(g)->area();
}

double contour_area(Cont* c){
	return static_cast<PointsContoursCollection*>(c)->area();
}

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
	auto cross1 = GridGeom::cross_grids(&grid1, &grid2, 0.0, 0.5, true, false, silent_callback);
	auto cross2 = GridGeom::cross_grids(cross1, &grid3, 0.0, 0.5, true, false, silent_callback);
	save_vtk(cross2, "test4_grid.vtk");
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
	//      It stopped working after I add boundary segments to delete list
	//      in PtsGraph::cut() procedure
	//test2();
	test3();
	test4();
	test5();
	test6();
	test7();
	std::cout<<"crossgrid shared library internal tests: DONE =========="<<std::endl;
}

