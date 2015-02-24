#include "crossgrid.h"
#include "grid.h"
#include "fileproc.h"
#include "wireframegrid.h"
#include "gmsh/Gmsh.h"

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
		std::cout<<e.what()<<std::endl;
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
	delete *pts; *pts = 0;
	delete *ed_pt; *ed_pt = 0;
	delete *cls_dims; *cls_dims = 0;
	delete *cls_eds; *cls_eds = 0;
}


void grid_save_vtk(Grid* g, const char* fn){
	save_vtk(static_cast<GridGeom*>(g), fn);
}

void grid_free(Grid* g){
	delete g;
}

Grid* cross_grids(Grid* gbase, Grid* gsecondary, double buffer_size, double density, int preserve_bp){
	return cross_grids_wcb(gbase, gsecondary, buffer_size, density, preserve_bp, global_callback);
}

Grid* cross_grids_wcb(Grid* gbase, Grid* gsecondary, double buffer_size, 
		double density, int preserve_bp, crossgrid_callback cb_fun){
	try{
		auto ret = GridGeom::cross_grids(
				static_cast<GridGeom*>(gbase),
				static_cast<GridGeom*>(gsecondary),
				buffer_size, density, (preserve_bp==1), cb_fun);
		return ret;
	} catch (const std::exception &e) {
		std::cout<<e.what()<<std::endl;
		return 0;
	}
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
	PtsGraph newgraph = PtsGraph::cut(gr, cc, 1);
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
	auto cross1 = GridGeom::cross_grids(&grid1, &grid2, 0.0, 0.5, true, silent_callback);
	auto cross2 = GridGeom::cross_grids(cross1, &grid3, 0.0, 0.5, true, silent_callback);
	save_vtk(cross2, "test4_grid.vtk");
	add_check(cross2->n_points()==362 && cross2->n_cells()==300, "combined grid topology");
	auto div = cross2->subdivide();
	add_check(div.size()==2, "number of single connected grids");
	add_check(div[0]->n_points()==121 && div[0]->n_cells()==100, "grid1 topology");
	add_check(div[1]->n_points()==242 && div[1]->n_cells()==200, "grid2 topology");
	delete cross1;
	delete cross2;
}

}

void crossgrid_internal_tests(){
	std::cout<<"crossgrid shared library internal tests ================"<<std::endl;
	test1();
	test2();
	test3();
	test4();
	std::cout<<"crossgrid shared library internal tests: DONE =========="<<std::endl;
}

