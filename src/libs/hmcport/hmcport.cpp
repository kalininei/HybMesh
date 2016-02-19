#include "hmcport.h"
#include "fileproc.h"
#include "grid.h"
#include "hmblay.hpp"

namespace{

// crossgrid callback
CrossGridCallback::func global_callback;

}//namespace

void crossgrid_cout_callback(){
	global_callback = CrossGridCallback::to_cout();
}

void crossgrid_silent_callback(){
	global_callback = CrossGridCallback::silent();
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
	try{
		return new PointsContoursCollection(p, e);
	} catch (const std::exception &e){
		std::cout<<"crossgrid error: "<<e.what()<<std::endl;
		return 0;
	}
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

void* create_ecollection_container(int Npts, double* pts, int Nedgs, int* edges){
	try{
		auto ret = new HMCont2D::Container<HMCont2D::ECollection>();
		//points
		for (int i=0; i<Npts; ++i){
			ret->pdata.add_value(Point(pts[2*i], pts[2*i+1]));
		}
		//edges
		for (int i=0; i<Nedgs; ++i){
			ret->add_value(HMCont2D::Edge(ret->point(edges[2*i]), ret->point(edges[2*i+1])));
		}
		return ret;
	} catch (const std::exception &e){
		std::cout<<"hybmesh_contours2d error: "<<e.what()<<std::endl;
		return 0;
	}
}
void free_ecollection_container(void* ecol){
	delete static_cast<HMCont2D::Container<HMCont2D::ECollection>*>(ecol);
}

int set_ecollection_bc(void* src, void* tar, int def, int* vsrc, int* vtar){
	try{
		auto esrc = static_cast<HMCont2D::ECollection*>(src);
		auto etar = static_cast<HMCont2D::ECollection*>(tar);
		for (int i=0; i<etar->size(); ++i){
			Point cpoint = etar->edge(i)->center();
			auto ce = HMCont2D::ECollection::FindClosestEdge(*esrc, cpoint);
			if (ISZERO(std::get<1>(ce))){
				vtar[i] = vsrc[std::get<3>(ce)];
			} else {
				vtar[i] = def;
			}
		}
		return 0;
	} catch (const std::exception &e){
		std::cout<<"set_ecollection_bc error: "<<e.what()<<std::endl;
		return 1;
	}
}


double ecollection_area(void* ecol){
	try{
		auto c = static_cast<HMCont2D::ECollection*>(ecol);
		auto tree = HMCont2D::ExtendedTree::Assemble(*c);
		return HMCont2D::Area(tree);
	} catch (const std::exception &e){
		std::cout<<"domain area calculation error: "<<e.what()<<std::endl;
		return 0;
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

template <class C>
C string_option(std::string opt,
		std::initializer_list<std::string> opt_val,
		std::initializer_list<C> ret_val){
	std::transform(opt.begin(), opt.end(), opt.begin(), ::toupper);
	auto ival = opt_val.begin();
	auto rval = ret_val.begin();
	while (ival != opt_val.end() || rval != ret_val.end()){
		auto s = *ival;
		std::transform(s.begin(), s.end(), s.begin(), ::toupper);
		if (*ival == opt) return *rval;
		++ival;
		++rval;
	}
	throw std::runtime_error("impossible to parse option " + opt);
}

Grid* boundary_layer_grid_wcb(int N, BoundaryLayerGridOption* popt, 
		crossgrid_callback cb_fun){
	try{
		vector<HMBlay::Input> vinp(N);
		for (int i=0; i<N; ++i){
			auto& opt = popt[i];
			auto& inp = vinp[i];
			inp.edges = static_cast<HMCont2D::ECollection*>(opt.cont);
			inp.direction = HMBlay::DirectionFromString(opt.tp);
			inp.bnd_step_method = HMBlay::MethFromString(opt.mesh_cont);
			inp.bnd_step = opt.mesh_cont_step;
			inp.acute_angle = opt.angle_range[0];
			inp.right_angle = opt.angle_range[1];
			inp.straight_angle = opt.angle_range[2];
			inp.reentrant_angle = opt.angle_range[3];
			inp.start = Point(opt.start[0], opt.start[1]);
			inp.end = Point(opt.end[0], opt.end[1]);
			inp.partition = vector<double>(opt.part, opt.part + opt.Npart);
			inp.force_conformal = (opt.force_conformal == 1);
		}
		GridGeom* ret = new GridGeom(HMBlay::BuildBLayerGrid(vinp));
		return ret;
	} catch (const std::exception &e){
		std::cout<<"Boundary layer builder error: "<<e.what()<<std::endl;
		return 0;
	}
}
