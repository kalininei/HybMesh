#include "hmcport.h"
#include "grid.h"
#include "hmblay.hpp"
#include "procgrid.h"
#include "hmmapping.hpp"
#include "hmtesting.hpp"
#include "hmcallback.hpp"
#include "debug_grid2d.h"
#include "hmxmlreader.hpp"
namespace{
int silent2_function(const char*, const char*, double, double){
	return HMCallback::OK;
}
}

void free_int_array(int* a){ if (a != NULL) delete[] a; }
void free_double_array(double* a){ if (a!= NULL) delete[] a; }
void free_char_array(char* a){ if (a!=NULL) delete[] a; }
void free_voidp_array(void** a){ if (a!=NULL) delete[] a; }

size_t get_ascii_file_hash(const char* fn){
	return HMTesting::calculate_file_hash(fn);
}
BoundaryNamesStruct* set_boundary_names(int n, const char** nm, int* vals){
	BoundaryNamesStruct* ret = new BoundaryNamesStruct();
	ret->n = n;
	ret->values = new int[n];
	ret->names = new const char*[n];
	for (int i=0; i<n; ++i){
		ret->values[i] = vals[i];
		ret->names[i] = nm[i];
	}
	return ret;
}
void free_boundary_names(BoundaryNamesStruct* s){
	delete[] s->values;
	delete[] s->names;
	delete s;
}

Grid* grid_construct(int Npts, int Ncells, double* pts, int* cells){
	try{
		return new GridGeom(Npts, Ncells, pts, cells);
	} catch (const std::exception &e){
		std::cout<<"crossgrid error: "<<e.what()<<std::endl;
		return NULL;
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


void grid_free(Grid* g){
	if (g!=NULL) delete static_cast<GridGeom*>(g);
}

Grid* cross_grids(Grid* gbase, Grid* gsecondary, double buffer_size, int preserve_bp, int eh, double angle0, int algo){
	return cross_grids_wcb(gbase, gsecondary, buffer_size, preserve_bp, eh, angle0, algo, silent2_function);
}

Grid* cross_grids_wcb(Grid* gbase, Grid* gsecondary, double buffer_size,
		int preserve_bp, int empty_holes, double angle0, int algo, hmcport_callback cb_fun){
	try{
		if (gbase == NULL || gsecondary == NULL)
			throw std::runtime_error("nullptr grid data");
		auto ret = GridGeom::cross_grids(
				static_cast<GridGeom*>(gbase),
				static_cast<GridGeom*>(gsecondary),
				buffer_size, 0.5, (preserve_bp==1),
				(empty_holes==1), angle0, algo, cb_fun);
		return ret;
	} catch (const std::exception &e) {
		std::cout<<"crossgrid error: "<<e.what()<<std::endl;
		return 0;
	}
}

//Cont* contour_construct(int Npts, int Ned, double* pts, int* edges){
//        auto p = vector<Point>();
//        auto e = vector<int>(edges, edges+2*Ned);
//        for (int i=0; i<Npts; ++i){
//                p.push_back(Point(pts[2*i], pts[2*i+1]));
//        }
//        try{
//                return new PointsContoursCollection(p, e);
//        } catch (const std::exception &e){
//                std::cout<<"crossgrid error: "<<e.what()<<std::endl;
//                return 0;
//        }
//}

//void cont_free(Cont* c){
//        delete c;
//}

Grid* grid_exclude_cont(Grid* grd, void* cont, int is_inner){
	return grid_exclude_cont_wcb(grd, cont, is_inner, silent2_function);
}

Grid* grid_exclude_cont_wcb(Grid* grd, void* cont, int is_inner,
		hmcport_callback cb_fun){
	try{
		HMCont2D::ECollection* pnt = static_cast<HMCont2D::ECollection*>(cont);
		PointsContoursCollection pcc(*pnt);
		auto ret = GridGeom::grid_minus_cont(
				static_cast<GridGeom*>(grd), &pcc,
				is_inner!=0, cb_fun);
		return ret;
	} catch (const std::exception &e){
		std::cout<<"grid_exclude_cont error: "<<e.what()<<std::endl;
		return 0;
	}
}

int simplify_grid_boundary(Grid* grd, double angle){
	GridGeom* g = static_cast<GridGeom*>(grd);
	ScaleBase sc = g->do_scale();
	int ret = 1;
	try{
		angle = angle/180*M_PI;
		if (angle<-geps || angle>M_PI+geps) throw std::runtime_error("invalid angle");
		GGeom::Modify::SimplifyBoundary(*g, angle);
		ret = 0;
	} catch (const std::exception &e){
		std::cout<<"simplify_grid_boundary error: "<<e.what()<<std::endl;
	}
	g->undo_scale(sc);
	return ret;
}

//void add_contour_bc(Cont* src, Cont* tar, int* vsrc, int* vtar, int def){
//        try{
//                auto src1 = static_cast<PointsContoursCollection*>(src);
//                auto tar1 = static_cast<PointsContoursCollection*>(tar);
//                vector<int> cor = PointsContoursCollection::edge_correlation(*src1, *tar1);
//                for (int i=0; i<tar1->n_edges(); ++i){
//                        vtar[i] =  (cor[i]>=0) ? vsrc[cor[i]] : def;
//                }
//        } catch (const std::exception &e){
//                std::cout<<"crossgrid error: "<<e.what()<<std::endl;
//        }
//}

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

void* domain_clip(void* c1, void* c2, int oper, int simplify){
	auto cont1 = static_cast<HMCont2D::ECollection*>(c1);
	auto cont2 = static_cast<HMCont2D::ECollection*>(c2);
	auto pnt1 = cont1->all_points();
	auto pnt2 = cont2->all_points();

	ScaleBase sc = ScaleBase::p_doscale(pnt1.begin(), pnt1.end()); 
	sc.p_scale(pnt2.begin(), pnt2.end());

	auto tree1 = HMCont2D::Assembler::ETree(*cont1);
	auto tree2 = HMCont2D::Assembler::ETree(*cont2);
	try{
		if (tree1.nodes.size() == 0) throw std::runtime_error("not a closed contour");
		if (tree2.nodes.size() == 0) throw std::runtime_error("not a closed contour");

		HMCont2D::Container<HMCont2D::ContourTree> res;
		switch (oper){
			case 1: res = HMCont2D::Clip::Union(tree1, tree2); break;
			case 2: res = HMCont2D::Clip::Difference(tree1, tree2); break;
			case 3: res = HMCont2D::Clip::Intersection(tree1, tree2); break;
			case 4: res = HMCont2D::Clip::XOR(tree1, tree2); break;
			default: throw std::runtime_error("unknown operation");
		}
		HMCont2D::Clip::Heal(res);
		if (res.nodes.size() !=0 && !simplify){
			vector<Point*> allpnt(pnt1); allpnt.insert(allpnt.end(), pnt2.begin(), pnt2.end());
			for (auto p: allpnt){
				auto fnd = HMCont2D::ECollection::FindClosestEdge(res, *p);
				if (std::get<0>(fnd) != 0 && std::get<1>(fnd)<geps && 
						std::get<2>(fnd) > geps && std::get<2>(fnd) < 1-geps){
					auto cont = res.get_contour(std::get<0>(fnd));
					cont->GuaranteePoint(*p, res.pdata);
					res.ReloadEdges();
				}
			}
		}
		HMCont2D::Container<HMCont2D::ECollection>* r = 0;
		if (res.nodes.size() > 0){
			r = new HMCont2D::Container<HMCont2D::ECollection>();
			r->Unite(res);
			sc.p_unscale(r->pdata.begin(), r->pdata.end());
		}
		sc.p_unscale(pnt1.begin(), pnt1.end());
		sc.p_unscale(pnt2.begin(), pnt2.end());
		return r;
	} catch (const std::exception &e){
		sc.p_unscale(pnt1.begin(), pnt1.end());
		sc.p_unscale(pnt2.begin(), pnt2.end());
		std::cout<<"domain clip error: "<<e.what()<<std::endl;
		return 0;
	}
}

void ecollection_edges_info(void* ecol, int* npts, int* neds, double** pts, int** eds){
	auto c = static_cast<HMCont2D::Container<HMCont2D::ECollection>*>(ecol);
	std::map<Point*, int> pind;
	for (auto p: c->all_points()) pind.emplace(p, pind.size());
	*npts = pind.size();
	*pts = new double[2 * *npts];
	for (auto s: pind){
		(*pts)[2*s.second] = s.first->x;
		(*pts)[2*s.second+1] = s.first->y;
	}
	*neds = c->size();
	*eds = new int[2 * *neds];
	for (int i=0; i<*neds; ++i){
		(*eds)[2*i] = pind[c->data[i]->pstart];
		(*eds)[2*i+1] = pind[c->data[i]->pend];
	}
}
void free_ecollection_edges_info(double* pts, int* eds){
	delete[] pts;
	delete[] eds;
}

void free_ecollection_container(void* ecol){
	delete static_cast<HMCont2D::Container<HMCont2D::ECollection>*>(ecol);
}

int set_ecollection_bc(void* src, void* tar, int def, int* vsrc, int* vtar){
	try{
		auto esrc = static_cast<HMCont2D::ECollection*>(src);
		auto etar = static_cast<HMCont2D::ECollection*>(tar);
		auto pntsrc = esrc->all_points();
		auto pnttar = etar->all_points();

		//scaling
		ScaleBase sc = ScaleBase::p_doscale(pntsrc.begin(), pntsrc.end()); 
		sc.p_scale(pnttar.begin(), pnttar.end());

		//searching
		for (int i=0; i<etar->size(); ++i){
			vtar[i] = def;
			Point cpoint = etar->edge(i)->center();
			auto ce = HMCont2D::ECollection::FindClosestEdge(*esrc, cpoint);
			if (fabs(std::get<1>(ce)) < 1e-3){
				int ied = std::get<3>(ce);
				double m1 = Point::meas_line(*etar->edge(i)->pstart,
						*esrc->edge(ied)->pstart, *esrc->edge(ied)->pend);
				double m2 = Point::meas_line(*etar->edge(i)->pend,
						*esrc->edge(ied)->pstart, *esrc->edge(ied)->pend);
				if (ISZERO(m1) && ISZERO(m2)) vtar[i] = vsrc[ied];
			}
		}

		//unscaling
		sc.p_unscale(pntsrc.begin(), pntsrc.end());
		sc.p_unscale(pnttar.begin(), pnttar.end());
		return 0;
	} catch (const std::exception &e){
		std::cout<<"set_ecollection_bc error: "<<e.what()<<std::endl;
		return 1;
	}
}

int set_ecollection_bc_force(void* src, void* tar, int* vsrc, int* vtar, int algo){
	HMCont2D::ECollection* src_col = static_cast<HMCont2D::ECollection*>(src);
	HMCont2D::ECollection* tar_col = static_cast<HMCont2D::ECollection*>(tar);
	ScaleBase sc = HMCont2D::Scale01(*src_col);
	HMCont2D::Scale(*tar_col, sc);
	double ret = 1;
	try{
		//target edge->source edge
		std::map<HMCont2D::Edge*, HMCont2D::Edge*> ans;
		//set of target edges
		std::set<HMCont2D::Edge*> unused;
		for (auto e: *tar_col) unused.insert(e.get());
		//assemble source contour tree
		auto stree = HMCont2D::Assembler::ETree(*src_col);
		//loop over all tree contours
		for (int i=0; i<stree.cont_count(); ++i){
			auto it = unused.begin();
			while (it!=unused.end()){
				HMCont2D::Edge* e = *it;
				HMCont2D::Contour& cont = *stree.get_contour(i);
				auto fnd1 = HMCont2D::ECollection::FindClosestEdge(cont, *e->pstart);
				auto fnd2 = HMCont2D::ECollection::FindClosestEdge(cont, *e->pend);
				auto fndc = HMCont2D::ECollection::FindClosestEdge(cont, e->center());
				//1) if all 3 edge points lie on contour return bc from center
				if (ISZERO(std::get<1>(fnd1)) && ISZERO(std::get<1>(fnd2)) && ISZERO(std::get<1>(fndc))){
					ans[e] = std::get<0>(fndc);
					it = unused.erase(it);
					continue;
				}
				//2) if start and end lie on contour get boundary of weighted center
				if (algo > 1 && ISZERO(std::get<1>(fnd1)) && ISZERO(std::get<1>(fnd2))){
					HMCont2D::Edge* r;
					if (std::get<0>(fnd1) == std::get<0>(fnd2)) r = std::get<0>(fnd1);
					else {
						double w1 = std::get<1>(cont.coord_at(*e->pstart));
						double w2 = std::get<1>(cont.coord_at(*e->pend));
						if (!cont.is_closed()){
							Point p = HMCont2D::Contour::WeightPoint(cont, (w1+w2)/2.0);
							r = std::get<0>(HMCont2D::ECollection::FindClosestEdge(cont, p));
						} else {
							//for closed contour we have to consider two points and
							//choose one that is closer to edge center
							double var1 = (w1 + w2)/2.0;
							double var2 = (1.0 + w1 + w2)/2.0;
							if (var2>1.0) var2 -= 1.0;
							Point pvar1 = HMCont2D::Contour::WeightPoint(cont, var1);
							Point pvar2 = HMCont2D::Contour::WeightPoint(cont, var2);
							double meas_var1 = Point::meas(pvar1, e->center());
							double meas_var2 = Point::meas(pvar2, e->center());
							if (meas_var1 < meas_var2){
								r = std::get<0>(HMCont2D::ECollection::FindClosestEdge(cont, pvar1));
							} else {
								r = std::get<0>(HMCont2D::ECollection::FindClosestEdge(cont, pvar2));
							}
						}
					}
					ans[e] = r;
					it = unused.erase(it);
					continue;
				}
				++it;
			}
		}
		//3) for all remaining cells get closest edge from source collection
		if (algo>2){
			for (auto e: unused){
				auto ce = HMCont2D::ECollection::FindClosestEdge(*src_col, e->center());
				ans[e] = std::get<0>(ce);
			}
		}
		//fill output arrays
		for (auto a: ans){
			vtar[tar_col->get_index(a.first)] = vsrc[src_col->get_index(a.second)];
		}
		ret = 0;
	} catch (const std::exception &e){
		std::cout<<"set_ecollection_bc_force error: "<<e.what()<<std::endl;
	}
	HMCont2D::Unscale(*src_col, sc);
	HMCont2D::Unscale(*tar_col, sc);
	return ret;
}


double ecollection_area(void* ecol){
	try{
		auto c = static_cast<HMCont2D::ECollection*>(ecol);
		auto tree = HMCont2D::Assembler::ETree(*c);
		return HMCont2D::Area(tree);
	} catch (const std::exception &e){
		std::cout<<"domain area calculation error: "<<e.what()<<std::endl;
		return 0;
	}
}

//void contour_get_info(Cont* c, int* Npnt, int* Neds,
//                double** pts,
//                int** eds){
//        PointsContoursCollection* cont = static_cast<PointsContoursCollection*>(c);
//        *Npnt = cont->n_total_points();
//        *Neds = cont->n_edges();
//        //points
//        *pts = new double[2*(*Npnt)];
//        double* pp = *pts;
//        for (int i=0; i<*Npnt; ++i){
//                auto p = cont->get_point(i);
//                *pp++ = p->x;
//                *pp++ = p->y;
//        }
//        //edges
//        *eds = new int[2*(*Neds)];
//        int *ee = *eds;
//        for (int i=0; i<*Neds; ++i){
//                auto e = cont->get_edge(i);
//                *ee++ = e.first;
//                *ee++ = e.second;
//        }
//}

//void contour_free_info(double** pts, int** eds){
//        delete[] (*pts);
//        delete[] (*eds);
//}

double grid_area(Grid* g){
	return static_cast<GridGeom*>(g)->area();
}

int report_skewness(void* grid, double threshold, double* max_skew, int* max_skew_cell,
		int* bad_cells_num, int* bad_indicies, double* bad_skew){
	try{
		vector<double> sc = GGeom::Info::Skewness(*static_cast<GridGeom*>(grid));
		*max_skew = 0.0;
		*max_skew_cell = -1;
		*bad_cells_num = 0;
		for (int i=0; i<sc.size(); ++i){
			if (sc[i]>*max_skew){
				*max_skew_cell = i;
				*max_skew = sc[i];
			}
			if (sc[i]>=threshold){
				*bad_cells_num = *bad_cells_num + 1;
				*bad_indicies++ = i;
				*bad_skew++ = sc[i];
			}
		}
		return 0;
	} catch (const std::exception &e){
		std::cout<<"Skewness calculation error: "<<e.what()<<std::endl;
		return 1;
	}
}

//double contour_area(Cont* c){
//        return static_cast<PointsContoursCollection*>(c)->area();
//}

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
		hmcport_callback cb_fun){
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

			if (inp.bnd_step_method == HMBlay::BndStepMethod::INCREMENTAL){
				if (inp.start == inp.end){
					throw std::runtime_error("Can not use incremental stepping "
							         "without divergent start/end points");
				}
				inp.bnd_step_basis.push_back(std::make_pair(Point(inp.start), opt.step_start));
				inp.bnd_step_basis.push_back(std::make_pair(Point(inp.end), opt.step_end));
			}
		}
		GridGeom* ret = new GridGeom(HMBlay::BuildBLayerGrid(vinp));
		return ret;
	} catch (const std::exception &e){
		std::cout<<"Boundary layer builder error: "<<e.what()<<std::endl;
		return 0;
	}
}

Grid* build_grid_mapping(void* base_grid, void* target_contour, int Npnt,
		double* pbase, double* ptarget, int snap_method, int algo,
		int reversed,
		int return_invalid, hmcport_callback cb){
	GridGeom* g = static_cast<GridGeom*>(base_grid);
	HMCont2D::ECollection* col = static_cast<HMCont2D::ECollection*>(target_contour);
	//scaling
	ScaleBase bscale = g->do_scale();
	ScaleBase cscale = HMCont2D::ECollection::Scale01(*col);
	Grid* ret = 0;
	try{
		std::vector<Point> p1, p2;
		for (int i=0; i<Npnt; ++i){
			p1.push_back(Point(pbase[2*i], pbase[2*i+1]));
			p2.push_back(Point(ptarget[2*i], ptarget[2*i+1]));
		}
		bscale.scale(p1.begin(), p1.end());
		cscale.scale(p2.begin(), p2.end());
		HMMap::Options opt;
		switch (snap_method){
			case 1: opt.snap = "NO"; break;
			case 2: opt.snap = "ADD_VERTICES"; break;
			case 3: opt.snap = "SHIFT_VERTICES"; break;
		}
		switch (algo){
			case 1: opt.algo = "direct-laplace"; break;
			case 2: opt.algo = "inverse-laplace"; break;
		}
		GridGeom ans = HMMap::MapGrid.WithCallback(cb, *g, *col, p1, p2, (reversed==1), opt);
		ans.undo_scale(cscale);
		ret = new GridGeom(std::move(ans));
	} catch (HMMap::EInvalidGrid &e){
		if (!return_invalid){
			std::cout<<e.what()<<std::endl;
			ret = 0;
		} else{
			std::cout<<"Ignored error: "<<e.what()<<std::endl;
			e.invalid_grid.undo_scale(cscale);
			ret = new GridGeom(std::move(e.invalid_grid));
		}
	} catch (const std::exception &e){
		std::cout<<e.what()<<std::endl;
		ret = 0;
	}
	//unscaling
	g->undo_scale(bscale);
	HMCont2D::ECollection::Unscale(*col, cscale);
	return ret;
}

void* new_writer(){
	try{
		HMXML::Reader* ret = HMXML::ReaderA::pcreate("HybMeshData");
		return ret;
	} catch (const std::exception &e){
		std::cout<<e.what()<<std::endl;
		return 0;
	}
}
void* new_reader(const char* fn){
	try{
		HMXML::Reader* ret = new HMXML::ReaderA(fn, "</HybMeshData>");
		return ret;
	} catch (const std::exception &e){
		std::cout<<e.what()<<std::endl;
		return 0;
	}
}
int finalize_writer(void* writer, const char* fn){
	try{
		HMXML::ReaderA* wr = static_cast<HMXML::ReaderA*>(writer);
		if (fn!=NULL) wr->write(fn);
		delete wr;
		return 1;
	} catch (const std::exception &e){
		std::cout<<e.what()<<std::endl;
		return 0;
	}
}

void free_hmxml_node(void* node){
	HMXML::Reader* n = static_cast<HMXML::Reader*>(node);
	if (n->isroot) n->Free();
	delete n;
}

int hmxml_query(void* node, const char* q, int* num, void** ans){
	try{
		HMXML::Reader* wr = static_cast<HMXML::Reader*>(node);
		vector<HMXML::Reader> fnd = wr->findall_by_path(q);
		*num = fnd.size();
		if (*num > 1000) throw std::runtime_error("query result entries number of too big (>1000)");
		for (int i=0; i<*num; ++i){
			ans[i] = new HMXML::Reader(fnd[i]);
		}
		return 1;
	} catch (const std::exception &e){
		std::cout<<e.what()<<std::endl;
		return 0;
	}
}

int hmxml_change_basenode(void* node, const char* q){
	try{
		HMXML::Reader* wr = static_cast<HMXML::Reader*>(node);
		HMXML::Reader fnd = wr->find_by_path(q, true);
		wr->_nd = fnd._nd;
		return 1;
	} catch (const std::exception &e){
		std::cout<<e.what()<<std::endl;
		return 0;
	}
}

char* hmxml_purged_string(void* node){
	try{
		HMXML::Reader* wr = static_cast<HMXML::Reader*>(node);
		HMXML::Reader copy = HMXML::Reader::copy_to_root(*wr);
		auto unused1 = copy.findall_by_path(".//GRID2D");
		auto unused2 = copy.findall_by_path(".//GRID3D");
		auto unused3 = copy.findall_by_path(".//CONTOUR2D");
		auto unused = unused1;
		std::copy(unused2.begin(), unused2.end(), std::back_inserter(unused));
		std::copy(unused3.begin(), unused3.end(), std::back_inserter(unused));
		for (auto& n: unused) n.unlink_node();
		std::string s = copy.tostring();
		char* out = new char[s.size()+10];
		strcpy(out, s.c_str());
		return out;
	} catch (const std::exception &e){
		std::cout<<e.what()<<std::endl;
		return 0;
	}
}
