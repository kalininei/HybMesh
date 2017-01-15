#include "hmcport.h"
#include "hmblay.hpp"
#include "hmmapping.hpp"
#include "hmtesting.hpp"
#include "hmcallback.hpp"
#include "hmxmlreader.hpp"
#include "contclipping.hpp"
#include "buildgrid.hpp"
#include "contabs2d.hpp"
#include "unite_grids.hpp"
#include "modgrid.hpp"
#include "infogrid.hpp"
#include "finder2d.hpp"

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
		return new HM2D::GridData(std::move(
			HM2D::Grid::Constructor::FromRaw(Npts, Ncells, pts, cells, -1)));
	} catch (const std::exception &e){
		std::cout<<"crossgrid error: "<<e.what()<<std::endl;
		return NULL;
	}
}

int grid_npoints(Grid* g){
	return static_cast<HM2D::GridData*>(g)->vvert.size();
}

int grid_ncells(Grid* g){
	return static_cast<HM2D::GridData*>(g)->vcells.size();
}
int grid_cellsdim(Grid* grid){
	int ret = 0;
	auto g = static_cast<HM2D::GridData*>(grid);
	for (int i=0; i<g->vcells.size(); ++i) ret += g->vcells[i]->edges.size();
	return ret;
}

void grid_get_points_cells(Grid* g, double* pts, int* cells){
	auto gg=static_cast<HM2D::GridData*>(g);
	//points
	if (pts!=0) for (int i=0; i<gg->vvert.size(); ++i){
		auto p = gg->vvert[i];
		*pts++ = p->x;
		*pts++ = p->y;
	}
	//cells
	aa::enumerate_ids_pvec(gg->vvert);
	if (cells!=0) for (int i=0; i<gg->vcells.size(); ++i){
		auto c = gg->vcells[i];
		*cells++ = c->edges.size();
		for (auto p: HM2D::Contour::OrderedPoints(c->edges)){
			*cells++ = p->id;
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
	HM2D::GridData* g = static_cast<HM2D::GridData*>(grd);
	//fill counts
	*Npnt = g->vvert.size();
	*Ncls = g->vcells.size();
	*Neds = g->vedges.size();
	//allocate arrays
	*pts = new double[2 * (*Npnt)];
	*ed_pt = new int[2 * (*Neds)];
	*cls_dims = new int[*Ncls];
	*cls_eds = new int[grid_cellsdim(grd)];
	//fill points
	double* p_pts = *pts;
	for (int i=0; i<g->vvert.size(); ++i){
		*p_pts++ = g->vvert[i]->x;
		*p_pts++ = g->vvert[i]->y;
	}
	//fill edges
	std::map<std::pair<int, int>, int> nd_eds;
	int* p_eds = *ed_pt;
	aa::enumerate_ids_pvec(g->vvert);
	for (size_t i=0; i<g->vedges.size(); ++i){
		auto& e = g->vedges[i];
		*p_eds++ = e->first()->id;
		*p_eds++ = e->last()->id;
		nd_eds.emplace(std::make_pair(*(p_eds-2), *(p_eds-1)), i);
	}
	//fill elements dimensions
	int* p_cdim = *cls_dims;
	for (int i=0; i<g->vcells.size(); ++i){
		*p_cdim++ = g->vcells[i]->edges.size();
	}
	//fill cell->edges array
	aa::enumerate_ids_pvec(g->vedges);
	int* p_ced = *cls_eds;
	for (auto c: g->vcells)
	for (auto e: c->edges){
		*p_ced++ = e->id;
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
	auto gg = static_cast<HM2D::GridData*>(g);
	gg->enumerate_all();
	auto& ed = gg->vedges;
	*Neds = ed.size();
	*ed_cell = new int[*Neds];
	int *ec = *ed_cell;
	int ret = 1;
	if (ed_pt == 0){
		for (auto& e: ed){
			if (e->has_left_cell()) *ec++ = e->left.lock()->id;
			else *ec++ = -1;
			if (e->has_right_cell()) *ec++ = e->right.lock()->id;
			else *ec++ = -1;
		}
	} else {
		auto vertedge = HM2D::Connectivity::VertexEdge(gg->vedges);
		int it=0;
		for (auto& ve: vertedge) ve.v->id = it++;
		for (int i=0; i<*Neds; ++i){
			int p1 = gg->vvert[ed_pt[2*i]]->id, p2 = gg->vvert[ed_pt[2*i+1]]->id;
			HM2D::Edge* fnd = 0;
			for (auto e: vertedge[p1].eind){
				if (gg->vedges[e]->first()->id == p1 && gg->vedges[e]->last()->id == p2){
					fnd = gg->vedges[e].get(); break;
				}
				if (gg->vedges[e]->last()->id == p1 && gg->vedges[e]->first()->id == p2){
					fnd = gg->vedges[e].get(); break;
				}
			}
			if (fnd != 0){
				if (p1 == fnd->first()->id){
					if (fnd->has_left_cell()) *ec++ = fnd->left.lock()->id;
					else *ec++ = -1;
					if (fnd->has_right_cell()) *ec++ = fnd->right.lock()->id;
					else *ec++ = -1;
				} else {
					if (fnd->has_right_cell()) *ec++ = fnd->right.lock()->id;
					else *ec++ = -1;
					if (fnd->has_left_cell()) *ec++ = fnd->left.lock()->id;
					else *ec++ = -1;
				}
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
	if (g!=NULL) delete static_cast<HM2D::GridData*>(g);
}

Grid* cross_grids(Grid* gbase, Grid* gsecondary, double buffer_size, int preserve_bp, int eh, double angle0, int algo){
	return cross_grids_wcb(gbase, gsecondary, buffer_size, preserve_bp, eh, angle0, algo, silent2_function);
}

Grid* cross_grids_wcb(Grid* gbase, Grid* gsecondary, double buffer_size,
		int preserve_bp, int empty_holes, double angle0, int algo, hmcport_callback cb_fun){
	HM2D::GridData* g0 = static_cast<HM2D::GridData*>(gbase);
	HM2D::GridData* g1 = static_cast<HM2D::GridData*>(gsecondary);
	ScaleBase sc = HM2D::Scale01(g0->vvert, 1.0 + sqrt(2.0)/100.0 + sqrt(3.0)/1000.0);
	HM2D::Scale(g1->vvert, sc);
	HM2D::GridData* ret;
	try{
		HM2D::Grid::Algos::OptUnite opt;
		opt.buffer_size = buffer_size/sc.L;
		opt.preserve_bp = (bool)preserve_bp;
		opt.empty_holes = (bool)empty_holes;
		opt.angle0 = angle0;
		opt.filler = algo;
		auto gret = HM2D::Grid::Algos::UniteGrids.WithCallback(cb_fun, *g0, *g1, opt);
		HM2D::Unscale(gret.vvert, sc);
		ret = new HM2D::GridData(std::move(gret));
	} catch (const std::exception &e) {
		std::cout<<"crossgrid error: "<<e.what()<<std::endl;
		ret = 0;
	}

	HM2D::Unscale(g0->vvert, sc);
	HM2D::Unscale(g1->vvert, sc);
	return ret;
}

Grid* grid_exclude_cont(Grid* grd, void* cont, int is_inner){
	return grid_exclude_cont_wcb(grd, cont, is_inner, silent2_function);
}

Grid* grid_exclude_cont_wcb(Grid* grd, void* cont, int is_inner,
		hmcport_callback cb_fun){
	HM2D::GridData* g0 = static_cast<HM2D::GridData*>(grd);
	HM2D::EdgeData* e0 = static_cast<HM2D::EdgeData*>(cont);
	//non g0 vertices
	HM2D::VertexData e0nog0;
	aa::constant_ids_pvec(HM2D::AllVertices(*e0), 1);
	aa::constant_ids_pvec(g0->vvert, 0);
	for (auto v: HM2D::AllVertices(*e0)) if (v->id == 1) e0nog0.push_back(v);
	//scaling
	ScaleBase sc = HM2D::Scale01(g0->vvert, 1.0 + sqrt(2.0)/100.0 + sqrt(3.0)/1000.0);
	HM2D::Scale(e0nog0, sc);
	HM2D::GridData* ret;
	try{
		auto tree = HM2D::Contour::Tree::Assemble(*e0);
		HM2D::GridData gret = HM2D::Grid::Algos::SubstractArea.WithCallback(
			cb_fun,
			*g0, tree, is_inner!=0);
		HM2D::Unscale(gret.vvert, sc);
		ret = new HM2D::GridData(std::move(gret));
	} catch (const std::exception &e){
		std::cout<<"grid_exclude_cont error: "<<e.what()<<std::endl;
		ret = 0;
	}
	//unscaling
	HM2D::Unscale(g0->vvert, sc);
	HM2D::Unscale(e0nog0, sc);
	return ret;
}

int simplify_grid_boundary(Grid* grd, double angle){
	HM2D::GridData* g = static_cast<HM2D::GridData*>(grd);
	ScaleBase sc = HM2D::Scale01(g->vvert);
	int ret = 1;
	try{
		if (angle<-geps || angle>180+geps) throw std::runtime_error("invalid angle");
		HM2D::Grid::Algos::SimplifyBoundary(*g, angle);
		ret = 0;
	} catch (const std::exception &e){
		std::cout<<"simplify_grid_boundary error: "<<e.what()<<std::endl;
		ret = 1;
	}
	HM2D::Unscale(g->vvert, sc);
	return ret;
}

void* create_ecollection_container(int Npts, double* pts, int Nedgs, int* edges){
	try{
		auto ret = new HM2D::EdgeData;
		HM2D::VertexData vd;
		for (int i=0; i<Npts; ++i){
			vd.push_back(std::make_shared<HM2D::Vertex>(pts[2*i], pts[2*i+1]));
		}
		//edges
		for (int i=0; i<Nedgs; ++i){
			ret->push_back(std::make_shared<HM2D::Edge>(
				vd[edges[2*i]], vd[edges[2*i+1]]));
		}
		return ret;
	} catch (const std::exception &e){
		std::cout<<"hybmesh_contours2d error: "<<e.what()<<std::endl;
		return 0;
	}
}

void* domain_clip(void* c1, void* c2, int oper, int simplify){
	auto cont1 = static_cast<HM2D::EdgeData*>(c1);
	auto cont2 = static_cast<HM2D::EdgeData*>(c2);
	auto pnt1 = HM2D::AllVertices(*cont1);
	auto pnt2 = HM2D::AllVertices(*cont2);

	ScaleBase sc = ScaleBase::p_doscale(pnt1.begin(), pnt1.end()); 
	sc.p_scale(pnt2.begin(), pnt2.end());

	auto tree1 = HM2D::Contour::Tree::Assemble(*cont1);
	auto tree2 = HM2D::Contour::Tree::Assemble(*cont2);
	try{
		if (tree1.nodes.size() == 0) throw std::runtime_error("not a closed contour");
		if (tree2.nodes.size() == 0) throw std::runtime_error("not a closed contour");

		HM2D::Contour::Tree res;
		switch (oper){
			case 1: res = HM2D::Contour::Clip::Union(tree1, tree2); break;
			case 2: res = HM2D::Contour::Clip::Difference(tree1, tree2); break;
			case 3: res = HM2D::Contour::Clip::Intersection(tree1, tree2); break;
			case 4: res = HM2D::Contour::Clip::XOR(tree1, tree2); break;
			default: throw std::runtime_error("unknown operation");
		}
		HM2D::Contour::Clip::Heal(res);
		if (res.nodes.size() !=0 && !simplify){
			HM2D::VertexData allpnt(pnt1);
			allpnt.insert(allpnt.end(), pnt2.begin(), pnt2.end());
			for (auto p: allpnt){
				auto fnd = HM2D::Finder::ClosestEdge(res.alledges(), *p);
				if (std::get<0>(fnd)>=0 && std::get<1>(fnd)<geps && 
						std::get<2>(fnd) > geps && std::get<2>(fnd) < 1-geps){
					auto e = res.alledges()[std::get<0>(fnd)].get();
					auto cont = res.find_node(e);
					HM2D::Contour::GuaranteePoint(cont->contour, *p);
				}
			}
		}
		HM2D::EdgeData* r = 0;
		if (res.nodes.size() > 0){
			r = new HM2D::EdgeData(res.alledges());
			auto av = HM2D::AllVertices(*r);
			sc.p_unscale(av.begin(), av.end());
		}
		sc.p_unscale(pnt1.begin(), pnt1.end());
		sc.p_unscale(pnt2.begin(), pnt2.end());
		return r;
	} catch (const std::exception &e){
		sc.p_unscale(pnt1.begin(), pnt1.end());
		sc.p_unscale(pnt2.begin(), pnt2.end());
		std::cout<<e.what()<<std::endl;
		return 0;
	}
}

void ecollection_edges_info(void* ecol, int* npts, int* neds, double** pts, int** eds){
	auto c = static_cast<HM2D::EdgeData*>(ecol);
	std::map<Point*, int> pind;
	for (auto p: HM2D::AllVertices(*c)) pind.emplace(p.get(), pind.size());
	*npts = pind.size();
	*pts = new double[2 * *npts];
	for (auto s: pind){
		(*pts)[2*s.second] = s.first->x;
		(*pts)[2*s.second+1] = s.first->y;
	}
	*neds = c->size();
	*eds = new int[2 * *neds];
	for (int i=0; i<*neds; ++i){
		(*eds)[2*i] = pind[(*c)[i]->first().get()];
		(*eds)[2*i+1] = pind[(*c)[i]->last().get()];
	}
}
void free_ecollection_edges_info(double* pts, int* eds){
	delete[] pts;
	delete[] eds;
}

void free_ecollection_container(void* ecol){
	delete static_cast<HM2D::EdgeData*>(ecol);
}

int set_ecollection_bc(void* src, void* tar, int def, int* vsrc, int* vtar){
	try{
		auto esrc = static_cast<HM2D::EdgeData*>(src);
		auto etar = static_cast<HM2D::EdgeData*>(tar);
		auto pntsrc = HM2D::AllVertices(*esrc);
		auto pnttar = HM2D::AllVertices(*etar);

		//scaling
		ScaleBase sc = ScaleBase::p_doscale(pntsrc.begin(), pntsrc.end()); 
		sc.p_scale(pnttar.begin(), pnttar.end());

		//searching
		for (int i=0; i<etar->size(); ++i){
			vtar[i] = def;
			Point cpoint = etar->at(i)->center();
			auto ce = HM2D::Finder::ClosestEdge(*esrc, cpoint);
			if (fabs(std::get<1>(ce)) < 1e-3){
				int ied = std::get<0>(ce);
				double m1 = Point::meas_line(*etar->at(i)->first(),
						*esrc->at(ied)->first(), *esrc->at(ied)->last());
				double m2 = Point::meas_line(*etar->at(i)->last(),
						*esrc->at(ied)->first(), *esrc->at(ied)->last());
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
	HM2D::EdgeData* src_col = static_cast<HM2D::EdgeData*>(src);
	HM2D::EdgeData* tar_col = static_cast<HM2D::EdgeData*>(tar);
	ScaleBase sc = HM2D::Scale01(*src_col);
	HM2D::Scale(*tar_col, sc);
	double ret = 1;
	try{
		//target edge->source edge
		std::map<HM2D::Edge*, HM2D::Edge*> ans;
		//set of target edges
		std::set<HM2D::Edge*> unused;
		for (auto e: *tar_col) unused.insert(e.get());
		//assemble source contour tree
		auto stree = HM2D::Contour::Tree::Assemble(*src_col);
		//loop over all tree contours
		for (int i=0; i<stree.nodes.size(); ++i){
			auto it = unused.begin();
			while (it!=unused.end()){
				HM2D::Edge* e = *it;
				HM2D::EdgeData& cont = stree.nodes[i]->contour;
				auto fnd1 = HM2D::Finder::ClosestEdge(cont, *e->first());
				auto fnd2 = HM2D::Finder::ClosestEdge(cont, *e->last());
				auto fndc = HM2D::Finder::ClosestEdge(cont, e->center());
				//1) if all 3 edge points lie on contour return bc from center
				if (ISZERO(std::get<1>(fnd1)) && ISZERO(std::get<1>(fnd2)) && ISZERO(std::get<1>(fndc))){
					ans[e] = cont[std::get<0>(fndc)].get();
					it = unused.erase(it);
					continue;
				}
				//2) if start and end lie on contour get boundary of weighted center
				if (algo > 1 && ISZERO(std::get<1>(fnd1)) && ISZERO(std::get<1>(fnd2))){
					HM2D::Edge* r;
					if (std::get<0>(fnd1) == std::get<0>(fnd2)) r = cont[std::get<0>(fnd1)].get();
					else {
						double w1 = std::get<1>(HM2D::Contour::CoordAt(cont, *e->first()));
						double w2 = std::get<1>(HM2D::Contour::CoordAt(cont, *e->last()));
						if (HM2D::Contour::IsOpen(cont)){
							Point p = HM2D::Contour::WeightPoint(cont, (w1+w2)/2.0);
							r = cont[std::get<0>(HM2D::Finder::ClosestEdge(cont, p))].get();
						} else {
							//for closed contour we have to consider two points and
							//choose one that is closer to edge center
							double var1 = (w1 + w2)/2.0;
							double var2 = (1.0 + w1 + w2)/2.0;
							if (var2>1.0) var2 -= 1.0;
							Point pvar1 = HM2D::Contour::WeightPoint(cont, var1);
							Point pvar2 = HM2D::Contour::WeightPoint(cont, var2);
							double meas_var1 = Point::meas(pvar1, e->center());
							double meas_var2 = Point::meas(pvar2, e->center());
							if (meas_var1 < meas_var2){
								r = cont[std::get<0>(HM2D::Finder::ClosestEdge(cont, pvar1))].get();
							} else {
								r = cont[std::get<0>(HM2D::Finder::ClosestEdge(cont, pvar2))].get();
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
				auto ce = HM2D::Finder::ClosestEdge(*src_col, e->center());
				ans[e] = (*src_col)[std::get<0>(ce)].get();
			}
		}
		//fill output arrays
		aa::constant_ids_pvec(*tar_col, -992997);
		for (int i=0; i<src_col->size(); ++i){
			(*src_col)[i]->id = vsrc[i];
		}
		for (auto a: ans){
			if (a.first->id == -992997) a.first->id = a.second->id;
		}
		for (int i=0; i<tar_col->size(); ++i){
			int b = (*tar_col)[i]->id;
			if (b != -992997) vtar[i] = b;
		}
		ret = 0;
	} catch (const std::exception &e){
		std::cout<<"set_ecollection_bc_force error: "<<e.what()<<std::endl;
	}
	HM2D::Unscale(*src_col, sc);
	HM2D::Unscale(*tar_col, sc);
	return ret;
}


double ecollection_area(void* ecol){
	try{
		auto c = static_cast<HM2D::EdgeData*>(ecol);
		auto tree = HM2D::Contour::Tree::Assemble(*c);
		return tree.area();
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
	return HM2D::Grid::Area(*static_cast<HM2D::GridData*>(g));
}

int report_skewness(void* grid, double threshold, double* max_skew, int* max_skew_cell,
		int* bad_cells_num, int* bad_indicies, double* bad_skew){
	try{
		vector<double> sc = HM2D::Grid::Skewness(*static_cast<HM2D::GridData*>(grid));
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
			inp.edges = static_cast<HM2D::EdgeData*>(opt.cont);
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
		HM2D::GridData* ret = new HM2D::GridData(HMBlay::BuildBLayerGrid(vinp));
		return ret;
	} catch (const std::exception &e){
		std::cout<<"Boundary layer builder error: "<<e.what()<<std::endl;
		return 0;
	}
}

Grid* build_grid_mapping(void* base_grid, void* target_contour,
		int Npnt, double* pbase, double* ptarget,
		int snap_method, int algo, int reversed,
		int return_invalid, hmcport_callback cb){
	HM2D::GridData* g = static_cast<HM2D::GridData*>(base_grid);
	HM2D::EdgeData* col = static_cast<HM2D::EdgeData*>(target_contour);
	//scaling
	ScaleBase bscale = HM2D::Scale01(g->vvert);
	ScaleBase cscale = HM2D::Scale01(*col);
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
		HM2D::GridData ans = HMMap::MapGrid.WithCallback(cb, *g, *col, p1, p2, (reversed==1), opt);
		HM2D::Unscale(ans.vvert, cscale);
		ret = new HM2D::GridData(std::move(ans));
	} catch (HMMap::EInvalidGrid &e){
		if (!return_invalid){
			std::cout<<e.what()<<std::endl;
			ret = 0;
		} else{
			std::cout<<"Ignored error: "<<e.what()<<std::endl;
			HM2D::Unscale(e.invalid_grid.vvert, cscale);
			ret = new HM2D::GridData(std::move(e.invalid_grid));
		}
	} catch (const std::exception &e){
		std::cout<<e.what()<<std::endl;
		ret = 0;
	}
	//unscaling
	HM2D::Unscale(g->vvert, bscale);
	HM2D::Unscale(*col, cscale);
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
		auto unused4 = copy.findall_by_path(".//SURFACE3D");
		auto unused = unused1;
		std::copy(unused2.begin(), unused2.end(), std::back_inserter(unused));
		std::copy(unused3.begin(), unused3.end(), std::back_inserter(unused));
		std::copy(unused4.begin(), unused4.end(), std::back_inserter(unused));
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
