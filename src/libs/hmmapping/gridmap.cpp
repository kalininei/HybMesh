#include "gridmap.hpp"
#include "femassembly.hpp"
#include "fileproc.h"
using namespace HMGMap::Impl;

void DoMapping::set_grid(const GridGeom& ig){
	GGeom::Modify::ClearAll(inpgrid);
	GGeom::Modify::ShallowAdd(&ig, &inpgrid);
}

void DoMapping::set_contour(const HMCont2D::ECollection& ecol){
	contdata = ecol;
}

void DoMapping::set_points(const vector<Point>& gridpnt, const vector<Point>& contpnt){
	gridpoints = gridpnt;
	contpoints = contpnt;
}

void DoMapping::prepare_mapped_contour(){
	//assemble all contours
	vector<HMCont2D::Contour> allconts = HMCont2D::Assembler::AllContours(contdata);
	//get contours which contains given points
	std::set<HMCont2D::Contour*> included_conts;
	for (auto p: contpoints){
		auto eres = HMCont2D::ECollection::FindClosestEdge(contdata, p);
		HMCont2D::Edge* e = std::get<0>(eres);
		for (auto& c: allconts){
			if (c.contains(e)){
				included_conts.insert(&c);
				break;
			}
		}
	}
	//assemble a tree from that contours
	mapped_outer.clear();
	for (auto c: included_conts){
		assert(c->is_closed());
		mapped_outer.AddContour(*c);
	}
}

void DoMapping::build_grid3(){
	//building a grid of opt.fem_nrec nodes
	//* preliminary data
	HMCont2D::ContourTree ct = GGeom::Info::Contour(inpgrid);
	ct = HMCont2D::Algos::Simplified(ct);
	double area = HMCont2D::Area(ct);
	auto elen = HMCont2D::ECollection::ELengths(ct);
	double longest_edge = *max_element(elen.begin(), elen.end());
	double shortest_edge = *min_element(elen.begin(), elen.end());
	{//* first try: uniform triangulation
		double l1=sqrt(area/opt.fem_nrec);
		if (l1 <= shortest_edge/opt.fem_nedge){
			g3 = HMFem::Grid43::Build3(ct, {}, l1);
			if (g3 != 0 && g3->n_points()<opt.fem_nmax) return;
		}
	}
	{//* TODO
	}
	throw MapException("failed to triangulate given area");
}

void DoMapping::build_mcc(){
	//mapping g3outer to mapped_outer
	int N = std::min(gridpoints.size(), contpoints.size());
	std::set<HMCont2D::Contour*> used_gcontours;
	for (int i=0; i<N; ++i){
		Point p1 = gridpoints[i];
		Point p2 = contpoints[i];
		HMCont2D::Edge* e1 = std::get<0>(HMCont2D::ECollection::FindClosestEdge(g3outer, p1));
		HMCont2D::Edge* e2 = std::get<0>(HMCont2D::ECollection::FindClosestEdge(mapped_outer, p2));
		HMCont2D::Contour* c1 = g3outer.get_contour(e1);
		HMCont2D::Contour* c2 = mapped_outer.get_contour(e2);
		used_gcontours.insert(c1);
		auto cfnd = mcol.insert(c1, c2);
		cfnd->add_connection(p1, p2);
	}
	if (g3outer.cont_count() != used_gcontours.size())
		throw HMGMap::MapException("All grid boundaries should contain at least one base point");
}

void DoMapping::prepare_grid(){
	build_grid3();
	g3outer = GGeom::Info::Contour(*g3);
	laplas_mat = HMFem::Assemble::PureLaplas(*g3);
	build_mcc();
}

vector<double> DoMapping::solve_u_problem(){
	auto laplas = HMFem::LaplasProblem(g3, laplas_mat);
	for (auto n: g3outer.nodes){
		HMCont2D::Contour* c = n.get();
		auto cmapping = mcol.find_by_base(c);
		auto dirfunc = [cmapping](const GridPoint* p)->double{
			return cmapping->map_from_base(*p).x;
		};
		laplas.SetDirichlet(*c, dirfunc);
	}
	vector<double> ret(g3->n_points(), 0.0);
	laplas.Solve(ret);
	return ret;
}

vector<double> DoMapping::solve_v_problem(){
	auto laplas = HMFem::LaplasProblem(g3, laplas_mat);
	for (auto n: g3outer.nodes){
		HMCont2D::Contour* c = n.get();
		auto cmapping = mcol.find_by_base(c);
		auto dirfunc = [cmapping](const GridPoint* p)->double{
			return cmapping->map_from_base(*p).y;
		};
		laplas.SetDirichlet(*c, dirfunc);
	}
	vector<double> ret(g3->n_points(), 0.0);
	laplas.Solve(ret);
	return ret;
}

GridGeom DoMapping::run(){
	//build solution area
	prepare_mapped_contour();
	//build a fem triangle grid
	prepare_grid();

	//solve problems
	vector<double> u = solve_u_problem();
	vector<double> v = solve_v_problem();
	//draw a resulting grid
	auto approx = g3->GetApprox();
	auto pmod = [&](GridPoint* p){
		double x = approx->Val(*p, u);
		double y = approx->Val(*p, v);
		p->set(x, y);
	};
	GridGeom ret = GGeom::Constructor::DeepCopy(inpgrid);
	GGeom::Modify::PointModify(ret, pmod);
	//snapping 
	for (auto p: ret.get_bnd_points()){
		Point p2 = HMCont2D::ECollection::ClosestPoint(mapped_outer, *p);
		const_cast<GridPoint*>(p)->set(p2.x, p2.y);
	}
	if (opt.snap == "ADD_VERTICES"){
		for (auto& n: mapped_outer.nodes)
			GGeom::Modify::SnapToContour(ret, *n, {});
	} else if (opt.snap == "SHIFT_VERTICES"){
		for (auto& n: mapped_outer.nodes)
			GGeom::Modify::ShiftToContour(ret, *n, {});
	} else if (opt.snap != "NO"){
		throw HMGMap::MapException(std::string("Unknown snapping option - ") + opt.snap);
	}
#ifndef NDEBUG
	if (!GGeom::Info::Check(ret)) save_vtk(ret, "_badmapped.vtk");
#endif
	if (!GGeom::Info::Check(ret)) throw HMGMap::MapException("Resulting grid is not valid");
	return ret;
}
