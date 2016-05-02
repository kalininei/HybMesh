#include "domapping.hpp"
#include "femassembly.hpp"
#include "hmtimer.hpp"

using namespace HMGMap::Impl;

// ========================== General DoMapping
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
GridGeom DoMapping::run(HMCallback::Caller2& cb){
	//build solution area
	cb.step_after(5, "Boundary mapping");
	prepare_mapped_contour();
	//build a fem triangle grid
	cb.step_after(50, "Triangulation");
	prepare_grid();
	//assemble fem problem
	cb.step_after(15, "FEM assembling");
	laplace.reset(new HMFem::LaplasProblem(g3));

	//solve direct problems
	//resulting u, v are functions defined in g3 at mapped_outer domain
	vector<double> u, v;
	cb.step_after(20, "SLAE solution");
	solve_uv_problems(u, v);

	//draw a resulting grid
	cb.step_after(5, "Grid building");
	auto approx = g3->GetApprox();
	auto pmod = [&](GridPoint* p){
		double x = approx->Val(*p, u);
		double y = approx->Val(*p, v);
		p->set(x, y);
	};
	GridGeom ret = GGeom::Constructor::DeepCopy(inpgrid);
	GGeom::Modify::PointModify(ret, pmod);

	cb.step_after(5, "Snapping");
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
	if (!GGeom::Info::Check(ret)) throw HMGMap::MapException("Resulting grid is not valid");
	return ret;
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
	mapped_outer = HMCont2D::Algos::Simplified(mapped_outer);

	//input grid contour
	inpgrid_outer = HMCont2D::Algos::Simplified(GGeom::Info::Contour(inpgrid));

	//build boundary mapping
	build_mcc();
}
void DoMapping::prepare_grid(){
	build_grid3();
	g3outer = GGeom::Info::Contour(*g3);
}
void DoMapping::build_mcc(){
	//mapping inpgrid_outer to mapped_outer
	int N = std::min(gridpoints.size(), contpoints.size());
	std::set<HMCont2D::Contour*> used_gcontours;
	for (int i=0; i<N; ++i){
		Point p1 = gridpoints[i];
		Point p2 = contpoints[i];
		HMCont2D::Edge* e1 = std::get<0>(HMCont2D::ECollection::FindClosestEdge(inpgrid_outer, p1));
		HMCont2D::Edge* e2 = std::get<0>(HMCont2D::ECollection::FindClosestEdge(mapped_outer, p2));
		HMCont2D::Contour* c1 = inpgrid_outer.get_contour(e1);
		HMCont2D::Contour* c2 = mapped_outer.get_contour(e2);
		used_gcontours.insert(c1);
		auto cfnd = mcol.insert(c1, c2);
		cfnd->add_connection(p1, p2);
	}
	if (inpgrid_outer.cont_count() != used_gcontours.size())
		throw HMGMap::MapException("All grid boundaries should contain at least one base point");
}

// =========================== DirectMapping
void DirectMapping::build_grid3(){
	g3.reset(new HMFem::Grid43(
		HMFem::AuxGrid3(inpgrid_outer, opt.fem_nrec, opt.fem_nmax)));
}
void DirectMapping::solve_uv_problems(vector<double>& u, vector<double>& v){
	u.resize(g3->n_points(), 0.0);
	v.resize(g3->n_points(), 0.0);

	// u - problem
	laplace->ClearBC();
	for (auto n: g3outer.nodes){
		HMCont2D::Edge *e = std::get<0>(HMCont2D::ECollection::FindClosestEdge(inpgrid_outer, *n->first()));
		HMCont2D::Contour* c = inpgrid_outer.get_contour(e);
		auto cmapping = mcol.find_by_base(c);
		assert(cmapping != nullptr);
		auto dirfunc = [cmapping](const GridPoint* p)->double{
			return cmapping->map_from_base(*p).x;
		};
		laplace->SetDirichlet(*n, dirfunc);
	}
	laplace->Solve(u);

	// v - problem
	laplace->ClearBC();
	for (auto n: g3outer.nodes){
		HMCont2D::Edge *e = std::get<0>(HMCont2D::ECollection::FindClosestEdge(inpgrid_outer, *n->first()));
		HMCont2D::Contour* c = inpgrid_outer.get_contour(e);
		auto cmapping = mcol.find_by_base(c);
		assert(cmapping != nullptr);
		auto dirfunc = [cmapping](const GridPoint* p)->double{
			return cmapping->map_from_base(*p).y;
		};
		laplace->SetDirichlet(*n, dirfunc);
	}
	laplace->QuickSolve_BC(v);
}

// =========================== InverseMapping
void InverseMapping::build_grid3(){
	//add corner points from inpgrid to mapped_outer
	HMCont2D::ContourTree mapped2;
	HMCont2D::ContourTree::DeepCopy(mapped_outer, mapped2);
	HMCont2D::PCollection pcol;
	for (int i=0; i<mcol.entry_num(); ++i){
		//find mapped contour
		auto cmapping = mcol.get(i);
		auto mapped_contour = cmapping->get_mapped();
		HMCont2D::Contour* copied_mapped_contour=0;
		for (int j=0; j<mapped_outer.nodes.size(); ++j){
			if (mapped_outer.nodes[j].get() == mapped_contour){
				copied_mapped_contour = mapped2.nodes[j].get();
				break;
			}
		}
		//for each point in contour
		for (auto p: cmapping->get_base()->all_points()){
			//find mapped point
			Point pmapped = cmapping->map_from_base(*p);
			copied_mapped_contour->GuaranteePoint(pmapped, pcol);
		}
	}
	
	//build grid
	g3.reset(new HMFem::Grid43(
		HMFem::AuxGrid3(mapped2, opt.fem_nrec, opt.fem_nmax)));
}

void InverseMapping::solve_uv_problems(vector<double>& u, vector<double>& v){
	u.resize(g3->n_points(), 0.0);
	v.resize(g3->n_points(), 0.0);

	//u problem
	for (auto n: g3outer.nodes){
		HMCont2D::Edge *e = std::get<0>(HMCont2D::ECollection::FindClosestEdge(mapped_outer, *n->first()));
		HMCont2D::Contour* c = mapped_outer.get_contour(e);
		auto cmapping = mcol.find_by_mapped(c);
		assert(cmapping != 0);
		auto dirfunc = [cmapping](const GridPoint* p)->double{
			return cmapping->map_from_mapped(*p).x;
		};
		laplace->SetDirichlet(*n, dirfunc);
	}
	laplace->Solve(u);

	//v problem
	laplace->ClearBC();
	for (auto n: g3outer.nodes){
		HMCont2D::Edge *e = std::get<0>(HMCont2D::ECollection::FindClosestEdge(mapped_outer, *n->first()));
		HMCont2D::Contour* c = mapped_outer.get_contour(e);
		auto cmapping = mcol.find_by_mapped(c);
		assert(cmapping != 0);
		auto dirfunc = [cmapping](const GridPoint* p)->double{
			return cmapping->map_from_mapped(*p).y;
		};
		laplace->SetDirichlet(*n, dirfunc);
	}
	laplace->QuickSolve_BC(v);

	//swap g3 coordinates and uv
	auto modfun = [&u, &v](GridPoint* p){
		std::swap(p->x, u[p->get_ind()]);
		std::swap(p->y, v[p->get_ind()]);
	};
	GGeom::Modify::PointModify(*g3, modfun);
}
