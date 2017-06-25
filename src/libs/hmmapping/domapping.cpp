#include "domapping.hpp"
#include "femassembly.hpp"
#include "hmtimer.hpp"
#include "assemble2d.hpp"
#include "modcont.hpp"
#include "treverter2d.hpp"
#include "healgrid.hpp"
#include "finder2d.hpp"
#include "snap_grid2cont.hpp"

using namespace HMMap::Impl;

// ========================== General DoMapping
void DoMapping::set_grid(const HM2D::GridData& ig){
	inpgrid = ig;
}
void DoMapping::set_contour(const HM2D::EdgeData& ecol){
	contdata = ecol;
}
void DoMapping::set_points(const vector<Point>& gridpnt, const vector<Point>& contpnt){
	gridpoints = gridpnt;
	contpoints = contpnt;
}
HM2D::GridData DoMapping::run(HMCallback::Caller2& cb){
	//build solution area
	cb.step_after(5, "Boundary mapping");
	prepare_mapped_contour();
	//build a fem triangle grid
	cb.step_after(45, "Triangulation");
	prepare_grid();
	//assemble fem problem
	cb.step_after(15, "FEM assembling");
	laplace.reset(new HMFem::LaplaceProblem(*g3));

	//solve direct problems
	//resulting u, v are functions defined in g3 at mapped_outer domain
	vector<double> u, v;
	cb.step_after(20, "SLAE solution");
	solve_uv_problems(u, v);

	//construct a resulting grid
	cb.step_after(5, "Grid building");
	HM2D::GridData ret;
	HM2D::DeepCopy(inpgrid, ret);
	//receive boundary points
	vector<bool> is_boundary(ret.vvert.size(), false);
	auto bp = AllVertices(HM2D::ECol::Assembler::GridBoundary(ret));
	aa::enumerate_ids_pvec(ret.vvert);
	for (auto p: bp) is_boundary[p->id] = true;
	auto approx = std::make_shared<HMFem::Grid43::Approximator>(g3.get());
	//do mappings
	for (int i=0; i<ret.vvert.size(); ++i){
		HM2D::Vertex* p = ret.vvert[i].get();
		if (is_boundary[i]){
			p->set(mcol.map_from_base(*p));
		} else {
			double x = approx->Val(*p, u);
			double y = approx->Val(*p, v);
			p->set(x, y);
		}
	}

	cb.step_after(5, "Check grid");
	if (mcol.is_reversed()){
		for (auto e: ret.vedges) std::swap(e->left, e->right);
		for (auto c: ret.vcells) HM2D::Contour::Algos::Reverse(c->edges);
	}
	if (!HM2D::Grid::Algos::Check(ret)) throw HMMap::EInvalidGrid(std::move(ret));

	cb.step_after(5, "Snapping");
	//snapping 
	if (opt.snap == "ADD_VERTICES"){
		for (auto& n: mapped_outer.nodes)
			HM2D::Grid::Algos::SnapToContour(ret, n->contour, {});
	} else if (opt.snap == "SHIFT_VERTICES"){
		for (auto& n: mapped_outer.nodes)
			HM2D::Grid::Algos::ShiftToContour(ret, n->contour, {});
	} else if (opt.snap != "NO"){
		throw HMMap::MapException(std::string("Unknown snapping option - ") + opt.snap);
	}

	//boundary types: now they are inherited from input grid as
	//a result of deepcopy procedure. if contour btype is needed use force assignment
	if (opt.btypes_from_contour){
		auto rc = HM2D::ECol::Assembler::GridBoundary(ret);
		HM2D::ECol::Algos::AssignBTypes(contdata, rc);
	}
	
	return ret;
}
void DoMapping::prepare_mapped_contour(){
	//assemble all contours
	vector<HM2D::EdgeData> allconts = HM2D::Contour::Assembler::AllContours(contdata);
	//get contours which contains given points
	std::set<HM2D::EdgeData*> included_conts;
	for (auto p: contpoints){
		auto eres = HM2D::Finder::ClosestEdge(contdata, p);
		HM2D::Edge* e = contdata[std::get<0>(eres)].get();
		for (auto& c: allconts){
			if (HM2D::Finder::Contains(c, e)){
				included_conts.insert(&c);
				break;
			}
		}
	}
	//assemble a tree from that contours
	mapped_outer.nodes.clear();
	for (auto c: included_conts){
		assert(HM2D::Contour::IsClosed(*c));
		mapped_outer.add_contour(*c);
	}
	mapped_outer = HM2D::Contour::Algos::Simplified(mapped_outer);

	//input grid contour
	inpgrid_outer = HM2D::Contour::Algos::Simplified(
		HM2D::Contour::Tree::GridBoundary(inpgrid));

	//build boundary mapping
	build_mcc();
}
void DoMapping::prepare_grid(){
	build_grid3();
	g3outer = HM2D::Contour::Tree::GridBoundary(*g3);
}
void DoMapping::build_mcc(){
	//mapping inpgrid_outer to mapped_outer
	int N = std::min(gridpoints.size(), contpoints.size());
	std::set<HM2D::EdgeData*> used_gcontours;
	for (int i=0; i<N; ++i){
		Point p1 = gridpoints[i];
		Point p2 = contpoints[i];
		auto _e1 = inpgrid_outer.alledges();
		auto _e2 = mapped_outer.alledges();
		int ie1 = std::get<0>(HM2D::Finder::ClosestEdge(_e1, p1));
		int ie2 = std::get<0>(HM2D::Finder::ClosestEdge(_e2, p2));
		HM2D::Edge* e1 = _e1[ie1].get();
		HM2D::Edge* e2 = _e2[ie2].get();
		HM2D::EdgeData* c1 = &inpgrid_outer.find_node(e1)->contour;
		HM2D::EdgeData* c2 = &mapped_outer.find_node(e2)->contour;
		used_gcontours.insert(c1);
		auto cfnd = mcol.insert(c1, c2);
		cfnd->add_connection(p1, p2);
	}
	if (inpgrid_outer.nodes.size() != used_gcontours.size())
		throw HMMap::MapException("All grid boundaries should contain at least one base point");
}

// =========================== DirectMapping
void DirectMapping::build_grid3(){
	//add corner points from inpgrid to mapped_outer
	auto inpgrid2 = HM2D::Contour::Tree::DeepCopy(inpgrid_outer);
	for (int i=0; i<mcol.entry_num(); ++i){
		//find mapped contour
		auto cmapping = mcol.get(i);
		auto base_contour = cmapping->get_base();
		HM2D::EdgeData* copied_base_contour=0;
		for (int j=0; j<inpgrid_outer.nodes.size(); ++j){
			if (&inpgrid_outer.nodes[j]->contour == base_contour){
				copied_base_contour = &inpgrid2.nodes[j]->contour;
				break;
			}
		}
		//for each point in contour
		for (auto p: HM2D::AllVertices(*cmapping->get_mapped())){
			//find mapped point
			Point pbase = cmapping->map_from_mapped(*p);
			HM2D::Contour::Algos::GuaranteePoint(*copied_base_contour, pbase);
		}
	}
	g3 = std::make_shared<HM2D::GridData>(
		HMFem::AuxGrid3(inpgrid2, opt.fem_nrec, opt.fem_nmax));
}
void DirectMapping::solve_uv_problems(vector<double>& u, vector<double>& v){
	u.resize(g3->vvert.size(), 0.0);
	v.resize(g3->vvert.size(), 0.0);
	HM2D::EdgeData eouter = inpgrid_outer.alledges();

	// u - problem
	laplace->ClearBC();
	for (auto n: g3outer.nodes){
		int ei = std::get<0>(HM2D::Finder::ClosestEdge(eouter, *HM2D::Contour::First(n->contour)));
		HM2D::Edge* e = eouter[ei].get();
		HM2D::EdgeData* c = &inpgrid_outer.find_node(e)->contour;
		auto cmapping = mcol.find_by_base(c);
		assert(cmapping != nullptr);
		auto dirfunc = [cmapping](const HM2D::Vertex* p)->double{
			return cmapping->map_from_base(*p).x;
		};
		laplace->SetDirichlet(n->contour, dirfunc);
	}
	laplace->Solve(u);

	// v - problem
	laplace->ClearBC();
	for (auto n: g3outer.nodes){
		int ei = std::get<0>(HM2D::Finder::ClosestEdge(eouter, *HM2D::Contour::First(n->contour)));
		HM2D::Edge* e = eouter[ei].get();
		HM2D::EdgeData* c = &inpgrid_outer.find_node(e)->contour;
		auto cmapping = mcol.find_by_base(c);
		assert(cmapping != nullptr);
		auto dirfunc = [cmapping](const HM2D::Vertex* p)->double{
			return cmapping->map_from_base(*p).y;
		};
		laplace->SetDirichlet(n->contour, dirfunc);
	}
	laplace->QuickSolve_BC(v);

}

// =========================== InverseMapping
void InverseMapping::build_grid3(){
	//add corner points from inpgrid to mapped_outer
	auto mapped2 = HM2D::Contour::Tree::DeepCopy(mapped_outer);
	for (int i=0; i<mcol.entry_num(); ++i){
		//find mapped contour
		auto cmapping = mcol.get(i);
		auto mapped_contour = cmapping->get_orig_mapped();
		HM2D::EdgeData* copied_mapped_contour=0;
		for (int j=0; j<mapped_outer.nodes.size(); ++j){
			if (&mapped_outer.nodes[j]->contour == mapped_contour){
				copied_mapped_contour = &mapped2.nodes[j]->contour;
				break;
			}
		}
		assert(copied_mapped_contour != 0);
		//for each point in contour
		for (auto p: HM2D::AllVertices(*cmapping->get_base())){
			//find mapped point
			Point pmapped = cmapping->map_from_base(*p);
			HM2D::Contour::Algos::GuaranteePoint(*copied_mapped_contour, pmapped);
		}
	}
	
	//build grid
	g3 = std::make_shared<HM2D::GridData>(
		HMFem::AuxGrid3(mapped2, opt.fem_nrec, opt.fem_nmax));
}

void InverseMapping::solve_uv_problems(vector<double>& u, vector<double>& v){
	u.resize(g3->vvert.size(), 0.0);
	v.resize(g3->vvert.size(), 0.0);
	HM2D::EdgeData eouter = mapped_outer.alledges();

	//u problem
	for (auto n: g3outer.nodes){
		int ei = std::get<0>(HM2D::Finder::ClosestEdge(eouter, *HM2D::Contour::First(n->contour)));
		HM2D::Edge* e = eouter[ei].get();
		HM2D::EdgeData* c = &mapped_outer.find_node(e)->contour;
		auto cmapping = mcol.find_by_mapped(c);
		assert(cmapping != 0);
		auto dirfunc = [cmapping](const HM2D::Vertex* p)->double{
			return cmapping->map_from_mapped(*p).x;
		};
		laplace->SetDirichlet(n->contour, dirfunc);
	}
	laplace->Solve(u);

	//v problem
	laplace->ClearBC();
	for (auto n: g3outer.nodes){
		int ei = std::get<0>(HM2D::Finder::ClosestEdge(eouter, *HM2D::Contour::First(n->contour)));
		HM2D::Edge* e = eouter[ei].get();
		HM2D::EdgeData* c = &mapped_outer.find_node(e)->contour;
		auto cmapping = mcol.find_by_mapped(c);
		assert(cmapping != 0);
		auto dirfunc = [cmapping](const HM2D::Vertex* p)->double{
			return cmapping->map_from_mapped(*p).y;
		};
		laplace->SetDirichlet(n->contour, dirfunc);
	}
	laplace->QuickSolve_BC(v);

	
	//swap g3 coordinates and uv
	for (int i=0; i<g3->vvert.size(); ++i){
		std::swap(g3->vvert[i]->x, u[i]);
		std::swap(g3->vvert[i]->y, v[i]);
	}
}
