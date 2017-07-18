#include "unite_grids.hpp"
#include "contour_tree.hpp"
#include "treverter2d.hpp"
#include "modgrid.hpp"
#include "buffergrid.hpp"
#include "modcont.hpp"
#include "clipdomain.hpp"
#include "finder2d.hpp"
#include "assemble2d.hpp"
#include "nodes_compare.h"
#include "inscribe_grid.hpp"
#include "debug_grid2d.hpp"
#include "wireframegrid.hpp"

using namespace HM2D;
using namespace HM2D::Grid;

HMCallback::FunctionWithCallback<Algos::TUniteGrids> Algos::UniteGrids;
HMCallback::FunctionWithCallback<Algos::TCombineGrids> Algos::CombineGrids;

namespace{

struct BoundaryShifter{
	VertexData tarpoints;
	vector<Point> oldcoords;
	const CellData* cells2;
	const Contour::Tree* tree2;
	BoundaryShifter(const Contour::Tree& c1, const Contour::Tree& c2, const CellData& grid2, double angle0){
		cells2 = &grid2;
		tree2 = &c2;
		for (auto& n1: c1.nodes)
		for (auto& n2: c2.nodes){
			process_pair(n1->contour, n2->contour, angle0);
		}
	}
	~BoundaryShifter(){
		for (int i=0; i<tarpoints.size(); ++i){
			tarpoints[i]->set(oldcoords[i]);
		}
	}
private:
	void process_pair(const EdgeData& c1, const EdgeData& c2, double angle0){
		//calculate crosses
		auto crosses = HM2D::Contour::Finder::CrossAll(c1, c2);
		//mark all c2 points which already equal c1 vertices
		auto c1v = AllVertices(c1), c2v = AllVertices(c2);
		Finder::VertexMatch vm(c1v);
		aa::constant_ids_pvec(c2v, 0);
		for (auto& v: c2v){
			if (vm.find(*v)) v->id = 1;
		}
		//loop over all crosses
		for (auto& c: crosses){
			Point p = std::get<1>(c);
			auto crd = HM2D::Contour::CoordAt(c2, p);
			double ksi = std::get<3>(crd);
			//if c2 has point at intersection
			if (!ISIN_NN(ksi, 0, 1)) continue;
			Edge* ied = c2[std::get<2>(crd)].get();
			shared_ptr<Vertex> p1 = ied->first();
			shared_ptr<Vertex> p2 = ied->last();

			//do not move points which already equal c1 vertices
			if (p1->id == 1) p1 = nullptr;
			if (p2->id == 1) p2 = nullptr;
			if (!p1 && !p2) continue;

			if (ksi > 0.5) std::swap(p1, p2);
			aa::RestoreIds<VertexData> rr(c2v);
			shift(p1, p2, p, angle0);
		}
	}
	//no need for complicated searchers since these functions
	//are called very rarely.
	std::pair<Point*, Point*> sibling_points(shared_ptr<Vertex> p1){
		std::pair<Point*, Point*> ret(0, 0);
		for (auto n: tree2->nodes){
			for (auto ed: n->contour)
			if (ed->first() == p1 || ed->last() == p1){
				if (ret.first == 0){
					ret.first = ed->sibling(p1.get()).get();
				} else if (ret.second == 0){
					ret.second = ed->sibling(p1.get()).get();
					return ret;
				}
			}
			assert(ret.first == 0);
		}
		assert(false);
		return ret;
	}
	vector<Cell*> sibling_cells(shared_ptr<Vertex> p1){
		vector<Cell*> ret;
		for (auto c: *cells2)
		for (auto v: AllVertices(c->edges)){
			if (v == p1){
				ret.push_back(c.get());
				break;
			}
		}
		assert(ret.size()>0);
		return ret;
	}

	bool iscorner(shared_ptr<Vertex> p1, double angle0){
		auto sibs = sibling_points(p1);
		double a = Angle(*sibs.first, *p1, *sibs.second)/M_PI*180;
		return !ISIN_EE(a, 180-angle0, 180+angle0);
	}

	bool try_to_shift(shared_ptr<Vertex> p1, Point newpoint){
		vector<Cell*> sibs = sibling_cells(p1);
		bool ret = true;
		Point old(*p1);
		p1->set(newpoint);
		for (auto c: sibs){
			if (std::get<0>(Contour::Finder::SelfCross(c->edges))){
				ret = false;
				break;
			}
		}
		if (ret == true){
			tarpoints.push_back(p1);
			oldcoords.push_back(old);
		} else {
			p1->set(old);
		}
		return ret;
	}

	void shift(shared_ptr<Vertex> p1, shared_ptr<Vertex> p2,
			Point newpoint, double a0){
		if (p1 && iscorner(p1, a0)) p1 = nullptr;
		if (p2 && iscorner(p2, a0)) p2 = nullptr;
		if (!p1 && !p2) return;
		if (p1 && try_to_shift(p1, newpoint)) return;
		if (p2) try_to_shift(p2, newpoint);
	}
};

int tree_highest_level(const Contour::Tree& tree){
	int ret = -1;
	for (auto& n: tree.nodes){
		if (n->level>ret) ret = n->level;
	}
	return ret;
}

Contour::Tree root_nodes_tree(const Contour::Tree& tree){
	Contour::Tree ret;
	for (auto& n: tree.nodes) if (n->level == 0){
		ret.add_detached_contour(n->contour);
	}
	for (auto& n: ret.nodes){
		n->level = 0;
	}
	return ret;
}

}

GridData Algos::TUniteGrids::_run(const GridData& base, const GridData& sec, const OptUnite& opt){
	callback->step_after(10, "Boundary analyzing", 5, 2);
	//----- contours assembling
	Contour::Tree contbase = Contour::Tree::GridBoundary(base);
	Contour::Tree contsec = Contour::Tree::GridBoundary(sec);
	Contour::R::RevertTree rev1(contbase);
	Contour::R::RevertTree rev2(contsec);
	GridData ret = base;
	Contour::Tree contret = contbase;

	//---- find contours intersection points and place gsec nodes there
	callback->subprocess_step_after(3);
	std::unique_ptr<BoundaryShifter> bsh;
	if (!opt.preserve_bp){
		bsh.reset(new BoundaryShifter(contbase, contsec, sec.vcells, opt.angle0));
	}
	callback->subprocess_fin();

	//---- combine grids without using buffers
	auto cb2 = callback->bottom_line_subrange(50);
	ret = CombineGrids.UseCallback(cb2, ret, sec, opt.empty_holes);
	if (opt.buffer_size < geps) return ret;
	vector<GridData> sg = SplitData(ret);

	//----- fill buffer
	callback->silent_step_after(30, "Filling buffer", sg.size()*contsec.bound_contours().size());
	for (int i=0; i<sg.size(); ++i)
	for (auto n: contsec.bound_contours()){
		callback->subprocess_step_after(1);
		EdgeData& csec = n->contour;
		Impl::BufferGrid bg(sg[i], csec, opt.buffer_size, opt.preserve_bp, opt.angle0);
		bg.update_original(opt.filler);
	}
	callback->subprocess_fin();

	//----- merge
	callback->step_after(10, "Final merging");
	ret = std::move(sg[0]);
	for (int i=1; i<sg.size(); ++i){
		Grid::Algos::MergeTo(sg[i], ret);
	}

	return ret;
}

GridData Algos::TCombineGrids::_run(const GridData& g1, const GridData& g2, bool keep_g2_holes){
	//1) build secondary contours
	callback->step_after(10, "Assemble boundary");
	Contour::Tree c2 = Contour::Tree::GridBoundary(g2);
	//if keep g2 holes, then leave only root nodes
	if (keep_g2_holes) c2 = root_nodes_tree(c2);
	//2) substract
	auto cb1 = callback->bottom_line_subrange(50);
	GridData ret = SubstractArea.UseCallback(cb1, g1, c2, true);
	//3) merge
	callback->step_after(40, "Merge grids");
	Algos::MergeBoundaries(g2, ret);
	//5) get rid of hanging + non-significant + boundary nodes
	//   They appear if boundaries of g1 and g2 partly coincide
	Grid::Algos::SimplifyBoundary(ret, 0);
	return ret;
}
