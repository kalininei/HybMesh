#include "inscribe_grid.hpp"
#include "infogrid.hpp"
#include "healgrid.hpp"
#include "clipdomain.hpp"
#include "trigrid.hpp"
#include "modgrid.hpp"
#include "assemble2d.hpp"
#include "treverter2d.hpp"
#include "finder2d.hpp"
#include "modcont.hpp"
#include "hmtimer.hpp"
#include "buildcont.hpp"
#include "sizefun.hpp"
#include "wireframegrid.hpp"

using namespace HM2D;
using namespace HM2D::Grid;

HMCallback::FunctionWithCallback<Algos::TSubstractCells> Algos::SubstractCells;
HMCallback::FunctionWithCallback<Algos::TInscribeGrid> Algos::InscribeGrid;
HMCallback::FunctionWithCallback<Algos::TInsertConstraints> Algos::InsertConstraints;
HMCallback::FunctionWithCallback<Algos::TSubstractArea> Algos::SubstractArea;

GridData Algos::TSubstractCells::_run(const GridData& base, const Contour::Tree& cont, SubstractCellsAlgo algo){
	auto invert_extraction = [&](const CellData& extr)->CellData{
		aa::constant_ids_pvec(base.vcells, 0);
		aa::constant_ids_pvec(extr, 1);
		CellData ret;
		std::copy_if(base.vcells.begin(), base.vcells.end(), std::back_inserter(ret),
				[](shared_ptr<Cell> c){ return c->id == 0; });
		return ret;
	};

	CellData goodcells;

	auto cb1 = callback->bottom_line_subrange(95);
	switch (algo){
	case SubstractCellsAlgo::FULLY_INSIDE:
	{
		CellData badcells = ExtractCells.WithCallback(cb1, base, cont, INSIDE);
		goodcells = invert_extraction(badcells);
		break;
	}
	case SubstractCellsAlgo::FULLY_OUTSIDE:
	{
		CellData badcells = ExtractCells.WithCallback(cb1, base, cont, OUTSIDE);
		goodcells = invert_extraction(badcells);
		break;
	}
	case SubstractCellsAlgo::PARTLY_INSIDE:
	{
		goodcells = ExtractCells.WithCallback(cb1, base, cont, OUTSIDE);
		break;
	}
	case SubstractCellsAlgo::PARTLY_OUTSIDE:
	{
		goodcells = ExtractCells.WithCallback(cb1, base, cont, INSIDE);
		break;
	}
	case SubstractCellsAlgo::CROSS:
	{
		goodcells = ExtractCells.WithCallback(cb1, base, cont, BOUND);
		break;
	}
	case SubstractCellsAlgo::NO_CROSS:
	{
		CellData badcells = ExtractCells.WithCallback(cb1, base, cont, BOUND);
		goodcells = invert_extraction(badcells);
		break;
	}
	}

	callback->step_after(5, "Assembling grid");
	GridData ret;
	DeepCopy(goodcells, ret.vcells, 2);
	RestoreFromCells(ret);

	return ret;
}

namespace{

#if 0
void id1_to_significant_prims(Contour::Tree& triarea, const EdgeData& keep_edges,
		const VertexData& keep_vert, double angle0){
	auto ae = triarea.alledges();
	auto av = HM2D::AllVertices(ae);

	//Significant edges
	EdgeData sig_edges;
	auto bb = HM2D::BBox(HM2D::AllVertices(keep_edges));
	Finder::RasterizeEdges rast(keep_edges, bb, bb.maxlen()/30);
	aa::constant_ids_pvec(av, 0);
	double ksi;
	for (auto v: av){
		for (int isusp: rast.bbfinder().suspects(*v)){
			auto esusp = keep_edges[isusp];
			if (isOnSection(*v, *esusp->pfirst(), *esusp->plast(), ksi)){
				v->id = 1;
				break;
			}
		}
	}
	for (auto e: ae) if (e->pfirst()->id == 1 && e->plast()->id == 1){
		Point pmid = Point::Weigh(*e->pfirst(), *e->plast(), 0.5);
		for (int isusp: rast.bbfinder().suspects(pmid)){
			auto esusp = keep_edges[isusp];
			if (isOnSection(pmid, *esusp->pfirst(), *esusp->plast(), ksi)){
				sig_edges.push_back(e);
				break;
			}
		}
	}

	//Significant vertices
	VertexData sig_vertices;
	if (angle0 >= 0){
		EdgeData nonsig_edges = ae;
		aa::constant_ids_pvec(nonsig_edges, 0);
		aa::constant_ids_pvec(sig_edges, 1);
		aa::keep_by_id(nonsig_edges, 0);
		for (auto c: Contour::Assembler::SimpleContours(nonsig_edges)){
			for (auto v: AllVertices(ECol::Algos::Simplified(c, angle0, false, keep_vert))){
				sig_vertices.push_back(v);
			}
		}
	} else sig_vertices = av;

	//do segmentation using id=1 to keep needed primitives
	aa::constant_ids_pvec(av, 0);
	aa::constant_ids_pvec(ae, 0);
	aa::constant_ids_pvec(sig_vertices, 1);
	aa::constant_ids_pvec(sig_edges, 1);
}

#else
struct _Id1ToSigPrimsND{
	HM2D::EdgeData* cont;
	shared_ptr<HM2D::Contour::R::ReallyDirect> rev;
	shared_ptr<BoundingBoxFinder> bbfinder;
	vector<double> accum_len;
	void set_data(HM2D::EdgeData& _cont){
		cont = &_cont;
		//!!!! during all struct operation *cont is a really directed contour.
		rev.reset(new HM2D::Contour::R::ReallyDirect(_cont));
		BoundingBox bb = HM2D::BBox(_cont);
		bbfinder.reset(new BoundingBoxFinder(bb, bb.maxlen()/30));
		for (auto& e: _cont){
			bbfinder->raw_addentry(
				bbfinder->sqrs_by_segment(*e->pfirst(), *e->plast()));
		}
		auto lens = HM2D::Contour::ELengths(*cont);
		accum_len.resize(lens.size()+1, 0);
		std::partial_sum(lens.begin(), lens.end(), accum_len.begin()+1);
	}
	double coords(const Point& p){
		double ksi;
		for (auto esus: bbfinder->suspects(p)){
			auto e = (*cont)[esus];
			if (isOnSection(p, *e->pfirst(), *e->plast(), ksi)){
				if (ISZERO(ksi)) return accum_len[esus]/contlen();
				assert(esus+1 < accum_len.size());
				if (ISZERO(ksi-1)) return accum_len[esus+1]/contlen();
				return (accum_len[esus] + ksi*elen(esus))/contlen();
			}
		}
		return -1;
	}
	bool get_coords(const Point& p1, const Point& p2, double& c1, double& c2){
		c1 = coords(p1);
		if (c1 == -1) return false;
		c2 = coords(p2);
		if (c2 == -1) return false;
	
		//if edge length equals distance along the contour between
		//edge first and last vertices.
		double len1 = Point::dist(p1, p2);
		double wlen;
		if (HM2D::Contour::IsClosed(*cont)){
			wlen = (c2>c1) ? std::min(c2-c1, 1+c1-c2)
		                       : std::min(c1-c2, 1+c2-c1);
		} else {
			wlen = fabs(c2-c1);
		}
		double len2 = contlen()*wlen;
		if (!ISEQ(len1, len2)) return false;
		return true;
	}

	double contlen(){ return accum_len.back(); }
	double elen(int ie){ return accum_len[ie+1]-accum_len[ie]; }
	shared_ptr<HM2D::Edge> guarantee_edge(double w1, double w2){
		if (ISEQ(w1, 0) && w2>0.5) w1 = 1;
		if (ISEQ(w2, 0) && w1>0.5) w2 = 1;
		double len1 = w1*contlen(), len2 = w2*contlen();
		if (len2<len1) std::swap(len1, len2);
		//get edge index
		auto it1 = std::lower_bound(accum_len.begin(), accum_len.end(), len1);
		auto it2 = std::lower_bound(accum_len.begin(), accum_len.end(), len2);
		if (!ISEQ(*it1, len1)) {--it1;}
		if (ISEQ(*(--it2), len2)) {--it2;}
		int iedge = it1 - accum_len.begin();
		int iedge2 = it2 - accum_len.begin();
		assert(iedge==iedge2);
		auto e = (*cont)[iedge];
		//get local edge coordinate. All edges are correctly directed.
		double loc1 = (len1-accum_len[iedge])/elen(iedge);
		double loc2 = (len2-accum_len[iedge])/elen(iedge);
		//if [w1-w2] is an existing edge return it.
		if (ISEQ(loc1, 0) && ISEQ(loc2, 1)) return e;
		//!!! Here we split an existing edge.
		//    accum_len and bbfinder will be adjusted
		if (ISEQ(loc1, 0)){
			HM2D::Contour::Algos::SplitEdge(*cont, iedge,
				{Point::Weigh(*e->pfirst(), *e->plast(), loc2)});
			accum_len.insert(accum_len.begin()+iedge+1,
				accum_len[iedge] + loc2*elen(iedge));
			bbfinder->raw_insertentry(bbfinder->sqrs_by_segment(
					*(*cont)[iedge+1]->pfirst(), *(*cont)[iedge+1]->plast()), iedge+1);
			return (*cont)[iedge];
		} else if (ISEQ(loc2, 1)){
			HM2D::Contour::Algos::SplitEdge(*cont, iedge,
				{Point::Weigh(*e->pfirst(), *e->plast(), loc1)});
			accum_len.insert(accum_len.begin()+iedge+1,
				accum_len[iedge] + loc1*elen(iedge));
			bbfinder->raw_insertentry(bbfinder->sqrs_by_segment(
					*(*cont)[iedge+1]->pfirst(), *(*cont)[iedge+1]->plast()), iedge+1);
			return (*cont)[iedge+1];
		} else{
			HM2D::Contour::Algos::SplitEdge(*cont, iedge,
				{Point::Weigh(*e->pfirst(), *e->plast(), loc1),
				 Point::Weigh(*e->pfirst(), *e->plast(), loc2)});
			double v1 = accum_len[iedge] + loc1*elen(iedge);
			double v2 = accum_len[iedge] + loc2*elen(iedge);
			accum_len.insert(accum_len.begin()+iedge+1, v2);
			accum_len.insert(accum_len.begin()+iedge+1, v1);
			bbfinder->raw_insertentry(bbfinder->sqrs_by_segment(
					*(*cont)[iedge+1]->pfirst(), *(*cont)[iedge+1]->plast()), iedge+1);
			bbfinder->raw_insertentry(bbfinder->sqrs_by_segment(
					*(*cont)[iedge+2]->pfirst(), *(*cont)[iedge+2]->plast()), iedge+2);
			return (*cont)[iedge+1];
		}
	}
};
struct _Id1ToSigPrims{
	Contour::Tree* src;
	vector<_Id1ToSigPrimsND> nds;
	_Id1ToSigPrims(Contour::Tree& _src){
		src = &_src;
		nds.resize(src->nodes.size());
		for (int i=0; i<src->nodes.size(); ++i){
			nds[i].set_data(src->nodes[i]->contour);
		}
	}
	double nodelen(int inode){ return nds[inode].contlen(); }

	shared_ptr<HM2D::Edge> guarantee_edge(int inode, double w1, double w2){
		return nds[inode].guarantee_edge(w1, w2);
	}
	int get_contour(const Point& p1, const Point& p2, double& c1, double& c2){
		for (int i=0; i<nds.size(); ++i){
			if (nds[i].get_coords(p1, p2, c1, c2)){
				return i;
			}
		}
		return -1;
	}
};

void id1_to_significant_prims(Contour::Tree& triarea, const EdgeData& keep_edges,
		const VertexData& keep_vert, double angle0){
	EdgeData ae = triarea.alledges();
	VertexData av = HM2D::AllVertices(ae);
	EdgeData sig_edges;
	VertexData sig_vertices;
	//-- get needed edges
	_Id1ToSigPrims helper(triarea);
	//auto keep_edges_v = HM2D::AllVertices(keep_edges);
	for (auto e: keep_edges){
		double coord1, coord2;
		int inode = helper.get_contour(*e->pfirst(), *e->plast(), coord1, coord2);
		if (inode>=0) sig_edges.push_back(helper.guarantee_edge(inode, coord1, coord2));
	}

	//-- get needed vertices
	if (angle0 >= 0){
		EdgeData nonsig_edges = ae;
		aa::constant_ids_pvec(nonsig_edges, 0);
		aa::constant_ids_pvec(sig_edges, 1);
		aa::keep_by_id(nonsig_edges, 0);
		for (auto c: Contour::Assembler::SimpleContours(nonsig_edges)){
			for (auto v: AllVertices(ECol::Algos::Simplified(c, angle0, false, keep_vert))){
				sig_vertices.push_back(v);
			}
		}
	} else sig_vertices = av;

	//do segmentation using id=1 to keep needed primitives
	aa::constant_ids_pvec(av, 0);
	aa::constant_ids_pvec(ae, 0);
	aa::constant_ids_pvec(sig_vertices, 1);
	aa::constant_ids_pvec(sig_edges, 1);
}
#endif

VertexData add_tree_crosses(HM2D::Contour::Tree& tree){
	VertexData ret;
	int n = tree.nodes.size();
	for (int i=0; i<n; ++i)
	for (int j=i+1; j<n; ++j){
		if (tree.nodes[i]->isbound() && tree.nodes[j]->isbound()) continue;
		auto& c1 = tree.nodes[i]->contour;
		auto& c2 = tree.nodes[j]->contour;
		auto crosses = HM2D::Contour::Finder::CrossAll(c1, c2);
		for (auto& c: crosses){
			auto res1 = HM2D::Contour::Algos::GuaranteePoint(c1, std::get<1>(c));
			auto res2 = HM2D::Contour::Algos::GuaranteePoint(c2, std::get<1>(c));
			ret.push_back(std::get<1>(res1));
			ret.push_back(std::get<1>(res2));
		}
	}
	return ret;
}

}

GridData Algos::TInscribeGrid::_run(const GridData& base, const Contour::Tree& cont2,
		OptInscribe opt){
	Contour::Tree cont = cont2; cont.remove_detached();
	//1) Throw away cells which are outside cont
	auto cb1 = callback->bottom_line_subrange(20);
	auto g1algo = opt.inside ? Algos::SubstractCellsAlgo::PARTLY_OUTSIDE
	                         : Algos::SubstractCellsAlgo::PARTLY_INSIDE;
	GridData g1 = SubstractCells.WithCallback(cb1, base, cont, g1algo);

	//2) offsetting
	//delta >> geps to avoid geometry errors
	callback->step_after(10, "Offset");
	opt.buffer_size = std::max(1e3*geps, opt.buffer_size);
	Contour::Tree bzone;
	for (auto& c: cont.bound_contours()){
		Contour::R::Clockwise rc(c->contour, false);
		double delta = opt.buffer_size;
		if (opt.inside) delta*=-1;
		if (c->isinner()) delta*=-1;
		Contour::Tree t1 = Contour::Algos::Offset(
			c->contour,
			delta,
			Contour::Algos::OffsetTp::RC_CLOSED_POLY);
		t1.add_contour(c->contour);
		bzone = Contour::Clip::Union(bzone, t1);
	}

	//3) Throw away cells within the buffer
	auto cb3 = callback->bottom_line_subrange(20);
	GridData g2 = SubstractCells.WithCallback(cb3, g1, bzone, Algos::SubstractCellsAlgo::PARTLY_INSIDE);
	if (opt.fillalgo == 99) return g2;

	//4) Assemble triangulation area
	callback->step_after(10, "Triangulation area");
	Contour::Tree triarea = Contour::Tree::DeepCopy(cont);
	for (auto c: Contour::Assembler::GridBoundary(g2)){
		triarea.add_contour(c);
	}
	if (!opt.inside){
		triarea = Contour::Clip::Difference(Contour::Tree::GridBoundary(base), triarea);
	}

	//5) 1D segmentation of triangulation area
	EdgeData keep_edges = ECol::Assembler::GridBoundary(g2);
	if (opt.keep_cont){
		for (auto e: cont.alledges()) keep_edges.push_back(e);
	}
	//sources
	id1_to_significant_prims(triarea, keep_edges, {}, opt.angle0);
	ApplySizeFunction(triarea);
	
	//6) Triangulation
	auto cb6 = callback->bottom_line_subrange(20);
	GridData g3;
	if (opt.fillalgo == 0) g3 = Mesher::UnstructuredTriangle.WithCallback(cb6, triarea);
	else if (opt.fillalgo == 1) g3 = Mesher::UnstructuredTriangleRecomb.WithCallback(cb6, triarea);

	//7) Merging
	callback->step_after(10, "Merging");
	Grid::Algos::MergeTo(g2, g3);
	
	return g3;
}

GridData Algos::TInsertConstraints::_run(const GridData& base,
		const vector<EdgeData>& cont,
		const vector<std::pair<Point, double>>& pnt,
		Algos::OptInsertConstraints opt){
	//offsetting from polylines
	callback->step_after(10, "Offset");
	opt.buffer_size = std::max(1e3*geps, opt.buffer_size);
	Contour::Tree bzone;
	for (auto& c: cont){
		assert(Contour::IsContour(c));
		Contour::R::Clockwise rc(c, false);
		Contour::Tree t1 = Contour::Algos::Offset(
			c, opt.buffer_size,
			Contour::Algos::OffsetTp::RC_OPEN_ROUND);
		bzone = Contour::Clip::Union(t1, bzone);
	}
	//offsetting from points
	for (auto& p: pnt){
		auto crc = Contour::Constructor::Circle(32, opt.buffer_size, p.first);
		Contour::Tree t1; t1.add_contour(crc);
		bzone = Contour::Clip::Union(t1, bzone);
	}

	//remove cells under the offset
	auto cb1 = callback->bottom_line_subrange(20);
	CellData needed_cells = ExtractCells.WithCallback(cb1, base, bzone, OUTSIDE);
	aa::constant_ids_pvec(base.vcells, 0);
	aa::constant_ids_pvec(needed_cells, 1);
	CellData not_needed_cells;
	not_needed_cells.reserve(base.vcells.size() - needed_cells.size());
	for (auto& c: base.vcells) if (c->id==0) not_needed_cells.push_back(c);
	GridData g2, g3;
	DeepCopy(needed_cells, g2.vcells);
	RestoreFromCells(g2);
	if (opt.fillalgo == 99) return g2;
	DeepCopy(not_needed_cells, g3.vcells);
	RestoreFromCells(g3);

	//triangulation area
	callback->step_after(35, "Assemble buffer area", 3, 1);
	Contour::Tree triarea = Contour::Tree::GridBoundary(g3);
	for (auto& c: cont) triarea.add_detached_contour(c);
	//paste crosses
	VertexData cr_points = add_tree_crosses(triarea);
	//id=1 to edges and vertices which should be kept
	EdgeData keep_edges = ECol::Assembler::GridBoundary(g2);
	if (opt.keep_cont) for (auto& n: triarea.detached_contours()){
		keep_edges.insert(keep_edges.end(),
			n->contour.begin(), n->contour.end());
	}
	id1_to_significant_prims(triarea, keep_edges, cr_points, opt.angle0);
	//build point constraints and check if we will need explicit size function in future.
	bool force_sizefun = false;
	vector<std::pair<Point, double>> src_tri;
	for (auto& p: pnt) {
		if (p.second > 0) src_tri.push_back(p);
		else { force_sizefun = true; }
	}
	//generate 1D mesh
	callback->subprocess_step_after(2);
	auto sfun = ApplySizeFunction(triarea, src_tri, force_sizefun);

	//calculate point constraint sizes if they were not given using sfun.
	CoordinateMap2D<double> pnt2;
	for (auto& p: pnt){
		double v = p.second;
		if (v <= 0) v = sfun->sz(p.first);
		pnt2.add(p.first, v);
	}

	//Triangulation
	auto cb2 = callback->bottom_line_subrange(30);
	GridData g4;
	if (opt.fillalgo == 0) g4 = Mesher::UnstructuredTriangle.WithCallback(cb2, triarea, pnt2);
	else if (opt.fillalgo == 1) g4 = Mesher::UnstructuredTriangleRecomb.WithCallback(cb2, triarea, pnt2);
	
	//Merge
	callback->step_after(5, "Merge");
	Grid::Algos::MergeTo(g2, g4);

	return  g4;
}

GridData Algos::TInsertConstraints::_run(const GridData& base,
		const vector<EdgeData>& cont,
		Algos::OptInsertConstraints opt){
	return _run(base, cont, vector<std::pair<Point, double>>(), opt);
}

GridData Algos::TInsertConstraints::_run(const GridData& base,
		const vector<std::pair<Point, double>>& pnt,
		Algos::OptInsertConstraints opt){
	return _run(base, {}, pnt, opt);
}

GridData Algos::TSubstractArea::_run(const GridData& ginp, const Contour::Tree& area, bool is_inner){
	GridData gg;
	callback->step_after(20, "Initialize");
	int goodpos = is_inner == true ? OUTSIDE : INSIDE;
	int badpos = is_inner == true ? INSIDE : OUTSIDE;
	//building the new area
	Contour::Tree ar = Contour::Tree::DeepCopy(area);
	ar.remove_detached();
	Finder::RasterFinder finder(ar, 100);

	//cells which are not crossed by area lines
	auto cb1 = callback->bottom_line_subrange(25);
	CellData not_crossed = ExtractCells.UseCallback(cb1, ginp, area, BOUND);
	aa::constant_ids_pvec(ginp.vcells, 0);
	aa::constant_ids_pvec(not_crossed, 1);
	CellData ambigious = aa::copy_by_id(ginp.vcells, 0);
	not_crossed = finder.extract_cells(not_crossed, goodpos);

	//create aux grids
	callback->step_after(10, "Ambigious area");
	GridData gamb;
	HM2D::DeepCopy(ambigious, gamb.vcells, 2);
	RestoreFromCells(gamb);
	GridData gnc;
	HM2D::DeepCopy(not_crossed, gnc.vcells, 2);
	RestoreFromCells(gnc);

	//extract needed primitives from crossed
	//building a graph
	callback->step_after(25, "Cross graph", 5, 1);
	Impl::PtsGraph graph(gamb);
	//add contour edges
	callback->subprocess_step_after(1);
	for (auto& n: ar.nodes){ graph.add_edges(n->contour); }

	//exclude edges
	//Contour::Tree grid_ar = Contour::Tree::GridBoundary(ginp);
	//Contour::Tree excl_ar;
	//if (is_inner == true) excl_ar = Contour::Clip::Intersection(grid_ar, ar);
	//else excl_ar = Contour::Clip::Difference(grid_ar, ar);
	//if (excl_ar.nodes.size() == 0){
	//        HM2D::DeepCopy(ginp, gg);
	//        return gg;
	//}
	//graph.exclude_area(ar, badpos);

	//build grid
	callback->subprocess_step_after(1);
	Contour::Tree grapharea;
	if (is_inner){
		grapharea = HM2D::Contour::Clip::Difference(
			HM2D::Contour::Tree::GridBoundary(gamb), ar);
	} else {
		grapharea = HM2D::Contour::Clip::Intersection(
			HM2D::Contour::Tree::GridBoundary(gamb), ar);
	}
	callback->subprocess_step_after(1);
	graph.purge_endpoints();
	callback->subprocess_step_after(1);
	gg = graph.togrid(grapharea);
	callback->subprocess_fin();

	//remove elements which can present in the grid due
	//to outer contour exclusion
	callback->step_after(10, "Merging");
	Algos::MergeBoundaries(gnc, gg);

	//place area boundary before grid boundary to
	//guarantee higher priority of area edges.
	callback->step_after(10, "Boundary types");
	EdgeData cd = area.alledges_bound();
	EdgeData gd = HM2D::ECol::Assembler::GridBoundary(ginp);
	cd.insert(cd.end(), gd.begin(), gd.end());
	EdgeData rd = HM2D::ECol::Assembler::GridBoundary(gg);
	HM2D::ECol::Algos::AssignBTypes(cd, rd);

	return gg;
}
