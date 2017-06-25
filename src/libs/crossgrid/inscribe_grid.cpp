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

using namespace HM2D;
using namespace HM2D::Grid;

HMCallback::FunctionWithCallback<Algos::TSubstractCells> Algos::SubstractCells;
HMCallback::FunctionWithCallback<Algos::TInscribeGrid> Algos::InscribeGrid;
HMCallback::FunctionWithCallback<Algos::TInsertConstraints> Algos::InsertConstraints;

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
		CellData badcells = ExtractCells.WithCallback(cb1, base, cont, BOUND);
		goodcells = invert_extraction(badcells);
		break;
	}
	case SubstractCellsAlgo::NO_CROSS:
	{
		goodcells = ExtractCells.WithCallback(cb1, base, cont, BOUND);
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

void segment_triarea(Contour::Tree& triarea, const EdgeData& keep_edges,
		const std::vector<std::pair<Point, double>>& src,
		double angle0){
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
			for (auto v: AllVertices(ECol::Algos::Simplified(c, angle0))){
				sig_vertices.push_back(v);
			}
		}
	}

	//do segmentation using id=1 to keep needed primitives
	aa::constant_ids_pvec(av, 0);
	aa::constant_ids_pvec(ae, 0);
	aa::constant_ids_pvec(sig_vertices, 1);
	aa::constant_ids_pvec(sig_edges, 1);
	for (auto c: triarea.bound_contours()){
		c->contour = Mesher::RepartSourceById(c->contour, src,
			c->isinner() ? 2 : 1);
	}
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
	std::vector<std::pair<Point, double>> src;
	for (auto e: keep_edges){
		src.emplace_back(e->center(), e->length());
	}
	segment_triarea(triarea, keep_edges, src, opt.angle0);
	
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
	opt.buffer_size = std::max(1e3*geps, opt.buffer_size);
	Contour::Tree bzone;
	for (auto& c: cont){
		Contour::R::Clockwise rc(c, false);
		Contour::Tree t1 = Contour::Algos::Offset(
			c, opt.buffer_size,
			Contour::Algos::OffsetTp::RC_OPEN_ROUND);
		bzone = Contour::Clip::Union(t1, bzone);
	}
	//offsetting from points
	for (auto& p: pnt){
		auto cp = Contour::Constructor::Circle(32, 1e3*geps, p.first);
		Contour::Tree t1 = Contour::Algos::Offset(cp,
			std::max(1e3*geps, opt.buffer_size-1e3*geps),
			Contour::Algos::OffsetTp::RC_CLOSED_POLY);
		bzone = Contour::Clip::Union(t1, bzone);
	}

	//remove cells under the offset
	//TODO: consider moving this block to InscribeProcedure
	CellData needed_cells = ExtractCells(base, bzone, OUTSIDE);
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
	Contour::Tree triarea = Contour::Tree::GridBoundary(g3);
	EdgeData keep_edges = ECol::Assembler::GridBoundary(g2);
	std::vector<std::pair<Point, double>> src_tri;
	for (auto& p: pnt) if (p.second > 0){
		src_tri.push_back(p);
	}
	if (opt.keep_cont) for (auto& c: cont){
		for (auto& e: c){
			src_tri.emplace_back(e->center(), e->length());
		}
	}
	segment_triarea(triarea, keep_edges, src_tri, opt.angle0);

	//part constraints
	vector<EdgeData> cont2 = cont;
	if (!opt.keep_cont){
		//TODO
	}

	//calculate point constraint sizes
	CoordinateMap2D<double> pnt2;
	for (auto& p: pnt){
		double v = p.second;
		if (v <= 0){
			//TODO
		}
		pnt2.add(p.first, v);
	}

	//Triangulation
	for (auto& c: cont2) triarea.add_detached_contour(c);
	GridData g4;
	if (opt.fillalgo == 0) g4 = Mesher::UnstructuredTriangle(triarea, pnt2);
	else if (opt.fillalgo == 1) g4 = Mesher::UnstructuredTriangleRecomb(triarea, pnt2);
	
	//Merge
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
