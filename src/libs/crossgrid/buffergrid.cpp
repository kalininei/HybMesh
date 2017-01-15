#include "buffergrid.hpp"
#include "trigrid.hpp"
#include "modgrid.hpp"
#include "healgrid.hpp"
#include "buildgrid.hpp"
#include "cont_assembler.hpp"
#include "algos.hpp"
#include "contabs2d.hpp"

using namespace HM2D;
using namespace HM2D::Grid;
using namespace HM2D::Grid::Impl;

BufferGrid::BufferGrid(GridData& main, const EdgeData& source, double buffer_size, bool preserve_bp, double angle0):
		angle0(angle0), preserve_bp(preserve_bp), source(&source), orig(&main){
	Contour::Tree buffer_zone = Contour::Algos::Offset(source, -buffer_size,
			Contour::Algos::OffsetTp::RC_CLOSED_POLY);
	CellData::operator=(Grid::Algos::ExtractCells(main, buffer_zone, INSIDE));
}

const BoundingBoxFinder& BufferGrid::source_finder() const{
	if (_sfinder == nullptr){
		auto bbox = BBox(*source);
		_sfinder.reset(new BoundingBoxFinder(bbox, bbox.maxlen()/40));
		for (auto& e: (*source)){
			_sfinder->addentry(BoundingBox(*e->first(), *e->last()));
		}
	}
	return *_sfinder;
}

int BufferGrid::is_on_source(const Point& p) const{ //0 - no, 1 - vertex, 2 - on the middle of the edge
	double ksi;
	vector<int> susp = source_finder().suspects(p);
	for (int i: susp){
		auto ed = (*source)[i];
		if (p == *ed->first() || p == *ed->last()) return 1;
	}
	for (int i: susp){
		auto ed = (*source)[i];
		if (isOnSection(*ed->first(), *ed->last(), p, ksi)) return 2;
	}
	return 0;
}

VertexData BufferGrid::significant_boundary_points(const EdgeData& bedges) const{
	if (preserve_bp || ISLOWER(angle0, 0)){
		return AllVertices(bedges);
	} else {
		return AllVertices(ECol::Algos::Simplified(bedges, angle0));
	}
}

EdgeData BufferGrid::define_source_edges(const EdgeData& cont) const{
	EdgeData ret;

	auto av = AllVertices(cont);
	aa::constant_ids_pvec(av, 0);
	for (auto& v: av){
		if (is_on_source(*v)>0) v->id = 1;
	}

	for (auto& e: cont) if (e->first()->id == 1 && e->last()->id == 1){
		ret.push_back(e);
	}

	return ret;
}

void BufferGrid::rebuild_source_edges(EdgeData& cont, EdgeData& sed) const{
	vector<int> remove_vert;
	VertexData badpnt;
	for (auto v: AllVertices(sed)){
		for (int ei: source_finder().suspects(*v)){
			if (*v == *(*source)[ei]->first() ||
			    *v == *(*source)[ei]->last()){
				goto GOODPOINT;
			}
		}
		badpnt.push_back(v);
	GOODPOINT:
		continue;
	}
	EdgeData oldcont = cont;
	Contour::Algos::RemovePoints(cont, remove_vert);

	//remove those sed edges which are no longer in cont
	aa::constant_ids_pvec(sed, 0);
	aa::constant_ids_pvec(cont, 1);
	aa::remove_by_id(sed, 0);

	//all newly added to cont edges are source edges, so add them to sed
	aa::constant_ids_pvec(oldcont, 0);
	for (auto& e: cont) if (e->id != 0)
		sed.push_back(e);
}

void BufferGrid::split_edges(EdgeData& cont, EdgeData& bs, EdgeData& bg, EdgeData& bb) const{
	//define source edges
	bs = define_source_edges(cont);

	//define grid edges
	aa::constant_ids_pvec(cont, 0);
	aa::constant_ids_pvec(bs, 1);
	for (auto e: cont) if (e->id == 0 && !e->is_boundary()){
		e->id = 2;
		bg.push_back(e);
	}

	//define boundary edges
	for (auto e: cont) if (e->id == 0) bb.push_back(e);

	//rebuild source edges throwing away superfluous nodes
	rebuild_source_edges(cont, bs);
}

void BufferGrid::contour_segmentation(EdgeData& cont) const{
	//split given contour by edges lying on source, original grid edges and
	//original boundary edges.
	//Recomputes source edges throwing away points which are not source vertices.
	EdgeData bsource, bgrid, bboundary;
	split_edges(cont, bsource, bgrid, bboundary);

	//get boundary points which should be kept using preserve_bp and angle0 option
	VertexData keep_points = significant_boundary_points(bboundary);

	//do segmentation using keep_id1=true to keep needed primitives
	aa::constant_ids_pvec(bsource, 1);
	aa::constant_ids_pvec(bgrid, 1);
	aa::constant_ids_pvec(bboundary, 0);
	aa::constant_ids_pvec(keep_points, 1);
	cont = Mesher::RepartSourceById(cont);
}

Contour::Tree BufferGrid::triangulation_boundary() const{
	//get buffer cells boundary
	auto ginv = std::make_shared<Constructor::InvokeGrid>(*this);
	vector<EdgeData> conts = Contour::Assembler::GridBoundary(ginv->grid);
	ginv.reset();

	//partition
	for (auto& c: conts){ contour_segmentation(c); }

	//assemble tree
	Contour::Tree ret;
	for (auto& c: conts) { ret.add_contour(std::move(c)); }

	return ret;
}

//0-triangle, 1-recombined triangle
void BufferGrid::update_original(int filler) const{
	Contour::Tree tree = triangulation_boundary();
	GridData g3 = (filler == 0) ? Mesher::UnstructuredTriangle(tree)
	                            : Mesher::UnstructuredTriangleRecomb(tree);
	aa::constant_ids_pvec(orig->vcells, 0);
	aa::constant_ids_pvec(*this, 1);
	aa::remove_by_id(orig->vcells, 1);
	Algos::RestoreFromCells(*orig);
	Algos::MergeTo(g3, *orig);
}
