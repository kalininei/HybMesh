#include "buffergrid.hpp"
#include "trigrid.hpp"
#include "modgrid.hpp"
#include "healgrid.hpp"
#include "buildgrid.hpp"
#include "assemble2d.hpp"
#include "modcont.hpp"
#include "contabs2d.hpp"
#include "infogrid.hpp"
#include "treverter2d.hpp"
#include "clipdomain.hpp"
#include "finder2d.hpp"
#include "debug_grid2d.hpp"

using namespace HM2D;
using namespace HM2D::Grid;
using namespace HM2D::Grid::Impl;

vector<EdgeData> BufferGrid::build_contacts(const EdgeData& bedges) const{
	EdgeData contact_edges;
	for (auto e: bedges){
		if (is_on_source(e->center()) != 0){
			contact_edges.push_back(e);
		}
	}
	if (contact_edges.size() == 0) return {};
	vector<EdgeData> contacts = Contour::Assembler::SimpleContours(contact_edges);
	//bedges have direction according to buffergrid contour tree.
	//we want contact lines to be directed according to it to specify offset direction.
	for (auto& c: contacts){
		if (!Contour::CorrectlyDirectedEdge(c, 0)) Contour::Algos::Reverse(c);
	}
	return contacts;
}

Contour::Tree BufferGrid::offset_closed_contact(const Contour::Tree& outer_bnd,
		const EdgeData& contact, double bsize) const{
	Contour::Tree tree = Contour::Algos::Offset(contact, -bsize, 
			Contour::Algos::OffsetTp::RC_CLOSED_POLY);
	tree.add_contour(contact);
	return Contour::Clip::Intersection(tree, outer_bnd);
}

Contour::Tree BufferGrid::offset_open_contact(const Contour::Tree& outer_bnd,
		const EdgeData& contact, double bsize) const{
	Contour::Tree tree = Contour::Algos::Offset(contact, bsize, 
			Contour::Algos::OffsetTp::RC_OPEN_ROUND);
	return Contour::Clip::Intersection(tree, outer_bnd);
}

Contour::Tree BufferGrid::offset_contact(const Contour::Tree& outer_bnd,
		const EdgeData& contact, double bsize) const{
	Contour::Tree ret = (Contour::IsClosed(contact)) ? offset_closed_contact(outer_bnd, contact, bsize)
	                                                 : offset_open_contact(outer_bnd, contact, bsize);
	//remove subtrees which have no source segment
	if (ret.nodes.size()<2) return ret;
	auto tt = Contour::Tree::CropLevel01(ret);
	if (tt.size() < 2) return ret;

	for (auto& t: tt){
		bool isgood=false;
		for (auto& e: t.alledges()){
			if (is_on_source(e->center()) != 0){
				isgood = true;
				break;
			}
		}
		if (isgood) continue;

		//remove t nodes from ret
		for (auto& n: t.nodes){
			auto e0 = n->contour[0];
			auto fnd = ret.find_node(e0.get());
			assert(fnd != nullptr);
			ret.remove_contour(fnd.get());
		}
	}

	return ret;
}

Contour::Tree BufferGrid::build_buffer_zone(const Contour::Tree& outer_bnd, double buffer_size) const{
	// get outstree edges which belong to source
	// assemble simple contours out of these edges
	vector<EdgeData> contacts = build_contacts(outer_bnd.alledges());

	// offset all contacts, leave only area lying within the original grid.
	// each contact line should give only one bounding contour.
	vector<Contour::Tree> contact_zones;
	for (auto& c: contacts) contact_zones.push_back(offset_contact(outer_bnd, c, buffer_size));

	// Unite all buffer zones
	if (contact_zones.size() == 0) return Contour::Tree();
	Contour::Tree buffer_zone(std::move(contact_zones[0]));
	for (int i=1; i<contact_zones.size(); ++i){
		buffer_zone = Contour::Clip::Union(buffer_zone, contact_zones[i]);
	}
	return buffer_zone;
}
BufferGrid::BufferGrid(GridData& main, const EdgeData& source, double buffer_size, bool preserve_bp, double angle0):
		angle0(angle0), preserve_bp(preserve_bp), source(&source), orig(&main){
	//1) throw away cells lying within the source
	CellData outs = Grid::ExtractCells(main, source, (Contour::Area(source)>0) ? OUTSIDE : INSIDE);
	Grid::Constructor::InvokeGrid inv(outs);

	//2) build buffer grid
	Contour::Tree outsbnd = Contour::Tree::GridBoundary(inv.grid);
	Contour::R::RevertTree bndrev(outsbnd);
	Contour::Tree buffer_zone = build_buffer_zone(outsbnd, buffer_size);
	if (buffer_zone.nodes.size() == 0) return;

	//3) get cells lying inside buffer
	auto outs2 = Grid::ExtractCells(inv.grid, buffer_zone, OUTSIDE);
	aa::constant_ids_pvec(inv.grid.vcells, 0);
	aa::constant_ids_pvec(outs2, 1);
	for (auto c: inv.grid.vcells) if (c->id != 1) push_back(c);
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
		if (isOnSection(p, *ed->first(), *ed->last(), ksi)) return 2;
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

namespace{
void remove_vertex(Vertex* v, GridData& from, HM2D::EdgeData& bnded){
	//do not remove superfluous edges but removes its points connectivity
	Edge *e1 = 0, *e2 = 0;
	auto it = bnded.begin();
	while (it != bnded.end()){
		if ((*it)->pfirst() == v || (*it)->plast() == v){
			e1 = it->get();
			break;
		}
		++it;
	}
	++it;
	while (it != bnded.end()){
		if ((*it)->pfirst() == v || (*it)->plast() == v){
			e2 = it->get();
			break;
		}
		++it;
	}
	assert(e1 != 0 && e2 != 0);
	//find cells
	//e1 and e2 should have same adjacent cells
	if (e1->no_left_cell()) e1->reverse();
	if (e2->no_left_cell()) e2->reverse();
	Cell *c1 = e1->left.lock().get();
	assert(e2->left.lock().get() == c1);

	//build new edge
	shared_ptr<Edge> newe(new Edge(*e1));
	if (e1->plast() == e2->pfirst()) newe->vertices[1] = e2->last();
	else if (e1->plast() == e2->plast()) newe->vertices[1] = e2->first();
	else if (e1->pfirst() == e2->pfirst()) newe->vertices[0] = e2->last();
	else newe->vertices[0] = e2->first();

	//adjust input data so that they could not be found in next
	//function call
	bnded.push_back(newe);
	e1->vertices[0].reset(); e1->vertices[1].reset();
	e2->vertices[0].reset(); e2->vertices[1].reset();

	//cell->edge connectivity
	int i0 = aa::shpvec_ifind(c1->edges, e1);
	int i1 = aa::shpvec_ifind(c1->edges, e2);
	c1->edges[i0] = newe;
	c1->edges.erase(c1->edges.begin()+i1);
}

}

//namespace {
//void rr2(Vertex* v, GridData* orig){
//        auto bnded = ECol::Assembler::GridBoundary(*orig);
//        aa::enumerate_ids_pvec(orig->vedges);
//        //find edges
//        Edge *e1 = 0, *e2 = 0;
//        auto it = bnded.begin();
//        while (it != bnded.end()){
//                if ((*it)->pfirst() == v || (*it)->plast() == v){
//                        e1 = it->get();
//                        break;
//                }
//                ++it;
//        }
//        ++it;
//        while (it != bnded.end()){
//                if ((*it)->pfirst() == v || (*it)->plast() == v){
//                        e2 = it->get();
//                        break;
//                }
//                ++it;
//        }
//        assert(e1 != 0 && e2 != 0);
//        //find cells
//        //e1 and e2 should have same adjacent cells
//        if (e1->no_left_cell()) e1->reverse();
//        if (e2->no_left_cell()) e2->reverse();
//        Cell *c1 = e1->left.lock().get();
//        assert(e2->left.lock().get() == c1);

//        //build new edge
//        shared_ptr<Edge> newe(new Edge(*e1));
//        if (e1->plast() == e2->pfirst()) newe->vertices[1] = e2->last();
//        else if (e1->plast() == e2->plast()) newe->vertices[1] = e2->first();
//        else if (e1->pfirst() == e2->pfirst()) newe->vertices[0] = e2->last();
//        else newe->vertices[0] = e2->first();

//        //cell->edge connectivity
//        int i0 = aa::shpvec_ifind(c1->edges, e1);
//        int i1 = aa::shpvec_ifind(c1->edges, e2);
//        c1->edges[i0] = newe;
//        c1->edges.erase(c1->edges.begin()+i1);
	
//        //remove edges from vedges
//        orig->vedges[e1->id] = newe;
//        orig->vedges.erase(orig->vedges.begin()+e2->id);
//}
//}

void BufferGrid::remove_vertices_from_orig(const VertexData& vd){
	auto bnded = ECol::Assembler::GridBoundary(*orig);
	int bndedlen = bnded.size();
	for (auto& v: vd) remove_vertex(v.get(), *orig, bnded);

	//adjust orit->vedges by removing edges with no point connectivity
	//and by addition of newly created edges passed to bnded list.
	orig->vedges.insert(orig->vedges.end(),
		bnded.begin() + bndedlen, bnded.end());
	auto rr = std::remove_if(orig->vedges.begin(), orig->vedges.end(),
		[](const shared_ptr<HM2D::Edge>& e){ return e->vertices[0] == nullptr; });
	orig->vedges.resize(rr - orig->vedges.begin());
}
void BufferGrid::rebuild_source_edges(EdgeData& cont, EdgeData& sed) {
	VertexData badpnt;
	vector<EdgeData> simpcont = Contour::Assembler::SimpleContours(sed);
	for (auto& sc: simpcont){
		auto op = Contour::OrderedPoints(sc);
		//ignore first last points for open contours;
		int imin=0, imax=op.size()-1;
		if (op.back() != op[0]) imin = 1;
		for (int i=imin; i<imax; ++i){
			if (is_on_source(*op[i]) != 1) badpnt.push_back(op[i]);
		}
	}
	//remove from cont by index
	aa::enumerate_ids_pvec(Contour::OrderedPoints1(cont));
	vector<int> remove_vert = aa::get_ids(badpnt);
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

	//save bad points to remove it from original grid prior to MergeTo execution
	badpoints.insert(badpoints.end(), badpnt.begin(), badpnt.end());
}

BufferGrid::splitR BufferGrid::split_edges(EdgeData& cont){
	splitR ret;
	//define source edges
	ret.bsource = define_source_edges(cont);

	//define grid edges
	aa::constant_ids_pvec(cont, 0);
	aa::constant_ids_pvec(ret.bsource, 1);
	for (auto e: cont) if (e->id == 0 && !e->is_boundary()){
		e->id = 2;
		ret.bgrid.push_back(e);
	}

	//define boundary edges
	for (auto e: cont) if (e->id == 0) ret.bboundary.push_back(e);

	//rebuild source edges throwing away superfluous nodes
	rebuild_source_edges(cont, ret.bsource);

	//keep boundary
	ret.keep_points = significant_boundary_points(ret.bboundary);

	return ret;
}

void BufferGrid::contour_segmentation(EdgeData& cont, const splitR& split,
			const vector<std::pair<Point, double>>& src){
	//do segmentation using id=1 to keep needed primitives
	aa::constant_ids_pvec(AllVertices(cont), 0);
	aa::constant_ids_pvec(cont, 0);
	aa::constant_ids_pvec(split.bsource, 1);
	aa::constant_ids_pvec(split.bgrid, 1);
	aa::constant_ids_pvec(split.bboundary, 0);
	aa::constant_ids_pvec(split.keep_points, 1);
	cont = Mesher::RepartSourceById(cont, src);
}

void BufferGrid::build_size_sources(const splitR& s, vector<std::pair<Point, double>>& addto){
	for (auto& e: s.bsource){
		addto.emplace_back(e->center(), e->length());
	}
	for (auto& e: s.bgrid){
		addto.emplace_back(e->center(), e->length());
	}
}

Contour::Tree BufferGrid::triangulation_boundary(){
	//get buffer cells boundary
	auto ginv = std::make_shared<Constructor::InvokeGrid>(*this);
	vector<EdgeData> conts = Contour::Assembler::GridBoundary(ginv->grid);
	ginv.reset();

	//split contour
	vector<splitR> splitted(conts.size());
	for (int i=0; i<conts.size(); ++i)
		splitted[i] = split_edges(conts[i]);

	//build size sources
	vector<std::pair<Point, double>> size_sources;
	for (int i=0; i<splitted.size(); ++i){
		build_size_sources(splitted[i], size_sources);
	}

	//partition
	for (int i=0; i<conts.size(); ++i)
		contour_segmentation(conts[i], splitted[i], size_sources);

	//assemble tree
	Contour::Tree ret;
	for (auto& c: conts) { ret.add_contour(std::move(c)); }

	return ret;
}

//0-triangle, 1-recombined triangle
void BufferGrid::update_original(int filler){
	if (size()==0) return;
	Contour::Tree tree = triangulation_boundary();
	GridData g3 = (filler == 0) ? Mesher::UnstructuredTriangle(tree)
	                            : Mesher::UnstructuredTriangleRecomb(tree);
	aa::constant_ids_pvec(orig->vcells, 0);
	aa::constant_ids_pvec(*this, 1);
	Algos::RemoveCellsById(*orig, 1);
	remove_vertices_from_orig(badpoints);
	Algos::MergeTo(g3, *orig);
}
