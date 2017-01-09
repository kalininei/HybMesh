#include "gpc_core.hpp"
#include "string.h"
#include "algos.hpp"
#include "cont_assembler.hpp"

using namespace HM2D::Impl;
using namespace HM2D;

namespace{

void gpc_fill_vertex_list(gpc_vertex_list& lst, const EdgeData& inp){
	assert(Contour::IsClosed(inp));
	auto pts = Contour::OrderedPoints(inp);
	lst.vertex = (gpc_vertex*)malloc(inp.size()*sizeof(gpc_vertex));
	for (int i=0; i<inp.size(); ++i){
		lst.vertex[i].x = pts[i]->x;
		lst.vertex[i].y = pts[i]->y;
	}
	lst.num_vertices = inp.size();
}

void gpc_copy_polygon(gpc_polygon& to, const gpc_polygon& from){
	int nn = from.num_contours;
	to.num_contours = nn;
	//copy holes
	to.hole = (int*)malloc(sizeof(int)*nn);
	memcpy(to.hole, from.hole, sizeof(int)*nn);
	//copy vertex_lists
	to.contour = (gpc_vertex_list*)malloc(sizeof(gpc_vertex_list)*nn);
	memcpy(to.contour, from.contour, sizeof(gpc_vertex_list)*nn);
	//copy verticies
	for (int i=0; i<nn; ++i){
		int nc = from.contour[i].num_vertices;
		to.contour[i].vertex = (gpc_vertex*)malloc(nc*sizeof(gpc_vertex));
		memcpy(to.contour[i].vertex, from.contour[i].vertex, sizeof(gpc_vertex)*nc);
	}
}

}
GpcTree::GpcTree(const EdgeData& inp): poly{0, 0, 0}{
	assert(Contour::IsClosed(inp));
	if (inp.size()>2){
		poly.num_contours = 1;
		poly.hole = (int*)malloc(sizeof(int));
		poly.hole[0] = 0;
		poly.contour = (gpc_vertex_list*)malloc(sizeof(gpc_vertex_list));
		gpc_fill_vertex_list(poly.contour[0], inp);
	}
}

GpcTree::GpcTree(const Contour::Tree& inp): poly{0, 0, 0}{
	int nc = inp.nodes.size();
	poly.num_contours = nc;
	poly.hole = (int*)malloc(sizeof(int)*nc);
	poly.contour = (gpc_vertex_list*)malloc(sizeof(gpc_vertex_list)*nc);
	int i=0;
	for (auto node: inp.nodes){
		if (node->contour.size()<3) continue;
		//holes
		bool is_inner = true;
		auto upper = node.get();
		while (!upper->parent.expired()){
			upper = upper->parent.lock().get();
			is_inner = !is_inner;
		}
		poly.hole[i] = is_inner?0:1;
		//contours
		gpc_fill_vertex_list(poly.contour[i++], node->contour);
	}
}

GpcTree::GpcTree(const GpcTree& other): poly{0, 0, 0}{
	gpc_copy_polygon(poly, other.poly);
}

GpcTree::GpcTree(GpcTree&& other) noexcept{
	poly = other.poly;
	other.poly = {0, 0, 0};
}

GpcTree& GpcTree::operator=(const GpcTree& other){
	if (&other == this) return *this;
	gpc_free_polygon(&poly);
	gpc_copy_polygon(poly, other.poly);
	return *this;
}

GpcTree::~GpcTree(){
	gpc_free_polygon(&poly);
}


Contour::Tree GpcTree::ToContourTree() const{
	Contour::Tree res;
	for (int i=0; i<poly.num_contours; ++i){
		auto& cont = poly.contour[i];
		int jmax = cont.num_vertices;
		if (jmax<3) continue;
		VertexData pdata;
		//build vector of points
		for (int j=0; j<jmax; ++j){
			pdata.emplace_back(new Vertex(cont.vertex[j].x, cont.vertex[j].y));
		}
		//add edges to contour
		EdgeData ecol = Contour::Assembler::Contour1(pdata, true);
		int nump_before = pdata.size();
		ECol::Algos::MergePoints(ecol);
		if (ecol.size() < 3) continue; 
		if (AllVertices(ecol).size() != nump_before){
			ecol = ECol::Algos::NoCrosses(ecol);
			vector<EdgeData> ac = Contour::Assembler::AllContours(ecol);
			vector<EdgeData*> closed_ones;
			for (auto& c: ac) if (Contour::IsClosed(c)){
				closed_ones.push_back(&c);
			}
			ecol.clear();
			for (auto cc: closed_ones){
				ecol.insert(ecol.end(), cc->begin(), cc->end());
			}
			for (auto cc:closed_ones) res.add_contour(*cc);
		} else {
			assert(ecol.size()>0 && Contour::IsClosed(ecol));
			//add contour to result
			res.add_contour(ecol);
		}
	}
	return res;
}

GpcTree GpcTree::Union(const GpcTree& c1, const GpcTree& c2){
	GpcTree res;
	gpc_polygon* sub1 = const_cast<gpc_polygon*>(&c1.poly);
	gpc_polygon* sub2 = const_cast<gpc_polygon*>(&c2.poly);
	gpc_polygon_clip(GPC_UNION, sub1, sub2, &res.poly);
	return res;
}

GpcTree GpcTree::Intersect(const GpcTree& c1, const GpcTree& c2){
	GpcTree res;
	gpc_polygon* sub1 = const_cast<gpc_polygon*>(&c1.poly);
	gpc_polygon* sub2 = const_cast<gpc_polygon*>(&c2.poly);
	gpc_polygon_clip(GPC_INT, sub1, sub2, &res.poly);
	return res;
}

GpcTree GpcTree::Substract(const GpcTree& c1, const GpcTree& c2){
	GpcTree res;
	gpc_polygon* sub1 = const_cast<gpc_polygon*>(&c1.poly);
	gpc_polygon* sub2 = const_cast<gpc_polygon*>(&c2.poly);
	gpc_polygon_clip(GPC_DIFF, sub1, sub2, &res.poly);
	return res;
}

GpcTree GpcTree::Xor(const GpcTree& c1, const GpcTree& c2){
	GpcTree res;
	gpc_polygon* sub1 = const_cast<gpc_polygon*>(&c1.poly);
	gpc_polygon* sub2 = const_cast<gpc_polygon*>(&c2.poly);
	gpc_polygon_clip(GPC_XOR, sub1, sub2, &res.poly);
	return res;
}
