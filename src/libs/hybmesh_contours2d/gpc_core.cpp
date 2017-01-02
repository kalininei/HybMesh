#include "gpc_core.hpp"
#include "string.h"

using namespace HMCont2D::Impl;

namespace{

void gpc_fill_vertex_list(gpc_vertex_list& lst, const HMCont2D::Contour& inp){
	assert(inp.is_closed());
	auto pts = inp.ordered_points();
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
GpcTree::GpcTree(const HMCont2D::Contour& inp): poly{0, 0, 0}{
	assert(inp.is_closed());
	if (inp.size()>2){
		poly.num_contours = 1;
		poly.hole = (int*)malloc(sizeof(int));
		poly.hole[0] = 0;
		poly.contour = (gpc_vertex_list*)malloc(sizeof(gpc_vertex_list));
		gpc_fill_vertex_list(poly.contour[0], inp);
	}
}

GpcTree::GpcTree(const HMCont2D::ContourTree& inp): poly{0, 0, 0}{
	int nc = inp.nodes.size();
	poly.num_contours = nc;
	poly.hole = (int*)malloc(sizeof(int)*nc);
	poly.contour = (gpc_vertex_list*)malloc(sizeof(gpc_vertex_list)*nc);
	int i=0;
	for (auto node: inp.nodes){
		if (node->size()<3) continue;
		//holes
		bool is_inner = true;
		auto upper = node.get();
		while (upper->parent != 0){
			upper = upper->parent;
			is_inner = !is_inner;
		}
		poly.hole[i] = is_inner?0:1;
		//contours
		gpc_fill_vertex_list(poly.contour[i++], *node);
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


HMCont2D::Container<HMCont2D::ContourTree> GpcTree::ToContourTree() const{
	HMCont2D::Container<HMCont2D::ContourTree> res;
	for (int i=0; i<poly.num_contours; ++i){
		auto& cont = poly.contour[i];
		int jmax = cont.num_vertices;
		if (jmax<3) continue;
		/*
		auto hmcont = shared_ptr<HMCont2D::Contour>(new HMCont2D::Contour);
		//build vector of points
		ShpVector<Point> shpp; shpp.reserve(poly.contour[i].num_vertices);
		for (int j=0; j<jmax; ++j){
			aa::add_shared(shpp, Point(cont.vertex[j].x, cont.vertex[j].y));
		}
		//add points
		res.pdata.add_values(shpp);
		//add edges to contour
		for (int j=0; j<jmax-1; ++j){
			hmcont->add_value(shared_ptr<Edge>(new Edge(shpp[j].get(), shpp[j+1].get())));
		}
		hmcont->add_value(shared_ptr<Edge>(new Edge(shpp.back().get(), shpp[0].get())));
		//add contour to result
		res.AddContour(hmcont);
		*/
		HMCont2D::Container<HMCont2D::ECollection> ecol;
		//build vector of points
		for (int j=0; j<jmax; ++j){
			ecol.pdata.add_value(Point(cont.vertex[j].x, cont.vertex[j].y));
		}
		//add edges to contour
		for (int j=0; j<jmax-1; ++j){
			ecol.add_value(shared_ptr<Edge>(
				new Edge(ecol.pdata.pvalue(j), ecol.pdata.pvalue(j+1))));
		}
		ecol.add_value(shared_ptr<Edge>(new Edge(ecol.pdata.data.back().get(), ecol.pdata.data[0].get())));
		int nump_before = ecol.pdata.size();
		HMCont2D::Algos::MergePoints(ecol);
		HMCont2D::Algos::DeleteUnusedPoints(ecol);
		if (ecol.size() < 3) continue; 
		if (ecol.pdata.size() != nump_before){
			ecol = HMCont2D::Algos::NoCrosses(ecol);
			vector<HMCont2D::Contour> ac = HMCont2D::Assembler::AllContours(ecol);
			vector<HMCont2D::Contour*> closed_ones;
			for (auto& c: ac) if (c.is_closed()){
				closed_ones.push_back(&c);
			}
			ecol.data.clear();
			for (auto cc: closed_ones){
				ecol.data.insert(ecol.data.end(), cc->data.begin(), cc->data.end());
			}
			HMCont2D::Algos::DeleteUnusedPoints(ecol);
			res.pdata.add_values(ecol.pdata.data);
			for (auto cc:closed_ones) res.AddContour(*cc);
		} else {
			auto hmcont = shared_ptr<HMCont2D::Contour>(new HMCont2D::Contour);
			hmcont->data=std::move(ecol.data);
			res.pdata.add_values(ecol.pdata.data);
			assert(hmcont->size()>0 && hmcont->is_closed());
			//add contour to result
			res.AddContour(hmcont);
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
