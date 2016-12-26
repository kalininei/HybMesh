#include "pyramid_layer.hpp"
#include "debug_grid3d.hpp"

using namespace HMGrid3D;
namespace{

CellData single_face_cells(const FaceData& fc){
	CellData ret;
	for (auto& f: fc){
		ret.emplace_back(new Cell());
		ret.back()->faces.push_back(f);
	}
	return ret;
}

struct PyrConstructor{
	PyrConstructor(const FaceData& limits){
	}

	void build_pyramid(shared_ptr<Cell>& c){
		auto fc = c->faces[0];
		auto ap = fc->sorted_vertices();

		Point3 cnt(0, 0, 0);
		double area = 0.;
		Vect3 rnrm = right_normal(*ap[0], *ap[1], *ap[2]);
		for (int i=1; i<ap.size() - 1; ++i){
			Point3 c1 = (*ap[0] + *ap[i] + *ap[i+1]) / 3.;
			double a1 = 0.5 * vecDot(vecCross(*ap[i] - *ap[0], *ap[i+1] - *ap[0]), rnrm);
			area += a1;
			cnt += c1 * a1;
		}
		assert(area > 0);
		cnt /= area;
		//double len = 1./sqrt(2) * sqrt(area);
		double len = 0.5 * sqrt(area);
		rnrm *= (len);

		shared_ptr<Vertex> vert(new Vertex(cnt - rnrm));
		EdgeData edges;
		for (auto& v: ap) edges.emplace_back(new Edge(v, vert));
		//build pyramid faces
		for (int n=0; n<edges.size(); ++n){
			int nnext = (n == edges.size() - 1) ? 0 : n+1;
			c->faces.emplace_back(new Face());
			c->faces.back()->edges.resize(3);
			c->faces.back()->edges[0] = edges[nnext];
			c->faces.back()->edges[1] = edges[n];
			c->faces.back()->edges[2] = c->faces[0]->edges[n];
		}
	
		//connectivity
		c->faces[0]->left = c;
		for (int i=1; i<c->faces.size(); ++i) c->faces[i]->right=c;
	}
private:
	void no_cross_check(Cell& c){
		_DUMMY_FUN_;
	}
};


vector<std::pair<int, int>> edge_cells_table(CellData& cells){
	vector<std::pair<int, int>> ret;

	enumerate_ids_pvec(cells);
	for (auto c: cells)
	for (auto e: c->faces[0]->edges) e->id = -1;

	for (auto c: cells)
	for (auto e: c->faces[0]->edges){
		if (e->id == -1){
			e->id = ret.size();
			ret.push_back(std::make_pair(c->id, -1));
		} else {
			ret[e->id].second = c->id;
		}
	}

	return ret;
}

bool reentrant_base(const Face& f1, const Face& f2){
	auto cp1 = f1.sorted_vertices(), cp2 = f2.sorted_vertices();
	//right normal
	Vect3 plane = vecCross(*cp1[1] - *cp1[0], *cp1[2] - *cp1[0]);
	double d = -vecDot(plane, *cp1[0]);
	for (auto p: cp2){
		auto fnd = std::find(cp1.begin(), cp1.end(), p);
		if (fnd != cp1.end()) continue;
		double s = vecDot(plane, *p) + d;
		return ISEQGREATER(s, 0.);
	}
	//faces have equal nodes
	//angle is 2pi
	return true;
}

bool has_pyramid(const Cell& c1){
	return c1.faces.size() > 3;
}

shared_ptr<Vertex> pyramid_vertex(const Cell& c){
	return c.faces[1]->edges[0]->last();
}

std::pair<int, int> find_eq_edges(const Face& f1, const Face& f2){
	for (int i=0; i<f1.edges.size(); ++i)
	for (int j=0; j<f2.edges.size(); ++j){
		if (f1.edges[i] == f2.edges[j]) return std::make_pair(i, j);
	}
	assert(false);
}

double faces_angle(const Face& f1, const Face& f2){
	//if f2 lies below f1 hence we've got self cross situaltion
	//and have to apply merging procedure.
	//We guarantee it by returning zero angle
	if (reentrant_base(f1, f2)) return 0;
	auto cp1 = f1.sorted_vertices(), cp2 = f2.sorted_vertices();
	Vect3 n1 = left_normal(*cp1[0], *cp1[1], *cp1[2]);
	Vect3 n2 = left_normal(*cp2[0], *cp2[1], *cp2[2]);
	double a = acos(vecDot(n1, n2))/M_PI*180;
	return 180 - a;
}

bool need_merge(const Cell& c1, const Cell& c2, double merge_angle){
	auto eqe = find_eq_edges(*c1.faces[0], *c2.faces[0]);

	shared_ptr<Face> f1, f2;
	if (has_pyramid(c1) == false) f1 = c1.faces[0];
	else f1 = c1.faces[1 + eqe.first];
	if (has_pyramid(c2) == false) f2 = c2.faces[0];
	else f2 = c2.faces[1 + eqe.second];

	double angle = faces_angle(*f1, *f2);

	return angle < merge_angle;
}

struct PyrVertices{
	CellData& cd;
	VertexData verts;
	vector<vector<Point3>> verts_hist;

	PyrVertices(CellData& cd): cd(cd){
		verts_hist.resize(cd.size());
		verts.resize(cd.size());

		for (int i=0; i<cd.size(); ++i) if (has_pyramid(*cd[i])){
			verts[i] = pyramid_vertex(*cd[i]);
		}
	}

	void merge(int c1, int c2, PyrConstructor& pc){
		if (has_pyramid(*cd[c1]) == false) tripyramid(c1, c2, pc);
		if (has_pyramid(*cd[c2]) == false) tripyramid(c2, c1, pc);

		std::set<Point3> avp;
		if (verts_hist[c1].size() == 0) avp.insert(*verts[c1]);
		else avp.insert(verts_hist[c1].begin(), verts_hist[c1].end());
		if (verts_hist[c2].size() == 0) avp.insert(*verts[c2]);
		else avp.insert(verts_hist[c2].begin(), verts_hist[c2].end());
		verts_hist[c1] = vector<Point3>(avp.begin(), avp.end());
		verts_hist[c2] = verts_hist[c1];

		Point3 p(0,0,0);
		for (auto& it: avp) p += it;
		p/=avp.size();

		verts[c1]->set(p);
		verts[c2]->set(p);
		
		for (int i=1; i<cd[c1]->faces.size(); ++i)
		for (int j=1; j<cd[c2]->faces.size(); ++j){
			auto f1 = cd[c1]->faces[i];
			auto f2 = cd[c2]->faces[j];
			auto e1 = f1->edges[2];
			auto e2 = f2->edges[2];
			if (e1 == e2){
				f2->reverse();
				f1->left = f2->left;
				cd[c2]->change_face(f2, f1);
				verts[c2] = verts[c1];
				return;
			}
		}
		assert(false);
	}
private:
	void tripyramid(int c, int ivert, PyrConstructor& pc){
		pc.build_pyramid(cd[c]);
		verts[c] = pyramid_vertex(*cd[c]);
	}
};

void merge_pyramids(CellData& cells, double merge_angle, PyrConstructor& pc){
	PyrVertices pv(cells);

	//cell-edge-cell connections
	vector<std::pair<int, int>> cec = edge_cells_table(cells);

	//loop over all cell-edge-cell connections
	for (auto ec: cec){
		Cell& c1 = *cells[ec.first];
		Cell& c2 = *cells[ec.second];
		if (has_pyramid(c1) == false && has_pyramid(c2) == false) continue;

		Face& base1 = *c1.faces[0];
		Face& base2 = *c2.faces[0];
		if (reentrant_base(base1, base2)) continue;

		if (need_merge(c1, c2, merge_angle))
			pv.merge(ec.first, ec.second, pc);
	}
}

}


GridData HMGrid3D::BuildPyramidLayer(const FaceData& faces, bool non3only, double merge_angle){
	//initial commit: add single face cells
	CellData ac = single_face_cells(faces);

	//Pyramid constructor which will build pyramids
	//which will not cross input surface
	PyrConstructor constructor(faces);

	//build initial pyramids
	for (auto c: ac){
		if (non3only && c->faces[0]->edges.size() < 4) continue;
		constructor.build_pyramid(c);
	}

	//merge pyramids
	merge_pyramids(ac, merge_angle, constructor);

	//assemble grid
	auto prim = AllPrimitives(ac);
	GridData ret;
	ret.vcells = std::move(ac);
	ret.vfaces = std::move(std::get<2>(prim));
	ret.vedges = std::move(std::get<1>(prim));
	ret.vvert = std::move(std::get<0>(prim));

	return ret;
}
