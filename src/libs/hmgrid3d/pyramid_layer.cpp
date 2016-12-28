#include "pyramid_layer.hpp"
#include "debug_grid3d.hpp"
#include <unordered_map>

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
		limvertices = AllVertices(limits);
		BoundingBox3D bball(limvertices);
		double L = bball.maxlen() / 30.0;
		bbfinder.reset(new BoundingBox3DFinder(bball, L));
		enumerate_ids_pvec(limvertices);
		limfaces.resize(limits.size());
		for (int i=0; i<limits.size(); ++i){
			auto f = limits[i].get();
			auto av = f->sorted_vertices();
			limfaces[i].resize(av.size());
			for (int j=0; j<av.size(); ++j)
				limfaces[i][j] = av[j]->id;
			BoundingBox3D avbb(av);
			bbfinder->addentry(avbb);
			faceind[f] = i;
		}
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
		double len = 0.5 * sqrt(area);
		rnrm *= (len);
		shared_ptr<Vertex> vert(new Vertex(cnt - rnrm));
		no_cross_check(c->faces[0].get(), cnt, *vert);

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
	static constexpr double CROSSLIMIT = 3;
	shared_ptr<BoundingBox3DFinder> bbfinder;
	VertexData limvertices;
	vector<vector<int>> limfaces;
	std::unordered_map<const Face*, int> faceind;

	void no_cross_check(const Face* srcface, const Point3& start, Point3& end){
		Point3 plus = Point3::Weigh(start, end, CROSSLIMIT);
		BoundingBox3D bbox(start, plus);
		vector<int> cind = bbfinder->suspects(bbox);
		double ksimin = 1;
		double xke[3];
		int iface = faceind[srcface];
		for (int i: cind) if (i!=iface){
			vector<int>& nds = limfaces[i];
			auto& p0 = *limvertices[nds[0]];
			for (int j=1; j<nds.size()-1; ++j){
				bool tricross = segment_triangle_cross3d(
					start, plus, 
					p0, *limvertices[nds[j]], *limvertices[nds[j+1]],
					xke);
				if (tricross){
					//found cross
					ksimin = std::min(ksimin, xke[0]);
					break;
				} else if (!ISIN_NN(xke[0], 0, 1)){
					//start-end doesn't cross face plane
					break;
				}
			}
		}
		if (ISLOWER(ksimin, 1)){
			end.set(Point3::Weigh(start, plus, ksimin/CROSSLIMIT));
		}
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
	struct pyr_group{
		pyr_group(int c1, int c2, shared_ptr<Vertex> v1, shared_ptr<Vertex> v2){
			icells.insert(c1); icells.insert(c2);
			origverts.insert(*v1); origverts.insert(*v2);
			pyrverts.insert(v1); pyrverts.insert(v2);
			assign_vertex();
		}
		std::set<int> icells;
		std::set<Point3> origverts;
		std::set<shared_ptr<Vertex>> pyrverts;
		Point3 vertex;
		bool has_cell(int icell){
			return icells.find(icell) != icells.end();
		}
		void add_cell(int c1, shared_ptr<Vertex> v1){
			icells.insert(c1);
			origverts.insert(*v1);
			pyrverts.insert(v1);
			assign_vertex();
		}
		void merge_with(const pyr_group& g){
			icells.insert(g.icells.begin(), g.icells.end());
			origverts.insert(g.origverts.begin(), g.origverts.end());
			pyrverts.insert(g.pyrverts.begin(), g.pyrverts.end());
			assign_vertex();
		}
		void assign_vertex(){
			vertex.set(0, 0, 0);
			for (auto& ov: origverts) vertex += ov;
			vertex /= origverts.size();
			for (auto& pv: pyrverts) pv->set(vertex);
		}
	};
	std::list<pyr_group> groups;
	void add_to_group(int c1, int c2, shared_ptr<Vertex> v1, shared_ptr<Vertex> v2){
		std::list<pyr_group>::iterator fnd1=groups.end();
		std::list<pyr_group>::iterator fnd2=groups.end();
		for (auto it=groups.begin(); it!=groups.end(); ++it){
			if (it->has_cell(c1)) fnd1 = it;
			if (it->has_cell(c2)) fnd2 = it;
		}
		if (fnd1 == groups.end() && fnd2 == groups.end()){
			groups.emplace_back(c1, c2, v1, v2);
		} else if (fnd1 != groups.end() && fnd2 == groups.end()){
			fnd1->add_cell(c2, v2);
		} else if (fnd1 == groups.end() && fnd2 != groups.end()){
			fnd2->add_cell(c1, v1);
		} else if (fnd1 != fnd2){
			fnd1->merge_with(*fnd2);
			groups.erase(fnd2);
		}
	}

	PyrVertices(CellData& cd): cd(cd){}

	void merge(int c1, int c2, PyrConstructor& pc){
		//build pyramids if they are absent
		if (has_pyramid(*cd[c1]) == false) tripyramid(c1, c2, pc);
		if (has_pyramid(*cd[c2]) == false) tripyramid(c2, c1, pc);
		//find common faces and quit if there are none
		shared_ptr<Face> fc1, fc2;
		for (int i=1; i<cd[c1]->faces.size(); ++i)
		for (int j=1; j<cd[c2]->faces.size(); ++j){
			auto f1 = cd[c1]->faces[i];
			auto f2 = cd[c2]->faces[j];
			auto e1 = f1->edges[2];
			auto e2 = f2->edges[2];
			if (e1 == e2){
				fc1 = f1;
				fc2 = f2;
				break;
			}
		}
		if (fc1 == nullptr) return;

		//make vertex nodes equal. Do not split primitives.
		//It will be done in supplementary_merge procedure.
		add_to_group(c1, c2, pyramid_vertex(*cd[c1]), pyramid_vertex(*cd[c2]));
	}

	void supplementary_merge(PyrConstructor& pc){
		//guarantee that all links within the group were merged
		for (auto g: groups) merge_group(g.icells, g.vertex);
	}

private:
	void tripyramid(int c, int ivert, PyrConstructor& pc){
		pc.build_pyramid(cd[c]);
	}
	void merge_group(const std::set<int>& g, Point3 vertex){
		FaceData lfaces;
		for (auto i: g) lfaces.push_back(cd[i]->faces[0]);
		auto prim = AllPrimitives(lfaces);
		VertexData& lvert(std::get<0>(prim));
		EdgeData& ledges(std::get<1>(prim));
		enumerate_ids_pvec(lvert);
		enumerate_ids_pvec(ledges);
		shared_ptr<Vertex> pv(new Vertex(vertex));
		//vertical edges
		EdgeData vedges;
		for (auto v: lvert)
			vedges.emplace_back(new Edge(v, pv));
		//vertical faces
		FaceData vfaces;
		for (auto e: ledges){
			vfaces.emplace_back(new Face());
			auto v1 = e->first();
			auto v2 = e->last();
			vfaces.back()->edges.push_back(vedges[v1->id]);
			vfaces.back()->edges.push_back(vedges[v2->id]);
			vfaces.back()->edges.push_back(e);
		}
		//construct cells
		for (auto i: g){
			auto cell = cd[i];
			cell->faces.resize(1);
			for (int i=0; i<cell->faces[0]->edges.size(); ++i){
				auto e = cell->faces[0]->edges[i];
				cell->faces.push_back(vfaces[e->id]);
				if (cell->faces[0]->is_positive_edge(i))
					vfaces[e->id]->left = cell;
				else
					vfaces[e->id]->right = cell;
			}
		}
		//guarantee that all bnd faces have cell to its right
		for (auto f: vfaces) if (f->is_boundary()){
			if (f->has_left_cell()) f->reverse();
		}
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
	pv.supplementary_merge(pc);
}

}


GridData HMGrid3D::BuildPyramidLayer(const FaceData& faces, bool non3only, double merge_angle){
	//initial commit: add single face cells
	CellData ac = single_face_cells(faces);
	bool need_pyramids = !non3only;
	for (auto f: faces) if (f->edges.size() > 3){
		need_pyramids = true;
		break;
	}

	if (need_pyramids){
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
	}

	//assemble grid
	auto prim = AllPrimitives(ac);
	GridData ret;
	ret.vcells = std::move(ac);
	ret.vfaces = std::move(std::get<2>(prim));
	ret.vedges = std::move(std::get<1>(prim));
	ret.vvert = std::move(std::get<0>(prim));

	return ret;
}
