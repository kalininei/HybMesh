#include "merge3d.hpp"
#include "debug_grid3d.hpp"

using namespace HMGrid3D;

namespace{

bool ptrsort(Point3* a, Point3* b){
	if (ISLOWER(a->x, b->x)) return true;
	else if (ISLOWER(b->x, a->x)) return false;
	else if (ISLOWER(a->y, b->y)) return true;
	else if (ISLOWER(b->y, a->y)) return false;
	else return ISLOWER(a->z, b->z);
};

vector<int> boundary_edges(const GridData& g){
	FaceData bfaces;
	for (auto f: g.vfaces) if (f->is_boundary()){
		bfaces.push_back(f);
	}
	EdgeData bedges = AllEdges(bfaces);
	enumerate_ids_pvec(g.vedges);
	vector<int> ret(bedges.size());
	for (int i=0; i<ret.size(); ++i) ret[i] = bedges[i]->id;
	return ret;
}

void assemble_duplicate_edges(EdgeData& efrom, EdgeData& eto,
		const vector<int>& used_efrom, const vector<int>& used_eto,
		vector<int>& from, vector<int>& to){
	for (int i=0; i<used_efrom.size(); ++i){
		auto e = efrom[used_efrom[i]];
		if (e->first()->id != -1 && e->last()->id != -1){
			from.push_back(used_efrom[i]);
		}
	}
	for (int i=0; i<used_eto.size(); ++i){
		auto e = eto[used_eto[i]];
		if (e->first()->id != -1 && e->last()->id != -1){
			to.push_back(used_eto[i]);
		}
	}

	std::vector<std::pair<int, int>> vpairs1(efrom.size());
	for (int i=0; i<from.size(); ++i){
		vpairs1[from[i]].first = efrom[from[i]]->first()->id;
		vpairs1[from[i]].second = efrom[from[i]]->last()->id;
		if (vpairs1[from[i]].first > vpairs1[from[i]].second)
			std::swap(vpairs1[from[i]].first, vpairs1[from[i]].second);
	}
	std::sort(from.begin(), from.end(),
		[&vpairs1](int a, int b)->bool{
			return vpairs1[a] < vpairs1[b];
		});

	std::vector<std::pair<int, int>> vpairs2(eto.size());
	for (int i=0; i<to.size(); ++i){
		vpairs2[to[i]].first = eto[to[i]]->first()->id;
		vpairs2[to[i]].second = eto[to[i]]->last()->id;
		if (vpairs2[to[i]].first > vpairs2[to[i]].second)
			std::swap(vpairs2[to[i]].first, vpairs2[to[i]].second);
	}
	std::sort(to.begin(), to.end(),
		[&vpairs2](int a, int b)->bool{
			return vpairs2[a] < vpairs2[b];
		});

	vector<int> to2, from2;
	vector<int>::iterator itto = to.begin(), itfrom = from.begin();
	while (itto!=to.end() && itfrom!=from.end()){
		if (vpairs1[*itfrom] < vpairs2[*itto]) ++itfrom;
		else if (vpairs1[*itfrom] > vpairs2[*itto]) ++itto;
		else{
			to2.push_back(*itto++);
			from2.push_back(*itfrom++);
		}
	}
	std::swap(to, to2);
	std::swap(from, from2);
}

void assemble_duplicate_faces(FaceData& ffrom, FaceData& fto, vector<int>& from, vector<int>& to){
	vector<int> from_candidates, to_candidates;
	for (int i=0; i<ffrom.size(); ++i) if (ffrom[i]->is_boundary()){
		if (all_of(ffrom[i]->edges.begin(), ffrom[i]->edges.end(),
				[](shared_ptr<Edge> e){ return e->id != -1; })){
			from_candidates.push_back(i);
		}
	}
	for (int i=0; i<fto.size(); ++i) if (fto[i]->is_boundary()){
		if (all_of(fto[i]->edges.begin(), fto[i]->edges.end(),
				[](shared_ptr<Edge> e){ return e->id != -1; })){
			to_candidates.push_back(i);
		}
	}
	//Candidate faces contain only duplicate edges.
	//However check for edges is not enough to filter out needed faces,
	//we have to compare all faces explicitly.
	//We do it using comparison of the sum of all face points.
	vector<Point3> from_candidates_cnt(from_candidates.size(), Point3(0, 0, 0));
	vector<Point3> to_candidates_cnt(to_candidates.size(), Point3(0, 0, 0));
	vector<Point3*> fromptr(from_candidates.size());
	vector<Point3*> toptr(to_candidates.size());
	for (int i=0; i<from_candidates.size(); ++i){
		fromptr[i] = &from_candidates_cnt[i];
		auto f = ffrom[from_candidates[i]];
		for (auto e: f->edges){
			from_candidates_cnt[i] += *e->first();
			from_candidates_cnt[i] += *e->last();
		}
	}
	for (int i=0; i<to_candidates.size(); ++i){
		toptr[i] = &to_candidates_cnt[i];
		auto f = fto[to_candidates[i]];
		for (auto e: f->edges){
			to_candidates_cnt[i] += *e->first();
			to_candidates_cnt[i] += *e->last();
		}
	}
	std::sort(fromptr.begin(), fromptr.end(), ptrsort);
	std::sort(toptr.begin(), toptr.end(), ptrsort);

	auto itfrom = fromptr.begin(), itto = toptr.begin();
	while (itfrom != fromptr.end() && itto != toptr.end()){
		if (ptrsort(*itfrom, *itto)) ++itfrom;
		else if (ptrsort(*itto, *itfrom)) ++itto;
		else{
			int ind1 = *itfrom - &from_candidates_cnt[0];
			from.push_back(from_candidates[ind1]);
			int ind2 = *itto - &to_candidates_cnt[0];
			to.push_back(to_candidates[ind2]);
			++itto; ++itfrom;
		}
	}
}

}

void HMGrid3D::MergeGrid(GridData& from, GridData& to,
		const vector<int>& from_vert, const vector<int>& to_vert){
	assert(
		//check points equality
		[&]()->bool{
			if (from_vert.size() != to_vert.size()) return false;
			for (int i=0; i<from_vert.size(); ++i){
				auto p1 = from.vvert[from_vert[i]];
				auto p2 = to.vvert[to_vert[i]];
				if (*p1 != *p2) return false;
			}
			return true;
		}()
	);
	//if (from_vert.size() == 0) return;
	vector<int> bnd_from_edges = boundary_edges(from);
	vector<int> bnd_to_edges = boundary_edges(to);
	//primitive id stores index of duplicate primitive
	//in another grid or -1 if there is no duplicate one.
	constant_ids_pvec(to.vvert, -1);
	constant_ids_pvec(to.vedges, -1);
	constant_ids_pvec(to.vfaces, -1);
	constant_ids_pvec(from.vvert, -1);
	constant_ids_pvec(from.vedges, -1);
	constant_ids_pvec(from.vfaces, -1);

	for (int i=0; i<from_vert.size(); ++i){
		to.vvert[to_vert[i]]->id = from_vert[i];
		from.vvert[from_vert[i]]->id = to_vert[i];
	}

	//1. change edges vertices
	for (auto& e: to.vedges){
		if (e->first()->id != -1) e->vertices[0] = from.vvert[e->first()->id];
		if (e->last()->id != -1) e->vertices.back() = from.vvert[e->last()->id];
	}
	vector<int> from_edges, to_edges;
	assemble_duplicate_edges(from.vedges, to.vedges,
			bnd_from_edges, bnd_to_edges,
			from_edges, to_edges);
	for (int i=0; i<from_edges.size(); ++i){
		to.vedges[to_edges[i]]->id = from_edges[i];
		from.vedges[from_edges[i]]->id = to_edges[i];
	}

	//2. change faces edges
	for (auto& f: to.vfaces){
		for (int i=0; i<f->edges.size(); ++i){
			if (f->edges[i]->id != -1){
				f->edges[i] = from.vedges[f->edges[i]->id];
			}
		}
	}
	vector<int> from_faces, to_faces;
	assemble_duplicate_faces(from.vfaces, to.vfaces, from_faces, to_faces);
	for (int i=0; i<from_faces.size(); ++i){
		to.vfaces[to_faces[i]]->id = from_faces[i];
		from.vfaces[from_faces[i]]->id = to_faces[i];
	}

	//3. change face->cells connectivity
	for (int i=0; i<from_faces.size(); ++i){
		auto ffrom = from.vfaces[from_faces[i]];
		auto fto = to.vfaces[to_faces[i]];
		assert(ffrom->is_boundary());
		assert(fto->is_boundary());
		auto cto = (fto->has_left_cell()) ? fto->left.lock() : fto->right.lock();
		ffrom->has_left_cell() ? ffrom->right = cto : ffrom->left = cto;
		auto ind = std::find(cto->faces.begin(), cto->faces.end(), fto) - cto->faces.begin();
		assert(ind < cto->faces.size());
		cto->faces[ind] = ffrom;
	}

	//4. change collection entries
	//cells
	for (auto c: from.vcells) if (c->faces.size() > 1){
		to.vcells.push_back(c);
	}
	//faces
	for (auto f: from.vfaces){
		if (f->id == -1) to.vfaces.push_back(f);
		else to.vfaces[f->id] = f;
	}
	//edges
	for (auto e: from.vedges){
		if (e->id == -1) to.vedges.push_back(e);
		else to.vedges[e->id] = e;
	}
	//vertices
	for (auto v: from.vvert){
		if (v->id == -1) to.vvert.push_back(v);
		else to.vvert[v->id] = v;
	}
}

GridData HMGrid3D::MergeGrids(const GridData& _g1, const GridData& _g2){
	//make deep copies of g1, g2
	GridData g1, g2;
	DeepCopy(_g1, g1);
	DeepCopy(_g2, g2);

	//assemble coincident points
	auto srf1 = Surface::GridSurface(g1);
	auto srf2 = Surface::GridSurface(g2);
	auto pc1 = AllVertices(srf1.faces);
	auto pc2 = AllVertices(srf2.faces);

	enumerate_ids_pvec(g1.vvert);
	enumerate_ids_pvec(g2.vvert);

	vector<Point3*> pp1(pc1.size());
	for (int i=0; i<pc1.size(); ++i) pp1[i] = pc1[i].get();
	std::sort(pp1.begin(), pp1.end(), ptrsort);

	vector<int> vfrom, vto;
	for (auto& p: pc2){
		auto fnd = std::lower_bound(pp1.begin(), pp1.end(), p.get(), ptrsort);
		if (fnd!=pp1.end() && !ptrsort(p.get(), *fnd)){
			vfrom.push_back(static_cast<Vertex*>(*fnd)->id);
			vto.push_back(p->id);
		}
	}

	MergeGrid(g1, g2, vfrom, vto);
	return g2;
}
