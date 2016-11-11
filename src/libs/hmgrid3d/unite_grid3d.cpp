#include "unite_grid3d.hpp"
#include "surface_grid3d.hpp"
#include "tetrahedral.hpp"
#include "hmgrid3dgeom.hpp"
#include "debug_grid3d.hpp"
using namespace HMGrid3D;

namespace{

GridData FilterCells(GridData& base, SurfaceTree& tree, double buf){
	if (tree.nodes.size()!=1){
		//implementation works only for first found root
		if (tree.nodes[0]->level != 0) _THROW_NOT_IMP_;
	}
	std::list<int> in;
	for (int i=0; i<base.vcells.size(); ++i) in.push_back(i);
	_UNCOMPLETE_FUN_("filter cells using bounding boxes");

	vector<BoundingBox3D> cellsbb;
	for (auto c: base.vcells) cellsbb.push_back(BoundingBox3D(*c));
	BoundingBox3D treebb = BoundingBox3D(*tree.nodes[0]);

	//1) check against enhanced bounding box
	BoundingBox3D treebb_plus = treebb.widen(buf);
	for (auto it=in.begin(); it!=in.end();){
		int ind = *it;
		if (treebb_plus.position(cellsbb[ind]) == OUTSIDE){
			it = in.erase(it);
		} else ++it;
	}

	//2) check against true bounding box
	for (auto it=in.begin(); it!=in.end();){
		int ind = *it;
		if (treebb.position(cellsbb[ind]) == OUTSIDE &&
				ISGREATER(treebb.meas(cellsbb[ind]), buf*buf)){
			it = in.erase(it);
		} else ++it;
	}

	//3) build a sphere inside tree and check if all vertices lie there.
	BoundingSphere3D treesph(*tree.nodes[0]);
	for (auto it=in.begin(); it!=in.end();){
		int ind = *it;
		if (treesph.position(cellsbb[ind]) == OUTSIDE &&
				ISGREATER(treesph.meas(cellsbb[ind]), buf*buf)){
			it = in.erase(it);
		} else ++it;
	}

	//TODO: proceed with check or write good algorithm

	//filter if cellsbb 
	GridData ret;
	for (auto i: in){
		ret.vcells.push_back(base.vcells[i]);
	}
	ret.fill_from_cells();
	return ret;
}

struct FaceCellAdopt{
	ShpVector<HMGrid3D::Face> used_faces;
	WpVector<HMGrid3D::Cell> old_cells;
	FaceCellAdopt(GridData& actual, GridData& old){
		constant_ids_pvec(old.vcells, -1);
		constant_ids_pvec(actual.vcells, 1);
		for (auto f: actual.vfaces){
			if (f->left.lock()->id != f->right.lock()->id){
				used_faces.push_back(f);
			}
		}
		for (auto f: used_faces){
			if (f->left.lock()->id == -1){
				old_cells.push_back(f->left);
				f->left.reset();
			} else {
				old_cells.push_back(f->right);
				f->right.reset();
			}
		}
	}
	void restore(){
		auto it = old_cells.begin();
		for (auto f: used_faces){
			if (f->left.expired()) f->left = *it++;
			else f->right = *it++;
		}
	}
};

void merge_all(GridData& ret){
	_UNCOMPLETE_FUN_("this should be done using gmsh data");
	ShpVector<Face> allbndfaces;
	for (auto f: ret.vfaces) if (f->is_boundary()) allbndfaces.push_back(f);
	vector<ShpVector<Vertex>> doubled_vertices;
	Surface s; s.faces = allbndfaces;
	auto _abv = s.allvertices();
	std::list<shared_ptr<Vertex>> allbndvert(_abv.begin(), _abv.end());
	for (auto it=allbndvert.begin(); it!=allbndvert.end(); ++it){
		ShpVector<Vertex> forit {*it};
		for (auto it2 = std::next(it); it2!=allbndvert.end();){
			if (ISEQ((*it)->x, (*it2)->x) && ISEQ((*it)->y, (*it2)->y) && ISEQ((*it)->z, (*it2)->z)){
				forit.push_back(*it2);
				it2 = allbndvert.erase(it2);
			} else ++it2;
		}
		if (forit.size()>1) doubled_vertices.push_back(forit);
	}
	constant_ids_pvec(ret.vvert, -1);
	for (int i=0; i<doubled_vertices.size(); ++i){
		for (int j=0; j<doubled_vertices[i].size(); ++j){
			doubled_vertices[i][j]->id = i;
		}
	}

	std::list<shared_ptr<Face>> suspects;
	auto is_single_face = [](shared_ptr<Face> f){
		return f->edges[0]->vertices[0]->id == -1;
	};
	auto is_doubled = [&doubled_vertices](shared_ptr<Face> f1, shared_ptr<Face> f2)->bool{
		auto sv1 = f1->sorted_vertices();
		auto sv2 = f2->sorted_vertices();
		for (auto p: sv1){
			int tp = p->id;
			assert(tp>-1);
			for (auto p2: sv2){
				if (p2->id == tp) { tp = -1; break; }
			}
			if (tp != -1) return false;
		}
		return true;
	};

	for (int i=0; i<allbndfaces.size(); ++i){
		if (!is_single_face(allbndfaces[i])) suspects.push_back(allbndfaces[i]);
	}
	//fill doubled faces
	vector<ShpVector<Face>> doubled_faces;
	while (suspects.size() > 0){
		auto f1 = suspects.back();
		suspects.pop_back();
		doubled_faces.push_back({f1});
		for (auto it = suspects.begin(); it!=suspects.end(); ++it){
			if (is_doubled(f1, *it)){
				doubled_faces.back().push_back(*it);
				suspects.erase(it);
				break;
			}
		}
	}
	for (auto& dv: doubled_vertices) dv[0]->id = -1;
	for (auto& v: doubled_faces){
		if (v[0]->edges[0]->first()->id != -1) std::swap(v[0], v[1]);
	}

	constant_ids_pvec(ret.vfaces, -1);
	for (int i=0; i<doubled_faces.size(); ++i){
		for (int j=1; j<doubled_faces[i].size(); ++j){
			doubled_faces[i][j]->id = i;
		}
	}
	//remove faces
	for (auto& c: ret.vcells){
		for (auto& f: c->faces){
			if (f->id>=0){
				auto f0 = doubled_faces[f->id][0];
				f = f0;
				if (f0->left.expired()) f0->left = c;
				else if (f0->right.expired()) f0->right = c;
				else assert(false);
			}
		}
	}
	//remove edges
	s.faces.clear(); 
	for (auto& it: doubled_faces) s.faces.push_back(it[0]);
	ShpVector<Edge> good_edges = s.alledges();


	auto find_good_edge = [&](shared_ptr<Edge> be)->shared_ptr<Edge>{
		auto shouldhave1 = doubled_vertices[be->first()->id][0];
		auto shouldhave2 = doubled_vertices[be->last()->id][0];

		for (auto ge: good_edges){
			if (ge->first() == shouldhave1 && ge->last() == shouldhave2) return ge;
			if (ge->last() == shouldhave1 && ge->first() == shouldhave2) return ge;
		}

		return shared_ptr<Edge>();
	};

	for (auto c: ret.vcells)
	for (auto f: c->faces)
	for (auto& e: f->edges)if (e->first()->id>=0 && e->last()->id>=0){
		auto en = find_good_edge(e);
		if (en) e = en;
	}

	//remove vertices
	for (auto c: ret.vcells)
	for (auto f: c->faces)
	for (auto e: f->edges)
	for (auto& v: e->vertices){
		if (v->id >= 0) v = doubled_vertices[v->id][0];
	}
	
	ret.fill_from_cells();
}

}


GridData Unite::Superimpose(GridData& base, GridData& imp, double buf){
	//=== assemble base cells which have to be deleted
	SurfaceTree imptree = Surface::BuildTree(imp, -1);
	GridData del = FilterCells(base, imptree, buf);
	
	//=== build buffer zone
	//temporary remove links to base cells from del grid to assemble deltree
	//and then restore
	FaceCellAdopt fca(del, base);

	SurfaceTree deltree = Surface::BuildTree(del, -1);
	fca.restore();
	SurfaceTree bzone = SurfaceTree::Substract(deltree, imptree);
	
	//=== unstructured fill
	GridData bgrid = Mesher::UnstructuredTetrahedral(bzone);
	
	//=== gather
	GridData ret;
	//base: copy cells and adopt faces for merging procedure
	constant_ids_pvec(base.vcells, 0);
	constant_ids_pvec(del.vcells, 1);
	for (int i=0; i<base.vcells.size(); ++i){
		if (base.vcells[i]->id == 0){
			ret.vcells.push_back(base.vcells[i]);
			for (auto f: base.vcells[i]->faces){
				if (!f->left.expired() && f->left.lock()->id == 1) f->left.reset();
				if (!f->right.expired() && f->right.lock()->id == 1) f->right.reset();
			}
		}
	}
	//imp
	std::copy(imp.vcells.begin(), imp.vcells.end(), std::back_inserter(ret.vcells));
	//buffer
	std::copy(bgrid.vcells.begin(), bgrid.vcells.end(), std::back_inserter(ret.vcells));

	//assemble
	ret.fill_from_cells();

	//merge doubled entities: very slow and is needed only if
	//resulting grid is not final
	if (HMDebug::get().i == 1) merge_all(ret);

	return ret;
}
