#include "tetramesh_preproc.hpp"
#include "vtk_export_grid3d.hpp"
#include "debug_grid3d.hpp"
#include "merge3d.hpp"
#include "pyramid_layer.hpp"
using namespace HMGrid3D::Mesher;

SurfacePreprocess::SurfacePreprocess(const HMGrid3D::SurfaceTree& tree, double split_angle){
	assert(
		//first node is root, all others are its children
		[&](){
			if (tree.nodes[0]->level != 0) return false;
			for (int i=1; i<tree.nodes.size(); ++i)
				if (tree.nodes[i]->level != 1) return false;
			return true;
		}()
	);
	allfaces1 = tree.allfaces();
	for (auto n: tree.nodes){
		decomposed_surfs1.push_back(HMGrid3D::Surface::ExtractSmooth(*n, split_angle));
	}
	assemble_bnd_grid();
	supplement_from_decomposed_surfs();
}

void SurfacePreprocess::assemble_bnd_grid(){
	//keep subsurface identifier in face id.
	//to restore split_bnd_grid vector
	int k = 0, ids=0, is=0;
	std::map<int, std::pair<int, int>> kmap;
	for (auto& ds: decomposed_surfs1){
		is = 0;
		for (auto& s: ds){
			kmap[k] = std::make_pair(ids, is);
			for (auto f: s.faces) f->id = k;
			++k;
			++is;
		}
		++ids;
	}
	
	//initial commit: add single face cells
	FaceData fd;
	DeepCopy(allfaces1, fd, 2);
	for (auto af: fd){ af->left.reset(); af->right.reset(); }
	shared_ptr<RestoreIds<FaceData>> resfd(new RestoreIds<FaceData>(fd));
	bnd_grid = HMGrid3D::BuildPyramidLayer(fd, true, 60);

	//split bnd grid
	resfd.reset();
	split_bnd_grid.resize(decomposed_surfs1.size());
	for (auto& c: bnd_grid.vcells){
		auto& ijc = kmap[c->faces[0]->id];
		auto& sgrid = split_bnd_grid[ijc.first];
		if (sgrid.size() <= ijc.second) sgrid.resize(ijc.second + 1);
		sgrid[ijc.second].vcells.push_back(c);
	}
	for (auto& it1: split_bnd_grid)
	for (auto& it2: it1){
		it2.vfaces = AllFaces(it2.vcells);
	}

	//assemble decomposed_surfs
	decomposed_surfs.resize(split_bnd_grid.size());
	for (int i=0; i<split_bnd_grid.size(); ++i){
		decomposed_surfs[i].resize(split_bnd_grid[i].size());
		for (int j=0; j<split_bnd_grid[i].size(); ++j){
			GridData& g = split_bnd_grid[i][j];
			for (auto& f: g.vfaces) if (!f->has_left_cell()){
				decomposed_surfs[i][j].faces.push_back(f);
			}
		}
	}
}

void SurfacePreprocess::supplement_from_decomposed_surfs(){
	//all edges and vertices
	for (auto& ds: decomposed_surfs)
	for (auto& s: ds) 
	for (auto& f: s.faces)
	for (auto& e: f->edges){
		e->id = 0;
		for (auto& v: e->vertices) v->id = 0;
	}
	for (auto& ds: decomposed_surfs)
	for (auto& s: ds) 
	for (auto& f: s.faces)
	for (auto& e: f->edges) if (e->id == 0){
		e->id = 1;
		ae.push_back(e);
		for (auto& v: e->vertices) if (v->id == 0){
			v->id = 1;
			av.push_back(v);
		}
	}
	//decomposed_edges
	for (auto& ds: decomposed_surfs){
		vector<vector<HMGrid3D::EdgeData>> de;
		for (auto& s: ds){
			vector<HMGrid3D::EdgeData> ds;
			Vect3 right_normal = s.faces[0]->left_normal()*(-1);
			auto ex = HMGrid3D::Surface::ExtractAllBoundaries(s, right_normal);
			assert(ex[0].size() == 1);
			ds.push_back(ex[0][0]);
			for (int i=0; i<ex[1].size(); ++i) ds.push_back(ex[1][i]);
			de.push_back(ds);
		}
		decomposed_edges.push_back(de);
	}
}

namespace{

bool vertsort(const Point3* a, const Point3* b){
	if (ISLOWER(a->x, b->x)) return true;
	if (ISGREATER(a->x, b->x)) return false;
	if (ISLOWER(a->y, b->y)) return true;
	if (ISGREATER(a->y, b->y)) return false;
	return ISLOWER(a->z, b->z);
}
std::pair<int, int> find_indicies(Point3* p1, Point3* p2,
		const vector<Point3>& v1, const vector<Point3>& v2){
	std::pair<int, int> ret(p1 - &v1[0], p2 - &v2[0]);
	if (ret.first <0 || ret.first >= v1.size()){
		ret.first = p2 - &v1[0];
		ret.second = p1 - &v2[0];
	}
	assert(ret.first >=0 && ret.second >=0 && ret.first < v1.size() && ret.second < v2.size());
	return ret;
}

}

void SurfacePreprocess::Restore(HMGrid3D::GridData& g, const VertexData& gvert, const VertexData& svert){
	enumerate_ids_pvec(g.vvert);
	enumerate_ids_pvec(bnd_grid.vvert);
	vector<int> gind(gvert.size()), sind(svert.size());
	for (int i=0; i<gvert.size(); ++i){
		gind[i] = gvert[i]->id;
		sind[i] = svert[i]->id;
		assert(gind[i]>=0 && gind[i]<=g.vvert.size());
		assert(sind[i]>=0 && sind[i]<=bnd_grid.vvert.size());
	}
	
	//merge with boundary grid with its btypes
	merge_with_bnd(g, gind, sind);
}

void SurfacePreprocess::merge_with_bnd(GridData& tar, const vector<int>& tar_points, const vector<int>& bnd_points){
	HMGrid3D::MergeGrid(bnd_grid, tar, bnd_points, tar_points);
}

