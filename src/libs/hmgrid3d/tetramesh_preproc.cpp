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

	//decomposed_surfs1 stores decomposition surfaces for all tree node surfaces
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
	
	//copy faces (is it really needed?) 
	FaceData fd;
	DeepCopy(allfaces1, fd, 2);
	for (auto af: fd){ af->left.reset(); af->right.reset(); }
	//Save faces id since they store decomposition information
	shared_ptr<RestoreIds<FaceData>> resfd(new RestoreIds<FaceData>(fd));
	//build pyramids
	bnd_grid = HMGrid3D::BuildPyramidLayer(fd, true, 60);

	//split bnd grid inner surfaces according to input surface decomposition
	resfd.reset();
	decomposed_surfs.resize(decomposed_surfs1.size());
	for (int i=0; i<decomposed_surfs1.size(); ++i)
		decomposed_surfs[i].resize(decomposed_surfs1[i].size());
	for (auto& c: bnd_grid.vcells){
		auto& ijc = kmap[c->faces[0]->id];
		auto& addto = decomposed_surfs[ijc.first][ijc.second];
		for (auto f: c->faces) if (!f->has_left_cell()){
			addto.faces.push_back(f);
		}
	}
	//remove zero length surfaces
	for (int i=0; i<decomposed_surfs.size(); ++i){
		auto rs = std::remove_if(decomposed_surfs[i].begin(), decomposed_surfs[i].end(),
			[](const Surface& s){ return s.faces.size() == 0; });
		decomposed_surfs[i].resize(rs - decomposed_surfs[i].begin());
	}
	auto rs = std::remove_if(decomposed_surfs.begin(), decomposed_surfs.end(),
			[](const vector<Surface>& sv){ return sv.size() == 0; });
	decomposed_surfs.resize(rs - decomposed_surfs.begin());

	//subdivide each decomposed surfs if necessary
	for (int i=0; i<decomposed_surfs.size(); ++i){
		int jmax = decomposed_surfs[i].size();
		for (int j=0; j<jmax; ++j){
			auto sd = Face::SubDivide(decomposed_surfs[i][j].faces);
			assert(sd.size() > 0);
			if (sd.size() < 2) continue;
			decomposed_surfs[i][j].faces = sd[0];
			for (int k=1; k<sd.size(); ++k){
				decomposed_surfs[i].emplace_back();
				decomposed_surfs[i].back().faces = sd[k];
			}
		}
	}

	//normals
	for (auto& ds: decomposed_surfs)
	for (auto& s: ds){
		assert(s.faces.size() > 0);
		Face* f;
		if (s.faces[0]->has_right_cell()) f = s.faces[0]->right.lock()->faces[0].get();
		else f = s.faces[0].get();
		surfs_rnormals[&s] = -1*f->left_normal();
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
		for (auto& s: ds) if (s.faces.size() > 0){
			vector<HMGrid3D::EdgeData> vds;
			assert(surfs_rnormals.find(&s) != surfs_rnormals.end());
			Vect3 right_normal = surfs_rnormals[&s];
			auto ex = HMGrid3D::Surface::ExtractAllBoundaries(s, right_normal);
			assert(ex[0].size() == 1);
			vds.push_back(ex[0][0]);
			for (int i=0; i<ex[1].size(); ++i) vds.push_back(ex[1][i]);
			de.push_back(vds);
		} else {
			assert(false);
			de.emplace_back();
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

