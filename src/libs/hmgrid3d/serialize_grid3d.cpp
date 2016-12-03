#include "serialize_grid3d.hpp"
#include "hmcallback.hpp"
#include "addalgo.hpp"

using namespace HMGrid3D;

namespace{
template<class SerClass>
vector<int> fvtab(const SerClass& s, int nface){
	vector<int> ret;
	const auto& fe = s.face_edge[nface];
	size_t lenf = fe.size();
	ret.reserve(lenf);

	//first vertex
	{
		int e1 = fe[0], e2 = fe[1];
		int p1 = s.edges[2*e1], p2 = s.edges[2*e1+1];
		int p3 = s.edges[2*e2], p4 = s.edges[2*e2+1];
		if (p1 == p3 || p1 == p4) std::swap(p1, p2);
		ret.push_back(p1); ret.push_back(p2);
	}
	//other vertices
	for (size_t k=1; k<lenf-1; ++k){
		int e1 = fe[k];
		int p1 = s.edges[2*e1], p2 = s.edges[2*e1+1];
		ret.push_back( (p1 == ret.back()) ? p2 : p1 );
	}
	return ret;
}
}


//================================================= Surfaces
void SimpleSerializeSurface::supplement(){
	//n_*
	n_vert = vert.size()/3;
	n_edges = edges.size()/2;
	n_faces = face_edge.size();
	cache.clear();
}

vector<vector<int>>& SimpleSerializeSurface::Cache::face_vertex(){
	if (!_face_vertex){
		_face_vertex.reset(new vector<vector<int>>(parent->n_faces));
		for (size_t i=0; i<parent->n_faces; ++i){
			(*_face_vertex)[i] = fvtab(*parent, i);
		}
	}
	return *_face_vertex;
}


SSurface::SSurface(SimpleSerializeSurface&& ss) noexcept:SimpleSerializeSurface(std::move(ss)), Surface(){
	actualize_data();
}
SSurface::SSurface(HMGrid3D::Surface&& srf) noexcept:SimpleSerializeSurface(), HMGrid3D::Surface(std::move(srf)){
	actualize_serial_data();
}

void SSurface::actualize_data() noexcept{
	ShpVector<Vertex> vvert;
	ShpVector<Edge> vedges;
	//vertices
	vvert.resize(n_vert);
	auto vit = vert.begin();
	for (int i=0; i<n_vert; ++i){
		vvert[i].reset(new Vertex(*vit, *(vit+1), *(vit+2)));
		vit+=3;
	}
	//edges
	vedges.resize(n_edges);
	auto eit = edges.begin();
	for (int i=0; i<n_edges; ++i){
		int p1 = *eit++;
		int p2 = *eit++;
		vedges[i].reset(new Edge(vvert[p1], vvert[p2]));
	}
	//faces
	faces.resize(n_faces);
	for (int i=0; i<n_faces; ++i){
		faces[i].reset(new Face());
		faces[i]->boundary_type = btypes[i];
		faces[i]->edges.resize(face_edge[i].size());
		for (size_t j=0; j<face_edge[i].size(); ++j){
			faces[i]->edges[j]=vedges[face_edge[i][j]];
		}
	}
}

void SSurface::actualize_serial_data() noexcept{
	cache.clear();
	//enumeration
	int kf=0;
	for (auto f: faces){
		f->id = kf++;
		for (auto e: f->edges){
			e->id = -1;
			for (auto v: e->vertices){
				v->id = -1;
			}
		}
	}
	//assembling
	ShpVector<Vertex> vvert;
	ShpVector<Edge> vedges;
	int kv = 0, ke = 0;
	for (int i=0; i<kf; ++i){
		auto f = faces[i];
		for (auto e: f->edges) if (e->id < 0){
			e->id = ke++;
			vedges.push_back(e);
			if (e->first()->id < 0){
				e->first()->id = kv++;
				vvert.push_back(e->first());
			}
			if (e->last()->id < 0){
				e->last()->id = kv++;
				vvert.push_back(e->last());
			}
		}
	}

	//vertices
	n_vert = vvert.size();
	vert.resize(3*n_vert);
	auto vit = vert.begin();
	for (int i=0; i<n_vert; ++i){
		*vit++ = vvert[i]->x;
		*vit++ = vvert[i]->y;
		*vit++ = vvert[i]->z;
	}
	//edges
	n_edges = vedges.size();
	edges.resize(2*n_edges);
	auto eit = edges.begin();
	for (int i=0; i<n_edges; ++i){
		*eit++ = vedges[i]->first()->id;
		*eit++ = vedges[i]->last()->id;
	}
	//faces
	n_faces = faces.size();
	face_edge.resize(n_faces);
	for (size_t i=0; i<n_faces; ++i){
		face_edge[i].resize(faces[i]->n_edges());
		for (size_t j=0; j<faces[i]->n_edges(); ++j){
			face_edge[i][j] = faces[i]->edges[j]->id;
		}
	}
	//boundary types
	btypes.resize(n_faces);
	for (size_t i=0; i<n_faces; ++i){
		btypes[i] = faces[i]->boundary_type;
	}
	//all other
	supplement();
}

//================================================= Grids
//fills everything from: vert, edges, faces, bnd.
void SimpleSerialize::supplement(){
	//n_*
	n_vert = vert.size()/3;
	n_edges = edges.size()/2;
	n_faces = face_edge.size();
	n_cells = *std::max_element(face_cell.begin(), face_cell.end())+1;
	cache.clear();
}

SimpleSerializeSurface SimpleSerialize::serialized_surface() const{
	SimpleSerializeSurface ret;
	//vertices
	vector<int> vindices(n_vert, -1);
	auto& bv = bvert();
	ret.vert.resize(bv.size()*3);
	int vit = 0;
	for (int i: bv){
		vindices[i] = vit;
		ret.vert[3*vit] = vert[3*i];
		ret.vert[3*vit+1] = vert[3*i+1];
		ret.vert[3*vit+2] = vert[3*i+2];
		++vit;
	}
	//edges
	vector<int> eindices(n_edges, -1);
	auto& be = bedges();
	ret.edges.resize(2*be.size());
	int eit = 0;
	for (int i: be){
		eindices[i] = eit;
		ret.edges[2*eit] = vindices[edges[2*i]];
		ret.edges[2*eit+1] = vindices[edges[2*i+1]];
		++eit;
	}
	//faces
	auto& bf = bfaces();
	ret.face_edge.resize(bf.size());
	ret.btypes.resize(bf.size());
	int fit = 0;
	for (int i: bf){
		ret.face_edge[fit].resize(face_edge[i].size());
		for (size_t j=0; j<ret.face_edge[fit].size(); ++j){
			ret.face_edge[fit][j] = eindices[face_edge[i][j]];
		}
		ret.btypes[fit] = btypes[i];
		++fit;
	}

	ret.supplement();
	return ret;
};

vector<vector<int>>& SimpleSerialize::Cache::face_vertex(){
	if (!_face_vertex){
		_face_vertex.reset(new vector<vector<int>>(parent->n_faces));
		for (size_t i=0; i<parent->n_faces; ++i){
			(*_face_vertex)[i] = fvtab(*parent, i);
		}
	}
	return *_face_vertex;
}
vector<int>& SimpleSerialize::Cache::bfaces(){
	if (!_bfaces){
		_bfaces.reset(new vector<int>());
		for (size_t i=0; i<parent->n_faces; ++i){
			if (parent->face_cell[2*i]<0 || parent->face_cell[2*i+1]<0){
				_bfaces->push_back((int)i);
			}
		}
	}
	return *_bfaces;
}
vector<int>& SimpleSerialize::Cache::bedges(){
	if (!_bedges){
		_bedges.reset(new vector<int>());
		vector<bool> used(parent->n_edges, false);
		for (auto bf: bfaces())
		for (auto e: parent->face_edge[bf])
			used[e]=true;
		for (size_t i=0; i<used.size(); ++i)
		if (used[i])
			_bedges->push_back(i);
	}
	return *_bedges;
}
vector<int>& SimpleSerialize::Cache::bvert(){
	if (!_bvert){
		_bvert.reset(new vector<int>());
		vector<bool> used(parent->n_vert, false);
		for (auto be: bedges()){
			used[parent->edges[2*be]]=true;
			used[parent->edges[2*be+1]]=true;
		}
		for (size_t i=0; i<used.size(); ++i)
		if (used[i])
			_bvert->push_back(i);
	}
	return *_bvert;
}
vector<vector<int>>& SimpleSerialize::Cache::cell_face(){
	if (!_cell_face){
		_cell_face.reset(new vector<vector<int>>(parent->n_cells));
		for (size_t i=0; i<parent->n_faces; ++i){
			int c1 = parent->face_cell[2*i];
			int c2 = parent->face_cell[2*i+1];
			if (c1>=0) (*_cell_face)[c1].push_back(i);
			if (c2>=0) (*_cell_face)[c2].push_back(i);
		}
	}
	return *_cell_face;
}
vector<vector<int>>& SimpleSerialize::Cache::cell_vertex(){
	if (!_cell_vertex){
		_cell_vertex.reset(new vector<vector<int>>(parent->n_cells));
		_THROW_NOT_IMP_;
	}
	return *_cell_vertex;
}

//fills SimpleSerialize on the basis of GridData
void SGrid::actualize_serial_data() noexcept{
	cache.clear();
	enumerate_all();
	//vertices
	n_vert = vvert.size();
	vert.resize(3*n_vert);
	auto vit = vert.begin();
	for (int i=0; i<n_vert; ++i){
		*vit++ = vvert[i]->x;
		*vit++ = vvert[i]->y;
		*vit++ = vvert[i]->z;
	}
	//edges
	n_edges = vedges.size();
	edges.resize(2*n_edges);
	auto eit = edges.begin();
	for (int i=0; i<n_edges; ++i){
		*eit++ = vedges[i]->first()->id;
		*eit++ = vedges[i]->last()->id;
	}
	//faces
	n_faces = vfaces.size();
	face_edge.resize(n_faces);
	for (size_t i=0; i<n_faces; ++i){
		face_edge[i].resize(vfaces[i]->n_edges());
		for (size_t j=0; j<vfaces[i]->n_edges(); ++j){
			face_edge[i][j] = vfaces[i]->edges[j]->id;
		}
	}
	face_cell.resize(2*n_faces);
	for (size_t i=0; i<n_faces; ++i){
		if (vfaces[i]->has_left_cell())
			face_cell[2*i] = vfaces[i]->left.lock()->id;
		else
			face_cell[2*i] = -1;
		if (vfaces[i]->has_right_cell())
			face_cell[2*i+1] = vfaces[i]->right.lock()->id;
		else
			face_cell[2*i+1] = -1;
	}
	//boundary types
	btypes.resize(n_faces);
	for (size_t i=0; i<n_faces; ++i){
		btypes[i] = vfaces[i]->boundary_type;
	}
	//all other
	supplement();
}

//fills GridData vectors from SimpleSerialize Data
void SGrid::actualize_data() noexcept{
	//vertices
	vvert.resize(n_vert);
	auto vit = vert.begin();
	for (int i=0; i<n_vert; ++i){
		vvert[i].reset(new Vertex(*vit, *(vit+1), *(vit+2)));
		vit+=3;
	}
	//edges
	vedges.resize(n_edges);
	auto eit = edges.begin();
	for (int i=0; i<n_edges; ++i){
		int p1 = *eit++;
		int p2 = *eit++;
		vedges[i].reset(new Edge(vvert[p1], vvert[p2]));
	}
	//faces
	vfaces.resize(n_faces);
	for (int i=0; i<n_faces; ++i){
		vfaces[i].reset(new Face());
		vfaces[i]->boundary_type = btypes[i];
		vfaces[i]->edges.resize(face_edge[i].size());
		for (size_t j=0; j<face_edge[i].size(); ++j){
			vfaces[i]->edges[j]=vedges[face_edge[i][j]];
		}
	}
	//cells
	vcells.resize(n_cells);
	for (int i=0; i<n_cells; ++i){ vcells[i].reset(new Cell()); }
	for (int i=0; i<n_faces; ++i){
		int lcell = face_cell[2*i];
		int rcell = face_cell[2*i+1];
		if (lcell >= 0){
			vcells[lcell]->faces.push_back(vfaces[i]);
			vfaces[i]->left = vcells[lcell];
		}
		if (rcell >= 0){
			vcells[rcell]->faces.push_back(vfaces[i]);
			vfaces[i]->right = vcells[rcell];
		}
	}
}

void SGrid::set_btype(std::function<int(Vertex, int)> func){
	Face::SetBoundaryTypes(vfaces, func);
	btypes.resize(n_faces);
	for (size_t i=0; i<n_faces; ++i){
		btypes[i] = vfaces[i]->boundary_type;
	}
}

void SGrid::renumber_by_cells(){
	constant_ids_pvec(vfaces, -10);
	constant_ids_pvec(vedges, -10);
	constant_ids_pvec(vvert, -10);

	int indv=0, indf=0, inde=0;
	for (auto c: vcells){
		for (auto f: c->faces){
			if (f->id<0) f->id = indf++;
			for (auto e: f->edges){
				if (e->id<0) e->id = inde++;
				for (auto v: e->vertices){
					if (v->id<0) v->id = indv++;
				}
			}
		}
	}
	assert(indv==n_vert && indf==n_faces && inde==n_edges);

	FaceData newfaces(vfaces.size());
	EdgeData newedges(vedges.size());
	VertexData newvert(vvert.size());

	for (int i=0; i<n_faces; ++i) newfaces[vfaces[i]->id] = vfaces[i];
	for (int i=0; i<n_edges; ++i) newedges[vedges[i]->id] = vedges[i];
	for (int i=0; i<n_vert; ++i) newvert[vvert[i]->id] = vvert[i];

	std::swap(vfaces, newfaces);
	std::swap(vedges, newedges);
	std::swap(vvert, newvert);

	actualize_serial_data();
}
