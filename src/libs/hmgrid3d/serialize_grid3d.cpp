#include "serialize_grid3d.hpp"
#include "hmcallback.hpp"
#include "addalgo.hpp"

using namespace HMGrid3D;

vector<int> SimpleSerialize::face_vertex(int nface) const{
	vector<int> ret;
	int kf = ifaces[nface];
	int lenf = faces[kf];
	ret.reserve(lenf);

	//first vertex
	int e1 = faces[kf+1], e2 = faces[kf+2];
	int p1 = edges[2*e1], p2 = edges[2*e1+1];
	int p3 = edges[2*e2], p4 = edges[2*e2+1];
	if (p1 == p3 || p1 == p4) std::swap(p1, p2);
	ret.push_back(p1); ret.push_back(p2);
	//other vertices
	for (int k=1; k<lenf-1; ++k){
		e1 = faces[kf + k];
		e2 = faces[kf + 1 + k];
		int p1 = edges[2*e1], p2 = edges[2*e1+1];
		int p3 = edges[2*e2], p4 = edges[2*e2+1];
		ret.push_back( (p3 == ret.back()) ? p4 : p3 );
	}
	return ret;
}

vector<vector<int>> SimpleSerialize::face_vertex() const{
	vector<vector<int>> ret; ret.reserve(n_faces);
	for (int iface=0; iface<n_faces; ++iface){
		ret.push_back(face_vertex(iface));
	}
	return ret;
}

namespace{

void btype_1to2(ShpVector<Face>& vdata, vector<int>& sdata){
	int n_faces = vdata.size();
	std::map<int, vector<int>> btmap;
	for (int i=0; i<n_faces; ++i) if (vdata[i]->is_boundary()){
		int bt = vdata[i]->boundary_type;
		auto fnd = btmap.find(bt);
		if (fnd == btmap.end()){
			auto emp = btmap.emplace(bt, vector<int>());
			fnd = emp.first;
		}
		fnd->second.push_back(i);
	}
	int bsz = 0;
	for (auto& it: btmap) bsz+=(2+it.second.size());
	sdata.resize(bsz);
	auto bit = sdata.begin();
	for (auto& v: btmap){
		*bit++ = v.first;
		*bit++ = v.second.size();
		bit=std::copy(v.second.begin(), v.second.end(), bit);
	}
}

}


//fills SimpleSerialize on the basis of GridData
void SGrid::actualize_serial_data(){
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
	int fsz = 0;
	for (int i=0; i<n_faces; ++i) fsz+=(3 + vfaces[i]->edges.size());
	faces.resize(fsz);
	ifaces.resize(n_faces+1);
	auto fit = faces.begin();
	for (int i=0; i<n_faces; ++i){
		ifaces[i] = fit - faces.begin();
		auto f = vfaces[i];
		*fit++ = f->edges.size();
		for (int j=0; j<f->edges.size(); ++j) *fit++ = f->edges[j]->id;
		*fit++ = (f->has_left_cell()) ? f->left.lock()->id : -1;
		*fit++ = (f->has_right_cell()) ? f->right.lock()->id : -1;
	}
	ifaces[n_faces] = fit - faces.begin();
	//cells
	n_cells = vcells.size();
	int csz = 0;
	for (int i=0; i<n_cells; ++i) csz+=(2+vcells[i]->faces.size());
	cells.resize(csz);
	icells.resize(n_cells+1);
	auto cit = cells.begin();
	for (int i=0; i<n_cells; ++i){
		icells[i] = cit - cells.begin();
		auto c = vcells[i];
		int n = c->faces.size();
		*cit++ = n;
		for (int j=0; j<n; ++j) *cit++ = c->faces[j]->id;
	}
	icells[n_cells] = cit-cells.begin();
	//boundary types
	btype_1to2(vfaces, bnd);
}

//fills GridData vectors from SimpleSerialize Data
void SGrid::actualize_data(){
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
	auto fit = faces.begin();
	for (int i=0; i<n_faces; ++i){
		vfaces[i].reset(new Face());
		int nf = *fit++;
		vfaces[i]->edges.resize(nf);
		for (int j=0; j<nf; ++j){
			vfaces[i]->edges[j] = vedges[*fit++];
		}
		fit+=2;
	}
	//cells
	vcells.resize(n_cells);
	auto cit = cells.begin();
	for (int i=0; i<n_cells; ++i){
		vcells[i].reset(new Cell());
		int nf = *cit++;
		vcells[i]->faces.resize(nf);
		for (int j=0; j<nf; ++j){
			vcells[i]->faces[j] = vfaces[*cit++];
		}
	}
	//face->cell connectivity
	fit = faces.begin();
	for (int i=0; i<n_faces; ++i){
		int n = *fit++;
		fit+=n;
		int cleft = *fit++;
		int cright = *fit++;
		if (cleft >= 0) vfaces[i]->left = vcells[cleft];
		if (cright >= 0) vfaces[i]->right = vcells[cright];
	}
	//boundary conditions
	auto bit = bnd.begin();
	while (bit != bnd.end()){
		int b = *bit++;
		int n = *bit++;
		for (int i=0; i<n; ++i){
			vfaces[*bit++]->boundary_type = b;
		}
	}
}

void SGrid::set_btype(std::function<int(Vertex, int)> func){
	Face::SetBoundaryTypes(vfaces, func);
	btype_1to2(vfaces, bnd);
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
