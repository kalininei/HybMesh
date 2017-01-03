#include "serialize3d.hpp"
#include "surface.hpp"

using namespace HM3D;
using namespace HM3D::Ser;

namespace{
template<class SerClass>
vector<int> fvtab(const SerClass& s, int nface){
	vector<int> ret;
	const auto& fe = s.face_edge()[nface];
	size_t lenf = fe.size();
	ret.reserve(lenf);

	auto& edge_vert = s.edge_vert();
	//first vertex
	{
		int e1 = fe[0], e2 = fe[1];
		int p1 = edge_vert[2*e1], p2 = edge_vert[2*e1+1];
		int p3 = edge_vert[2*e2], p4 = edge_vert[2*e2+1];
		if (p1 == p3 || p1 == p4) std::swap(p1, p2);
		ret.push_back(p1); ret.push_back(p2);
	}
	//other vertices
	for (size_t k=1; k<lenf-1; ++k){
		int e1 = fe[k];
		int p1 = edge_vert[2*e1], p2 = edge_vert[2*e1+1];
		ret.push_back( (p1 == ret.back()) ? p2 : p1 );
	}
	return ret;
}
}

// ======================================= Surface
struct Ser::Surface::Cache{
	const Ser::Surface* parent;
	Cache(const Ser::Surface& par): parent(&par){}

	//=============== data
	//main tables
	vector<double> _vert;
	vector<int> _btypes;
	vector<int> _edge_vert;
	vector<vector<int>> _face_edge;
	vector<vector<int>> _face_vertex;

	// ============== callers
	vector<double>& vert(){
		if (_vert.size() == 0){
			auto av = AllVertices(parent->surface);
			_vert.reserve(av.size()*3);
			for (auto& v: av){
				_vert.push_back(v->x);
				_vert.push_back(v->y);
				_vert.push_back(v->z);
			}
		}
		return _vert;
	}
	vector<int>& edge_vert(){
		if (_edge_vert.size() == 0){
			auto ave = AllPrimitives(parent->surface);
			auto& av = std::get<0>(ave);
			auto& ae = std::get<1>(ave);
			aa::enumerate_ids_pvec(av);
			_edge_vert.reserve(2*ae.size());
			for (auto e: ae){
				_edge_vert.push_back(e->first()->id);
				_edge_vert.push_back(e->last()->id);
			}
		}
		return _edge_vert;
	}
	vector<int>& btypes(){
		if (_btypes.size() == 0){
			for (auto f: parent->surface){
				_btypes.push_back(f->boundary_type);
			}
		}
		return _btypes;
	}
	vector<vector<int>>& face_edge(){
		if (_face_edge.size() == 0){
			auto ae = AllEdges(parent->surface);
			aa::enumerate_ids_pvec(ae);
			_face_edge.resize(parent->n_faces());
			for (int i=0; i<_face_edge.size(); ++i){
				for (auto e: parent->surface[i]->edges){
					_face_edge[i].push_back(e->id);
				}
			}
		}
		return _face_edge;
	}
	vector<vector<int>>& face_vertex(){
		if (_face_vertex.size() == 0){
			_face_vertex.resize(parent->n_faces());
			for (size_t i=0; i<parent->n_faces(); ++i){
				_face_vertex[i] = fvtab(*parent, i);
			}
		}
		return _face_vertex;
	}

	int n_edges() { return edge_vert().size()/2; }
	int n_vert() { return vert().size()/3; }
};


Ser::Surface::Surface(){
	cache.reset(new Cache(*this));
}
Ser::Surface::Surface(HM3D::FaceData&& srf) noexcept{
	std::swap(surface, srf);
	cache.reset(new Cache(*this));
}
Ser::Surface::Surface(const HM3D::FaceData& srf){
	surface = srf;
	cache.reset(new Cache(*this));
}
void Ser::Surface::empty_cache() const {
	cache.reset(new Cache(*this));
}
Ser::Surface::~Surface(){}

int Ser::Surface::n_edges() const { return cache->n_edges(); }
int Ser::Surface::n_vert() const { return cache->n_vert(); }
const vector<double>& Ser::Surface::vert() const { return cache->vert(); }
const vector<int>& Ser::Surface::edge_vert() const { return cache->edge_vert(); }
const vector<vector<int>>& Ser::Surface::face_edge() const { return cache->face_edge(); }
const vector<int>& Ser::Surface::btypes() const { return cache->btypes(); }
const vector<int>& Ser::Surface::face_vertex(int num_face) const{
	return cache->face_vertex()[num_face];
}
const vector<vector<int>>& Ser::Surface::face_vertex() const{
	return cache->face_vertex();
}
void Ser::Surface::fill_from_serial(const vector<double>& vert,
		const vector<int>& edgevert,
		const vector<vector<int>>& faceedge,
		const vector<int>& btypes){
	//fill cache
	empty_cache();
	cache->_vert = vert;
	cache->_edge_vert = edgevert;
	cache->_face_edge = faceedge;
	cache->_btypes = btypes;
	//fill grid
	VertexData vvert;
	EdgeData vedges;
	surface.clear();
	//vertices
	vvert.resize(vert.size()/3);
	for (int i=0; i<vert.size(); i+=3){
		vvert[i/3].reset(new Vertex(vert[i], vert[i+1], vert[i+2]));
	}
	//edges
	vedges.resize(edgevert.size()/2);
	for (int i=0; i<edgevert.size(); i+=2){
		vedges[i/2].reset(new Edge(vvert[edgevert[i]], vvert[edgevert[i+1]]));
	}
	//faces
	surface.resize(faceedge.size());
	for (int i=0; i<faceedge.size(); ++i){
		surface[i].reset(new Face());
		for (int j=0; j<faceedge[i].size(); ++j){
			surface[i]->edges.push_back(vedges[faceedge[i][j]]);
		}
		surface[i]->boundary_type = btypes[i];
	}
}


// ======================================= Grid
struct Ser::Grid::Cache{
	const Ser::Grid* parent;
	Cache(const Ser::Grid& par): parent(&par){}

	// ================ data
	//main tables
	vector<double> _vert;
	vector<int> _edge_vert;
	vector<vector<int>> _face_edge;
	vector<int> _face_cell;
	vector<int> _btypes;
	//aux tables
	vector<int> _bfaces;
	vector<int> _bedges;
	vector<int> _bvert;
	vector<vector<int>> _face_vertex;

	// ================ callers
	//main tables
	vector<double>& vert(){
		if (_vert.size() == 0){
			_vert.reserve(3*parent->n_vert());
			for (auto& v: parent->grid.vvert){
				_vert.push_back(v->x);
				_vert.push_back(v->y);
				_vert.push_back(v->z);
			}
		}
		return _vert;
	}
	vector<int>& edge_vert(){
		if (_edge_vert.size() == 0){
			_edge_vert.reserve(parent->n_edges()*2);
			aa::enumerate_ids_pvec(parent->grid.vvert);
			for (auto e: parent->grid.vedges){
				_edge_vert.push_back(e->first()->id);
				_edge_vert.push_back(e->last()->id);
			}
		}
		return _edge_vert;
	}
	vector<vector<int>>& face_edge(){
		if (_face_edge.size() == 0){
			aa::enumerate_ids_pvec(parent->grid.vedges);
			_face_edge.resize(parent->n_faces());
			for (int i=0; i<_face_edge.size(); ++i){
				for (auto e: parent->grid.vfaces[i]->edges){
					_face_edge[i].push_back(e->id);
				}
			}
		}
		return _face_edge;
	}
	vector<int>& face_cell(){
		if (_face_cell.size() == 0){
			_face_cell.resize(parent->n_faces()*2, -1);
			aa::enumerate_ids_pvec(parent->grid.vcells);
			for (int i=0; i<parent->n_faces(); ++i){
				auto f=parent->grid.vfaces[i];
				if (f->has_left_cell())
					_face_cell[2*i] = f->left.lock()->id;
				if (f->has_right_cell())
					_face_cell[2*i+1] = f->right.lock()->id;
			}
		}
		return _face_cell;
	}
	vector<int>& btypes(){
		if (_btypes.size()==0){
			_btypes.resize(parent->n_faces());
			for (int i=0; i<parent->n_faces(); ++i){
				_btypes[i] = parent->grid.vfaces[i]->boundary_type;
			}
		}
		return _btypes;
	}
	//aux tables
	vector<vector<int>>& face_vertex(){
		if (_face_vertex.size() == 0){
			_face_vertex.resize(parent->n_faces());
			for (size_t i=0; i<parent->n_faces(); ++i){
				_face_vertex[i] = fvtab(*parent, i);
			}
		}
		return _face_vertex;
	}

	vector<int>& bfaces(){
		if (_bfaces.size() == 0){
			for (int i=0; i<parent->n_faces(); ++i){
				if (parent->grid.vfaces[i]->is_boundary())
					_bfaces.push_back(i);
			}
		}
		return _bfaces;
	}

	vector<int>& bedges(){
		if (_bedges.size() == 0){
			vector<bool> used(parent->n_edges(), false);
			auto& fe = face_edge();
			for (auto bf: bfaces())
			for (auto e: fe[bf])
				used[e]=true;
			for (size_t i=0; i<used.size(); ++i)
			if (used[i]) _bedges.push_back(i);
		}
		return _bedges;
	}

	vector<int>& bvert(){
		if (_bvert.size() == 0){
			vector<bool> used(parent->n_vert(), false);
			auto& ed = edge_vert();
			for (auto be: bedges()){
				used[ed[2*be]]=true;
				used[ed[2*be+1]]=true;
			}
			for (size_t i=0; i<used.size(); ++i)
			if (used[i]) _bvert.push_back(i);
		}
		return _bvert;
	}
};

Ser::Grid::Grid(){
	cache.reset(new Cache(*this));
}
Ser::Grid::Grid(HM3D::GridData&& g) noexcept{
	cache.reset(new Cache(*this));
	std::swap(grid, g);
}
Ser::Grid::Grid(const HM3D::GridData& g){
	cache.reset(new Cache(*this));
	grid = g;
}
Ser::Grid::~Grid(){}
void Ser::Grid::empty_cache() const {
	cache.reset(new Cache(*this));
}

const vector<double>& Ser::Grid::vert() const { return cache->vert(); }
const vector<int>& Ser::Grid::edge_vert() const { return cache->edge_vert(); }
const vector<vector<int>>& Ser::Grid::face_edge() const { return cache->face_edge(); }
const vector<int>& Ser::Grid::face_cell() const { return cache->face_cell(); }
const vector<int>& Ser::Grid::bfaces() const { return cache->bfaces(); }
const vector<int>& Ser::Grid::bedges() const { return cache->bedges(); }
const vector<int>& Ser::Grid::bvert() const { return cache->bvert(); }
const vector<int>& Ser::Grid::btypes() const { return cache->btypes(); }
const vector<vector<int>>& Ser::Grid::face_vertex() const { return cache->face_vertex(); }
const vector<int>& Ser::Grid::face_vertex(int n) const { return cache->face_vertex()[n]; }

void Ser::Grid::set_btype(std::function<int(Vertex, int)> func){
	auto bsurf = HM3D::Surface::GridSurface(grid);
	HM3D::Surface::SetBoundaryTypes(bsurf, func);
	cache->_btypes.clear();
}

void Ser::Grid::renumber_by_cells(){
	aa::constant_ids_pvec(grid.vfaces, -10);
	aa::constant_ids_pvec(grid.vedges, -10);
	aa::constant_ids_pvec(grid.vvert, -10);

	int indv=0, indf=0, inde=0;
	for (auto c: grid.vcells){
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
	assert(indv==n_vert() && indf==n_faces() && inde==n_edges());

	FaceData newfaces(grid.vfaces.size());
	EdgeData newedges(grid.vedges.size());
	VertexData newvert(grid.vvert.size());

	for (int i=0; i<n_faces(); ++i) newfaces[grid.vfaces[i]->id] = grid.vfaces[i];
	for (int i=0; i<n_edges(); ++i) newedges[grid.vedges[i]->id] = grid.vedges[i];
	for (int i=0; i<n_vert(); ++i) newvert[grid.vvert[i]->id] = grid.vvert[i];

	std::swap(grid.vfaces, newfaces);
	std::swap(grid.vedges, newedges);
	std::swap(grid.vvert, newvert);

	reset_geometry();
}

void Ser::Grid::fill_from_serial(const vector<double>& vert,
		const vector<int>& edgevert,
		const vector<vector<int>>& faceedge,
		const vector<int>& facecell,
		const vector<int>& btypes){
	//fill cache
	empty_cache();
	cache->_vert = vert;
	cache->_edge_vert = edgevert;
	cache->_face_edge = faceedge;
	cache->_face_cell = facecell;
	cache->_btypes = btypes;
	//fill grid
	grid.clear();
	//vertices
	grid.vvert.resize(vert.size()/3);
	for (int i=0; i<vert.size(); i+=3){
		grid.vvert[i/3].reset(new Vertex(vert[i], vert[i+1], vert[i+2]));
	}
	//edges
	grid.vedges.resize(edgevert.size()/2);
	for (int i=0; i<edgevert.size(); i+=2){
		grid.vedges[i/2].reset(new Edge(grid.vvert[edgevert[i]], grid.vvert[edgevert[i+1]]));
	}
	//faces
	grid.vfaces.resize(faceedge.size());
	for (int i=0; i<faceedge.size(); ++i){
		grid.vfaces[i].reset(new Face());
		for (int j=0; j<faceedge[i].size(); ++j){
			grid.vfaces[i]->edges.push_back(grid.vedges[faceedge[i][j]]);
		}
		grid.vfaces[i]->boundary_type = btypes[i];
	}
	//cells
	int cmax = *std::max_element(facecell.begin(), facecell.end());
	grid.vcells.resize(cmax + 1);
	for (auto& c: grid.vcells) c.reset(new Cell());
	for (int i=0; i<facecell.size(); i+=2){
		auto f = grid.vfaces[i/2];
		if (facecell[i]>=0){
			f->left = grid.vcells[facecell[i]];
			grid.vcells[facecell[i]]->faces.push_back(f);
		}
		if (facecell[i+1]>=0){
			f->right = grid.vcells[facecell[i+1]];
			grid.vcells[facecell[i+1]]->faces.push_back(f);
		}
	}
}
