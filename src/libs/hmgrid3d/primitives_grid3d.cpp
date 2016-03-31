#include "primitives_grid3d.hpp"
#include "addalgo.hpp"

using namespace HMGrid3D;
// ================= Edge
double Edge::measure() const{
	double ret = 0;
	for (int i = 0; i<(int)vertices.size()-1; ++i){
		double xd = vertices[i]->x - vertices[i+1]->x;
		double yd = vertices[i]->y - vertices[i+1]->y;
		double zd = vertices[i]->z - vertices[i+1]->z;
		ret += xd*xd + yd*yd + zd*zd;
	}
	return ret;
}
double Edge::length() const{
	return sqrt(measure());
}

// ================= Face
ShpVector<Vertex> Face::sorted_vertices() const{
	assert(edges.size() > 1);
	ShpVector<Vertex> ret;
	auto it = edges.begin();
	auto it2 = std::next(it);
	while (it2 != edges.end()){
		if ((*it)->vertices.back() == (*it2)->vertices[0] ||
		    (*it)->vertices.back() == (*it2)->vertices.back()){
			ret.insert(ret.end(), (*it)->vertices.begin(), (*it)->vertices.end()-1);
		} else {
			ret.insert(ret.end(), (*it)->vertices.rbegin(), (*it)->vertices.rend()-1);
		}
		++it; ++it2;
	}
	it2 = edges.begin();
	if ((*it)->vertices.back() == (*it2)->vertices[0] ||
	    (*it)->vertices.back() == (*it2)->vertices.back()){
		ret.insert(ret.end(), (*it)->vertices.begin(), (*it)->vertices.end()-1);
	} else {
		ret.insert(ret.end(), (*it)->vertices.rbegin(), (*it)->vertices.rend()-1);
	}
	return ret;
}

ShpVector<Vertex> Face::allvertices() const{
	ShpVector<Vertex> ret;
	for (auto e: edges){
		ret.insert(ret.end(), e->vertices.begin(), e->vertices.end());
	}
	return aa::no_dublicates(ret);
}

// ================== Cell
int Cell::n_faces() const{
	return faces.size();
}

int Cell::n_edges() const{
	return alledges().size();
}
int Cell::n_vertices() const{
	return allvertices().size();
}

ShpVector<Vertex> Cell::allvertices() const{
	ShpVector<Vertex> ret;
	for (auto f: faces){
		auto dt = f->allvertices();
		ret.insert(ret.end(), dt.begin(), dt.end());
	}
	return aa::no_dublicates(ret);
}

ShpVector<Face> Cell::allfaces() const{ return faces; }

ShpVector<Edge> Cell::alledges() const{
	ShpVector<Edge> ret;
	for (auto c: faces){
		auto dt = c->edges;
		ret.insert(ret.end(), dt.begin(), dt.end());
	}
	return aa::no_dublicates(ret);
}

// ================== Grid
int Grid::n_cells() const{
	return cells.size();
}
int Grid::n_faces() const{
	return allfaces().size();
}
int Grid::n_edges() const{
	return alledges().size();
}
int Grid::n_vertices() const{
	return allvertices().size();
}
ShpVector<Vertex> Grid::allvertices() const{
	ShpVector<Vertex> ret;
	for (auto c: cells){
		auto dt = c->allvertices();
		ret.insert(ret.end(), dt.begin(), dt.end());
	}
	return aa::no_dublicates(ret);
}

ShpVector<Edge> Grid::alledges() const{
	ShpVector<Edge> ret;
	for (auto c: cells){
		auto dt = c->alledges();
		ret.insert(ret.end(), dt.begin(), dt.end());
	}
	return aa::no_dublicates(ret);
}

ShpVector<Face> Grid::allfaces() const{
	ShpVector<Face> ret;
	for (auto c: cells) ret.insert(ret.end(), c->faces.begin(), c->faces.end());
	return aa::no_dublicates(ret);
}

ShpVector<Cell> Grid::allcells() const{ return cells; }
