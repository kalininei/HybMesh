#include "buildgrid.hpp"
#include "modgrid.hpp"
#include <unordered_map>
#include "healgrid.hpp"

using namespace HM2D;
namespace hgc = HM2D::Grid::Constructor;

GridData hgc::RectGrid01(int Nx, int Ny){
	vector<double> x(Nx+1), y(Ny+1);
	for (int i=0; i<x.size(); ++i){ x[i] = (double)i/(Nx); }
	for (int i=0; i<y.size(); ++i){ y[i] = (double)i/(Ny); }
	return hgc::RectGrid(x, y);
}

GridData hgc::RectGrid(Point p0, Point p1, int Nx, int Ny){
	vector<double> x(Nx+1), y(Ny+1);
	for (int i=0; i<x.size(); ++i){ x[i] = (p1.x-p0.x)*(double)i/(Nx)+p0.x; }
	for (int i=0; i<y.size(); ++i){ y[i] = (p1.y-p0.y)*(double)i/(Ny)+p0.y; }
	return hgc::RectGrid(x, y);
}

GridData hgc::RectGrid(const vector<double>& part_x, const vector<double>& part_y){
	int Npts = part_x.size()*part_y.size();
	int Ncls = (part_x.size()-1)*(part_y.size()-1);
	vector<double> points;
	for (int j=0; j<part_y.size(); ++j){
		for (int i=0; i<part_x.size(); ++i){
			points.push_back(part_x[i]);
			points.push_back(part_y[j]);
		}
	}
	vector<int> cells;
	for (int j=0; j<part_y.size() - 1; ++j){
		for (int i=0; i<part_x.size() - 1; ++i){
			int i0 = j*part_x.size() + i;
			int i1 = j*part_x.size() + i + 1;
			int i2 = (j + 1) * part_x.size() + i + 1;
			int i3 = (j + 1) * part_x.size() + i;
			cells.push_back(i0);
			cells.push_back(i1);
			cells.push_back(i2);
			cells.push_back(i3);
		}
	}
	return hgc::FromRaw(Npts, Ncls, &points[0], &cells[0], 4);
}

GridData hgc::Ring(Point p0, double rad1, double rad2, int narc, int nrad){
	assert(rad1>rad2);
	vector<double> pts;
	vector<int> cells;
	//1) build points
	for (int i=0; i<nrad+1; ++i){
		double r = rad1 + (rad2 - rad1)/nrad*i;
		for (int j=0; j<narc; ++j){
			double phi = 2*M_PI/narc*j;
			pts.push_back(r*cos(phi) + p0.x);
			pts.push_back(r*sin(phi) + p0.y);
		}
	}
	//2) build cells
	for (int i=0; i<nrad; ++i){
		for (int j=0; j<narc; ++j){
			int jnext = (j==narc-1)?0:j+1;
			cells.push_back(narc*i+j);
			cells.push_back(narc*i+jnext);
			cells.push_back(narc*(i+1)+jnext);
			cells.push_back(narc*(i+1)+j);
		}
	}
	//3) build grid
	auto ret = FromRaw(pts.size()/2, cells.size()/4, &pts[0], &cells[0], 4);
	return ret;
}

GridData hgc::Circle(Point p0, double rad, int narc, int nrad, bool tri_center){
	GridData ret;
	if (nrad > 1){
		double rad2 = rad*(1.0)/nrad;
		ret = Ring(p0, rad, rad2, narc, nrad-1);
	} else {
		for (int i=0; i<narc; ++i){
			double phi = 2.0*M_PI/narc*i;
			ret.vvert.emplace_back(
				new Vertex(rad*cos(phi) + p0.x, rad*sin(phi) + p0.y));
		}
	}
	VertexData tmp(ret.vvert);
	vector<vector<int>> cellvert;
	if (!tri_center){
		cellvert.emplace_back();
		for (int i=0; i<narc; ++i){
			int ind = ret.vvert.size() - narc + i;
			cellvert.back().push_back(ind);
		}
	} else {
		tmp.emplace_back(new Vertex(p0));
		for (int i=0; i<narc; ++i){
			int ind1 = ret.vvert.size() - narc + i;
			int ind2 = (i == narc - 1) ? ret.vvert.size() - narc
			                           : ind1 + 1;
			cellvert.emplace_back();
			cellvert.back().push_back(ind1);
			cellvert.back().push_back(ind2);
			cellvert.back().push_back(tmp.size()-1);
		}
	}
	GridData g3 = FromTab(std::move(tmp), cellvert);
	Algos::MergeTo(g3, ret);
	return ret;
}

GridData hgc::FromRaw(int npnt, int ncls, double* pnt, int* cls, int dim){
	VertexData vv(npnt);
	for (auto i=0; i<npnt; ++i){
		vv[i] = std::make_shared<Vertex>(pnt[2*i], pnt[2*i+1]);
	}
	vector<vector<int>> cell_vert(ncls);
	int* itc = cls;
	int ic = 0;
	while (ic<ncls){
		int d = (dim < 2) ? *itc++ : dim;
		for (int i=0; i<d; ++i){
			cell_vert[ic].push_back(*itc++);
		}
		++ic;
	}
	return FromTab(std::move(vv), cell_vert);
}

//uses only those vert which present in vert_cell tabs
GridData hgc::FromTab(const VertexData& vert, const vector<vector<int>>& cell_vert){
	return FromTab(VertexData(vert), cell_vert);
}

GridData hgc::FromTab(VertexData&& vert, const vector<vector<int>>& cell_vert){
	GridData r;
	r.vvert = std::move(vert);

	std::vector<std::map<int, int>> used_edges(r.vvert.size());

	vector<int> edge_p1p2;
	vector<int> edge_leftright;
	vector<vector<int>> celledge(cell_vert.size());
	int neds=0;

	for (int i=0; i<cell_vert.size(); ++i){
		for (int j=0; j<cell_vert[i].size(); ++j){
			int i1 = cell_vert[i][j];
			int i2 = cell_vert[i][(j==cell_vert[i].size()-1) ? 0 : j+1];
			bool rev = false;
			if (i1 > i2) {std::swap(i1, i2); rev = true; }
			auto er = used_edges[i1].emplace(i2, neds);
			if (er.second){
				edge_p1p2.push_back(i1);
				edge_p1p2.push_back(i2);
				edge_leftright.push_back(-1);
				edge_leftright.push_back(-1);
				neds++;
			}
			if (rev){
				edge_leftright[2*er.first->second+1] = i;
			} else {
				edge_leftright[2*er.first->second] = i;
			}
			celledge[i].push_back(er.first->second);
		}
	}

	//edges
	vector<int>::iterator it=edge_p1p2.begin();
	while (it != edge_p1p2.end()){
		auto p1 = r.vvert[*it++];
		auto p2 = r.vvert[*it++];
		r.vedges.push_back(std::make_shared<Edge>(p1, p2));
	}
	//cells
	for (int i=0; i<celledge.size(); ++i){
		r.vcells.push_back(std::make_shared<Cell>());
		for (int j=0; j<celledge[i].size(); ++j){
			r.vcells.back()->edges.push_back(r.vedges[celledge[i][j]]);
		}
	}
	//edges->cells
	for (int i=0; i<r.vedges.size(); ++i){
		int lc = edge_leftright[2*i];
		int rc = edge_leftright[2*i+1];
		if (lc>=0) r.vedges[i]->left = r.vcells[lc];
		if (rc>=0) r.vedges[i]->right = r.vcells[rc];
	}

	//get rid of unused vertices
	for (int i=0; i<r.vvert.size(); ++i){
		if (r.vvert[i] != nullptr) r.vvert[i]->id = 0;
	}
	for (int i: edge_p1p2){
		assert(r.vvert[i]!=nullptr);
		r.vvert[i]->id = 1;
	}
	r.vvert.resize(std::remove_if(r.vvert.begin(), r.vvert.end(),
	                              [](const shared_ptr<Vertex>& v){return v == nullptr;}) -
	               r.vvert.begin());
	aa::remove_by_id(r.vvert, 0);

	return r;
}

hgc::InvokeGrid::InvokeGrid(const CellData& data){
	grid.vcells.reserve(data.size());
	std::copy_if(data.begin(), data.end(),
			std::back_inserter(grid.vcells),
			[](shared_ptr<Cell> c){ return c!=nullptr; });
	grid.vedges = AllEdges(grid.vcells);
	grid.vvert = AllVertices(grid.vedges);

	for (auto e: grid.vedges){
		if (e->has_left_cell()) e->left.lock()->id = 1;
		if (e->has_right_cell()) e->right.lock()->id = 1;
	}
	aa::constant_ids_pvec(data, 0);
	for (int i=0; i<grid.vedges.size(); ++i){
		auto& e = grid.vedges[i];
		if (e->has_left_cell() && e->left.lock()->id == 1){
			oldleft.emplace(i, e->left);
			e->left.reset();
		}
		if (e->has_right_cell() && e->right.lock()->id == 1){
			oldright.emplace(i, e->right);
			e->right.reset();
		}
	}
};

void hgc::InvokeGrid::make_permanent(){
	oldright.clear();
	oldleft.clear();
}

hgc::InvokeGrid::~InvokeGrid(){
	for (auto& it: oldright){
		grid.vedges[it.first]->right = it.second;
	}
	for (auto& it: oldleft){
		grid.vedges[it.first]->left = it.second;
	}
}
