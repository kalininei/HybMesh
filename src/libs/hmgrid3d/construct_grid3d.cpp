#include "construct_grid3d.hpp"
#include "debug_grid3d.hpp"

using namespace HMGrid3D;

namespace cns = Constructor;

HMGrid3D::Grid cns::Cuboid(HMGrid3D::Vertex leftp, double lx, double ly, double lz, int nx, int ny, int nz){
	HMGrid3D::Grid ret;
	//Vertices
	ShpVector<Vertex> vert; vert.reserve(nx*ny*nz);
	double hx = lx/(nx), hy = ly/(ny), hz = lz/(nz);
	vector<double> x(nx+1), y(ny+1), z(nz+1);
	for (int i=0; i<nx+1; ++i) x[i] = leftp.x + i*hx;
	for (int j=0; j<ny+1; ++j) y[j] = leftp.y + j*hy;
	for (int k=0; k<nz+1; ++k) z[k] = leftp.z + k*hz;
	for (int k=0; k<nz+1; ++k){
		for (int j=0; j<ny+1; ++j){
			for (int i=0; i<nx+1; ++i){
				vert.emplace_back(new Vertex(x[i], y[j], z[k]));
			}
		}
	}
	int sx=1, sy=nx+1, sz=(nx+1)*(ny+1), N=(nx+1)*(ny+1)*(nz+1);
	//Edges
	ShpVector<Edge> xedges(N), yedges(N), zedges(N);
	int gi=0;
	for (int k=0; k<nz+1; ++k){
		for (int j=0; j<ny+1; ++j){
			for (int i=0; i<nx+1; ++i){
				auto p1 = vert[gi];
				//x edge
				if (i<nx){
					auto p2 = vert[gi+sx];
					xedges[gi].reset(new Edge(p1, p2));
				}
				//y edge
				if (j<ny){
					auto p2 = vert[gi+sy];
					yedges[gi].reset(new Edge(p1, p2));
				}
				//z edge
				if (k<nz){
					auto p2 = vert[gi+sz];
					zedges[gi].reset(new Edge(p1, p2));
				}
				++gi;
			}
		}
	}
	//Faces
	ShpVector<Face> xfaces(N), yfaces(N), zfaces(N);
	gi=0;
	for (int k=0; k<nz+1; ++k){
		for (int j=0; j<ny+1; ++j){
			for (int i=0; i<nx+1; ++i){
				//z face
				if (i<nx && j<ny){
					auto e1=yedges[gi], e2=xedges[gi+sy], e3=yedges[gi+sx], e4=xedges[gi];
					zfaces[gi].reset(new Face({e1, e2, e3, e4}));
				}
				//y face
				if (i<nx && k<nz){
					auto e1=xedges[gi], e2=zedges[gi+sx], e3=xedges[gi+sz], e4=zedges[gi];
					yfaces[gi].reset(new Face({e1, e2, e3, e4}));
				};
				//x face
				if (j<ny && k<nz){
					auto e1=zedges[gi], e2=yedges[gi+sz], e3=zedges[gi+sy], e4=yedges[gi];
					xfaces[gi].reset(new Face({e1, e2, e3, e4}));
				}
				++gi;
			}
		}
	}
	//Cells
	ret.cells.reserve(nx*ny*nz);
	for (int k=0; k<nz; ++k){
		for (int j=0; j<ny; ++j){
			for (int i=0; i<nx; ++i){
				gi = i*sx + j*sy + k*sz;
				ret.cells.emplace_back(new Cell);
				auto c = ret.cells.back();
				auto xleft = xfaces[gi]; xleft->left = c;
				auto xright = xfaces[gi+sx]; xright->right = c;
				auto yleft =  yfaces[gi]; yleft->left = c;
				auto yright = yfaces[gi+sy]; yright->right = c;
				auto zleft = zfaces[gi]; zleft->left = c;
				auto zright = zfaces[gi+sz]; zright->right = c;
				c->faces = {xleft, xright, yleft, yright, zleft, zright};
			}
		}
	}
	//Boundary Types
	auto setbt = [](shared_ptr<Face> f, int nor, int nol){
		if (f){
			if (!f->right) f->boundary_type = nor;
			if (!f->left) f->boundary_type = nol;
		}
	};
	for (auto f: xfaces) setbt(f, 1, 2);
	for (auto f: yfaces) setbt(f, 3, 4);
	for (auto f: zfaces) setbt(f, 5, 6);

	return ret;
}

HMGrid3D::Grid cns::SweepGrid2D(const GridGeom& g, const vector<double>& zcoords){
	HMGrid3D::Grid ret;
	int n2 = g.n_points(), nz = zcoords.size();
	//Vertices
	ShpVector<Vertex> vert; vert.reserve(n2*nz);
	for (int i=0; i<nz; ++i){
		auto p = g.get_point(i);
		for (auto z: zcoords) vert.emplace_back(new Vertex(p->x, p->y, z));
	}

	//Edges
	auto e2d = g.get_edges();
	ShpVector<Edge> xyedges(e2d.size()*nz), zedges((nz-1)*n2);
	//xy edges
	for (int i=0; i<nz; ++i){
		int j=0;
		for (auto& e: e2d){
			int i1 = i*e2d.size() + e.p1;
			int i2 = i*e2d.size() + e.p2;
			xyedges[i*nz+j].reset(new Edge(vert[i1], vert[i2]));
			++j;
		}
	}
	//z edges
	for (int i=0; i<n2; ++i){
		for (int j=0;
	}

	/*
	//Faces
	ShpVector<Face> xfaces(N), yfaces(N), zfaces(N);
	gi=0;
	for (int k=0; k<nz+1; ++k){
		for (int j=0; j<ny+1; ++j){
			for (int i=0; i<nx+1; ++i){
				//z face
				if (i<nx && j<ny){
					auto e1=yedges[gi], e2=xedges[gi+sy], e3=yedges[gi+sx], e4=xedges[gi];
					zfaces[gi].reset(new Face({e1, e2, e3, e4}));
				}
				//y face
				if (i<nx && k<nz){
					auto e1=xedges[gi], e2=zedges[gi+sx], e3=xedges[gi+sz], e4=zedges[gi];
					yfaces[gi].reset(new Face({e1, e2, e3, e4}));
				};
				//x face
				if (j<ny && k<nz){
					auto e1=zedges[gi], e2=yedges[gi+sz], e3=zedges[gi+sy], e4=yedges[gi];
					xfaces[gi].reset(new Face({e1, e2, e3, e4}));
				}
				++gi;
			}
		}
	}
	//Cells
	ret.cells.reserve(nx*ny*nz);
	for (int k=0; k<nz; ++k){
		for (int j=0; j<ny; ++j){
			for (int i=0; i<nx; ++i){
				gi = i*sx + j*sy + k*sz;
				ret.cells.emplace_back(new Cell);
				auto c = ret.cells.back();
				auto xleft = xfaces[gi]; xleft->left = c;
				auto xright = xfaces[gi+sx]; xright->right = c;
				auto yleft =  yfaces[gi]; yleft->left = c;
				auto yright = yfaces[gi+sy]; yright->right = c;
				auto zleft = zfaces[gi]; zleft->left = c;
				auto zright = zfaces[gi+sz]; zright->right = c;
				c->faces = {xleft, xright, yleft, yright, zleft, zright};
			}
		}
	}
	//Boundary Types
	auto setbt = [](shared_ptr<Face> f, int nor, int nol){
		if (f){
			if (!f->right) f->boundary_type = nor;
			if (!f->left) f->boundary_type = nol;
		}
	};
	for (auto f: xfaces) setbt(f, 1, 2);
	for (auto f: yfaces) setbt(f, 3, 4);
	for (auto f: zfaces) setbt(f, 5, 6);
	*/

	return ret;
}
