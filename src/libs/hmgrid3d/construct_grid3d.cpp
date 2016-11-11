#include "construct_grid3d.hpp"
#include "debug_grid2d.h"
#include "debug_grid3d.hpp"
#include "procgrid.h"
#include "revolve_grid3d.hpp"
#include "vtk_export_grid3d.hpp"
using namespace HMGrid3D;

namespace cns = Constructor;

HMGrid3D::SGrid cns::Cuboid(HMGrid3D::Vertex leftp, double lx, double ly, double lz, int nx, int ny, int nz){
	HMGrid3D::SGrid ret;
	//Vertices
	VertexData vert; vert.reserve(nx*ny*nz);
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
	ret.vvert = vert;
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
	std::copy_if(xedges.begin(), xedges.end(), std::back_inserter(ret.vedges),
			[](shared_ptr<Edge>& e)->bool{ return e != nullptr; });
	std::copy_if(yedges.begin(), yedges.end(), std::back_inserter(ret.vedges),
			[](shared_ptr<Edge>& e)->bool{ return e != nullptr; });
	std::copy_if(zedges.begin(), zedges.end(), std::back_inserter(ret.vedges),
			[](shared_ptr<Edge>& e)->bool{ return e != nullptr; });
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
	std::copy_if(xfaces.begin(), xfaces.end(), std::back_inserter(ret.vfaces),
			[](shared_ptr<Face>& f){ return f != nullptr; });
	std::copy_if(yfaces.begin(), yfaces.end(), std::back_inserter(ret.vfaces),
			[](shared_ptr<Face>& f){ return f != nullptr; });
	std::copy_if(zfaces.begin(), zfaces.end(), std::back_inserter(ret.vfaces),
			[](shared_ptr<Face>& f){ return f != nullptr; });
	//Cells
	ret.vcells.reserve(nx*ny*nz);
	for (int k=0; k<nz; ++k){
		for (int j=0; j<ny; ++j){
			for (int i=0; i<nx; ++i){
				gi = i*sx + j*sy + k*sz;
				ret.vcells.emplace_back(new Cell);
				auto c = ret.vcells.back();
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
			if (!f->has_right_cell()) f->boundary_type = nor;
			if (!f->has_left_cell()) f->boundary_type = nol;
		}
	};
	for (auto f: xfaces) setbt(f, 1, 2);
	for (auto f: yfaces) setbt(f, 3, 4);
	for (auto f: zfaces) setbt(f, 5, 6);

	ret.actualize_serial_data();
	return ret;
}

//HMGrid3D::SGrid cns::SphericalShell(HMGrid3D::Vertex center, double rinner, double router, double hr, double harc){
//        assert(router>rinner);
//        int nr = (router - rinner)/hr + 1;
//        int na = 2*M_PI*router/harc;
//        na += (4 - na%4);
//        GridGeom g2 = GGeom::Constructor::Ring(Point(center.x, center.y), router, rinner, na, nr);
//        vector<const ::Cell*> rmc;
//        for (int i=0; i<g2.n_cells(); ++i){
//                const ::Cell* c = g2.get_cell(i);
//                for (int ip=0; ip<c->dim(); ++ip){
//                        const Point* p = c->get_point(ip);
//                        if (ISLOWER(p->x, center.x)){
//                                rmc.push_back(c);
//                                break;
//                        }
//                }
//        }
//        GGeom::Modify::RemoveCells(g2, rmc);
//        vector<double> phi(na+1, 0);
//        for (int i=0; i<na; ++i) phi[i+1] = phi[i] + 360.0/na;

//        std::set<::Edge> edset = g2.get_edges();
//        vector<::Edge> ed(edset.begin(), edset.end());
//        auto side_bt = [&ed, &na](int i)->int{
//                if (ed[i].p1 < na) return 2;
//                else return 1;
//        };
//        HMGrid3D::SGrid ret = RevolveGrid2D(g2, phi,
//                        Point(center.x, center.y+router),
//                        Point(center.x, center.y-router),
//                        true, side_bt);
//        for (auto p: ret.vvert) p->z += center.z;
//        return ret;
//}

//HMGrid3D::SGrid cns::SphericalShell(HMGrid3D::Vertex center, double rinner, double router, double hr, double harc){
//        int nr_out = (int)(router/hr) + 1;
//        int nr_in = (int)(rinner/hr) + 1;
//        nr_out = 6;
//        nr_in = 4;
//        Vertex leftp(-nr_out, -nr_out, -nr_out);
//        auto cube = Cuboid(leftp, 2*nr_out, 2*nr_out, 2*nr_out, 2*nr_out, 2*nr_out, 2*nr_out);
//        //purge inner cells
//        ShpVector<HMGrid3D::Cell> shell_cells;
//        int g=0;
//        for (int k=0; k<2*nr_out; ++k)
//        for (int j=0; j<2*nr_out; ++j)
//        for (int i=0; i<2*nr_out; ++i){
//                if ((k>=nr_out-nr_in && k<nr_out+nr_in) &&
//                    (j>=nr_out-nr_in && j<nr_out+nr_in) &&
//                    (i>=nr_out-nr_in && i<nr_out+nr_in)){
//                } else shell_cells.push_back(cube.vcells[g]);
//                ++g;
//        }
//        std::cout<<shell_cells.size()<<std::endl;
//        SGrid ret;
//        ret.vcells = shell_cells;
//        ret.fill_from_cells();
//        //translate to sphere
//        for (auto p: ret.vvert){
//                double m = std::max(fabs(p->x), std::max(fabs(p->y), fabs(p->z)));
//                int layer = std::round(m);
//                double r = double(m - nr_in)/(nr_out-nr_in)*(router-rinner) + rinner;
//                double truerad = sqrt(p->x*p->x + p->y*p->y + p->z*p->z);
//                double coef = r/truerad;
//                p->x*=coef; p->y*=coef; p->z*=coef;
//        }
//        ret.actualize_serial_data();
//        return ret;
//}

namespace{

ShpVector<HMGrid3D::Edge> rv(const ShpVector<HMGrid3D::Edge>& e){
	ShpVector<HMGrid3D::Edge> ret(e);
	std::reverse(ret.begin(), ret.end());
	return ret;
};
ShpVector<HMGrid3D::Vertex> sorted_points(const ShpVector<HMGrid3D::Edge>& e){
	auto p0 = e[0]->first(), p1 = e[0]->last();
	if (p0 == e[1]->first() || p0 == e[1]->last()) std::swap(p0, p1);
	ShpVector<HMGrid3D::Vertex> ret{p0, p1};
	for (int i=1; i<e.size(); ++i){
		auto pnext = e[i]->first();
		if (pnext == ret.back()) pnext = e[i]->last();
		ret.push_back(pnext);
	}

	return ret;
}


ShpVector<HMGrid3D::Edge> build_edges(
		shared_ptr<HMGrid3D::Vertex> p0,
		shared_ptr<HMGrid3D::Vertex> p1, int n){
	ShpVector<HMGrid3D::Edge> ret;
	shared_ptr<HMGrid3D::Vertex> last = p0;
	for (int i=0; i<n; ++i){
		double w = (i+1.)/n;
		shared_ptr<HMGrid3D::Vertex> next(new HMGrid3D::Vertex());
		next->x = (1.-w)*p0->x + w*p1->x;
		next->y = (1.-w)*p0->y + w*p1->y;
		next->z = (1.-w)*p0->z + w*p1->z;
		ret.emplace_back(new HMGrid3D::Edge(last, next));
		last = next;
	}
	ret.back()->vertices.back() = p1;
	return ret;
}
ShpVector<HMGrid3D::Face> build_faces(
		const ShpVector<HMGrid3D::Edge>& e0,
		const ShpVector<HMGrid3D::Edge>& e1,
		const ShpVector<HMGrid3D::Edge>& e2,
		const ShpVector<HMGrid3D::Edge>& e3){
	int n=e0.size();
	ShpVector<HMGrid3D::Vertex> p0 = sorted_points(e0);
	ShpVector<HMGrid3D::Vertex> p1 = sorted_points(e1);
	ShpVector<HMGrid3D::Vertex> p2 = sorted_points(rv(e2));
	ShpVector<HMGrid3D::Vertex> p3 = sorted_points(rv(e3));

	vector<ShpVector<HMGrid3D::Vertex>> points(n+1);
	points[0] = p3; points[n] = p1;
	for (int i=1; i<n; ++i) points[i] = sorted_points(build_edges(p0[i], p2[i], n));

	vector<ShpVector<HMGrid3D::Edge>> xedges(n+1), yedges(n+1);
	for (int i=0; i<n+1; ++i) {xedges[i].resize(n+1); yedges[i].resize(n+1); }
	for (int j=0; j<n+1; ++j)
	for (int i=0; i<n+1; ++i){
		if (i!=n){
			if (j==0) xedges[i][j] = e0[i];
			else if (j==n) xedges[i][j] = e2[n-1-i];
			else xedges[i][j].reset(new HMGrid3D::Edge(points[i][j], points[i+1][j]));
		}

		if (j!=n){
			if (i==0) yedges[i][j] = e3[n-1-j];
			else if (i==n) yedges[i][j] = e1[j];
			else yedges[i][j].reset(new HMGrid3D::Edge(points[i][j], points[i][j+1]));
		}
	}

	ShpVector<HMGrid3D::Face> ret;
	for (int j=0; j<n; ++j)
	for (int i=0; i<n; ++i){
		auto f = aa::add_shared(ret, HMGrid3D::Face());
		f->edges = { xedges[i][j], yedges[i+1][j], xedges[i][j+1], yedges[i][j] };
	}
	return ret;
}
HMGrid3D::Surface build_cubic_surface(int n){
	//corner points
	shared_ptr<HMGrid3D::Vertex> corner000(new Vertex(0, 0, 0));
	shared_ptr<HMGrid3D::Vertex> corner001(new Vertex(0, 0, 1));
	shared_ptr<HMGrid3D::Vertex> corner010(new Vertex(0, 1, 0));
	shared_ptr<HMGrid3D::Vertex> corner011(new Vertex(0, 1, 1));
	shared_ptr<HMGrid3D::Vertex> corner100(new Vertex(1, 0, 0));
	shared_ptr<HMGrid3D::Vertex> corner101(new Vertex(1, 0, 1));
	shared_ptr<HMGrid3D::Vertex> corner110(new Vertex(1, 1, 0));
	shared_ptr<HMGrid3D::Vertex> corner111(new Vertex(1, 1, 1));
	
	//edges
	ShpVector<HMGrid3D::Edge> yz00 = build_edges(corner000, corner100, n);
	ShpVector<HMGrid3D::Edge> yz01 = build_edges(corner001, corner101, n);
	ShpVector<HMGrid3D::Edge> yz10 = build_edges(corner010, corner110, n);
	ShpVector<HMGrid3D::Edge> yz11 = build_edges(corner011, corner111, n);

	ShpVector<HMGrid3D::Edge> zx00 = build_edges(corner000, corner010, n);
	ShpVector<HMGrid3D::Edge> zx01 = build_edges(corner100, corner110, n);
	ShpVector<HMGrid3D::Edge> zx10 = build_edges(corner001, corner011, n);
	ShpVector<HMGrid3D::Edge> zx11 = build_edges(corner101, corner111, n);

	ShpVector<HMGrid3D::Edge> xy00 = build_edges(corner000, corner001, n);
	ShpVector<HMGrid3D::Edge> xy01 = build_edges(corner010, corner011, n);
	ShpVector<HMGrid3D::Edge> xy10 = build_edges(corner100, corner101, n);
	ShpVector<HMGrid3D::Edge> xy11 = build_edges(corner110, corner111, n);

	//faces
	ShpVector<Face> xy_face_0 = build_faces(zx00, yz10, rv(zx01), rv(yz00));
	ShpVector<Face> xy_face_1 = build_faces(yz01, zx11, rv(yz11), rv(zx10)); 
	ShpVector<Face> yz_face_0 = build_faces(xy00, zx10, rv(xy01), rv(zx00));
	ShpVector<Face> yz_face_1 = build_faces(zx01, xy11, rv(zx11), rv(xy10));
	ShpVector<Face> zx_face_0 = build_faces(yz00, xy10, rv(yz01), rv(xy00));
	ShpVector<Face> zx_face_1 = build_faces(xy01, yz11, rv(xy11), rv(yz10));

	//return
	Surface ret;
	std::copy(xy_face_0.begin(), xy_face_0.end(), std::back_inserter(ret.faces));
	std::copy(xy_face_1.begin(), xy_face_1.end(), std::back_inserter(ret.faces));
	std::copy(yz_face_0.begin(), yz_face_0.end(), std::back_inserter(ret.faces));
	std::copy(yz_face_1.begin(), yz_face_1.end(), std::back_inserter(ret.faces));
	std::copy(zx_face_0.begin(), zx_face_0.end(), std::back_inserter(ret.faces));
	std::copy(zx_face_1.begin(), zx_face_1.end(), std::back_inserter(ret.faces));

	return ret;
}
}

HMGrid3D::SGrid cns::SphericalShell(HMGrid3D::Vertex center, double rinner, double router, double hr, double harc){
	//create cubic surface
	int n = (int)(2*M_PI*router/harc/4) + 1;
	HMGrid3D::Surface cubic = build_cubic_surface(n);
	//project to sphere
	for (auto v: cubic.allvertices()){
		v->x -= 0.5; v->y -= 0.5; v->z -=0.5;
		double rad = sqrt(v->x*v->x + v->y*v->y + v->z*v->z);
		double coef = rinner/rad;
		v->x*=coef; v->y*=coef; v->z*=coef;
	}

	//build cells by normal
	auto calc_normal = [](const HMGrid3D::Vertex& v){return v;};
	int nr = (router-rinner)/hr + 1;
	vector<double> rads(nr+1, 0); rads.back() = router-rinner;
	for (int i=1; i<nr; ++i) rads[i] = (double)i*rads.back()/nr;
	GridData ret = NormalToSurface(cubic, rads, calc_normal);
	for (auto p: ret.vvert){
		p->x += center.x; p->y += center.y; p->z += center.z;
	}

	return SGrid(ret);
}

HMGrid3D::SGrid cns::SweepGrid2D(const GridGeom& g, const vector<double>& zcoords){
	return SweepGrid2D(g, zcoords,
			[](int){ return 1; },
			[](int){ return 2; },
			[](int){ return 3; });
}

HMGrid3D::SGrid cns::SweepGrid2D(const GridGeom& g, const vector<double>& zcoords,
		std::function<int(int)> bottom_bt,
		std::function<int(int)> top_bt,
		std::function<int(int)> side_bt){

	HMGrid3D::SGrid ret;

	//Needed Data
	auto e2d = g.get_edges();
	vector<::Edge> e2dvector(e2d.begin(), e2d.end());
	int n2p = g.n_points(), n2c = g.n_cells(), n2e = e2d.size(), nz = zcoords.size();
	vector<vector<int>> cell_edges2d(n2c);
	for (int i=0; i<n2c; ++i){
		auto c2d = g.get_cell(i);
		for (int j=0; j<c2d->dim(); ++j){
			int p1 = c2d->get_point(j)->get_ind();
			int p2 = c2d->get_point(j+1)->get_ind();
			::Edge e(p1, p2);
			int fnd = std::find(e2dvector.begin(), e2dvector.end(), e) - e2dvector.begin();
			cell_edges2d[i].push_back(fnd);
		}
	}

	//Vertices
	ShpVector<Vertex>& vert = ret.vvert; vert.reserve(n2p*nz);
	{
		for (int i=0; i<nz; ++i){
			double z = zcoords[i];
			for (int j=0; j<n2p; ++j){
				auto p = g.get_point(j);
				vert.emplace_back(new Vertex(p->x, p->y, z));
			}
		}
	}

	//Edges
	ShpVector<Edge> xyedges(n2e*nz), zedges((nz-1)*n2p);
	{
		//xy edges
		int j=0;
		for (int i=0; i<nz; ++i){
			for (auto& e: e2dvector){
				int i1 = i*n2p + e.p1;
				int i2 = i*n2p + e.p2;
				xyedges[j++].reset(new Edge(vert[i1], vert[i2]));
			}
		}
		//z edges
		int k=0;
		for (int i=0; i<nz-1; ++i){
			for (int j=0; j<n2p; ++j){
				int i1 = i*n2p + j;
				int i2 = (i+1)*n2p + j;
				zedges[k++].reset(new Edge(vert[i1], vert[i2]));
			}
		}
	}
	ret.vedges.insert(ret.vedges.end(), xyedges.begin(), xyedges.end());
	ret.vedges.insert(ret.vedges.end(), zedges.begin(), zedges.end());

	//Faces
	ShpVector<Face> xyfaces; xyfaces.reserve(n2c*nz);
	ShpVector<Face> zfaces; zfaces.reserve(n2e*(nz-1));
	{
		//xyfaces
		for (int i=0; i<nz; ++i){
			for (int j=0; j<n2c; ++j){
				Face* f = new Face();
				for (auto k: cell_edges2d[j]){
					f->edges.push_back(xyedges[i*n2e+k]);
				}
				xyfaces.emplace_back(f);
			}
		}
		//zfaces
		for (int i=0; i<nz-1; ++i){
			for (int j=0; j<n2e; ++j){
				int i12d = e2dvector[j].p1;
				int i22d = e2dvector[j].p2;
				Face* f = new Face();
				auto e1 = xyedges[j + i*n2e];
				auto e2 = zedges[i22d + i*n2p];
				auto e3 = xyedges[j+ (i+1)*n2e];
				auto e4 = zedges[i12d + i*n2p];
				f->edges = {e1, e2, e3, e4};
				zfaces.emplace_back(f);
			}
		}
	}
	ret.vfaces.insert(ret.vfaces.end(), xyfaces.begin(), xyfaces.end());
	ret.vfaces.insert(ret.vfaces.end(), zfaces.begin(), zfaces.end());

	//Cells
	{
		ret.vcells.reserve(n2c*(nz-1));
		for (int i=0; i<nz-1; ++i){
			for (int j=0; j<n2c; ++j){
				ret.vcells.emplace_back(new Cell());
				auto& c = ret.vcells.back();
				int ifc = i*n2c + j;
				auto bot = xyfaces[ifc]; bot->right = c;
				auto top = xyfaces[ifc + n2c]; top->left = c;
				c->faces = {bot, top};
				for (int k: cell_edges2d[j]){
					auto f1 = zfaces[k + n2e*i];
					c->faces.push_back(f1);
					bool isleft = (e2dvector[k].cell_left == j);
					if (isleft) f1->left = c;
					else f1->right = c;
				}
			}
		}
	}

	//Boundary Types
	{
		//top, bottom
		int j = n2c*(nz-1);
		for (int i=0; i<n2c; ++i){
			xyfaces[i]->boundary_type = bottom_bt(i);
			xyfaces[j++]->boundary_type = top_bt(i);
		}
		//sides
		for (int i=0; i<n2e; ++i) if (e2dvector[i].is_boundary()){
			int bt = side_bt(i);
			for (int k=0; k<nz-1; ++k){
				zfaces[i + k*n2e]->boundary_type = bt;
			}
		}
	}

	ret.actualize_serial_data();
	return ret;
}



HMGrid3D::SGrid cns::Copy::ShallowVertices(const SGrid& g){
	HMGrid3D::SGrid ret(g);
	ret.actualize_data();
	enumerate_ids_pvec(ret.vvert);
	for (auto& e: ret.vedges){
		for (auto& n: e->vertices) n = g.vvert[n->id];
	}
	ret.vvert = g.vvert;
	return ret;
}

HMGrid3D::GridData cns::NormalToSurface(HMGrid3D::Surface& src, const vector<double>& part,
		std::function<HMGrid3D::Vertex(const HMGrid3D::Vertex&)> normal){
	ShpVector<HMGrid3D::Edge> srcedges = src.alledges();
	ShpVector<HMGrid3D::Vertex> srcpoints = src.allvertices();
	//build surfaces, no sweep
	ShpVector<HMGrid3D::Surface> _sdata(part.size()-1);
	vector<HMGrid3D::Surface*> srfs(part.size());
	srfs[0] = &src;
	for (int i=0; i<part.size()-1; ++i){
		_sdata[i].reset(new Surface(src));
		srfs[i+1] = _sdata[i].get();
		ReallocateAll(srfs[i+1]->faces);
	}
	//assemble unique edges and points
	vector<ShpVector<HMGrid3D::Edge>> alledges(srfs.size());
	vector<ShpVector<HMGrid3D::Vertex>> allpoints(srfs.size());
	alledges[0] = srcedges; allpoints[0] = srcpoints;
	for (int i=1; i<srfs.size(); ++i){
		alledges[i] = srfs[i]->alledges();
		allpoints[i] = srfs[i]->allvertices();
	}
	//sweep nodes
	for (int j=0; j<srcpoints.size(); ++j){
		HMGrid3D::Vertex norm = normal(*srcpoints[j]);
		double nlen = norm.x*norm.x + norm.y*norm.y + norm.z*norm.z;
		nlen = sqrt(nlen);
		for (int i=0; i<srfs.size()-1; ++i){
			allpoints[i+1][j]->x += part[i+1]/nlen*norm.x;
			allpoints[i+1][j]->y += part[i+1]/nlen*norm.y;
			allpoints[i+1][j]->z += part[i+1]/nlen*norm.z;
		}
	}

	//build vertical edges
	vector<ShpVector<HMGrid3D::Edge>> vert_edges(srfs.size()-1);
	for (int i=0; i<srfs.size()-1; ++i){
		vert_edges[i].resize(srcpoints.size());
		for (int j=0; j<srcpoints.size(); ++j){
			vert_edges[i][j].reset(new HMGrid3D::Edge(
				allpoints[i][j], allpoints[i+1][j]));
		}
	}

	//ordering of edges within faces
	vector<vector<bool>> face_edge_ordering(src.n_faces());
	for (int i=0; i<src.n_faces(); ++i){
		face_edge_ordering[i].resize(src.faces[i]->n_edges());
		auto sv = src.faces[i]->sorted_vertices();
		std::rotate(sv.begin(), sv.begin()+1, sv.end()); //to start with 1st edge 
		for (int j=0; j<src.faces[i]->n_edges(); ++j){
			face_edge_ordering[i][j] = (sv[j] == src.faces[i]->edges[j]->first());
		}
	}

	//build cells
	enumerate_ids_pvec(srcedges);
	enumerate_ids_pvec(srcpoints);
	HMGrid3D::GridData ret;
	for (int i=0; i<srfs.size()-1; ++i){
		HMGrid3D::Surface* s0 = srfs[i];
		HMGrid3D::Surface* s1 = srfs[i+1];
		ShpVector<HMGrid3D::Face> vertfaces(alledges[i].size(), 0);
		for (int j=0; j<src.n_faces(); ++j){
			ret.vcells.emplace_back(new HMGrid3D::Cell());
			auto c = ret.vcells.back();
			//bottom, top
			c->faces.push_back(s0->faces[j]); s0->faces[j]->right = c;
			c->faces.push_back(s1->faces[j]); s1->faces[j]->left = c;
			//side 
			for (int k=0; k<src.faces[i]->n_edges(); ++k){
				int eid = src.faces[j]->edges[k]->id;
				shared_ptr<HMGrid3D::Face> vf;
				if (vertfaces[eid]){
					vf = vertfaces[eid];
					vf->right = c;
				} else {
					vf.reset(new Face());
					int pstart = src.faces[j]->edges[k]->first()->id;
					int pend = src.faces[j]->edges[k]->last()->id;
					if (face_edge_ordering[j][k]) std::swap(pstart, pend);
					shared_ptr<HMGrid3D::Edge> vedge1 = vert_edges[i][pend];
					shared_ptr<HMGrid3D::Edge> vedge2 = vert_edges[i][pstart];
					vf->edges = {s0->faces[j]->edges[k], vedge1, s1->faces[j]->edges[k], vedge2};
					vf->left = c;
					vertfaces[eid] = vf;
				}
				c->faces.push_back(vf);
			}
		}
	}
	ret.fill_from_cells();
	return ret;
}
