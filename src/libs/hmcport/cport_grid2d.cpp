#include "hmcport.h"
#include "cport_grid2d.h"
#include "primitives2d.hpp"
#include "infogrid.hpp"
#include "cont_assembler.hpp"
#include "c2cpp_helper.hpp"
#include "buildgrid.hpp"
#include "finder2d.hpp"
#include "import2d_hm.hpp"
#include "pebi.hpp"
#include "tscaler.hpp"
#include "trigrid.hpp"
#include "gridmap.hpp"
#include "rectangle_grid_builder.hpp"
#include "circrect.hpp"
#include "hmblay.hpp"
#include "cport_cont2d.h"
#include "healgrid.hpp"
#include "modgrid.hpp"
#include "unite_grids.hpp"
#include "export2d_fluent.hpp"
#include "export2d_tecplot.hpp"
#include "export2d_hm.hpp"


int g2_dims(void* obj, int* ret){
	try{
		auto g = static_cast<HM2D::GridData*>(obj);
		ret[0] = g->vvert.size();
		ret[1] = g->vedges.size();
		ret[2] = g->vcells.size();
		return HMSUCCESS;
	} catch (std::exception& e){
		std::cout<<e.what()<<std::endl;
		return HMERROR;
	}
}
int g2_area(void* obj, double* ret){
	try{
		auto g = static_cast<HM2D::GridData*>(obj);
		*ret = HM2D::Grid::Area(*g);
		return HMSUCCESS;
	} catch (std::exception& e){
		std::cout<<e.what()<<std::endl;
		return HMERROR;
	}
}
int g2_bnd_dims(void* obj, int* ret){
	try{
		auto g = static_cast<HM2D::GridData*>(obj);
		auto c = HM2D::ECol::Assembler::GridBoundary(*g);
		return c2_dims(&c, ret);
	} catch (std::exception& e){
		std::cout<<e.what()<<std::endl;
		return HMERROR;
	}
}
int g2_bnd_length(void* obj, double* ret){
	try{
		auto g = static_cast<HM2D::GridData*>(obj);
		auto c = HM2D::ECol::Assembler::GridBoundary(*g);
		return c2_length(&c, ret);
	} catch (std::exception& e){
		std::cout<<e.what()<<std::endl;
		return HMERROR;
	}
}
int g2_skewness(void* obj, double threshold, double* maxskew, int* maxskewindex,
		int* badnum, int** badindex, double** badvals){
	try{
		vector<double> sc = HM2D::Grid::Skewness(*static_cast<HM2D::GridData*>(obj));
		*maxskew = 0.0;
		*maxskewindex = -1;
		*badnum = 0;
		vector<int> badindex_;
		vector<double> badvals_;
		for (int i=0; i<sc.size(); ++i){
			if (sc[i]>*maxskew){
				*maxskewindex = i;
				*maxskew = sc[i];
			}
			if (sc[i]>=threshold){
				badindex_.push_back(i);
				badvals_.push_back(sc[i]);
			}
		}
		*badnum = badindex_.size();
		*badindex = new int[*badnum];
		*badvals = new double[*badnum];
		std::copy(badindex_.begin(), badindex_.end(), *badindex);
		std::copy(badvals_.begin(), badvals_.end(), *badvals);
		return HMSUCCESS;
	} catch (std::exception& e){
		std::cout<<e.what()<<std::endl;
		return HMERROR;
	}
}
int g2_deepcopy(void* obj, void** ret){
	try{
		HM2D::GridData ret_;
		HM2D::DeepCopy(*static_cast<HM2D::GridData*>(obj), ret_);
		c2cpp::to_pp(ret_, ret);
		return HMSUCCESS;
	} catch (std::exception& e){
		std::cout<<e.what()<<std::endl;
		return HMERROR;
	}
}
int g2_free(void* obj){
	try{
		delete static_cast<HM2D::GridData*>(obj);
		return HMSUCCESS;
	} catch (std::exception& e){
		std::cout<<e.what()<<std::endl;
		return HMERROR;
	}
}
int g2_concatenate(int nobjs, void** objs, void** ret){
	try{
		auto gg = c2cpp::to_pvec<HM2D::GridData>(nobjs, objs);
		HM2D::GridData ret_;
		for (auto g: gg){
			ret_.vvert.insert(ret_.vvert.end(), g->vvert.begin(), g->vvert.end());
			ret_.vedges.insert(ret_.vedges.end(), g->vedges.begin(), g->vedges.end());
			ret_.vcells.insert(ret_.vcells.end(), g->vcells.begin(), g->vcells.end());
		}
		c2cpp::to_pp(ret_, ret);
		return HMSUCCESS;
	} catch (std::exception& e){
		std::cout<<e.what()<<std::endl;
		return HMERROR;
	}
}
int g2_move(void* obj, double* dx){
	try{
		for (auto& v: static_cast<HM2D::GridData*>(obj)->vvert){
			v->x += dx[0];
			v->y += dx[1];
		}
		return HMSUCCESS;
	} catch (std::exception& e){
		std::cout<<e.what()<<std::endl;
		return HMERROR;
	}
}
int g2_scale(void* obj, double* pc, double* p0){
	try{
		auto g = static_cast<HM2D::GridData*>(obj);
		double xc = pc[0] / 100., yc = pc[1] / 100.;
		for (auto& p: g->vvert){
			p->x -= p0[0]; p->x *= xc;
			p->y -= p0[1]; p->y *= yc;
		}
		return HMSUCCESS;
	} catch (std::exception& e){
		std::cout<<e.what()<<std::endl;
		return HMERROR;
	}
}
int g2_reflect(void* obj, double* v0, double* v1){
	try{
		auto g = static_cast<HM2D::GridData*>(obj);
		double lx = v1[0] - v0[0], ly = v1[1] - v0[1];
		double r2 = sqrt(lx*lx + ly*ly);
		lx /= r2; ly /= r2;
		double A[4] = {lx*lx - ly*ly, 2*lx*ly,
		               2*lx*ly, ly*ly-lx*lx};
		for (auto& p: g->vvert){
			double a = p->x - v0[0];
			double b = p->y - v0[1];
			p->x = A[0]*a + A[1]*b + v0[0];
			p->y = A[2]*a + A[3]*b + v0[1];
		}
		for (auto& c: g->vcells){
			std::reverse(c->edges.begin(), c->edges.end());
		}
		for (auto& e: g->vedges){
			std::swap(e->left, e->right);
		}
		return HMSUCCESS;
	} catch (std::exception& e){
		std::cout<<e.what()<<std::endl;
		return HMERROR;
	}
}
int g2_rotate(void* obj, double* p0, double a){
	try{
		auto g = static_cast<HM2D::GridData*>(obj);
		double cs = cos(a / 180. * M_PI), sn = sin(a / 180. * M_PI);
		for (auto& p: g->vvert){
			double a = p->x - p0[0];
			double b = p->y - p0[1];
			p->x = a*cs - b*sn + p0[0];
			p->y = a*sn + b*cs + p0[1];
		}
		return HMSUCCESS;
	} catch (std::exception& e){
		std::cout<<e.what()<<std::endl;
		return HMERROR;
	}
}
int g2_from_points_edges(int npoints, double* points, int neds, int* eds, void** ret){
	try{
		vector<double> pts(npoints*2);
		std::copy(points, points+2*npoints, pts.begin());
		vector<int> edgevert, edgecell, bnds;
		int* it = eds;
		for (int i=0; i<neds; ++i){
			edgevert.push_back(*it++);
			edgevert.push_back(*it++);
			edgecell.push_back(*it++);
			edgecell.push_back(*it++);
			bnds.push_back(*it++);
		}
		HM2D::GridData ret_ = HM2D::Import::GridFromTabs(pts, edgevert, edgecell);
		for (int i=0; i<neds; ++i) ret_.vedges[i]->boundary_type = bnds[i];
		c2cpp::to_pp(ret_, ret);
		return HMSUCCESS;
	} catch (std::exception& e){
		std::cout<<e.what()<<std::endl;
		return HMERROR;
	}
}
int g2_from_points_cells(int npoints, double* points, int ncells, int* cellsizes, int* cellvert,
		int nbedges, int* bedges, void** ret){
	try{
		//build grid
		vector<int> cv;
		int* it = cellvert;
		for (int i=0; i<ncells; ++i){
			cv.push_back(cellsizes[i]);
			for (int j=0; j<cellsizes[i]; ++j){
				cv.push_back(*it++);
			}
		}
		HM2D::GridData ret_ = HM2D::Grid::Constructor::FromRaw(
			npoints, ncells, points, &cv[0], -1);
		//boundary types
		auto be = HM2D::ECol::Assembler::GridBoundary(ret_);
		auto efinder = HM2D::Finder::EdgeFinder(be);
		for (int i=0; i<nbedges; ++i){
			auto v1 = ret_.vvert[bedges[3*i]].get();
			auto v2 = ret_.vvert[bedges[3*i+1]].get();
			auto fnd = efinder.find(v1, v2);
			if (std::get<0>(fnd)){
				std::get<0>(fnd)->boundary_type = bedges[3*i+2];
			}
		}
		c2cpp::to_pp(ret_, ret);
		return HMSUCCESS;
	} catch (std::exception& e){
		std::cout<<e.what()<<std::endl;
		return HMERROR;
	}
}
int g2_rect_grid(int nx, double* xdata, int ny, double* ydata, int* bnds, void** ret){
	try{
		//build grid
		vector<double> xdata_(nx), ydata_(ny);
		std::copy(xdata, xdata+nx, xdata_.begin());
		std::copy(ydata, ydata+ny, ydata_.begin());
		HM2D::GridData ret_ = HM2D::Grid::Constructor::RectGrid(xdata_, ydata_);
		//boundary types
		for (auto e: HM2D::ECol::Assembler::GridBoundary(ret_)){
			Point cnt = e->center();
			std::array<double, 4> d {
				fabs(cnt.y - ydata_[0]),
				fabs(cnt.x - xdata_.back()),
				fabs(cnt.y - ydata_.back()),
				fabs(cnt.x - xdata_[0])};
			switch (std::min_element(d.begin(), d.end()) - d.begin()){
				case 0: e->boundary_type = bnds[0]; break;
				case 1: e->boundary_type = bnds[1]; break;
				case 2: e->boundary_type = bnds[2]; break;
				case 3: e->boundary_type = bnds[3]; break;
			};
		}
		c2cpp::to_pp(ret_, ret);
		return HMSUCCESS;
	} catch (std::exception& e){
		std::cout<<e.what()<<std::endl;
		return HMERROR;
	}
}
int g2_circ_grid(double* p0, int nr, double* rdata, int na, double* adata,
		int istrian, int bnd, void** ret){
	try{
		//build unit grid
		HM2D::GridData ret_ = HM2D::Grid::Constructor::Circle(
			Point(0, 0), 1., na-1, nr-1, istrian);
		//adopt coordinates
		int k = 0;
		for (int ir = nr-1; ir>0; --ir)
		for (int ia = 0; ia<na-1; ++ia){
			HM2D::Vertex* v = ret_.vvert[k++].get();
			v->x = rdata[ir]*cos(adata[ia]/180.*M_PI);
			v->y = rdata[ir]*sin(adata[ia]/180.*M_PI);
		}
		for (auto v: ret_.vvert){
			v->x += p0[0];
			v->y += p0[1];
		}
		//boundary
		for (auto e: HM2D::ECol::Assembler::GridBoundary(ret_)){
			e->boundary_type = bnd;
		}
		c2cpp::to_pp(ret_, ret);
		return HMSUCCESS;
	} catch (std::exception& e){
		std::cout<<e.what()<<std::endl;
		return HMERROR;
	}
}
int g2_ring_grid(double* p0, int nr, double* rdata, int na, double* adata, int* bnds, void** ret){
	try{
		Point p(p0[0], p0[1]);
		//build unit grid
		HM2D::GridData ret_ = HM2D::Grid::Constructor::Ring(
			Point(0, 0), 1., 0.1, na-1, nr-1);
		//adopt coordinates
		int k = 0;
		for (int ir = nr-1; ir>=0; --ir)
		for (int ia = 0; ia<na-1; ++ia){
			HM2D::Vertex* v = ret_.vvert[k++].get();
			v->x = rdata[ir]*cos(adata[ia]/180.*M_PI) + p.x;
			v->y = rdata[ir]*sin(adata[ia]/180.*M_PI) + p.y;
		}
		//boundary
		for (auto e: HM2D::ECol::Assembler::GridBoundary(ret_)){
			double dc = Point::dist(*e->pfirst(), p);
			double d1 = fabs(dc-rdata[0]);
			double d2 = fabs(dc-rdata[nr-1]);
			if (d1 < d2) e->boundary_type = bnds[0];
			else e->boundary_type = bnds[1];
		}
		c2cpp::to_pp(ret_, ret);
		return HMSUCCESS;
	} catch (std::exception& e){
		std::cout<<e.what()<<std::endl;
		return HMERROR;
	}
}
int g2_tri_grid(double* verts, int ne, int* bnds, void** ret){
	try{
		Point p0(verts[0], verts[1]);
		Point p1(verts[2], verts[3]);
		Point p2(verts[4], verts[5]);
		//check points order
		double x0 = p0.x - p2.x, y0 = p0.y - p2.y;
		double x1 = p1.x - p2.x, y1 = p1.y - p2.y;
		if (x0 * y1 - y0 * x1 > 0){
			std::swap(p0, p2);
			std::swap(bnds[0], bnds[2]);
		}

		//build points
		vector<HM2D::VertexData> plines;
		for (int i=0; i<ne+1; ++i){
			double w1 = double (i) / ne;
			double x1 = w1 * p2.x + (1 - w1) * p0.x;
			double y1 = w1 * p2.y + (1 - w1) * p0.y;
			double x2 = w1 * p2.x + (1 - w1) * p1.x;
			double y2 = w1 * p2.y + (1 - w1) * p1.y;
			plines.emplace_back();
			for (int j=0; j<ne+1-i; ++j){
				double w2 = (ne != i) ? double(j)/(ne-i) : 1;
				double x3 = (1 - w2) * x1 + w2 * x2;
				double y3 = (1 - w2) * y1 + w2 * y2;
				plines.back().push_back(std::make_shared<HM2D::Vertex>(x3, y3));
			}
		}
		HM2D::VertexData allpoints;
		for (auto& it1: plines) allpoints.insert(allpoints.end(), it1.begin(), it1.end());
		//build cells
		aa::enumerate_ids_pvec(allpoints);
		vector<vector<int>> allcells;
		for (int i=0; i<ne; ++i){
			int p1 = plines[i][0]->id;
			int p2 = plines[i + 1][0]->id;
			int p3 = plines[i][1]->id;
			allcells.emplace_back(vector<int>{p1, p2, p3});
			for (int j=0; j<ne-i-1; ++j){
				int p1 = plines[i][j + 1]->id;
				int p2 = plines[i + 1][j]->id;
				int p3 = plines[i + 1][j + 1]->id;
				int p4 = plines[i][j + 2]->id;
				allcells.emplace_back(vector<int>{p1, p2, p3, p4});
			}
		}
		//assemble grid
		HM2D::GridData ret_ = HM2D::Grid::Constructor::FromTab(
			std::move(allpoints), allcells);
		//boundaries
		for (auto e: HM2D::ECol::Assembler::GridBoundary(ret_)){
			double d0 = Point::meas_line(e->center(), p0, p1);
			double d1 = Point::meas_line(e->center(), p1, p2);
			double d2 = Point::meas_line(e->center(), p2, p0);
			if (d0 < d1 && d0 < d2){
				e->boundary_type = bnds[0];
			} else if (d1 < d0 && d1 < d2){
				e->boundary_type = bnds[1];
			} else {
				e->boundary_type = bnds[2];
			}
		}
		c2cpp::to_pp(ret_, ret);
		return HMSUCCESS;
	} catch (std::exception& e){
		std::cout<<e.what()<<std::endl;
		return HMERROR;
	}
}
int g2_hex_grid(const char* areatype, double* area, double crad, int strict, void** ret){
	try{
		ScaleBase sc;
		HM2D::GridData ret_;
		switch (std::map<std::string, int>{
				{"rect", 1},
				{"hex", 2} }[areatype]){
		case 1:
			sc = BoundingBox(area[0], area[1], area[2], area[3]).to_scale();
			ret_ = HM2D::Grid::Constructor::RegularHexagonal(
				Point(0, 0),
				Point((area[2]-area[0])/sc.L, (area[3]-area[1])/sc.L),
				crad/sc.L, (bool)strict);
			break;
		case 2:
			sc = ScaleBase(area[0], area[1], area[2]);
			ret_ = HM2D::Grid::Constructor::RegularHexagonal(
				Point(0, 0),
				area[2]/sc.L,
				crad/sc.L,
				(bool)strict);
			break;
		default:
			throw std::runtime_error("unknown algo");
		}
		HM2D::Unscale(ret_.vvert, sc);
		c2cpp::to_pp(ret_, ret);
		return HMSUCCESS;
	} catch (std::exception& e){
		std::cout<<e.what()<<std::endl;
		return HMERROR;
	}
}
int g2_extract_contour(void* obj, void** ret){
	try{
		auto cont = HM2D::ECol::Assembler::GridBoundary(
			*static_cast<HM2D::GridData*>(obj));
		return c2_deepcopy(&cont, ret);
	} catch (std::exception& e){
		std::cout<<e.what()<<std::endl;
		return HMERROR;
	}
}
int g2_set_btypes(void* obj, int* bndlist, int for_all_edges){
	try{
		auto g = static_cast<HM2D::GridData*>(obj);
		if (for_all_edges){
			return c2_set_btypes(&g->vedges, bndlist);
		} else {
			auto cont = HM2D::ECol::Assembler::GridBoundary(*g);
			return c2_set_btypes(&cont, bndlist);
		}
	} catch (std::exception& e){
		std::cout<<e.what()<<std::endl;
		return HMERROR;
	}
}
int g2_unstructured_fill(void* domain, void* constraint,
		int nembpts, double* embpts, const char* filler, void** ret){
	using namespace HM2D::Mesher;
	try{
		HM2D::EdgeData* dom = static_cast<HM2D::EdgeData*>(domain);
		HM2D::EdgeData* con = static_cast<HM2D::EdgeData*>(constraint);

		Autoscale::D2 sc(dom);
		if (con!=0) sc.add_data(con);
		std::map<Point, double> ep;
		for (int i=0; i<nembpts; ++i){
			Point p(embpts[3*i], embpts[3*i+1]);
			double sz = embpts[3*i+2];
			sc.scale(p);
			sc.scale(sz);
			ep.emplace(p, sz);
		}

		auto etree = HM2D::Contour::Tree::Assemble(*dom);
		if (etree.roots().size() == 0)
			throw std::runtime_error("Failed to find bounding contour");

		vector<HM2D::EdgeData> cas;
		if (con!=0) cas = HM2D::Contour::Assembler::AllContours(*con);
		for (auto& c: cas){ etree.add_detached_contour(c); }

		HM2D::GridData ret_;
		switch (std::map<std::string, int>{
				{"3", 3},
				{"4", 4},
				{"pebi", 6} }[filler]){
			case 3: ret_ = UnstructuredTriangle(etree, ep); break;
			case 4: ret_ = UnstructuredTriangleRecomb(etree, ep); break;
			case 6: ret_ = HM2D::Mesher::UnstructuredPebi(etree, ep); break;
			default: throw std::runtime_error("unknown unstructured fill algorithm");
		};
		sc.unscale(&ret_);
		c2cpp::to_pp(ret_, ret);
		return HMSUCCESS;
	} catch (std::exception& e){
		std::cout<<e.what()<<std::endl;
		return HMERROR;
	}
}
int g2_custom_rect_grid(const char* algo, void* left, void* bottom, void* right, void* top,
		double* herw, int rinvalid, void** ret, hmcport_callback cb){
	ScaleBase sc_backup;
	try{
		//gather data
		HM2D::EdgeData* left1 = static_cast<HM2D::EdgeData*>(left);
		HM2D::EdgeData* bot1 = static_cast<HM2D::EdgeData*>(bottom);
		HM2D::EdgeData* right1 = static_cast<HM2D::EdgeData*>(right);
		HM2D::EdgeData* top1 = static_cast<HM2D::EdgeData*>(top);
		//scale
		Autoscale::D2 sc(vector<HM2D::EdgeData*>{left1, bot1, right1, top1});
		sc_backup = sc.get_scale(); //copy scale for InvalidReturn treatment
		//assemble contours
		HM2D::EdgeData left2 = HM2D::Contour::Assembler::SimpleContours(*left1)[0];
		HM2D::EdgeData bot2 = HM2D::Contour::Assembler::SimpleContours(*bot1)[0];
		HM2D::EdgeData right2 = HM2D::Contour::Assembler::SimpleContours(*right1)[0];
		HM2D::EdgeData top2 = HM2D::Contour::Assembler::SimpleContours(*top1)[0];
		//deep copy because builder procedures shift nodes
		HM2D::EdgeData left3, bot3, right3, top3;
		HM2D::DeepCopy(left2, left3);
		HM2D::DeepCopy(bot2, bot3);
		HM2D::DeepCopy(right2, right3);
		HM2D::DeepCopy(top2, top3);
		//assemble grid
		HM2D::GridData ret_;
		switch (std::map<std::string, int>{
			{"linear", 1},
			{"inverse_laplace", 2},
			{"direct_laplace", 3},
			{"orthogonal", 4},
			{"linear_tfi", 5},
			{"hermite_tfi", 6} }[algo]){
		case 1: 
			ret_ = HMMap::LinearRectGrid(left3, bot3, right3, top3);
			break;
		case 2: 
			ret_ = HMMap::LaplaceRectGrid.WithCallback(
				cb, left3, bot3, right3, top3, "inverse-laplace");
			break;
		case 3:
			ret_ = HMMap::LaplaceRectGrid.WithCallback(
				cb, left3, bot3, right3, top3, "direct-laplace");
			break;
		case 4:
			ret_ = HMMap::OrthogonalRectGrid.WithCallback(
				cb, left3, bot3, right3, top3);
			break;
		case 5:
			ret_ = HMMap::LinearTFIRectGrid(
				left3, bot3, right3, top3);
			break;
		case 6:
			ret_ = HMMap::CubicTFIRectGrid(
				left3, bot3, right3, top3,
				{herw[0], herw[1], herw[2], herw[3]});
			break;
		default:
			throw std::runtime_error("unknown algorithm");
		}
		sc.unscale(&ret_);
		c2cpp::to_pp(ret_, ret);
		return HMSUCCESS;
	} catch (HMMap::EInvalidGrid &e){
		if (rinvalid){
			std::cout<<"Ignored error: "<<e.what()<<std::endl;
			HM2D::Unscale(e.invalid_grid.vvert, sc_backup);
			c2cpp::to_pp(e.invalid_grid, ret);
			return HMSUCCESS;
		} 
		std::cout<<e.what()<<std::endl;
		return HMERROR;
	} catch (std::runtime_error &e){
		std::cout<<e.what()<<std::endl;
		return HMERROR;
	}
}

int g2_circ4grid(const char* algo, double* p0, double rad, double step, double sqrside, double rcoef, void** ret){
	try{
		double n = 2*M_PI*rad/step;
		int n1 = round(n/8.0);
		std::string stralgo;
		switch (std::map<std::string, int>{
				{"linear", 1},
				{"laplace", 2},
				{"orthogonal_circ", 3},
				{"orthogonal_rect", 4} }[algo]){
			case 1: stralgo="linear"; break;
			case 2: stralgo="laplace"; break;
			case 3: stralgo="orthogonal-circ"; break;
			case 4: stralgo="orthogonal-rect"; break;
			default: throw std::runtime_error("unknown algorithm");
		};
		HM2D::GridData ret_ = HMMap::Circ4Prototype(Point(0, 0), 1.0, 8*n1,
				stralgo, sqrside, rcoef);
		//unscale
		ScaleBase sc(p0[0], p0[1], rad);
		HM2D::Unscale(ret_.vvert, sc);
		c2cpp::to_pp(ret_, ret);
		return HMSUCCESS;
	} catch (std::exception& e){
		std::cout<<e.what()<<std::endl;
		return HMERROR;
	}
}

int g2_stripe_grid(void* obj, int npart, double* part, const char* tipalgo, int* bnd, void** ret,
		hmcport_callback cb){
	try{
		auto cont = static_cast<HM2D::EdgeData*>(obj);
		Autoscale::D2 sc(cont);

		//scaling
		HM2D::EdgeData cbot, cright, cleft, ctop;
		vector<double> spart(part, part+npart);
		sc.scale(spart);

		//main procedure
		Point bl, br, tr, tl;
		HM2D::EdgeData ic = HM2D::Contour::Assembler::Contour1(*cont);
		int ta = 0;
		if (c2cpp::eqstring(tipalgo, "radial")) ta = 1;
		HM2D::GridData ret_ = HMBlay::BuildStripeGrid.WithCallback(cb, ic, spart, ta,
				bl, br, tr, tl);

		//boundary types
		auto gc = HM2D::Contour::Tree::GridBoundary(ret_);
		auto av = HM2D::AllVertices(gc.alledges());
		auto bl2 = av[std::get<0>(HM2D::Finder::ClosestPoint(av, bl))].get();
		auto br2 = av[std::get<0>(HM2D::Finder::ClosestPoint(av, br))].get();
		auto tl2 = av[std::get<0>(HM2D::Finder::ClosestPoint(av, tl))].get();
		auto tr2 = av[std::get<0>(HM2D::Finder::ClosestPoint(av, tr))].get();
		if (bl != br){
			//if open
			auto& cont = gc.nodes[0]->contour;
			cbot = HM2D::Contour::Assembler::Contour1(cont, bl2, br2);
			ctop = HM2D::Contour::Assembler::Contour1(cont, tr2, tl2);
			cleft = HM2D::Contour::Assembler::Contour1(cont, tl2, bl2);
			cright = HM2D::Contour::Assembler::Contour1(cont, br2, tr2);
			for (auto& e: cbot) e->boundary_type = bnd[0];
			for (auto& e: cright) e->boundary_type = bnd[1];
			for (auto& e: ctop) e->boundary_type = bnd[2];
			for (auto& e: cleft) e->boundary_type = bnd[3];
		} else {
			//if closed
			auto& gout = gc.roots()[0]->contour;
			auto& gin = gc.roots()[0]->children[0].lock()->contour;
			cbot = HM2D::Contour::Assembler::Contour1(gin, bl2, br2);
			ctop = HM2D::Contour::Assembler::Contour1(gout, tr2, tl2);
			for (auto& e: cbot) e->boundary_type = bnd[0];
			for (auto& e: ctop) e->boundary_type = bnd[2];
		}

		sc.unscale(&ret_);
		c2cpp::to_pp(ret_, ret);
		return HMSUCCESS;
	} catch (std::exception& e){
		std::cout<<e.what()<<std::endl;
		return HMERROR;
	}
}
int g2_map_grid(void* base_obj, void* target_obj,
		int npoints, double* base_points, double* target_points,
		const char* snap, int bt_from_contour, const char* algo,
		int is_reversed, int rinvalid, void** ret, hmcport_callback cb){
	ScaleBase sc_backup;
	try{
		HM2D::GridData* g = static_cast<HM2D::GridData*>(base_obj);
		HM2D::EdgeData* cont = static_cast<HM2D::EdgeData*>(target_obj);
		//scaling
		Autoscale::D2 gscale(g);
		Autoscale::D2 cscale(cont);
		sc_backup = cscale.get_scale();
		vector<Point> p1 = c2cpp::to_points2(npoints, base_points);
		vector<Point> p2 = c2cpp::to_points2(npoints, target_points);
		gscale.scale(p1);
		cscale.scale(p2);

		HMMap::Options opt;
		opt.btypes_from_contour = bt_from_contour;

		switch (std::map<std::string, int>{
				{"no", 1},
				{"add_vertices", 2},
				{"shift_vertices", 3}}[snap]){
			case 1: opt.snap = "NO"; break;
			case 2: opt.snap = "ADD_VERTICES"; break;
			case 3: opt.snap = "SHIFT_VERTICES"; break;
			default: throw std::runtime_error("unknown snap method");
		}

		switch (std::map<std::string, int>{
				{"direct_laplace", 1},
				{"inverse_laplace", 2}}[algo]){
			case 1: opt.algo = "direct-laplace"; break;
			case 2: opt.algo = "inverse-laplace"; break;
			default: throw std::runtime_error("unknown algorithm");
		}
		
		HM2D::GridData ret_ = HMMap::MapGrid.WithCallback(cb, *g, *cont,
				p1, p2, is_reversed, opt);
		cscale.unscale(&ret_);
		c2cpp::to_pp(ret_, ret);
		return HMSUCCESS;
	} catch (HMMap::EInvalidGrid &e){
		if (rinvalid){
			std::cout<<"Ignored error: "<<e.what()<<std::endl;
			HM2D::Unscale(e.invalid_grid.vvert, sc_backup);
			c2cpp::to_pp(e.invalid_grid, ret);
			return HMSUCCESS;
		}
		std::cout<<e.what()<<std::endl;
		return HMERROR;
	} catch (std::exception& e){
		std::cout<<e.what()<<std::endl;
		return HMERROR;
	}
}

int g2_closest_points(void* obj, int npts, double* pts, const char* proj, double* ret){
	try{
		auto g = static_cast<HM2D::GridData*>(obj);
		return c2_closest_points(&g->vedges, npts, pts, proj, ret);
	} catch (std::exception& e){
		std::cout<<e.what()<<std::endl;
		return HMERROR;
	}
}
int g2_point_at(void* obj, int index, double* ret){
	try{
		auto g = static_cast<HM2D::GridData*>(obj);
		if (index >= g->vvert.size()) throw std::runtime_error("index is out of range");
		ret[0] = g->vvert[index]->x;
		ret[1] = g->vvert[index]->y;
		return HMSUCCESS;
	} catch (std::exception& e){
		std::cout<<e.what()<<std::endl;
		return HMERROR;
	}
}

int g2_tab_btypes(void* obj, int* ret){
	try{
		auto g = static_cast<HM2D::GridData*>(obj);
		return c2_tab_btypes(&g->vedges, ret);
	} catch (std::exception& e){
		std::cout<<e.what()<<std::endl;
		return HMERROR;
	}
}
int g2_tab_vertices(void* obj, double* ret){
	try{
		auto g = static_cast<HM2D::GridData*>(obj);
		for (auto& v: g->vvert){
			*ret++ = v->x;
			*ret++ = v->y;
		}
		return HMSUCCESS;
	} catch (std::exception& e){
		std::cout<<e.what()<<std::endl;
		return HMERROR;
	}
}
int g2_tab_edgevert(void* obj, int* ret){
	try{
		auto g = static_cast<HM2D::GridData*>(obj);
		aa::enumerate_ids_pvec(g->vvert);
		for (auto& e: g->vedges){
			*ret++ = e->pfirst()->id;
			*ret++ = e->plast()->id;
		}
		return HMSUCCESS;
	} catch (std::exception& e){
		std::cout<<e.what()<<std::endl;
		return HMERROR;
	}
}
int g2_tab_cellsizes(void* obj, int* ret){
	try{
		auto g = static_cast<HM2D::GridData*>(obj);
		for (auto& c: g->vcells){
			*ret++ = c->edges.size();
		}
		return HMSUCCESS;
	} catch (std::exception& e){
		std::cout<<e.what()<<std::endl;
		return HMERROR;
	}
}
int g2_tab_cellvert(void* obj, int* nret, int** ret){
	try{
		auto g = static_cast<HM2D::GridData*>(obj);
		vector<int> ret_;
		aa::enumerate_ids_pvec(g->vvert);
		for (auto& c: g->vcells){
			for(auto& p: HM2D::Contour::OrderedPoints1(c->edges)){
				ret_.push_back(p->id);
			}
		}
		*nret = ret_.size();
		*ret = new int[*nret];
		std::copy(ret_.begin(), ret_.end(), *ret);
		return HMSUCCESS;
	} catch (std::exception& e){
		std::cout<<e.what()<<std::endl;
		return HMERROR;
	}
}
int g2_tab_celledge(void* obj, int* nret, int** ret){
	try{
		auto g = static_cast<HM2D::GridData*>(obj);
		vector<int> ret_;
		aa::enumerate_ids_pvec(g->vedges);
		for (auto& c: g->vcells){
			for(auto& e: c->edges){
				ret_.push_back(e->id);
			}
		}
		*nret = ret_.size();
		*ret = new int[*nret];
		std::copy(ret_.begin(), ret_.end(), *ret);
		return HMSUCCESS;
	} catch (std::exception& e){
		std::cout<<e.what()<<std::endl;
		return HMERROR;
	}
}
int g2_tab_centers(void* obj, double* ret){
	try{
		auto g = static_cast<HM2D::GridData*>(obj);
		for (auto& c: g->vcells){
			Point p(0, 0);
			for (auto& v: HM2D::AllVertices(c->edges)){
				p+=*v;
			}
			p/=c->edges.size();
			*ret++ = p.x;
			*ret++ = p.y;
		}
		return HMSUCCESS;
	} catch (std::exception& e){
		std::cout<<e.what()<<std::endl;
		return HMERROR;
	}
}
int g2_tab_bedges(void* obj, int* nret, int** ret){
	try{
		auto g = static_cast<HM2D::GridData*>(obj);
		auto cont = HM2D::ECol::Assembler::GridBoundary(*g);
		*nret = cont.size();
		*ret = new int[*nret];
		aa::enumerate_ids_pvec(g->vedges);
		for (int i=0; i<*nret; ++i) (*ret)[i] = cont[i]->id;
		return HMSUCCESS;
	} catch (std::exception& e){
		std::cout<<e.what()<<std::endl;
		return HMERROR;
	}
}
int g2_simplify_bnd(void* obj, double angle, void** ret){
	try{
		if (angle<-geps || angle>180+geps) throw std::runtime_error("invalid angle");
		HM2D::GridData* g = static_cast<HM2D::GridData*>(obj);
		Autoscale::D2 sc(g);
		HM2D::GridData ret_;
		HM2D::DeepCopy(*g, ret_);
		HM2D::Grid::Algos::SimplifyBoundary(ret_, angle);
		sc.unscale(&ret_);
		c2cpp::to_pp(ret_, ret);
		return HMSUCCESS;
	} catch (std::exception& e){
		std::cout<<e.what()<<std::endl;
		return HMERROR;
	}
}
int g2_convex_cells(void* obj, double angle, void** ret){
	try{
		if (angle<-geps || angle>180+geps) throw std::runtime_error("invalid angle");
		HM2D::GridData* g = static_cast<HM2D::GridData*>(obj);
		Autoscale::D2 sc(g);
		HM2D::GridData ret_;
		HM2D::DeepCopy(*g, ret_);
		HM2D::Grid::Algos::NoConcaveCells(ret_, angle);
		sc.unscale(&ret_);
		c2cpp::to_pp(ret_, ret);
		return HMSUCCESS;
	} catch (std::exception& e){
		std::cout<<e.what()<<std::endl;
		return HMERROR;
	}
}
int g2_exclude_cont(void* obj, void* cont, int isinner, void** ret, hmcport_callback cb){
	try{
		HM2D::GridData* g0 = static_cast<HM2D::GridData*>(obj);
		HM2D::EdgeData* e0 = static_cast<HM2D::EdgeData*>(cont);

		//scaling
		//using non-unity scaling to minimize risk of almost doubling points
		//when using uniform rectangles
		double unity = 1.0 + sqrt(2.0)/100.0 + sqrt(3.0)/1000.0;
		Autoscale::D2 sc(g0, unity);
		sc.add_data(e0);
		//building contour tree
		auto tree = HM2D::Contour::Tree::Assemble(*e0);
		//processing
		HM2D::GridData ret_ = HM2D::Grid::Algos::SubstractArea.WithCallback(
			cb, *g0, tree, isinner);
		//finalizing
		sc.unscale(&ret_);
		c2cpp::to_pp(ret_, ret);
		return HMSUCCESS;
	} catch (std::exception& e){
		std::cout<<e.what()<<std::endl;
		return HMERROR;
	}
}
int g2_unite_grids(void* obj1, void* obj2, double buf, int fixbnd,
		int emptyholes, double angle0, const char* filler,
		void** ret, hmcport_callback cb){
	try{
		HM2D::GridData* g0 = static_cast<HM2D::GridData*>(obj1);
		HM2D::GridData* g1 = static_cast<HM2D::GridData*>(obj2);
		//using non-unity scaling to minimize risk of almost doubling points
		//when using uniform rectangles
		double unity = 1.0 + sqrt(2.0)/100.0 + sqrt(3.0)/1000.0;
		Autoscale::D2 sc(g0, unity);
		sc.add_data(g1);
		sc.scale(buf);

		HM2D::Grid::Algos::OptUnite opt;
		opt.buffer_size = buf;
		opt.preserve_bp = (bool)fixbnd;
		opt.empty_holes = (bool)emptyholes;
		opt.angle0 = angle0;
		opt.filler = (c2cpp::eqstring(filler, "3")) ? 0 : 1;
		auto ret_ = HM2D::Grid::Algos::UniteGrids.WithCallback(cb, *g0, *g1, opt);
		sc.unscale(&ret_);
		c2cpp::to_pp(ret_, ret);
		return HMSUCCESS;
	} catch (std::exception& e){
		std::cout<<e.what()<<std::endl;
		return HMERROR;
	}
}

namespace{
HM2D::Export::BNamesFun construct_bnames(const BoundaryNamesStruct& bnames){
	std::map<int, std::string> bnames_map;
	for (int i=0; i<bnames.n; ++i){
		bnames_map[bnames.index[i]] = std::string(bnames.name[i]);
	}
	auto fnames = [bnames_map](int i)->std::string{
			auto fnd = bnames_map.find(i);
			if (fnd == bnames_map.end()) return "boundary" + std::to_string(i);
			else return fnd->second;
		};
	return fnames;
}
}

int g2_to_msh(void* obj, const char* fname, BoundaryNamesStruct btypes, int n_per_data, int* per_data){
	try{
		auto g = static_cast<HM2D::GridData*>(obj);
		auto fnames = construct_bnames(btypes);
		HM2D::Export::PeriodicData pd;
		for (int i=0; i<n_per_data; ++i){
			int b1 = *per_data++;
			int b2 = *per_data++;
			bool is_rev = (bool)(*per_data++);
			pd.add_data(b1, b2, is_rev);
		}
		HM2D::Export::GridMSH(*g, fname, fnames, pd);
		return HMSUCCESS;
	} catch (std::exception& e){
		std::cout<<e.what()<<std::endl;
		return HMERROR;
	}
}
int g2_to_tecplot(void* obj, const char* fname, BoundaryNamesStruct btypes){
	try{
		auto g = static_cast<HM2D::GridData*>(obj);
		auto fnames = construct_bnames(btypes);
		HM2D::Export::GridTecplot(*g, fname, fnames);
		return HMSUCCESS;
	} catch (std::exception& e){
		std::cout<<e.what()<<std::endl;
		return HMERROR;
	}
}
int g2_to_hm(void* doc, void* node, void* obj, const char* name, const char* fmt, int naf, const char** af){
	try{
		HMXML::ReaderA* wr = static_cast<HMXML::ReaderA*>(doc);
		HMXML::Reader* sn = static_cast<HMXML::Reader*>(node);
		HM2D::GridData* g = static_cast<HM2D::GridData*>(obj);
		HM2D::Export::GridWriter gw(*g, wr, sn, name, fmt);

		//additional fields
		for (int i=0; i<naf; ++i){
			switch (std::map<std::string, int>{
					{"cell_edges", 1},
					{"cell-edges", 2},
					{"cell_vertices", 3},
					{"cell-vertices", 4}}[af[i]]){
			case 1: case 2:
				gw.AddCellEdgeConnectivity();
				break;
			case 3: case 4:
				gw.AddCellVertexConnectivity();
				break;
			default:
				throw std::runtime_error("unknown field "+std::string(af[i]));
			}
		}
		return HMSUCCESS;
	} catch (std::exception& e){
		std::cout<<e.what()<<std::endl;
		return HMERROR;
	}
}
int g2_boundary_layer(int nopt, BoundaryLayerGridOption* opt, void** ret, hmcport_callback cb){
	try{
		vector<HMBlay::Input> vinp(nopt);
		for (int i=0; i<nopt; ++i){
			auto& opt_ = opt[i];
			auto& inp = vinp[i];
			inp.edges = static_cast<HM2D::EdgeData*>(opt_.cont);
			inp.direction = HMBlay::DirectionFromString(opt_.tp);
			inp.bnd_step_method = HMBlay::MethFromString(opt_.mesh_cont);
			inp.bnd_step = opt_.mesh_cont_step;
			inp.acute_angle = opt_.angle_range[0];
			inp.right_angle = opt_.angle_range[1];
			inp.straight_angle = opt_.angle_range[2];
			inp.reentrant_angle = opt_.angle_range[3];
			inp.start = Point(opt_.start[0], opt_.start[1]);
			inp.end = Point(opt_.end[0], opt_.end[1]);
			inp.partition = vector<double>(opt_.part, opt_.part + opt_.Npart);
			inp.force_conformal = (opt_.force_conformal == 1);

			if (inp.bnd_step_method == HMBlay::BndStepMethod::INCREMENTAL){
				if (inp.start == inp.end){
					throw std::runtime_error("Can not use incremental stepping "
							         "without divergent start/end points");
				}
				inp.bnd_step_basis.push_back(std::make_pair(Point(inp.start), opt_.step_start));
				inp.bnd_step_basis.push_back(std::make_pair(Point(inp.end), opt_.step_end));
			}
		}
		*ret = new HM2D::GridData(HMBlay::BuildBLayerGrid(vinp));
		return HMSUCCESS;
	} catch (std::exception& e){
		std::cout<<e.what()<<std::endl;
		return HMERROR;
	}
}
