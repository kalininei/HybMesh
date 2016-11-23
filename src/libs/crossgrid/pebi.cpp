#include "pebi.h"
#include "trigrid.h"
#include "procgrid.h"

namespace{
struct PebiBndPoint{
	int ind;
	shared_ptr<GridPoint> prev, next;
	shared_ptr<GridPoint> hprev, hnext;
};

std::vector<PebiBndPoint> assemble_bnd(const TriGrid& g){
	std::vector<PebiBndPoint> ret;
	auto cont = GGeom::Info::Contour(g);
	auto gpoints = GGeom::Info::SharePoints(g);
	for (auto c: cont.nodes){
		auto op = c->ordered_points();
		for (int i=0; i<op.size()-1; ++i){
			int im = (i==0)?op.size()-2:i-1;
			int cur = static_cast<const GridPoint*>(op[i])->get_ind();
			int prev = static_cast<const GridPoint*>(op[im])->get_ind();
			int next = static_cast<const GridPoint*>(op[i+1])->get_ind();
			ret.push_back(PebiBndPoint {cur, gpoints[prev], gpoints[next]});
			ret.back().hnext.reset(new GridPoint((*gpoints[cur]+*gpoints[next])/2.0));
			if (i!=0) ret.back().hprev = ret.end()[-2].hnext;
			if (i==op.size()-2) ret[ret.size()-1-i].hprev = ret.back().hnext;
		}
	}
	return ret;
}

std::vector<std::vector<int>> ordered_points_cells(const TriGrid& g, std::vector<PebiBndPoint>& bnd){
	std::vector<std::vector<int>> ret(g.n_points());
	auto edges=g.get_edges();
	std::vector<std::list<Edge>> point_edges(g.n_points());
	for (auto& e: edges){
		point_edges[e.p1].push_back(e);
		point_edges[e.p2].push_back(e);
	}
	for (int i=0; i<g.n_points(); ++i){
		auto& pe = point_edges[i];
		Edge estart(-1, -1);
		//looking for boundary edge
		for (auto it = pe.begin(); it!=pe.end(); ++it){
			int right_cell = (it->p1 == i)?it->cell_right:it->cell_left;
			if (right_cell < 0){ estart = *it; pe.erase(it); break; }
		}
		//if no boundary edges start from first
		if (estart.p1 == -1) {estart = pe.front(); pe.pop_front();}
		Edge ecur = estart;
		while (1){
			int cur_cell = (ecur.p1 == i)?ecur.cell_left:ecur.cell_right;
			if (cur_cell<0) break;
			ret[i].push_back(cur_cell);
			//find cur_edge as edge of cur_cell which is not ecur but contains i;
			if (pe.size() == 0) break;
			else for (auto it=pe.begin(); it!=pe.end(); ++it){
				//find edge which has a link to cur_cell
				if (it->cell_left == cur_cell || it->cell_right == cur_cell){
					ecur = *it; pe.erase(it); break;
				}
			}
		}
	}

	return ret;
}

Point calc_pebi(Point p1, Point p2, Point p3){
	auto line2p = [](double x1, double y1, double x2, double y2){
		double xc = (x1+x2)/2.0, yc = (y1+y2)/2.0;
		return std::array<double, 3> { x2-x1, y2-y1, -(x2-x1)*xc-(y2-y1)*yc };
	};
	auto l1 = line2p(p1.x, p1.y, p2.x, p2.y);
	auto l2 = line2p(p2.x, p2.y, p3.x, p3.y);
	double A[] = {l1[0], l1[1], l2[0], l2[1]};
	double det = A[0]*A[3]-A[1]*A[2];
	if (ISZERO(det)) throw std::runtime_error("pebi error");
	double B[] = {A[3]/det, -A[1]/det, -A[2]/det, A[0]/det};
	return Point(-B[0]*l1[2]-B[1]*l2[2], -B[2]*l1[2]-B[3]*l2[2]);
}

ShpVector<GridPoint> build_pebi_pts(const TriGrid& g){
	ShpVector<GridPoint> ret;
	vector<vector<int>> cc = g.cell_cell();
	auto within = [&](Point p, int cn)->bool{
		const Cell* c=g.get_cell(cn);
		auto &x1 = c->points[0]->x, &x2 = c->points[1]->x, &x3 = c->points[2]->x;
		auto &y1 = c->points[0]->y, &y2 = c->points[1]->y, &y3 = c->points[2]->y;
		double j11 = x2 - x1, j21 = x3 - x1;
		double j12 = y2 - y1, j22 = y3 - y1;
		double modj = (j22*j11 - j21*j12);
		Point ksieta;
		ksieta.x = ( j22*(p.x - x1) - j21*(p.y - y1))/modj;
		ksieta.y = (-j12*(p.x - x1) + j11*(p.y - y1))/modj;
		return (ksieta.x >= 0 && ksieta.x <= 1.0 && ksieta.y>=0 && ksieta.y<=1.0-ksieta.x);
	};
	for (int i=0; i<g.n_cells(); ++i){
		Point p1 = *g.get_cell(i)->points[0];
		Point p2 = *g.get_cell(i)->points[1];
		Point p3 = *g.get_cell(i)->points[2];
		Point p = calc_pebi(p1, p2, p3);

		//bad point is point which lies far away outside parent triangle.
		//Usually even out of area. Hence we snap it to parant triangle edge.
		bool bad_point=true;
		//within itself
		if (within(p, i)) bad_point = false;
		//within neighbours
		if (bad_point) for (int j=0; j<cc[i].size(); ++j){
			if (within(p, cc[i][j])) { bad_point = false; break; }
		}
		//within neighbours of neighbours
		if (bad_point) for (int j=0; j<cc[i].size(); ++j){
			int in = cc[i][j];
			for (int k=0; k<cc[in].size(); ++k)
				if (within(p, cc[in][k])) { bad_point = false; break; }
			if (!bad_point) break;
		}

		if (bad_point){
			double d1 = Point::meas_line(p, p1, p2);
			double d2 = Point::meas_line(p, p1, p3);
			double d3 = Point::meas_line(p, p2, p3);
			if (d1 < d2 && d1 < d3) p = (p1 + p2)/2.0;
			else if (d2<d1 && d2<d3) p = (p1+p3)/2.0;
			else p = (p2+p3)/2.0;
		}
		aa::add_shared(ret, GridPoint(p));
	}
	return ret;
}


}

GridGeom TriToPebi(const TriGrid& g){
	//calculate pebi points
	ShpVector<GridPoint> pebi_points = build_pebi_pts(g);

	//boundary information
	std::vector<PebiBndPoint> bnd_left_right = assemble_bnd(g);

	//ordered points->cells table
	std::vector<std::vector<int>> points_cells = ordered_points_cells(g, bnd_left_right);

	//asseble points
	ShpVector<GridPoint> rp(pebi_points.begin(), pebi_points.end());
	for (int i=0; i<bnd_left_right.size(); ++i){
		rp.push_back(bnd_left_right[i].hnext);
	}
	//assemble cells
	ShpVector<Cell> rc;
	//internal cells
	for (int i=0; i<g.n_points(); ++i){
		Cell* c=aa::add_shared(rc, Cell());
		for (auto& x: points_cells[i]) c->points.push_back(pebi_points[x].get());
	}
	//boundary segments
	for (auto& b: bnd_left_right){
		Cell* c = rc[b.ind].get();
		c->points.insert(c->points.begin(), b.hnext.get());
		c->points.push_back(b.hprev.get());
		double ksi;
		if (!isOnSection(*g.get_point(b.ind), *b.hnext, *b.hprev, ksi)){
			c->points.push_back(aa::add_shared(rp, *g.get_point(b.ind)));
		}
	}
	
	auto ret = GGeom::Constructor::FromData(rp, rc);
	//post processing: collapse reversed edges which provoke сell self intersection if possible
	std::set<const GridPoint*> bpts = g.get_bnd_points();
	for (int i=0; i<rc.size(); ++i){
		const Cell* c = rc[i].get();
		for (int k=0; k<c->dim(); ++k){
			GridPoint* p0 = const_cast<GridPoint*>(c->get_point(k));
			GridPoint* p1 = const_cast<GridPoint*>(c->get_point(k+1));
			int w = LinePointWhereIs(*g.get_point(i), *p0, *p1);
			if (w != 2) continue;
			if (bpts.find(p0) != bpts.end() || bpts.find(p1) != bpts.end()) continue;
			GridPoint* prev = const_cast<GridPoint*>(c->get_point(c->dim()+k-1));
			GridPoint* next = const_cast<GridPoint*>(c->get_point(k+2));
			if (*prev == *p0 || *next == *p1) continue;
			double ksieta[2];
			if (SectCross(*prev, *p0, *p1, *next, ksieta)){
				p0->set((*p0 + *p1)/2.0);
				p1->set(*p0);
			}
		}
	}
	//remove short edges
	GGeom::Repair::RemoveShortEdges(ret, 0.1);
	return ret;
}

namespace{

std::array<Point, 6> build_regular_hex(Point c, double rad, bool ox){
	static const double y = sqrt(3)/2.0;
	static const std::array<Point, 6> retox {
		Point(1, 0), Point(0.5, y),
		Point(-0.5, y), Point(-1, 0),
		Point(-0.5, -y), Point(0.5, -y)};
	static const std::array<Point, 6> retoy {
		Point(y, 0.5), Point(0, 1),
		Point(-y, 0.5), Point(-y, -0.5),
		Point(0, -1), Point(y, -0.5)};

	std::array<Point, 6> ret(ox?retox:retoy);
	for (auto& p: ret) p = p*rad + c;
	return ret;
}
void build_regular_hex(Point cnt, double rad, ShpVector<GridPoint>& pts, ShpVector<Cell>& cls){
	auto ap = build_regular_hex(cnt, rad, true);
	Cell* c = aa::add_shared(cls, Cell());
	for (int i=0; i<6; ++i){
		c->points.push_back(aa::add_shared(pts, GridPoint(ap[i])));
	}
}
}

GridGeom RegularHexagonal(Point cnt, double area_rad, double cr, bool strict_area){
	double hy = sqrt(3)/2*cr;
	int nmax;
	if (!strict_area){
		double R = 2*hy;
		nmax = 1;
		while (R < area_rad){ R += 2*hy; ++nmax; }
	} else {
		nmax = (int)std::round(area_rad/2/hy)-1;
		if (nmax < 0) nmax = 0;
		hy = area_rad/2/(nmax+1);
		cr = hy/sqrt(3)*2;
		++nmax;
	}
	//calculate hexagon centers
	vector<Point> cnts {cnt};
	for (int n=1; n<=nmax; ++n){
		double R = 2*hy*n;
		auto hx = build_regular_hex(cnt, R, false);
		for (int k=0; k<6; ++k){
			Point& pprev = hx[k];
			Point& pnext = hx[(k+1)%6];
			for (int i=0; i<n; ++i){
				cnts.push_back(Point::Weigh(pprev, pnext, (double)i/n));
			}
		}
	}

	//assemble grid
	ShpVector<GridPoint> pts;
	ShpVector<Cell> cls;
	for (auto& p: cnts) build_regular_hex(p, cr, pts, cls);
	GridGeom ret = GGeom::Constructor::FromData(pts, cls);
	GGeom::Repair::Heal(ret);
	return ret;
}
GridGeom RegularHexagonal(Point cnt1, Point cnt2, double cr, bool strict_area){
	double hy = sqrt(3)/2*cr;
	ShpVector<GridPoint> pts;
	ShpVector<Cell> cls;

	Point curp = cnt1, pmax = cnt1;
	int nx, ny;
	if (!strict_area){
		nx = (int)std::ceil( (cnt2.x - cnt1.x)/(1.5*cr) );
		ny = (int)std::ceil( (cnt2.y - cnt1.y)/(hy) );
	} else {
		nx = (int)std::round( (cnt2.x - cnt1.x)/(1.5*cr) );
		ny = (int)std::round( (cnt2.y - cnt1.y)/(hy) );
		nx = std::max(1, nx);
		ny = std::max(1, ny);
	}

	int mx = 1;
	for (int j=0; j<=ny; ++j){
		pmax.y = curp.y;
		for (int i=(mx==1)?0:1; i<=nx; i+=2){
			if (curp.x > pmax.x) pmax.x = curp.x;
			build_regular_hex(curp, cr, pts, cls);
			curp.x += 3*cr;
		}
		mx *= -1;
		curp.x = cnt1.x+0.75*cr*(1-mx);
		curp.y = curp.y + hy;
	}

	GridGeom ret = GGeom::Constructor::FromData(pts, cls);
	GGeom::Repair::Heal(ret);
	if (strict_area){
		double xcoef = (cnt2.x - cnt1.x)/(pmax.x - cnt1.x);
		double ycoef = (cnt2.y - cnt1.y)/(pmax.y - cnt1.y);
		GGeom::Modify::PointModify(ret, [&](GridPoint* p){
					p->x = (p->x - cnt1.x)*xcoef + cnt1.x;
					p->y = (p->y - cnt1.y)*ycoef + cnt1.y;
				});
	}
	return ret;
}
