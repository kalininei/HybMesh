#include "pebi.hpp"
#include "cont_assembler.hpp"
#include "contabs2d.hpp"
#include "buildgrid.hpp"
#include "healgrid.hpp"
#include "trigrid.hpp"
using namespace HM2D;
using namespace HM2D::Grid;

namespace{
struct PebiBndPoint{
	int ind;
	shared_ptr<Vertex> prev, next;
	shared_ptr<Vertex> hprev, hnext;
};

std::vector<PebiBndPoint> assemble_bnd(const GridData& g){
	std::vector<PebiBndPoint> ret;
	auto cont = Contour::Assembler::GridBoundary(g);
	auto& gpoints = g.vvert;
	aa::enumerate_ids_pvec(gpoints);

	for (auto& c: cont){
		auto op = HM2D::Contour::OrderedPoints(c);
		for (int i=0; i<op.size()-1; ++i){
			int im = (i==0)?op.size()-2:i-1;
			int cur = op[i]->id;
			int prev = op[im]->id;
			int next = op[i+1]->id;
			ret.push_back(PebiBndPoint {cur, gpoints[prev], gpoints[next]});
			ret.back().hnext.reset(new Vertex((*gpoints[cur]+*gpoints[next])/2.0));
			if (i!=0) ret.back().hprev = ret.end()[-2].hnext;
			if (i==op.size()-2) ret[ret.size()-1-i].hprev = ret.back().hnext;
		}
	}
	return ret;
}

std::vector<std::vector<int>> ordered_points_cells(const GridData& g, std::vector<PebiBndPoint>& bnd){
	std::vector<std::vector<int>> ret(g.vvert.size());
	auto& edges=g.vedges;
	auto ve = Connectivity::VertexEdge(edges, g.vvert);
	std::vector<std::list<int>> point_edges(g.vvert.size());
	int iv=0;
	for (auto& it: ve){
		auto& at = point_edges[iv++];
		for (int ie: it.eind) at.push_back(ie);
	}

	g.enumerate_all();
	for (int i=0; i<g.vvert.size(); ++i){
		auto& pe = point_edges[i];
		Edge* estart=0;
		//looking for boundary edge
		for (auto it = pe.begin(); it!=pe.end(); ++it){
			if (edges[*it]->is_boundary()){
				auto ed = edges[*it].get();
				if ((ed->first()->id == i && ed->has_left_cell()) ||
				    (ed->last()->id == i && ed->has_right_cell())){
					estart = ed;
					pe.erase(it);
					break;
				}
			}
		}
		//if no boundary edges start from first
		if (estart == 0) {estart = edges[pe.front()].get(); pe.pop_front();}
		Edge* ecur = estart;
		while (1){
			auto wcur_cell = (ecur->first()->id == i) ? ecur->left : ecur->right;
			if (wcur_cell.expired()) break;
			int cur_cell = wcur_cell.lock()->id;
			ret[i].push_back(cur_cell);
			size_t pesize = pe.size();
			if (pesize == 0) break;
			//find cur_edge as edge of cur_cell which is not ecur but contains i;
			for (auto it=pe.begin(); it!=pe.end(); ++it){
				auto eit = edges[*it];
				//find edge which has a link to cur_cell
				if ((eit->has_left_cell() && eit->left.lock()->id == cur_cell) ||
				    (eit->has_right_cell() && eit->right.lock()->id == cur_cell)){
					ecur = edges[*it].get();
					pe.erase(it);
					break;
				}
			}
			assert(pe.size() == pesize-1);
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
	if (det==0) throw std::runtime_error("pebi builder error");
	double B[] = {A[3]/det, -A[1]/det, -A[2]/det, A[0]/det};
	return Point(-B[0]*l1[2]-B[1]*l2[2], -B[2]*l1[2]-B[3]*l2[2]);
}

VertexData build_pebi_pts(const GridData& g){
	VertexData ret;
	vector<vector<int>> cc = Connectivity::CellCell(g.vcells);
	auto within = [&](Point p, const std::array<Point, 3>& ptri)->bool{
		auto &x1 = ptri[0].x, &x2 = ptri[1].x, &x3 = ptri[2].x;
		auto &y1 = ptri[0].y, &y2 = ptri[1].y, &y3 = ptri[2].y;
		double j11 = x2 - x1, j21 = x3 - x1;
		double j12 = y2 - y1, j22 = y3 - y1;
		double modj = (j22*j11 - j21*j12);
		Point ksieta;
		ksieta.x = ( j22*(p.x - x1) - j21*(p.y - y1))/modj;
		ksieta.y = (-j12*(p.x - x1) + j11*(p.y - y1))/modj;
		return (ksieta.x >= 0 && ksieta.x <= 1.0 && ksieta.y>=0 && ksieta.y<=1.0-ksieta.x);
	};
	vector<std::array<Point, 3>> cellvert(g.vcells.size());
	for (int i=0; i<g.vcells.size(); ++i){
		const Cell* c=g.vcells[i].get();
		auto op = Contour::OrderedPoints(c->edges);
		cellvert[i][0] = *op[0];
		cellvert[i][1] = *op[1];
		cellvert[i][2] = *op[2];
	}
	for (int i=0; i<g.vcells.size(); ++i){
		auto& cv = cellvert[i];
		Point p = calc_pebi(cv[0], cv[1], cv[2]);

		//bad point is point which lies far away outside parent triangle.
		//Usually even out of area. Hence we snap it to parant triangle edge.
		bool bad_point=true;
		//within itself
		if (within(p, cellvert[i])) bad_point = false;
		//within neighbours
		if (bad_point) for (int j=0; j<cc[i].size(); ++j){
			if (within(p, cellvert[cc[i][j]])) { bad_point = false; break; }
		}
		//within neighbours of neighbours
		if (bad_point) for (int j=0; j<cc[i].size(); ++j){
			int in = cc[i][j];
			for (int k=0; k<cc[in].size(); ++k)
				if (within(p, cellvert[cc[in][k]])) { bad_point = false; break; }
			if (!bad_point) break;
		}

		if (bad_point){
			double d1 = Point::meas_line(p, cv[0], cv[1]);
			double d2 = Point::meas_line(p, cv[0], cv[2]);
			double d3 = Point::meas_line(p, cv[1], cv[2]);
			if (d1 < d2 && d1 < d3) p = (cv[0] + cv[1])/2.0;
			else if (d2<d1 && d2<d3) p = (cv[0]+cv[2])/2.0;
			else p = (cv[1]+cv[2])/2.0;
		}
		aa::add_shared(ret, Vertex(p));
	}
	return ret;
}


}

GridData Constructor::TriToPebi(const GridData& g){
	//calculate pebi points
	ShpVector<Vertex> pebi_points = build_pebi_pts(g);

	//boundary information
	std::vector<PebiBndPoint> bnd_left_right = assemble_bnd(g);

	//ordered points->cells table
	std::vector<std::vector<int>> points_cells = ordered_points_cells(g, bnd_left_right);

	//asseble points
	VertexData rp = pebi_points;
	for (int i=0; i<bnd_left_right.size(); ++i){
		rp.push_back(bnd_left_right[i].hnext);
	}
	aa::enumerate_ids_pvec(rp);
	//assemble cells
	vector<vector<int>> cellvert;
	//internal cells
	for (int i=0; i<g.vvert.size(); ++i){
		cellvert.emplace_back();
		auto& cc = cellvert.back();
		for (auto& x: points_cells[i]) cc.push_back(pebi_points[x]->id);
	}
	//boundary segments
	for (auto& b: bnd_left_right){
		auto& cc = cellvert[b.ind];
		cc.insert(cc.begin(), b.hnext->id);
		cc.push_back(b.hprev->id);
		double ksi;
		if (!isOnSection(*g.vvert[b.ind], *b.hnext, *b.hprev, ksi)){
			auto np = std::make_shared<Vertex>(*g.vvert[b.ind]);
			np->id = rp.size();
			rp.push_back(np);
			cc.push_back(np->id);
		}
	}
	auto ret = Grid::Constructor::FromTab(rp, cellvert);

	//post processing: collapse reversed edges which provoke сell self intersection if possible
	VertexData bpts = AllVertices(ECol::Assembler::GridBoundary(ret));
	std::sort(bpts.begin(), bpts.end());
	for (int i=0; i<cellvert.size(); ++i){
		const Cell* c = ret.vcells[i].get();
		for (int k=0; k<cellvert[i].size(); ++k){
			int kn = (k == cellvert[i].size()-1) ? 0 : k+1;
			int km = (k == 0) ? cellvert[i].size() - 1 : k-1;
			int knn = (kn == cellvert[i].size()-1) ? 0 : kn+1;
			auto p0 = rp[cellvert[i][k]];
			auto p1 = rp[cellvert[i][kn]];
			auto prev = rp[cellvert[i][km]];
			auto next = rp[cellvert[i][knn]];
			int w = LinePointWhereIs(*rp[i], *p0, *p1);
			if (w != 2) continue;
			if (std::binary_search(bpts.begin(), bpts.end(), p0)) continue;
			if (std::binary_search(bpts.begin(), bpts.end(), p1)) continue;
			if (*prev == *p0 || *next == *p1) continue;
			double ksieta[2];
			if (SectCross(*prev, *p0, *p1, *next, ksieta)){
				p0->set((*p0 + *p1)/2.0);
				p1->set(*p0);
			}
		}
	}

	//remove short edges
	Grid::Algos::RemoveShortEdges(ret, 0.1);
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
void build_regular_hex(Point cnt, double rad, VertexData& pts, vector<vector<int>>& cls){
	auto ap = build_regular_hex(cnt, rad, true);
	cls.emplace_back();
	auto& cc = cls.back();
	for (int i=0; i<6; ++i){
		cc.push_back(pts.size());
		pts.push_back(std::make_shared<Vertex>(ap[i]));
	}
}
}

GridData Constructor::RegularHexagonal(Point cnt, double area_rad, double cr, bool strict_area){
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
	ShpVector<Vertex> pts;
	vector<vector<int>> cls;
	for (auto& p: cnts) build_regular_hex(p, cr, pts, cls);
	GridData ret = Grid::Constructor::FromTab(pts, cls);
	Grid::Algos::Heal(ret);
	return ret;
}

GridData Constructor::RegularHexagonal(Point cnt1, Point cnt2, double cr, bool strict_area){
	double hy = sqrt(3)/2*cr;
	ShpVector<Vertex> pts;
	vector<vector<int>> cls;

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

	GridData ret = Grid::Constructor::FromTab(pts, cls);
	Grid::Algos::Heal(ret);
	if (strict_area){
		double xcoef = (cnt2.x - cnt1.x)/(pmax.x - cnt1.x);
		double ycoef = (cnt2.y - cnt1.y)/(pmax.y - cnt1.y);
		for (auto& v: ret.vvert){
			v->x = (v->x - cnt1.x)*xcoef + cnt1.x;
			v->y = (v->y - cnt1.y)*ycoef + cnt1.y;
		}
	}
	return ret;
}

HMCallback::FunctionWithCallback<Mesher::TUnstructuredPebi> Mesher::UnstructuredPebi;

GridData Mesher::TUnstructuredPebi::_run(const Contour::Tree& source, const std::map<Point, double>& embedded){
	//triangulation
	auto cb = callback->subrange(100, 100);
	GridData tri = UnstructuredTriangle.UseCallback(cb, source, embedded);
	//build pebi
	callback->step_after(10, "assembling polygons");
	return Constructor::TriToPebi(tri);
}

GridData Mesher::TUnstructuredPebi::_run(const Contour::Tree& source){
	return _run(source, std::map<Point, double>());
}
