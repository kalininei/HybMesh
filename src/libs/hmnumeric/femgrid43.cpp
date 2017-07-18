#include "femgrid43.hpp"
#include "modcont.hpp"
#include "coarsencont.hpp"
#include "assemble2d.hpp"
#include "contabs2d.hpp"
#include "modgrid.hpp"
#include "trigrid.hpp"
#include "finder2d.hpp"
#include "buildcont.hpp"
#include "treverter2d.hpp"

using namespace HMFem;

void Grid43::AddSegments(HM2D::GridData& grid, const vector<vector<Point>>& pts){
	if (pts.size()==0) return;

	//find grid nodes
	HM2D::Finder::VertexMatch vfnd(grid.vvert);
	aa::enumerate_ids_pvec(grid.vvert);
	vector<std::tuple<HM2D::Vertex*,
	                  HM2D::Vertex*,
	                  vector<Point>>> ppsplit;
	for (int i=0; i<pts.size(); ++i){
		auto& pe = pts[i];
		shared_ptr<HM2D::Vertex> fnd1 = vfnd.find(pe[0]);
		shared_ptr<HM2D::Vertex> fnd2 = vfnd.find(pe.back());
		if (fnd1 == nullptr || fnd2 == nullptr) continue;
		ppsplit.emplace_back(fnd1.get(), fnd2.get(),
			vector<Point>(pts[i].begin()+1, pts[i].end()-1));
	}

	//find edges
	vector<std::pair<HM2D::Edge*,
	                 vector<Point>>> esplit;
	HM2D::Finder::EdgeFinder efnd(grid.vedges);
	for (auto& pp: ppsplit){
		auto fres = efnd.find(std::get<0>(pp), std::get<1>(pp));
		assert(std::get<0>(fres));
		esplit.emplace_back(std::get<0>(fres).get(), std::get<2>(pp));
		if (!std::get<1>(fres)) {
			auto& pv = std::get<1>(esplit.back());
			std::reverse(pv.begin(), pv.end());
		}
	}
	//split edges
	aa::enumerate_ids_pvec(grid.vedges);
	for (auto& es: esplit){
		int eind = std::get<0>(es)->id;
		bool r = HM2D::Grid::Algos::SplitEdge(grid, eind, std::get<1>(es));
		if (!r) throw std::runtime_error("refinement of coarse grid failed");
	}
	//grid back to triangle
	HM2D::Grid::Algos::CutCellDims(grid, 3);

	//all cells are triangle
	assert([grid](){
		for (auto c: grid.vcells){
			if (c->edges.size() != 3) return false;
		}
		return true;
	}());
}

Grid43::Approximator::Approximator(const HM2D::GridData* g, int n): grid(g){
	auto bbox = HM2D::BBox(g->vcells);
	double L = bbox.maxlen()/n;
	cfinder.reset(new BoundingBoxFinder(bbox, L));
	int nc = g->vcells.size();
	icellvert.resize(nc);
	cellvert.resize(nc);
	is3.resize(nc, true);
	for (int ic=0; ic<nc; ++ic){
		auto op = HM2D::Contour::OrderedPoints1(g->vcells[ic]->edges);
		auto bb = HM2D::BBox(op);
		cfinder->addentry(bb);
		if (op.size() < 3 || op.size() > 4) throw std::runtime_error(
			"invalid grid was passed to approximator");
		cellvert[ic][0] = op[0].get();
		cellvert[ic][1] = op[1].get();
		cellvert[ic][2] = op[2].get();
		if (op.size() == 4){
			cellvert[ic][3] = op[3].get();
			is3[ic] = false;
		}
	}
	aa::enumerate_ids_pvec(g->vvert);
	for (int ic=0; ic<nc; ++ic){
		icellvert[ic][0] = cellvert[ic][0]->id;
		icellvert[ic][1] = cellvert[ic][1]->id;
		icellvert[ic][2] = cellvert[ic][2]->id;
		if (!is3[ic]) icellvert[ic][3] = cellvert[ic][3]->id;
	}
	for (auto& e: HM2D::ECol::Assembler::GridBoundary(*g)){
		e->pfirst()->id = -1;
		e->plast()->id = -1;
	}
	isbnd.resize(grid->vvert.size(), false);
	for (int i=0; i<isbnd.size(); ++i){
		isbnd[i] = grid->vvert[i]->id == -1;
	}
};
Point Grid43::Approximator::LocalCoordinates(int c, Point p) const{
	if (is3[c]) return LocalCoordinates3(c, p);
	else return LocalCoordinates4(c, p);
}

Point Grid43::Approximator::LocalCoordinates3(int c, Point p) const{
	auto &x1 = cellvert[c][0]->x, &x2 = cellvert[c][1]->x, &x3 = cellvert[c][2]->x;
	auto &y1 = cellvert[c][0]->y, &y2 = cellvert[c][1]->y, &y3 = cellvert[c][2]->y;

	double j11 = x2 - x1, j21 = x3 - x1;
	double j12 = y2 - y1, j22 = y3 - y1;
	double modj = (j22*j11 - j21*j12);

	if (fabs(modj) < geps*geps){
		double ksi;
		if (isOnSection(p, *cellvert[c][0], *cellvert[c][1], ksi)){
			return Point(ksi, 0);
		} else if (isOnSection(p, *cellvert[c][1], *cellvert[c][2], ksi)){
			return Point(1-ksi, ksi);
		} else if (isOnSection(p, *cellvert[c][0], *cellvert[c][2], ksi)){
			return Point(0, ksi);
		} else {
			assert(false);
		}
	}

	Point ksieta;
	ksieta.x = ( j22*(p.x - x1) - j21*(p.y - y1))/modj;
	ksieta.y = (-j12*(p.x - x1) + j11*(p.y - y1))/modj;

	//
	//assert( ksieta.x > -0.1 && ksieta.x < 1.1 );
	//assert( ksieta.y > -0.1 && ksieta.y < 1-ksieta.x+0.1 );
	//

	return ksieta;
}

Point Grid43::Approximator::LocalCoordinates4(int c, Point p) const{
	_THROW_NOT_IMP_;
}

double Grid43::Approximator::Interpolate(int c, Point ksieta, const vector<double>& fun) const{
	if (is3[c]) return Interpolate3(c, ksieta, fun);
	else return Interpolate4(c, ksieta, fun);
}
double Grid43::Approximator::Interpolate3(int c, Point ksieta, const vector<double>& fun) const{
	int i1 = icellvert[c][0],
	    i2 = icellvert[c][1],
	    i3 = icellvert[c][2];
	return (1-ksieta.x-ksieta.y)*fun[i1] + (ksieta.x)*fun[i2] + (ksieta.y)*fun[i3];
}

double Grid43::Approximator::Interpolate4(int c, Point ksieta, const vector<double>& fun) const{
	_THROW_NOT_IMP_;
}

void Grid43::Approximator::FillJ3(std::array<double, 5>& J, int c) const{
	auto &x1 = cellvert[c][0]->x, &x2 = cellvert[c][1]->x, &x3 = cellvert[c][2]->x;
	auto &y1 = cellvert[c][0]->y, &y2 = cellvert[c][1]->y, &y3 = cellvert[c][2]->y;

	J[1] = x2 - x1; J[3] = x3 - x1;
	J[2] = y2 - y1, J[4] = y3 - y1;
	J[0] = (J[4]*J[1] - J[3]*J[2]);
}

void Grid43::Approximator::FillJ4(std::array<double, 5>& J, const Point& p, int c) const{
	_THROW_NOT_IMP_;
}

int Grid43::Approximator::FindPositive(const Point& p, Point& ksieta) const{
	auto candidates = cfinder->suspects(p);
	std::array<double, 5> J; //modj, j11, j12, j21, j22
	std::vector<Point> bad_ksieta(candidates.size());
	for (int i=0; i<candidates.size(); ++i){
		auto& c = candidates[i];
		auto& ke = bad_ksieta[i];
		if (is3[c]){FillJ3(J, c);}
		else { _THROW_NOT_IMP_;}
		if (fabs(J[0]) < geps*geps){
			double ksi;
			if (isOnSection(p, *cellvert[c][0], *cellvert[c][1], ksi)){
				ke = Point(ksi, 0);
			} else if (isOnSection(p, *cellvert[c][1], *cellvert[c][2], ksi)){
				ke = Point(1-ksi, ksi);
			} else if (isOnSection(p, *cellvert[c][0], *cellvert[c][2], ksi)){
				ke = Point(0, ksi);
			} else assert(false);
		} else {
			auto cp = cellvert[c][0];
			ke.x = ( J[4]*(p.x - cp->x) - J[3]*(p.y - cp->y))/J[0];
			ke.y = (-J[2]*(p.x - cp->x) + J[1]*(p.y - cp->y))/J[0];
		}
		if (J[0] > -geps*geps && ke.x>-geps && ke.x<1+geps && ke.y>-geps && ke.y<1-ke.x+geps){
			ksieta = ke;
			return c;
		}
	}
	//find best outer approximation
	int besti=-1;
	double bestdist = std::numeric_limits<double>::max();
	for (int i=0; i<candidates.size(); ++i){
		auto& ke = bad_ksieta[i];
		//if negative triangle contains point
		if (ke.x>-geps && ke.x<1+geps && ke.y>-geps && ke.y<1-ke.x+geps){
			besti = i;
			break;
		}
		//point lies outside a positive triangle: find measure to the basis triangle
		double d = std::min(Point::meas_section(ke, Point(0,0), Point(1,0)),
			std::min(Point::meas_section(ke, Point(1,0), Point(0,1)),
				Point::meas_section(ke, Point(0,1), Point(0,0))));
		if (d<bestdist) besti = i;
	}
	if (besti>-1){
		ksieta = bad_ksieta[besti];
		return candidates[besti];
	} else throw Grid43::Approximator::EOutOfArea();
}

double Grid43::Approximator::Val(Point p, const vector<double>& fun) const{
	Point ksieta;
	int c = FindPositive(p, ksieta);
	return Interpolate(c, ksieta, fun);
}


vector<double> Grid43::Approximator::Vals(Point p, const vector<const vector<double>*>& funs) const{
	Point ksieta;
	int c = FindPositive(p, ksieta);
	vector<double> ret;
	for (auto& fun: funs) ret.push_back(Interpolate(c, ksieta, *fun));
	return ret;
}

std::tuple<int, int, double> Grid43::Approximator::BndCoordinates(Point p) const{
	vector<int> susp = cfinder->suspects(p);
	int e1=-1, e2;
	double mindist = 1e32;
	double ksi;
	for (int c: susp){
		int mn = (is3[c]) ? 3 : 4;
		int i = 0;
		for (auto e: grid->vcells[c]->edges){
			if (e->is_boundary()){
				int in = (i==mn-1)?0:i+1;
				double ksi1;
				double dist = Point::meas_section(p, *cellvert[c][i], *cellvert[c][in], ksi1);
				if (dist < mindist){
					e1 = icellvert[c][i];
					e2 = icellvert[c][in];
					mindist = dist;
					ksi = ksi1;
				}
			}
			++i;
		}
	}
	if (e1<0) throw EOutOfArea();
	return std::make_tuple(e1, e2, ksi);
}

double Grid43::Approximator::BndVal(const Point& p, const vector<double>& fun) const{
	auto ap = BndCoordinates(p);
	int& e1 = std::get<0>(ap);
	int& e2 = std::get<1>(ap);
	double& ksi = std::get<2>(ap);
	return (1 - ksi)*fun[e1] + ksi*fun[e2];
}
vector<double> Grid43::Approximator::BndVals(const Point& p, const vector<const vector<double>*>& funs) const{
	auto ap = BndCoordinates(p);
	int& e1 = std::get<0>(ap);
	int& e2 = std::get<1>(ap);
	double& ksi = std::get<2>(ap);
	vector<double> ret;
	for (auto& fun: funs){
		double v = (1 - ksi)*(*fun)[e1] + ksi*(*fun)[e2];
		ret.push_back(v);
	}
	return ret;
}

Grid43::Approximator::NodePos Grid43::Approximator::FindNodePos(Point p) const{
	assert(std::all_of(is3.begin(), is3.end(), [](bool a){ return a; }));
	Grid43::Approximator::NodePos ret;
	try {
		double& k = ret.ksieta.x;
		double& e = ret.ksieta.y;
		ret.ncell = FindPositive(p, ret.ksieta);
		//have to check once more because FindPositive may return result
		//even for outer points.
		if (k<-geps || k>1+geps || e<-geps || e>1-k+geps) throw EOutOfArea();
		int inode = -1, iedge = -1;
		if (fabs(k)<geps){
			if (fabs(e)<geps){
				inode = 0;
			} else if (fabs(e-1)<geps){
				inode = 2;
			} else {
				iedge = 2;
			}
		} else if (fabs(e)<geps){
			if (fabs(k-1)<geps){
				inode = 1;
			} else {
				iedge = 0;
			}
		} else if (fabs(k+e-1.)<geps){
				iedge = 1;
		} else {
			ret.pos = ret.Internal;
			return ret;
		}
		if (inode >= 0){
			ret.nvert = icellvert[ret.ncell][inode];
			ret.pos = (isbnd[ret.nvert]) ? ret.BndVertex : ret.InternalVertex;
		} else {
			int next = (iedge==2) ? 0 : iedge+1;
			ret.nve1 = icellvert[ret.ncell][iedge];
			ret.nve2 = icellvert[ret.ncell][next];
			if (grid->vcells[ret.ncell]->edges[iedge]->is_boundary()){
				ret.pos = ret.BndEdge;
			} else ret.pos = ret.InternalEdge;
		}
	} catch (EOutOfArea){
		ret.pos = ret.Out;
	}
	return ret;
}

shared_ptr<HM2D::Vertex> Grid43::GuaranteePoint(HM2D::GridData& grid, Point pts,
		HMFem::Grid43::Approximator& approx){
	HMFem::Grid43::Approximator::NodePos pos;
	return GuaranteePoint(grid, pts, approx, pos);
}

shared_ptr<HM2D::Vertex> Grid43::GuaranteePoint(HM2D::GridData& grid, Point pts,
		HMFem::Grid43::Approximator& approx,
		HMFem::Grid43::Approximator::NodePos& pos){
	auto reassemble_cell = [&](int ic, shared_ptr<HM2D::Edge> e0, shared_ptr<HM2D::Edge> e1,
			shared_ptr<HM2D::Edge> e2, int v0, int v1, int v2){
		bool new_cell = false;
		if (grid.vcells.size() < ic+1){
			grid.vcells.resize(ic+1);
			grid.vcells[ic].reset(new HM2D::Cell());
			approx.cellvert.resize(ic+1);
			approx.icellvert.resize(ic+1);
			approx.is3.resize(ic+1, true);
			new_cell = true;
		}
		grid.vcells[ic]->edges.resize(3);
		grid.vcells[ic]->edges[0] = e0;
		grid.vcells[ic]->edges[1] = e1;
		grid.vcells[ic]->edges[2] = e2;
		approx.icellvert[ic][0] = v0;
		approx.icellvert[ic][1] = v1;
		approx.icellvert[ic][2] = v2;
		approx.cellvert[ic][0] = grid.vvert[v0].get();
		approx.cellvert[ic][1] = grid.vvert[v1].get();
		approx.cellvert[ic][2] = grid.vvert[v2].get();
		if (new_cell) approx.cfinder->addentry(HM2D::BBox(grid.vcells[ic]->edges));
	};

	auto direct_cell = [&](int ic, int ind){
		if (ind != approx.icellvert[ic][0]){
			int istart = (ind == approx.icellvert[ic][1]) ? 1 : 2;
			std::rotate(approx.cellvert[ic].begin(),
					approx.cellvert[ic].begin()+istart,
					approx.cellvert[ic].begin()+3);
			std::rotate(approx.icellvert[ic].begin(),
					approx.icellvert[ic].begin()+istart,
					approx.icellvert[ic].begin()+3);
			std::rotate(grid.vcells[ic]->edges.begin(),
					grid.vcells[ic]->edges.begin()+istart,
					grid.vcells[ic]->edges.end());
		}
		HM2D::Contour::R::ReallyDirect::Permanent(grid.vcells[ic]->edges);
	};
	auto set_neighbour = [&](int ie, int left, int right){
		if (left>=0) grid.vedges[ie]->left = grid.vcells[left];
		if (right>=0) grid.vedges[ie]->right = grid.vcells[right];
	};
	auto change_neighbour = [&](shared_ptr<HM2D::Edge> ed, int from, int to){
		if (ed->left.lock() == grid.vcells[from]) ed->left = grid.vcells[to];
		else ed->right = grid.vcells[to];
	};

	pos = approx.FindNodePos(pts);

	switch (pos.pos){
	case Approximator::NodePos::Out: return nullptr;
	case Approximator::NodePos::BndVertex: return grid.vvert[pos.nvert];
	case Approximator::NodePos::InternalVertex: return grid.vvert[pos.nvert];
	}

	int Nv = grid.vvert.size(), Ne = grid.vedges.size(), Nc = grid.vcells.size();
	auto ret = std::make_shared<HM2D::Vertex>(pts);
	approx.isbnd.push_back(pos.pos == pos.BndEdge);
	grid.vvert.push_back(ret);

	switch (pos.pos){
	case Approximator::NodePos::BndEdge:
	case Approximator::NodePos::InternalEdge:{
		direct_cell(pos.ncell, pos.nve1);
		HM2D::EdgeData& oe = grid.vcells[pos.ncell]->edges;
		auto icv = approx.icellvert[pos.ncell];
		grid.vcells[pos.ncell]->edges[0]->vertices[1] = ret;
		grid.vedges.push_back(std::make_shared<HM2D::Edge>(ret, grid.vvert[icv[1]])); //Ne
		grid.vedges.push_back(std::make_shared<HM2D::Edge>(ret, grid.vvert[icv[2]])); //Ne+1
		reassemble_cell(Nc, grid.vedges[Ne], oe[1], grid.vedges[Ne+1], //Nc
				Nv, icv[1], icv[2]);
		reassemble_cell(pos.ncell, oe[0], grid.vedges[Ne+1], oe[2],
				icv[0], Nv, icv[2]);

		set_neighbour(Ne, Nc, -1);
		set_neighbour(Ne+1, pos.ncell, Nc);
		change_neighbour(grid.vcells[Nc]->edges[1], pos.ncell, Nc);
		if (pos.pos == Approximator::NodePos::BndEdge) break;
		//second cell for non-boundary edges
		auto c2 = oe[0]->left.lock();
		if (c2 == grid.vcells[pos.ncell]) c2 = oe[0]->right.lock();
		int c2i = std::find(grid.vcells.begin(), grid.vcells.end(), c2) -
			grid.vcells.begin();
		direct_cell(c2i, pos.nve1);
		auto icv2 = approx.icellvert[c2i];
		HM2D::EdgeData& oe2 = grid.vcells[c2i]->edges;
		grid.vedges.push_back(std::make_shared<HM2D::Edge>(ret, grid.vvert[icv2[1]])); //Ne+2
		reassemble_cell(Nc+1, grid.vedges[Ne], grid.vedges[Ne+2], oe2[1], //Nc+1
				icv2[2], Nv, icv2[1]);
		reassemble_cell(c2i, oe2[0], grid.vedges[Ne+2], oe2[2],
				icv2[0], icv2[1], Nv);
		set_neighbour(Ne, Nc, Nc+1);
		set_neighbour(Ne+2, Nc+1, c2i);
		change_neighbour(grid.vcells[Nc+1]->edges[2], c2i, Nc+1);

		break;
	}
	case Approximator::NodePos::Internal:{
		HM2D::EdgeData& oe = grid.vcells[pos.ncell]->edges;
		auto icv = approx.icellvert[pos.ncell];

		grid.vedges.push_back(std::make_shared<HM2D::Edge>(ret, grid.vvert[icv[0]])); //Ne
		grid.vedges.push_back(std::make_shared<HM2D::Edge>(ret, grid.vvert[icv[1]])); //Ne+1
		grid.vedges.push_back(std::make_shared<HM2D::Edge>(ret, grid.vvert[icv[2]])); //Ne+2

		reassemble_cell(Nc,  oe[1], grid.vedges[Ne+2], grid.vedges[Ne+1],  //Nc
				icv[1], icv[2], Nv);
		reassemble_cell(Nc+1, oe[2], grid.vedges[Ne], grid.vedges[Ne+2], //Nc+1
				icv[2], icv[0], Nv);
		reassemble_cell(pos.ncell, oe[0], grid.vedges[Ne+1], grid.vedges[Ne],
				icv[0], icv[1], Nv);
		
		set_neighbour(Ne, pos.ncell, Nc+1);
		set_neighbour(Ne+1, Nc, pos.ncell);
		set_neighbour(Ne+2, Nc+1, Nc);
		change_neighbour(grid.vcells[Nc]->edges[0], pos.ncell, Nc);
		change_neighbour(grid.vcells[Nc+1]->edges[0], pos.ncell, Nc+1);
		break;
	}
	}
	return ret;
}


// ============================ Isoline Builder
Grid43::IsolineBuilder::IsolineBuilder(shared_ptr<HM2D::GridData>& grid43,
		shared_ptr<Grid43::Approximator> approx43){
	grid = grid43;
	if (approx43) approx = approx43;
	else approx.reset(new Approximator(grid.get()));
}

bool Grid43::IsolineBuilder::AddLine(int c,
		HM2D::EdgeData& ecol,
		double value, const vector<double>& fun) const{
	int i0 = approx->icellvert[c][0];
	int i1 = approx->icellvert[c][1];
	int i2 = approx->icellvert[c][2];
	double f0 = fun[i0], f1 = fun[i1], f2 = fun[i2];
	if (ISEQ(f0, f1) && ISEQ(f1, f2)) return false;
	//sort f0<f1<f2
	if (f1>f2) {std::swap(i1, i2); std::swap(f1, f2); };
	if (f0>f1) {std::swap(i0, i1); std::swap(f0, f1); };
	if (f1>f2) {std::swap(i1, i2); std::swap(f1, f2); };

	//check all possible positions
	if (value<f0 || value > f2) return false;
	else if (ISEQ(value, f0)){
		if (ISGREATER(f1, f0)) return false;
		else BuildLine(i0, i1, 0, i0, i1, 1, ecol);
	} else if (ISEQ(value, f1)){
		if (!ISEQ(f1,f2)) return false;
		else BuildLine(i1, i2, 0, i1, i2, 1, ecol);
	} else if (ISEQ(value, f2)){
		return false;
	} else if (value < f1){
		double w1 = (value - f0)/(f1 - f0);
		double w2 = (value - f0)/(f2 - f0);
		BuildLine(i0, i1, w1, i0, i2, w2, ecol);
	} else if (value < f2){
		double w1 = (value - f1)/(f2 - f1);
		double w2 = (value - f0)/(f2 - f0);
		BuildLine(i1, i2, w1, i0, i2, w2, ecol);
	} else return false;

	return true;
}

void Grid43::IsolineBuilder::BuildLine(int i0, int i1, double wi,
		int j0, int j1, double wj,
		HM2D::EdgeData& ecol) const{
	auto pi = std::make_shared<HM2D::Vertex>(
		Point::Weigh(*grid->vvert[i0], *grid->vvert[i1], wi));
	auto pj = std::make_shared<HM2D::Vertex>(
		Point::Weigh(*grid->vvert[j0], *grid->vvert[j1], wj));
	ecol.emplace_back(new HM2D::Edge(pi, pj));
}

HM2D::EdgeData Grid43::IsolineBuilder::FromPoint(Point pstart, const vector<double>& fun) const{
	//find value
	double value = approx->Val(pstart, fun);
	//assemble all edges of contour
	HM2D::EdgeData edges;
	for (int i=0; i<grid->vcells.size(); ++i)
		AddLine(i, edges, value, fun);
	//collect into contours
	auto et = HM2D::Contour::Tree::Assemble(edges);
	//choose contour which passes pstart
	HM2D::EdgeData* cnt;
	if (et.nodes.size() == 1) cnt = &et.nodes[0]->contour;
	else{
		for (int i=0; i<et.nodes.size(); ++i){
			cnt = &et.nodes[i]->contour;
			if (ISZERO(std::get<4>(HM2D::Contour::CoordAt(*cnt, pstart)))) break;
		}
	}
	
	HM2D::EdgeData ret;
	HM2D::DeepCopy(*cnt, ret);
	return ret;
}

HMCallback::FunctionWithCallback<HMFem::TAuxGrid3> HMFem::AuxGrid3;
double TAuxGrid3::step_estimate(const HM2D::Contour::Tree& tree, int nrec, double hrec){
	if (nrec <= 0) return hrec;
	double area = tree.area();
	double triarea = area/2/nrec;
	//from an area of equilateral triangle
	double hrec2 = sqrt(2.3094017677*triarea);
	if (hrec <= 0) return hrec2;
	return std::max(hrec, hrec2);
}

void TAuxGrid3::input(const HM2D::Contour::Tree& _tree, const vector<HM2D::EdgeData>& _constraints){
	clear();
	//deep copy edges input to internal structure
	tree = HM2D::Contour::Tree::DeepCopy(_tree);
	tree.remove_detached();
	constraints.reserve(_constraints.size());
	for (auto& c: _constraints){
		constraints.emplace_back(new HM2D::EdgeData());
		HM2D::DeepCopy(c, *constraints.back());
	}
	//fill mandatory points
	for (int i=0; i<constraints.size(); ++i){
		for (int j=0; j<constraints.size(); ++j) if (i!=j){
			mandatory_intersections(*constraints[i], *constraints[j]);
		}
		for (auto& c: tree.nodes) mandatory_intersections(*constraints[i], c->contour);
	}
}

void TAuxGrid3::mandatory_intersections(HM2D::EdgeData& c1, HM2D::EdgeData& c2){
	auto cres = HM2D::Contour::Finder::CrossAll(c1, c2);
	for (auto& it: cres){
		Point& p = std::get<1>(it);
		auto gp1 = HM2D::Contour::Algos::GuaranteePoint(c1, p);
		auto gp2 = HM2D::Contour::Algos::GuaranteePoint(c2, p);
		mandatory_points.insert(std::get<1>(gp1));
		mandatory_points.insert(std::get<1>(gp2));
	}
}

void TAuxGrid3::adopt_contour(HM2D::EdgeData& cont, double h,
		vector<vector<Point>>& lost){
	auto repart = HM2D::Contour::Algos::Coarsening(cont, mandatory_points, h,
			ZEROANGLE, 0.2);
	HM2D::VertexData ret;
	for (int i=0; i<repart.size(); ++i){
		auto& v = repart[i];
		ret.push_back(v[0]);
		if (v.size() > 1){
			lost.push_back(vector<Point>());
			for (int i=0; i<v.size(); ++i){
				lost.back().push_back(*v[i]);
			}
			lost.back().push_back(*repart[(i+1) % repart.size()][0]);
		}
	}
	cont = HM2D::Contour::Assembler::Contour1(ret, HM2D::Contour::IsClosed(cont));
}

void TAuxGrid3::adopt_boundary(HM2D::Contour::Tree& tree, double h,
		vector<vector<Point>>& lost){
	for (auto n: tree.nodes) adopt_contour(n->contour, h, lost);
}

namespace{
void adopt_lost(int ilost, vector<vector<Point>>& veclost, Point& ret){
	vector<Point>& lost = veclost[ilost];
	//find closest point amoung lost internals
	int iclosest = 1;
	double mclosest = Point::meas(lost[iclosest], ret);
	for (int i=2; i<lost.size()-1; ++i){
		double m1 = Point::meas(lost[i], ret);
		if (m1<mclosest){ iclosest = i; mclosest = m1; }
	}
	//lost index at which it will be divided
	int idiv=-1;
	double edge_meas = Point::meas(lost[0], lost.back());
	if (mclosest < 0.04 * edge_meas){
		//if distance between closest lost and ret is less than 0.2*edge_length
		//move ret point to lost points
		ret = lost[iclosest]; 
		idiv = iclosest;
	} else {
		//place ret to veclost points.
		//We hope they are sorted (is it guaranteed ???).
		HM2D::EdgeData contlost = HM2D::Contour::Constructor::FromPoints(lost);
		HM2D::Contour::Algos::GuaranteePoint(contlost, ret);
		lost.clear();
		for (auto p: HM2D::Contour::OrderedPoints1(contlost)) {
			lost.push_back(*p);
			if (*p == ret) idiv = lost.size();
		}
	}
	if (idiv<1 && idiv > lost.size()-2){
	} else if (idiv == 1){
		lost.erase(lost.begin());
	} else if (idiv == lost.size()-2){
		lost.resize(lost.size()-1);
	} else {
		vector<Point> newlost(lost.begin()+idiv, lost.end());
		lost = vector<Point>(lost.begin(), lost.begin() + idiv + 1);
		if (newlost.size()>2) veclost.push_back(newlost);
	}
	if (lost.size()<3){ veclost.erase(veclost.begin()+ilost); }
}

Point divide_edge(const HM2D::Edge* ed, vector<vector<Point>>& lost){
	//if this edge exists amoung lost points, restore it from those points
	Point p1 = *ed->first(), p2 = *ed->last();
	Point ret = ed->center();
	for (int i=0; i<lost.size(); ++i){
		if ((p1 == lost[i][0] && p2 == lost[i].back()) ||
		    (p2 == lost[i][0] && p1 == lost[i].back())){
			adopt_lost(i, lost, ret);
			break;
		}
	}
	return ret;
}
};

void TAuxGrid3::adopt_complicated_connections(HM2D::Contour::Tree& tree,
		vector<vector<Point>>& lost){
	//this searhes all complicated connections (more than 2 edges for vertex)
	//and checks for smooth section lengths transitions. If one edge is
	//3 or more times shorter than any sibling edge, divides this edge
	//All simple connections are guaranteed to be smooth as a result of adopt_contour procedure.
	HM2D::EdgeData cont = tree.alledges();
	
	//get complicated connections with bad size transitions
	vector<HM2D::Connectivity::VertexEdgeR> ae;
	vector<vector<double>> alens;
	for (auto& it: HM2D::Connectivity::VertexEdge(cont)){
		if (it.size()<3) continue;
		std::vector<double> lens;
		for (auto ei: it.eind) lens.push_back(cont[ei]->length());
		if (*max_element(lens.begin(), lens.end()) >
				*min_element(lens.begin(), lens.end())*3.5){
			ae.push_back(it);
			alens.push_back(lens);
		}
	}
	if (ae.size() == 0) return;

	std::set<int> div_edges;
	//find edges which has to be divided
	for (int i=0; i<ae.size(); ++i){
		double h = *min_element(alens[i].begin(), alens[i].end()) * 3.5;
		for (int j=0; j<alens[i].size(); ++j){
			if (alens[i][j]>h) div_edges.insert(ae[i].eind[j]);
		}
	}
	//divide
	for (auto de: div_edges) {
		auto node = tree.find_node(cont[de].get());
		Point pnew = divide_edge(cont[de].get(), lost);
		aa::enumerate_ids_pvec(node->contour);
		HM2D::Contour::Algos::SplitEdge(node->contour, cont[de]->id, {pnew});
	}
	//check one more time
	return adopt_complicated_connections(tree, lost);
}

HM2D::GridData HMFem::TAuxGrid3::_run(const HM2D::Contour::Tree& _tree,
		const vector<HM2D::EdgeData>& _constraints,
		int nrec, int nmax, double hrec){
	assert(nrec<nmax);
	//1)
	callback->step_after(5, "Estimate step");
	double hest = step_estimate(_tree, nrec, hrec);
	//2)
	callback->step_after(15, "Adopt boundary");
	input(_tree, _constraints);
	vector<vector<Point>> lost_points;
	adopt_boundary(tree, CORRECTION_FACTOR*hest, lost_points);
	for (auto c: constraints) adopt_contour(*c, CORRECTION_FACTOR*hest, lost_points);
	for (auto c: constraints) tree.add_detached_contour(std::move(*c));
	if (tree.detached_contours().size() > 0){
		auto ae = tree.alledges();
		HM2D::ECol::Algos::MergePoints(ae);
		assert(ae.size() == tree.alledges().size());
		adopt_complicated_connections(tree, lost_points);
	}
	//3) Triangulation
	auto subcaller = callback->bottom_line_subrange(60);
	auto ret = HM2D::Mesher::UnstructuredTriangle.UseCallback(subcaller, tree);
	if (ret.vvert.size() > nmax) { clear(); throw std::runtime_error("Failed to build auxiliary triangle grid");}
	if (lost_points.size() == 0) return ret;
	//4)
	callback->step_after(20, "Snapping");
	Grid43::AddSegments(ret, lost_points);
	return ret;
}
void HMFem::TAuxGrid3::clear(){
	//explicit clear because there is the reusable entry of TAuxGrid class
	tree = HM2D::Contour::Tree();
	constraints.clear();
	mandatory_points.clear();
}

HM2D::GridData HMFem::TAuxGrid3::_run(const HM2D::Contour::Tree& tree, int nrec, int nmax, double hrec){
	return _run(tree, {}, nrec, nmax, hrec);
}
HM2D::GridData HMFem::TAuxGrid3::_run(const HM2D::EdgeData& cont, int nrec, int nmax, double hrec){
	HM2D::Contour::Tree tree;
	tree.add_contour(cont);
	return _run(tree, nrec, nmax, hrec);
}
