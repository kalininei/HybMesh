#include "finder2d.hpp"
#include "contour.hpp"
#include "clipper_core.hpp"
using namespace HM2D;

// =============== contains procedures
shared_ptr<Vertex> Finder::Contains(const VertexData& data, const Point* pnt){
	auto fnd = std::find_if(data.begin(), data.end(),
		[&pnt](const shared_ptr<Vertex>& pd){ return pd.get() == pnt; });
	if (fnd != data.end()) return *fnd;
	else return nullptr;
}
shared_ptr<Vertex> Finder::Contains(const EdgeData& data, const Point* pnt){
	for (auto e: data){
		if (e->vertices[0].get() == pnt) return e->vertices[0];
		if (e->vertices[1].get() == pnt) return e->vertices[1];
	}
	return nullptr;
}
shared_ptr<Edge> Finder::Contains(const EdgeData& data, const Edge* ed){
	auto fnd = std::find_if(data.begin(), data.end(),
		[&ed](const shared_ptr<Edge>& pd){ return pd.get() == ed; });
	if (fnd != data.end()) return *fnd;
	else return nullptr;
}
shared_ptr<Vertex> Finder::Contains(const CellData& data, const Point* pnt){
	return Contains(AllVertices(data), pnt);
}
shared_ptr<Edge> Finder::Contains(const CellData& data, const Edge* ed){
	return Contains(AllEdges(data), ed);
}
shared_ptr<Cell> Finder::Contains(const CellData& data, const Cell* c){
	auto fnd = std::find_if(data.begin(), data.end(),
		[&c](const shared_ptr<Cell>& pd){ return pd.get() == c; });
	if (fnd != data.end()) return *fnd;
	else return nullptr;
}

// ================ find closest
std::tuple<int, double, double>
Finder::ClosestEdge(const EdgeData& dt, const Point& p){
	std::tuple<int, double, double> ret;
	int& ind = std::get<0>(ret);
	double& dist = std::get<1>(ret);
	double& ksi = std::get<2>(ret);
	dist = 1e99; ksi = 1e99; ind = -1;

	double k;
	for (int i=0; i<dt.size(); ++i){
		auto& e = dt[i];
		double dnew = Point::meas_section(p, *e->first(), *e->last(), k);
		//if (dnew < dist){
		if (dnew - dist < geps*geps){
			dist = dnew;
			ksi = k;
			ind = i;
			if (dist<geps*geps) break;
		}
	}

	dist = sqrt(dist);
	return ret;
}

Point Finder::ClosestEPoint(const EdgeData& dt, const Point& p){
	auto fec = ClosestEdge(dt, p);
	auto e = dt[std::get<0>(fec)];
	return Point::Weigh(*e->first(), *e->last(), std::get<2>(fec));
}

std::tuple<int, double> Finder::ClosestPoint(const VertexData& dt, const Point& p){
	if (dt.size() == 0) return std::tuple<int, double>(-1, -1);
	std::tuple<int, double> ret(0, 0);
	int& bind = std::get<0>(ret);
	double& bm = std::get<1>(ret);
	bm = Point::meas(*dt[0], p);
	for (int i=1; i<dt.size(); ++i){
		double m = Point::meas(*dt[i], p);
		if (m < bm){
			bm = m;
			bind = i;
		}
	}
	bm = sqrt(bm);
	return ret;
}

// =============== EdgeFinder
Finder::EdgeFinder::EdgeFinder(const EdgeData& d){
	ve = Connectivity::VertexEdge(d);
	data = &d;
	//sort it using pointer compare
	auto srtve = [](const Connectivity::VertexEdgeR& a, const Connectivity::VertexEdgeR& b)->bool{
		return a.v.get() < b.v.get();
	};
	std::sort(ve.begin(), ve.end(), srtve);
}

std::tuple<shared_ptr<Edge>, bool>              
Finder::EdgeFinder::find(Vertex* v1, Vertex* v2){
	std::tuple<shared_ptr<Edge>, bool>ret(0, false);
	auto fndve = [](const Connectivity::VertexEdgeR& a, Vertex* b)->bool{
		return a.v.get() < b;
	};
	auto fnd = std::lower_bound(ve.begin(), ve.end(), v1, fndve);

	//v1 was not found
	if (fnd == ve.end() || v1 != fnd->v.get()) return ret;

	for (auto ei: fnd->eind){
		shared_ptr<Edge> e = (*data)[ei];
		if (e->first().get() == v1 && e->last().get() == v2){
			std::get<0>(ret) = e;
			std::get<1>(ret) = true;
			break;
		}
		if (e->last().get() == v1 && e->first().get() == v2){
			std::get<0>(ret) = e;
			std::get<1>(ret) = false;
			break;
		}
	}

	return ret;
}

// ================== VertexFinder
Finder::VertexMatch::VertexMatch(const VertexData& vd): srt(vd){
	auto _less = [](const shared_ptr<Vertex>& p1,
			const shared_ptr<Vertex>& p2)->bool{
		return *p1 < *p2;
	};
	std::sort(srt.begin(), srt.end(), _less);
}
shared_ptr<Vertex> Finder::VertexMatch::find(const Point& p){
	auto _less = [](const shared_ptr<Vertex>& p1, const Point& p2)->bool{
		return *p1 < p2;
	};
	auto fnd = std::lower_bound(srt.begin(), srt.end(), p, _less);
	if (fnd != srt.end() && **fnd == p) return *fnd;
	else return nullptr;
}

VertexData Finder::VertexMatch::find(const vector<Point>& p){
	VertexData ret;
	for (auto& it: p) ret.push_back(find(it));
	return ret;
}


int Contour::Finder::WhereIs(const EdgeData& ed, const Point& p, BoundingBox* bb){
	//whereis for contour trees uses this routine
	//passing all tree edges as 'ed'.
	//Hence 'ed' is not always a closed contour.
	assert(!IsContour(ed) || IsClosed(ed));
	BoundingBox bbox;
	if (bb == nullptr){
		bbox = BBox(ed, 0);
		bb = &bbox;
	}
	if (bb->whereis(p) == OUTSIDE) return OUTSIDE;

	//Contour::Finderculate number of crosses between ed and [p, pout]
	//where pout is random outside point
	Point pout(bb->xmax + 1.56749, bb->ymax + 1.06574);

	for (int tries=0; tries<100; ++tries){
		double ksieta[2];
		int ncrosses = 0;
		for (auto& e: ed){
			SectCross(p, pout, *e->first(), *e->last(), ksieta);
			if (ISIN_NN(ksieta[0], 0, 1) && ISIN_NN(ksieta[1], 0, 1)){
				++ncrosses;
			} else if (ISEQ(ksieta[0], 0) && ISIN_EE(ksieta[1], 0, 1)){
			       return BOUND;
			} else if (ISEQ(ksieta[1], 0) || ISEQ(ksieta[1], 1)){
				//[p, pout] crosses one of ed vertices.
				//This is ambiguous so we change pout and try again
				goto NEXTTRY;
			}
		}
		if (ncrosses % 2 == 0) return OUTSIDE;
		else return INSIDE;
NEXTTRY:
		pout+=Point(-0.4483, 0.0342);
	}

	throw std::runtime_error("failed to detect point-contour relation");
}

// =============== Crosses
namespace{
vector<std::tuple<bool, Point, double, double>>
cross_core(const EdgeData& c1, const EdgeData& c2, bool is1){
	vector<std::tuple<bool, Point, double, double>> retv;
	auto bb1 = HM2D::BBox(c1), bb2 = HM2D::BBox(c2);
	if (!bb1.has_common_points(bb2)) return retv;

	VertexData op1 = Contour::OrderedPoints(c1);
	VertexData op2 = Contour::OrderedPoints(c2);

	auto lens1 = Contour::ELengths(c1), lens2 = Contour::ELengths(c2);
	double flen1 = std::accumulate(lens1.begin(), lens1.end(), 0.0);
	double flen2 = std::accumulate(lens2.begin(), lens2.end(), 0.0);

	auto addcross = [&](Point p, double w1, double w2){
		retv.push_back(std::make_tuple(true, p, w1, w2));
	};
	double ksieta[2];
	double L1=0.0, L2=0.0;
	for (int i=0; i<op1.size()-1; ++i){
		L2 = 0;
		for (int j=0; j<op2.size()-1; ++j){
			SectCross(*op1[i], *op1[i+1], *op2[j], *op2[j+1], ksieta);
			if (ksieta[0]>-geps && ksieta[0]<1+geps && ksieta[1]>-geps && ksieta[1]<1+geps){
				addcross(
					Point::Weigh(*op1[i], *op1[i+1], ksieta[0]),
					(L1 + lens1[i]*ksieta[0])/flen1,
					(L2 + lens2[j]*ksieta[1])/flen2
				);
				if (is1) return retv;
			}
			L2+=lens2[j];
		}
		L1 += lens1[i];
	}
	if (retv.size() < 2) return retv;

	//clear duplicates
	TCoordSet w1;
	vector<std::tuple<bool, Point, double, double>> ret;
	for (auto& v: retv){
		if (w1.find(std::get<2>(v)) == w1.end()){
			w1.insert(std::get<2>(v));
			ret.push_back(v);
		}
	}
	if (ret.size() == 1) return ret;

	return ret;
}
}

std::tuple<bool, Point, double, double>
Contour::Finder::Cross(const EdgeData& c1, const EdgeData& c2){
	auto retv = cross_core(c1, c2, true);
	if (retv.size() == 0) return std::make_tuple(false, Point(0,0), 0.0, 0.0);
	else return retv[0];
}

vector<std::tuple<bool, Point, double, double>>
Contour::Finder::CrossAll(const EdgeData& c1, const EdgeData& c2){
	return cross_core(c1, c2, false);
}
std::tuple<bool, Point, double, double>
Contour::Finder::SelfCross(const EdgeData& c1){
	double ksi[2];
	bool closed = IsClosed(c1);

	for (int i=0; i<c1.size(); ++i)
	for (int j=i+2; j<c1.size(); ++j){
		if (closed && i==0 && j==c1.size()-1) continue;
		SectCross(*c1[i]->first(), *c1[i]->last(),
		          *c1[j]->first(), *c1[j]->last(), ksi);
		if (ISIN_EE(ksi[0], 0, 1) && ISIN_EE(ksi[1], 0, 1))
			return std::make_tuple(
				true, Point::Weigh(*c1[i]->pfirst(),
				                   *c1[i]->plast(), ksi[0]),
				ksi[0], ksi[1]);
	}
	return std::make_tuple(false, Point(0, 0), 0., 0.);
}

Finder::RasterizeEdges::RasterizeEdges(const EdgeData& ed, const BoundingBox& bb, double step){
	bbf.reset(new BoundingBoxFinder(bb, step));
	black_squares.resize(bbf->nsqr(), false);
	for (auto& e: ed){
		vector<int> is = bbf->sqrs_by_segment(*e->pfirst(), *e->plast());
		bbf->raw_addentry(is);
		for (auto i: is){ black_squares[i] = true; }
	}
}

namespace{
vector<int> fill_group(int first, vector<bool>& used, int nx, int ny){
	std::array<int, 4> adj;
	auto fill_adj = [&nx, &ny, &adj](int i){
		int ix = i % nx;
		int iy = i / nx;
		adj[0] = (ix != 0) ? (iy*nx+ix-1) : -1;
		adj[1] = (ix != nx-1) ? (iy*nx+ix+1) : -1;
		adj[2] = (iy != 0) ? ((iy-1)*nx+ix) : -1;
		adj[3] = (iy != ny-1) ? ((iy+1)*nx+ix) : -1;
	};

	vector<int> ret(1, first); used[first] = true;
	int uu=0;
	while (uu<ret.size()){
		fill_adj(ret[uu]);
		for (auto a: adj) if (a != -1 && !used[a]){
			ret.push_back(a);
			used[a] = true;
		}
		++uu;
	}

	return ret;
}
}

vector<int> Finder::RasterizeEdges::colour_squares(bool use_groups) const{
	vector<int> ret(bbf->nsqr(), 1);
	for (int i=0; i<bbf->nsqr(); ++i)
		if (black_squares[i])
			ret[i] = 0;
	if (!use_groups) return ret;

	//grouping non zero using graph split
	vector<bool> used = black_squares;
	int curgroup = 1;
	while (1){
		int first = std::find(used.begin(), used.end(), false)-used.begin();
		if (first >= used.size()) break;

		for (auto i: fill_group(first, used, bbf->nx(), bbf->ny())){
			ret[i] = curgroup;
		}
		++curgroup;
	}

	return ret;
}

vector<int> Contour::Finder::SortOutPoints(const EdgeData& t1, const vector<Point>& pnt){
	Contour::Tree tree;
	tree.add_contour(t1);
	auto ret = SortOutPoints(tree, pnt);
	//if t1 is inner contour
	if (Contour::Area(t1) < 0){
		for (auto& r: ret){
			if (r == OUTSIDE) r = INSIDE;
			else if (r == INSIDE) r = OUTSIDE;
		}
	}
	return ret;
}

namespace{
//does sorting using clipper procedure.
//It could be used only if points are known not to lie on edges
vector<int> clipper_sort_out_points(const Contour::Tree& t1, const vector<Point>& pnt){
	auto ctree = Impl::ClipperTree::Build(t1);
	vector<int> ret = ctree.SortOutPoints(pnt);
	assert([&](){
		//Boundary points should be treated somewhere else
		for (auto i: ret) if (i == BOUND) return false;
		return true;
	}());
	return ret;
}
//does sort by calculating intersections
vector<int> raw_sort_out_points(const Contour::Tree& t1, const vector<Point>& pnt){
	//bounding boxes
	std::map<Contour::Tree::TNode*, BoundingBox> bb;
	for (auto& n: t1.bound_contours()){
		bb[n.get()] = HM2D::BBox(n->contour);
	}

	std::function<void(const ShpVector<Contour::Tree::TNode>&, Point&, int&)>
	lvwithin = [&lvwithin, &bb](const ShpVector<Contour::Tree::TNode>& lv, Point& p, int& a){
		int indwithin = -1;
		for (int i=0; i<lv.size(); ++i){
			int rs = Contour::Finder::WhereIs(lv[i]->contour, p, &bb[lv[i].get()]);
			if (rs == BOUND){ a = -2; return; }
			else if (rs == INSIDE) { indwithin = i; break; }
		}
		if (indwithin == -1) return;
		else {
			ShpVector<Contour::Tree::TNode> newlv;
			for (auto w: lv[indwithin]->children) newlv.push_back(w.lock());
			return lvwithin(newlv, p, ++a);
		}
	};

	vector<int> ret;
	auto roots = t1.roots();
	for (auto p: pnt){
		int maxlevel = -1;
		lvwithin(roots, p, maxlevel);
		if (maxlevel == -2) ret.push_back(BOUND);
		else if (maxlevel == -1) ret.push_back(OUTSIDE);
		else if (maxlevel % 2 == 0) ret.push_back(INSIDE);
		else ret.push_back(OUTSIDE);
	}
	return ret;
}

//does sort by discretizing input contour, grouping discretized squares
//and calling raw_sort_out_points for single point for each group
vector<int> qsort_out_points(const Contour::Tree& t1, const vector<Point>& pnt){
	vector<int> ret(pnt.size());
	auto pntbb = BoundingBox::Build(pnt.begin(), pnt.end());

	HM2D::EdgeData ae = t1.alledges();
	HM2D::Finder::RasterizeEdges discr(ae, pntbb, pntbb.maxlen()/50);
	vector<int> sqr_colours = discr.colour_squares(true);

	//colours of input points
	vector<int> pnt_colours(pnt.size());
	for (int i=0; i<pnt.size(); ++i){
		vector<int> sqrs = discr.bbfinder().sqrs_by_point(pnt[i]);
		pnt_colours[i] = sqr_colours[sqrs[0]];
		for (int j=1; j<sqrs.size(); ++j){
			pnt_colours[i] = std::min(pnt_colours[i], sqr_colours[sqrs[j]]);
		}
	}

	//build groups->points map
	std::map<int, vector<int>> used_colours = { {0, vector<int>()} };
	for (int i=0; i<pnt.size(); ++i) if (pnt_colours[i] != 0){
		auto fnd = used_colours.find(pnt_colours[i]);
		if (fnd == used_colours.end()){
			fnd = used_colours.emplace(pnt_colours[i], vector<int>{}).first;
		}
		fnd->second.push_back(i);
	} else {
		//check [0] (black) group if those points lie on edges
		ret[i] = INSIDE;
		for (int sus: discr.bbfinder().suspects(pnt[i])){
			double m = Point::meas_section(pnt[i], *ae[sus]->pfirst(), *ae[sus]->plast());
			if (m < geps*geps) {ret[i] = BOUND; break; }
		}
		//point is not on edge -> add it for future analysis
		if (ret[i] == INSIDE) used_colours[0].push_back(i);
	}

	//assembling point list to be passed to RawSortOut procedure
	//groups points, then points of undefined group
	vector<Point> pnt_to_raw;
	for (auto& it: used_colours) if (it.first > 0){
		//first points of non-black groups
		pnt_to_raw.push_back(pnt[it.second[0]]);
	} else {
		//all points of the black group
		for (auto& it2: it.second){
			pnt_to_raw.push_back(pnt[it2]);
		}
	}
	//call raw procedure
	//vector<int> raw_output = raw_sort_out_points(t1, pnt_to_raw);
	vector<int> raw_output = clipper_sort_out_points(t1, pnt_to_raw);

	//restore result
	auto rit = raw_output.begin();
	for (auto& it: used_colours) if (it.first > 0){
		int pos = *rit++;
		for (auto& ind: it.second){ ret[ind] = pos; }
	} else {
		//all points of the black group
		for (auto& it2: it.second){ ret[it2] = *rit++; }
	}
	return ret;
}
}

vector<int> Contour::Finder::SortOutPoints(const Tree& t1, const vector<Point>& pnt){
	if (pnt.size() < 100) return raw_sort_out_points(t1, pnt);
	else return qsort_out_points(t1, pnt);
}

