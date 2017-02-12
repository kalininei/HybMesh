#include "clipper_core.hpp"
#include <limits>
#include "assemble2d.hpp"
#include "treverter2d.hpp"

using namespace HM2D;
using namespace HM2D::Impl;

void ClipperObject::ApplyBoundingBox(const BoundingBox& newbbox){
	bbox = newbbox;
	//factor and p0, p1
	long double maxlen = (long double) std::max(bbox.lenx(), bbox.leny()) / 2.0;
	if (maxlen<geps){
		factor = (long double)1.0;
	} else {
		factor = (long double)CLIPPER_RESOLUTION / maxlen;
		factor = std::min(factor, (long double) std::numeric_limits<ClipperLib::cInt>::max()/100.0);
	}
	p0 = newbbox.center();
}

ClipperLib::IntPoint ClipperObject::ToIntGeom(const Point& p) const{
	long double x = (long double)(p.x - p0.x);
	long double y = (long double)(p.y - p0.y);
	x *= factor; y *= factor; 
	return ClipperLib::IntPoint((ClipperLib::cInt)round(x), (ClipperLib::cInt)round(y));
}

Point ClipperObject::ToRealGeom(const ClipperLib::IntPoint& p) const{
	long double x =(long double)p.X, y = (long double)p.Y;
	x /= factor; y /= factor;
	return Point((double)x+p0.x, (double)y+p0.y);
}


ClipperPath::ClipperPath(const EdgeData& path){
	auto pnt = Contour::OrderedPoints(path);
	ApplyBoundingBox(BoundingBox::PBuild(pnt.begin(), pnt.end()));
	std::transform(pnt.begin(), pnt.end(), std::back_inserter(data),
		[&](shared_ptr<Point> p){ return this->ToIntGeom(*p);
	});
	if (Contour::IsClosed(path)) data.resize(data.size()-1);
}

void ClipperPath::AddPointToEnd(const Point& p){
	data.push_back(ToIntGeom(p));
}

void ClipperPath::ApplyBoundingBox(const BoundingBox& newbbox){
	//data -> backup
	vector<Point> tmp; tmp.reserve(data.size());
	std::transform(data.begin(), data.end(), std::back_inserter(tmp),
		[&](ClipperLib::IntPoint& p){ return this->ToRealGeom(p);
	});

	//write data
	ClipperObject::ApplyBoundingBox(newbbox);

	//backup -> data
	std::transform(tmp.begin(), tmp.end(), data.begin(),
		[&](Point& p){ return this->ToIntGeom(p);
	});
}

int ClipperPath::WhereIs(Point p) const{
	return ClipperLib::PointInPolygon(ToIntGeom(p), data);
}

Contour::Tree ClipperPath::Offset(double delta, ClipperLib::EndType et) const{
	long double dd = (long double)delta * factor;
	double arctol = ClipperArcTolerance * factor;
	ClipperLib::ClipperOffset worker;
	
	worker.AddPath(this->data, ClipperLib::jtRound, et);
	worker.ArcTolerance = arctol;
	ClipperLib::PolyTree sol;
	worker.Execute(sol, dd);
	return ClipperTree::HMContainer(sol, bbox);
}

EdgeData ClipperPath::ToHMContainer(){
	return HMContainer(this->data, bbox);
}

namespace{

bool clipper_heal(ClipperLib::Path& path, const ClipperObject& bbox, bool is_closed){
	//Sometimes Union gives trees with singular contours with negative areas.
	//Here we ignore them.
	//if (is_closed && fabs(ClipperLib::Area(path)) <= 1.0) return false;
	//integer size which is neglectable in real geometry
	//long int eps = geps*bbox.factor;
	//long int eps = 0;
	////get rid of edges which turns at 2*pi
	////get rid of doubling nodes
	//std::set<int> not_needed;
	//auto eqnodes = [&](int i, int j){
	//        return (labs(path[i].X - path[j].X)<=eps &&
	//                        labs(path[i].Y - path[j].Y)<=eps);
	//};
	//auto is_collinear = [&](int i0, int i1, int i2){
	//        long int x0 = path[i1].X - path[i0].X;
	//        long int x1 = path[i2].X - path[i0].X;
	//        long int y0 = path[i1].Y - path[i0].Y;
	//        long int y1 = path[i2].Y - path[i0].Y;
	//        long int ar = x0*y1 - y0*x1;
	//        return ar<=eps*eps;
	//};
	//for (int i=0; i<path.size(); ++i){
	//        if (!is_closed && (i==0 || i==path.size()-1)) continue;
	//        int iprev = (i == 0) ? path.size() - 1 : i-1;
	//        int inext = (i == path.size() - 1) ? 0 : i+1;
	//        if (eqnodes(i,iprev)) not_needed.insert(i);
	//        if (eqnodes(iprev, inext)) {not_needed.insert(i); not_needed.insert(inext); }
	//        if (is_collinear(iprev, i, inext)) {not_needed.insert(i); }
	//}

	//if (not_needed.size() == 0) return true;
	//else if (path.size() - not_needed.size() < 2) return false;
	//else if (path.size() - not_needed.size() == 2) return !is_closed;
	//else{
	//        ClipperLib::Path np;
	//        for (int i=0; i<path.size(); ++i){
	//                if (not_needed.find(i) == not_needed.end()){
	//                        np.push_back(path[i]);
	//                }
	//        }
	//        path = np;
	//        return true;
	//}

	//If this is a closed area all points of which lie on the same
	//line return invalid status
	long int eps = geps*1e-2*sqr(bbox.factor);
	auto is_collinear = [&](int i0, int i1, int i2){
		long int x0 = path[i1].X - path[i0].X;
		long int x1 = path[i2].X - path[i0].X;
		long int y0 = path[i1].Y - path[i0].Y;
		long int y1 = path[i2].Y - path[i0].Y;
		long int ar = x0*y1 - y0*x1;
		return std::abs(ar)<=eps;
	};
	if (is_closed){
		if (path.size() < 3) return false;
		bool okstatus = false;
		for (int i1=0; i1<path.size(); ++i1){
			int i0 = (i1==0) ? path.size()-1 : i1-1;
			int i2 = (i1==path.size()-1) ? 0 : i1+1;
			if (!is_collinear(i0, i1, i2)) { okstatus = true; break; }
			
		}
		if (!okstatus) return false;
	}

	return true;
}

}

EdgeData ClipperPath::HMContainer(ClipperLib::Path& path, const BoundingBox& bbox, bool is_closed){
	EdgeData res;
	ClipperObject tob(bbox);
	//fixing some clipper bugs
	bool valid = clipper_heal(path, tob, is_closed);
	if (!valid) return res;
	
	if (path.size() < 2) return res;
	if (is_closed && path.size() < 3) return res;

	VertexData pdata;
	pdata.push_back(std::make_shared<Vertex>(tob.ToRealGeom(*path.begin())));
	auto it = path.begin() + 1;
	while (it != path.end()){
		Point pnext = tob.ToRealGeom(*it++);
		if (Point::meas(pnext, *pdata.back()) > 1e3*geps*geps){
			pdata.push_back(std::make_shared<Vertex>(pnext));
		}
	}
	if (is_closed){
		if (Point::meas(*pdata[0], *pdata.back()) <= 1e3*geps*geps){
			pdata.resize(pdata.size()-1);
		}
	}
	res = Contour::Assembler::Contour1(pdata, is_closed);
	return res;
}

Contour::Tree ClipperTree::HMContainer(const ClipperLib::PolyTree& tree, const BoundingBox& bbox){
	Contour::Tree ret;
	if (tree.GetFirst() == 0) return ret;
	//1) create data filling only parents fields
	std::function<void(ClipperLib::PolyNode*, shared_ptr<Contour::Tree::TNode>)>
	setcont = [&ret, &bbox, &setcont](ClipperLib::PolyNode* nd, shared_ptr<Contour::Tree::TNode> parent){
		if (nd->IsOpen()) return;
		//build contour
		EdgeData c = ClipperPath::HMContainer(nd->Contour, bbox);
		if (c.size() == 0) return;
		//Using raw push_backs to eschew ret tree rebuilding
		//add edges
		shared_ptr<Contour::Tree::TNode> newnode(new Contour::Tree::TNode);
		newnode->contour = c;
		newnode->parent = parent;
		ret.nodes.push_back(newnode);
		for (auto& child: nd->Childs) setcont(child, newnode);
	
	};
	for (int i=0; i<tree.ChildCount(); ++i) setcont(tree.Childs[i], nullptr);

	//2) restore children
	ret.update_topology();
	return ret;
}

void ClipperTree::ApplyBoundingBox(const BoundingBox& newbbox){
	ClipperObject::ApplyBoundingBox(newbbox);
	for (auto& p: data) p.ApplyBoundingBox(newbbox);
}

//from tree
ClipperTree ClipperTree::Build(const Contour::Tree& tree){
	HM2D::Contour::R::RevertTree rv(tree);
	ClipperTree ret;
	ret.ApplyBoundingBox(HM2D::BBox(tree.alledges()));
	for (auto n: tree.nodes){
		ret.isopen.push_back(false);
		ret.data.push_back(ClipperPath());
		ClipperPath& last = ret.data.back();
		last.ApplyBoundingBox(ret.bbox);
		auto sp = Contour::OrderedPoints(n->contour);
		for (int i=0; i<n->contour.size(); ++i){
			last.data.push_back(ret.ToIntGeom(*sp[i]));
		}
	}
	return ret;
}

namespace{
int whereis(const ClipperLib::Paths& paths, const ClipperLib::IntPoint& pnt){
	int level = 0;
	for (auto& path: paths){
		int r = ClipperLib::PointInPolygon(pnt, path);
		if (r == -1) return BOUND;
		if (r == 1) level +=1;
	}
	return level%2==0 ? OUTSIDE : INSIDE;
}
int whereis(const vector<ClipperPath>& paths, const ClipperLib::IntPoint& pnt){
	int level = 0;
	for (auto& path: paths){
		int r = ClipperLib::PointInPolygon(pnt, path.data);
		if (r == -1) return BOUND;
		if (r == 1) level +=1;
	}
	return level%2==0 ? OUTSIDE : INSIDE;
}
}

//sorting physical points: INSIDE/OUTSIDE/BOUND for each point
vector<int> ClipperTree::SortOutPoints(const vector<Point>& pts) const{
	vector<int> ret;
	//scaling points
	vector<ClipperLib::IntPoint> ipts;
	for (auto& p: pts) ipts.push_back(ToIntGeom(p));

	for (auto& p: ipts){
		ret.push_back(whereis(data, p));
	}

	//This was a code which tried to overcome difficulties
	//with BOUND feature. Now current procedure is used only for
	//non-bound nodes and it is not need. I hope.
	//
	////assembling offset procedure
	//ClipperLib::ClipperOffset os;
	//for (int i=0; i<data.size(); ++i) if (!isopen[i]){
		//os.AddPath(data[i].data, ClipperLib::jtSquare, ClipperLib::etClosedPolygon);
	//}
	////inner, outer
	//ClipperLib::Paths inner, outer;
	//double k = geps * factor;
	//os.Execute(inner, -k);
	//os.Execute(outer, k);

	////sorting points
	//for (auto& p: ipts){
		//if (whereis(inner, p) == INSIDE) ret.push_back(INSIDE);
		//else{
			//if (whereis(outer, p) == OUTSIDE) ret.push_back(OUTSIDE);
			//else ret.push_back(BOUND);
		//}
	//}
	
	return ret;
}


namespace{
BoundingBox SameBoundingBox(
		vector<ClipperPath>& pths1,
		vector<ClipperPath>& pths2){
	BoundingBox ret;
	bool uninit = true;
	for (int i=0; i<pths1.size(); ++i){
		if (pths1[i].data.size() == 0) continue;
		if (uninit) {ret = pths1[i].bbox; uninit = false; }
		else{
			ret.widen(pths1[i].bbox.pmin());
			ret.widen(pths1[i].bbox.pmax());
		}
	}
	for (int i=0; i<pths2.size(); ++i){
		if (pths2[i].data.size() == 0) continue;
		if (uninit) {ret = pths2[i].bbox; uninit = false; }
		else{
			ret.widen(pths2[i].bbox.pmin());
			ret.widen(pths2[i].bbox.pmax());
		}
	}
	assert(!uninit);
	for (int i=0; i<pths1.size(); ++i) pths1[i].ApplyBoundingBox(ret);
	for (int i=0; i<pths2.size(); ++i) pths2[i].ApplyBoundingBox(ret);

	return ret;
}

}

Contour::Tree ClipperPath::Intersect(
		vector<ClipperPath>& pths1,
		vector<ClipperPath>& pths2,
		bool embedded1,
		bool embedded2){
	auto bbox = SameBoundingBox(pths1, pths2);
	//collect pathes
	ClipperLib::Clipper clrp;
	clrp.StrictlySimple(true);
	for (auto& c: pths1) clrp.AddPath(c.data, ClipperLib::ptSubject, true);
	for (auto& c: pths2) clrp.AddPath(c.data, ClipperLib::ptClip, true);
	//execute
	ClipperLib::PolyTree res;
	auto t1 = (embedded1) ? ClipperLib::pftEvenOdd : ClipperLib::pftNonZero;
	auto t2 = (embedded2) ? ClipperLib::pftEvenOdd : ClipperLib::pftNonZero;
	clrp.Execute(ClipperLib::ctIntersection, res, t1, t2);
	return ClipperTree::HMContainer(res, bbox);
}

Contour::Tree ClipperPath::Union(
		vector<ClipperPath>& pths1,
		vector<ClipperPath>& pths2,
		bool embedded1,
		bool embedded2){
	auto bbox = SameBoundingBox(pths1, pths2);
	//collect pathes
	ClipperLib::Clipper clrp;
	clrp.StrictlySimple(true);
	for (auto& c: pths1) clrp.AddPath(c.data, ClipperLib::ptSubject, true);
	for (auto& c: pths2) clrp.AddPath(c.data, ClipperLib::ptClip, true);
	//execute
	ClipperLib::PolyTree res;
	auto t1 = (embedded1) ? ClipperLib::pftEvenOdd : ClipperLib::pftNonZero;
	auto t2 = (embedded2) ? ClipperLib::pftEvenOdd : ClipperLib::pftNonZero;
	bool isok = clrp.Execute(ClipperLib::ctUnion, res, t1, t2);
	assert(isok);
	return ClipperTree::HMContainer(res, bbox);
}

Contour::Tree ClipperPath::Substruct(
		vector<ClipperPath>& pths1,
		vector<ClipperPath>& pths2,
		bool embedded1,
		bool embedded2){
	auto bbox = SameBoundingBox(pths1, pths2);
	//collect pathes
	ClipperLib::Clipper clrp;
	clrp.StrictlySimple(true);
	for (auto& c: pths1) clrp.AddPath(c.data, ClipperLib::ptSubject, true);
	for (auto& c: pths2) clrp.AddPath(c.data, ClipperLib::ptClip, true);
	//execute
	ClipperLib::PolyTree res;
	auto t1 = (embedded1) ? ClipperLib::pftEvenOdd : ClipperLib::pftNonZero;
	auto t2 = (embedded2) ? ClipperLib::pftEvenOdd : ClipperLib::pftNonZero;
	clrp.Execute(ClipperLib::ctDifference, res, t1, t2);
	return ClipperTree::HMContainer(res, bbox);
}

Contour::Tree ClipperPath::CutLines(
		vector<ClipperPath>& area,
		vector<ClipperPath>& lines,
		bool embedded_area){
	_THROW_NOT_IMP_;
	auto bbox = SameBoundingBox(area, lines);
	ClipperLib::Clipper clrp;
	clrp.StrictlySimple(true);
	for (auto& c: area) clrp.AddPath(c.data, ClipperLib::ptClip, true);
	for (auto& c: lines) clrp.AddPath(c.data, ClipperLib::ptSubject, false);

	//execute
	ClipperLib::PolyTree res;
	auto t2 = (embedded_area) ? ClipperLib::pftEvenOdd : ClipperLib::pftNonZero;
	clrp.Execute(ClipperLib::ctIntersection, res);
	auto treecont = ClipperTree::HMContainer(res, bbox);
	Contour::Tree ret = treecont;
	ClipperLib::Paths op;
	ClipperLib::OpenPathsFromPolyTree(res, op);

	for (auto& path: op){
		auto hm = ClipperPath::HMContainer(path, bbox, false);
		//HMCont2D::Debug::geogebra_contour(hm);
	}

	return ret;
}
