#include "clipper_core.hpp"

using namespace HMCont2D;
using namespace HMCont2D::Impl;

void ClipperObject::ApplyBoundingBox(const BoundingBox& newbbox){
	bbox = newbbox;
	//factor and p0, p1
	long double maxlen = (long double) std::max(bbox.lenx(), bbox.leny()) / 2.0;
	factor = (long double)CLIPPER_RESOLUTION / maxlen;
	p0 = newbbox.Center();
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


ClipperPath::ClipperPath(const Contour& path){
	auto pnt = path.ordered_points();
	ApplyBoundingBox(BoundingBox::Build(pnt.begin(), pnt.end()));
	std::transform(pnt.begin(), pnt.end(), std::back_inserter(data),
		[&](Point* p){ return this->ToIntGeom(*p);
	});
	if (path.is_closed()) data.resize(data.size()-1);
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

Container<ContourTree> ClipperPath::Offset(double delta, HMCont2D::OffsetTp tp) const{
	long double dd = (long double)delta * factor;
	double arctol = ClipperArcTolerance * factor;
	ClipperLib::ClipperOffset worker;
	
	ClipperLib::EndType et;
	switch (tp){
		case HMCont2D::OffsetTp::CLOSED_POLY: et = ClipperLib::etClosedPolygon; break;
		case HMCont2D::OffsetTp::OPEN_ROUND: et = ClipperLib::etOpenRound; break;
		case HMCont2D::OffsetTp::OPEN_BUTT: et = ClipperLib::etOpenButt; break;
	}
	worker.AddPath(this->data, ClipperLib::jtRound, et);
	worker.ArcTolerance = arctol;
	ClipperLib::PolyTree sol;
	worker.Execute(sol, dd);
	return ClipperTree::HMContainer(sol, bbox);
}

Container<Contour> ClipperPath::ToHMContainer(){
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
		return fabs(ar)<=eps;
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

Container<Contour> ClipperPath::HMContainer(ClipperLib::Path& path, const BoundingBox& bbox, bool is_closed){
	Container<Contour> res;
	ClipperObject tob(bbox);
	//fixing some clipper bugs
	bool valid = clipper_heal(path, tob, is_closed);
	if (!valid) return res;
	
	//points
	std::transform(path.begin(), path.end(), std::back_inserter(res.pdata.data),
		[&tob](const ClipperLib::IntPoint& p){ return std::make_shared<Point>(tob.ToRealGeom(p)); });
	//edges
	for (int i=0; i<res.pdata.size()-1; ++i) res.add_value(Edge(res.pdata.point(i), res.pdata.point(i+1)));
	res.add_value(Edge(res.pdata.data.back().get(), res.pdata.data[0].get()));
	return res;
}

Container<ContourTree> ClipperTree::HMContainer(const ClipperLib::PolyTree& tree, const BoundingBox& bbox){
	Container<ContourTree> ret;
	if (tree.GetFirst() == 0) return ret;
	//1) create data filling only parents fields
	std::function<void(ClipperLib::PolyNode*, ContourTree::TreeNode*)>
	setcont = [&ret, &bbox, &setcont](ClipperLib::PolyNode* nd, ContourTree::TreeNode* parent){
		if (nd->IsOpen()) return;
		//build contour
		Container<Contour> c = ClipperPath::HMContainer(nd->Contour, bbox);
		if (c.size() == 0) return;
		//Using raw push_backs to eschew ret tree rebuilding
		//add points
		ret.pdata.Unite(c.pdata);
		//add edges
		std::copy(c.data.begin(), c.data.end(), std::back_inserter(ret.data));
		shared_ptr<ContourTree::TreeNode> newnode(new ContourTree::TreeNode);
		newnode->Unite(c); //add only edges
		newnode->parent = parent;
		ret.nodes.push_back(newnode);
		for (auto& child: nd->Childs) setcont(child, newnode.get());
	
	};
	for (int i=0; i<tree.ChildCount(); ++i) setcont(tree.Childs[i], nullptr);

	//2) restore children
	ret.UpdateTopology();
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
			ret.WidenWithPoint(pths1[i].bbox.BottomLeft());
			ret.WidenWithPoint(pths1[i].bbox.TopRight());
		}
	}
	for (int i=0; i<pths2.size(); ++i){
		if (pths2[i].data.size() == 0) continue;
		if (uninit) {ret = pths2[i].bbox; uninit = false; }
		else{
			ret.WidenWithPoint(pths2[i].bbox.BottomLeft());
			ret.WidenWithPoint(pths2[i].bbox.TopRight());
		}
	}
	assert(!uninit);
	for (int i=0; i<pths1.size(); ++i) pths1[i].ApplyBoundingBox(ret);
	for (int i=0; i<pths2.size(); ++i) pths2[i].ApplyBoundingBox(ret);

	return ret;
}

}

Container<ContourTree> ClipperPath::Intersect(
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

Container<ContourTree> ClipperPath::Union(
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

Container<ContourTree> ClipperPath::Substruct(
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

Container<ExtendedTree> ClipperPath::CutLines(
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
	Container<ExtendedTree> ret;
	ret.data = treecont.data;
	ret.nodes = treecont.nodes;
	ret.pdata = treecont.pdata;
	ClipperLib::Paths op;
	ClipperLib::OpenPathsFromPolyTree(res, op);

	for (auto& path: op){
		auto hm = ClipperPath::HMContainer(path, bbox, false);
		//HMCont2D::Debug::geogebra_contour(hm);
	}

	return ret;
}


