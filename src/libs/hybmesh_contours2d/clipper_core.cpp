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

Container<Contour> ClipperPath::ToHMContainer() const{
	return HMContainer(this->data, bbox);
}

Container<Contour> ClipperPath::HMContainer(const ClipperLib::Path& path, const BoundingBox& bbox){
	Container<Contour> res;
	ClipperObject tob(bbox);
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
		assert(!nd->IsOpen());
		Container<Contour> c = ClipperPath::HMContainer(nd->Contour, bbox);
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





