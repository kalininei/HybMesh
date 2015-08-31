#include "clipper_core.h"


using namespace HMCont2DCore;

Path::Path(){
	p0 = Point(0, 0);
	factor = 1.0;
}

ClipperLib::IntPoint Path::ToIntGeom(const Point& p){
	long double x = (long double)(p.x - p0.x);
	long double y = (long double)(p.y - p0.y);
	x *= factor; y *= factor; 
	return ClipperLib::IntPoint((ClipperLib::cInt)round(x), (ClipperLib::cInt)round(y));
}

Point Path::ToRealGeom(const ClipperLib::IntPoint& p){
	long double x =(long double)p.X, y = (long double)p.Y;
	x /= factor; y /= factor;
	return Point((double)x+p0.x, (double)y+p0.y);
}

int Path::ClosedDirection() const{
	return ClipperLib::Orientation(data) ? INSIDE : OUTSIDE;
}

bool Path::IsWithin(const Point& p){
	return ClipperLib::PointInPolygon(ToIntGeom(p), data) == 1;
}

void Path::ApplyBoundingBox(const BoundingBox& newbbox){
	//data -> backup
	vector<Point> tmp; tmp.reserve(data.size());
	for (auto& p: data) tmp.push_back(ToRealGeom(p));
	//new scaling factors
	bbox = newbbox;
	p0 = bbox.Center();
	long double maxlen = (long double) std::max(bbox.lenx(), bbox.leny()) / 2.0;
	factor = (long double)CLIPPER_RESOLUTION / maxlen;
	//backup -> data
	for (int i=0; i<data.size(); ++i) data[i] = ToIntGeom(tmp[i]);
}

void Path::SetPoints(const ShpVector<Point>& p){
	data.clear();
	ApplyBoundingBox(BoundingBox::Build(p.begin(), p.end()));
	for (auto& pnt: p) data.push_back(ToIntGeom(*pnt));
}

void Path::AddPointToEnd(const Point& p){
	if (data.size() == 0){
		BoundingBox newbbox(p, 1.0);
		Path::ApplyBoundingBox(newbbox);
	} else if (bbox.whereis(p) == OUTSIDE){
		BoundingBox newbbox(bbox); newbbox.WidenWithPoint(p);
		Path::ApplyBoundingBox(newbbox);
	}
	data.push_back(ToIntGeom(p));
}

double Path::Area() const{ return fabs(SignedArea()); }
double Path::SignedArea() const{
	return ClipperLib::Area(data)/factor/factor;
}
void Path::Reverse(){
	std::reverse(data.begin(), data.end());
}
