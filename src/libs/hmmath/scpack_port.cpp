#include "scpack_port.hpp"

using namespace HMMath::Conformal::Impl::SCPack;

extern "C"{

double scpack_init_(
		//input
		const int* n,
		const Pt wcoords[],
		const Pt* w0,
		const int corners[],
		const int* prec,
		//ouput
		Pt zcoords[],  //array[n]
		Pt* factor,
		Pt* factor2,
		Pt wcoords2[], //array[4]
		Pt* w02,
		double betam[],   //array[n]
		double qwork[],   //array[prec*(2*n+3)]
		double qwork2[]   //array[prec*(2*4+3)]
);

Pt scpack_forward_(
		const Pt* point,
		const int* n,
		const Pt wcoords[],
		const Pt* w0,
		const Pt wcoords2[],
		const Pt* w02,
		const Pt zcoords[], 
		const Pt zcoords2[],
		const Pt* factor,
		const Pt* factor2,
		const double betam[],
		const double qwork[],
		const double qwork2[],
		const int *prec
);

Pt scpack_backward_(
		const Pt* point,
		const int* n,
		const Pt wcoords[],
		const Pt* w0,
		const Pt wcoords2[],
		const Pt* w02,
		const Pt zcoords[], 
		const Pt zcoords2[],
		const Pt* factor,
		const Pt* factor2,
		const double betam[],
		const double qwork[],
		const double qwork2[],
		const int *prec
);

}

ToRect::ToRect(const vector<Pt>& dt, Pt p0, std::array<int, 4> i1, int _prec):
		wcoords(dt), w0(p0),
		corners({i1[0]+1, i1[1]+1, i1[2]+1, i1[3]+1}),
		betam(dt.size()), zcoords(dt.size()),
		qwork(_prec*(2*dt.size()+3)), qwork2(_prec*11), wcoords2(4),
		N(dt.size()), prec(_prec),
		_module(scpack_init_(
			&N,
			&wcoords[0],
			&w0,
			corners.data(),
			&prec,
			&zcoords[0],
			&factor,
			&factor2,
			&wcoords2[0],
			&w02,
			&betam[0],
			&qwork[0],
			&qwork2[0]
		)){
	zcoords2 = {zcoords[i1[0]], zcoords[i1[1]], zcoords[i1[2]], zcoords[i1[3]]};
}

shared_ptr<ToRect>
ToRect::Build(const vector<Point>& pnt, int i1, int i2, int i3){
	//find center point: rough procedure.
	Point p1(0,0), p2(0,0), p3(0,0), p4(0,0);
	for (int i=0; i<i1+1; ++i) p1+=pnt[i]; p1/=(i1+1);
	for (int i=i1; i<i2+1; ++i) p2+=pnt[i]; p2/=(i2-i1+1);
	for (int i=i2; i<i3+1; ++i) p3+=pnt[i]; p3/=(i3-i2+1);
	for (int i=i3; i<pnt.size(); ++i) p4+=pnt[i]; p4+=pnt[0]; p4/=(pnt.size()-i3+1);
	Pt pc = {(p1.x + p2.x + p3.x + p4.x)/4, (p1.y + p2.y + p3.y + p4.y)/4};
	vector<Pt> input(pnt.size());
	std::transform(pnt.begin(), pnt.end(), input.begin(),
			[](const Point& p){ return Pt{p.x, p.y}; }
	);
	
	auto ret =  new ToRect(input, pc, std::array<int, 4>{0, i1, i2, i3}, 12);
	if (ret->module() < 0) return 0;
	else return shared_ptr<ToRect>(ret);
}

shared_ptr<ToRect>
ToRect::Build(const HMCont2D::Contour& left, const HMCont2D::Contour& right,
		const HMCont2D::Contour& bot, const HMCont2D::Contour& top){
	auto prep = FactoryInput(left, right, bot, top);
	auto& vp = std::get<0>(prep);
	auto& a = std::get<1>(prep);
	return Build(vp, a[1], a[2], a[3]);
}

vector<Point> ToRect::MapToPolygon(const vector<Point>& input) const{
	vector<Point> ret; ret.reserve(input.size());
	for (auto& p: input){
		//find among corner points
		auto fnd = std::find_if(wcoords2.begin(), wcoords2.end(), [&](const Pt& x){
			return ISEQ(x.x, p.x) && ISEQ(x.y, p.y);
		});
		if (fnd!=wcoords2.end()){
			int i = fnd - wcoords2.begin();
			const Pt& pt = wcoords[corners[i] - 1];
			ret.push_back(Point(pt.x, pt.y));
		} else {
		//Do mapping
			Pt t{p.x, p.y};
			Pt t2 = scpack_backward_(
				&t, &N,
				&wcoords[0], &w0,
				&wcoords2[0], &w02,
				&zcoords[0], &zcoords2[0],
				&factor, &factor2,
				&betam[0], &qwork[0], &qwork2[0],
				&prec
			);
			ret.push_back(Point(t2.x, t2.y));
		}
	}
	return ret;
}


vector<Point> ToRect::MapToRectangle(const vector<Point>& input) const{
	vector<Point> ret; ret.reserve(input.size());
	for (auto& p: input){
		//find among corner points
		auto fnd = std::find_if(wcoords.begin(), wcoords.end(), [&](const Pt& x){
			return ISEQ(x.x, p.x) && ISEQ(x.y, p.y);
		});
		int i = fnd - wcoords.begin();
		if (i == corners[0] - 1){
			ret.push_back(Point(wcoords2[0].x, wcoords2[0].y));
		} else if (i == corners[1] - 1){
			ret.push_back(Point(wcoords2[1].x, wcoords2[1].y));
		} else if (i == corners[2] - 1){
			ret.push_back(Point(wcoords2[2].x, wcoords2[2].y));
		} else if (i == corners[3] - 1){
			ret.push_back(Point(wcoords2[3].x, wcoords2[3].y));
		} else {
		//Do mapping
			Pt t{p.x, p.y};
			Pt t2 = scpack_forward_(
				&t, &N,
				&wcoords[0], &w0,
				&wcoords2[0], &w02,
				&zcoords[0], &zcoords2[0],
				&factor, &factor2,
				&betam[0], &qwork[0], &qwork2[0],
				&prec
			);
			ret.push_back(Point(t2.x, t2.y));
		}
	}
	return ret;
}

vector<Point> ToRect::RectPoints() const{
	vector<Point> inp;
	std::transform(wcoords.begin(), wcoords.end(), std::back_inserter(inp),
			[](const Pt& p){ return Point(p.x, p.y); });
	return MapToRectangle(inp);
}








