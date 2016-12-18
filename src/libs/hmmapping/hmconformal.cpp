#include "hmconformal.hpp"
#include "scpack_port.hpp"
#include "dscpack_port.hpp"
#include "hmcompute.hpp"
#include "confrect_fem.hpp"

using namespace HMMap::Conformal;

namespace{

double get_len(const vector<Point>& path, int i0, int i2){
	double len = 0.0;
	for (int i=i0; i<i2; ++i){
		int inext = (i==path.size()-1) ? 0 : i+1;
		len+=Point::dist(path[i], path[inext]);
	}
	return len;
}

}

bool Options::CanUseRectApprox(const vector<Point>& path, int i1, int i2, int i3) const{
	if (!use_rect_approx) return false;
	//check length
	if (length_weight >= 1.0){
		auto check_length = [&](int st, int en)->bool{
			int i1 = (en == path.size()) ? 0 : en;
			return get_len(path, st, en) <
				length_weight*Point::dist(path[st], path[i1]);
		};
		if (check_length( 0, i1)==false) return false;
		if (check_length(i1, i2)==false) return false;
		if (check_length(i2, i3)==false) return false;
		if (check_length(i3, path.size())==false) return false;
	}
	//check angles
	auto check_ang = [&](int j0, int j1, int j2){
		double a = Angle(path[j0], path[j1], path[j2]);
		return  ISEQGREATER(a, M_PI/2.0-right_angle_eps) && 
			ISEQGREATER(M_PI/2.0+right_angle_eps, a);
	};
	if (check_ang( 0, i1, i2)==false) return false;
	if (check_ang(i1, i2, i3)==false) return false;
	if (check_ang(i2, i3, 0 )==false) return false;
	if (check_ang(i3, 0 , i1)==false) return false;

	return true;
}

bool Options::CanUseSCPACK(const vector<Point>& path, int i1, int i2, int i3) const{
	if (!use_scpack) return false;
	double lenmin = get_len(path, 0, i1);
	double lenmax = lenmin;
	double l2 = get_len(path, i1, i2);
	double l3 = get_len(path, i2, i3);
	double l4 = get_len(path, i3, path.size());
	lenmin = std::min(lenmin, l2);
	lenmin = std::min(lenmin, l3);
	lenmin = std::min(lenmin, l4);
	lenmax = std::max(lenmax, l2);
	lenmax = std::max(lenmax, l3);
	lenmax = std::max(lenmax, l4);
	double hap = std::max(lenmin/lenmax, lenmax/lenmin);
	return (hap < scpack_ratio_limit);
}
bool Options::CanUseDSCPACK() const{
	return use_scpack;
}

vector<Point> Rect::MapToPolygon(const vector<Point*>& input) const{
	vector<Point> pvec(input.size());
	for (int i=0; i<input.size(); ++i){pvec[i].set(input[i]->x, input[i]->y);}
	return MapToPolygon(pvec);
}
vector<Point> Rect::MapToRectangle(const vector<Point*>& input) const{
	vector<Point> pvec(input.size());
	for (int i=0; i<input.size(); ++i){pvec[i].set(input[i]->x, input[i]->y);}
	return MapToRectangle(pvec);
}

shared_ptr<Rect> Rect::Factory(
		const HMCont2D::Contour& left,
		const HMCont2D::Contour& right,
		const HMCont2D::Contour& bot,
		const HMCont2D::Contour& top,
		const Options& opt){
	auto input = FactoryInput(left, right, bot, top);
	return Factory(std::get<0>(input), std::get<1>(input), opt);
}

shared_ptr<Rect> Rect::Factory(
		const vector<Point>& _path,
		std::array<int, 4> corners,
		const Options& opt){
	//rebuild path so that corners[0] = 0;
	vector<Point> path(_path);
	std::rotate(path.begin(), path.begin() + corners[0], path.end());
	int i1 = corners[1] - corners[0]; if (i1<0) i1+=path.size();
	int i2 = corners[2] - corners[0]; if (i2<0) i2+=path.size();
	int i3 = corners[3] - corners[0]; if (i3<0) i3+=path.size();

	shared_ptr<Rect> ret;
	//check if we can use rectangle approximation
	//without building conformal mapping
	if (opt.CanUseRectApprox(path, i1, i2, i3)){
		ret = Impl::RectApprox::Build(path, i1, i2, i3);
	}
	if (ret) return ret;
	
	//So we can not. 
	//Lets check if polygon is not
	//too oblong and we can use SCPACK routines
	if (opt.CanUseSCPACK(path, i1, i2, i3)){
		ret = Impl::SCPack::ToRect::Build(path, i1, i2, i3);
	}
	if (ret && std::max(ret->module(), 1.0/ret->module()) < opt.scpack_ratio_limit) return ret;
	
/*
	//Maybe rectangle can be mirrored into a closed
	//doubly connected one and be mapped into a ring?
	_THROW_NOT_IMP_;
	
	//Here should be a check for the possibility of
	//using modified Schwarz-Christoffel algorithm
	//see "A modified Schwarz-Christoffel transformation for elongated regions (1990)"
	//     by Louis H , Howellt , Lloyd N. Trefethen 
	_THROW_NOT_IMP_;
*/
	
	//The only chance left is a direct fem solution
	//of the Laplace problem
	ret = std::make_shared<Impl::ConfFem::ToRect>(Impl::ConfFem::ToRect::Build(path, i1, i2, i3, opt));

	if (!ret) throw std::runtime_error("Failed to build mapping to rectangle");

	return ret;
}

std::tuple<vector<Point>, std::array<int, 4>>
Rect::FactoryInput(const HMCont2D::Contour& left, const HMCont2D::Contour& right,
		const HMCont2D::Contour& bot, const HMCont2D::Contour& top){
	std::tuple<vector<Point>, std::array<int, 4>> ret;
	auto& vp = std::get<0>(ret);
	auto& crn = std::get<1>(ret);

	assert(*left.first() == *bot.first());
	assert(*left.last() == *top.first());
	assert(*right.first() == *bot.last());
	assert(*right.last() == *top.last());
	vector<Point*> lp = left.ordered_points(); std::reverse(lp.begin(), lp.end());
	vector<Point*> rp = right.ordered_points();
	vector<Point*> tp = top.ordered_points(); std::reverse(tp.begin(), tp.end());
	vector<Point*> bp = bot.ordered_points();
	for (int i=0; i<lp.size()-1; ++i) vp.push_back(*lp[i]);
	for (int i=0; i<bp.size()-1; ++i) vp.push_back(*bp[i]);
	for (int i=0; i<rp.size()-1; ++i) vp.push_back(*rp[i]);
	for (int i=0; i<tp.size()-1; ++i) vp.push_back(*tp[i]);
	int i1 = lp.size()-1;
	int i2 = i1 + bp.size()-1;
	int i3 = i2 + rp.size()-1;
	crn = {0, i1, i2, i3};
	return ret;
}

std::tuple<vector<Point>, std::array<int, 4>>
Rect::FactoryInput(const HMCont2D::Contour& cont, const std::array<int,4>& a){
	assert(cont.is_closed());
	assert(a[0]<a[1] && a[1]<a[2] && a[2]<a[3]);
	assert(a[0]>=0 && a[3]<cont.size());

	std::tuple<vector<Point>, std::array<int, 4>> ret;
	auto& rcont = std::get<0>(ret);
	auto& ra = std::get<1>(ret);

	//copy points from contour
	for (auto p: cont.ordered_points()) rcont.push_back(*p);
	rcont.pop_back();
	ra = a;

	//rotating points to ensure ra[0] = 0
	if (ra[0] > 0){
		std::rotate(rcont.begin(), rcont.begin()+ra[0], rcont.end());
		ra[1]-=ra[0]; if (ra[1] < 0) ra[1] += cont.size();
		ra[2]-=ra[0]; if (ra[2] < 0) ra[2] += cont.size();
		ra[3]-=ra[0]; if (ra[3] < 0) ra[3] += cont.size();
		ra[0] = 0;
	}

	return ret;
}

HMCallback::FunctionWithCallback<HMMap::Conformal::TBuildRect> HMMap::Conformal::BuildRect;
shared_ptr<Rect> TBuildRect::_run(const HMCont2D::Contour& left,
		const HMCont2D::Contour& right,
		const HMCont2D::Contour& bot,
		const HMCont2D::Contour& top,
			const Options& opt){
	callback->step_after(5, "Sorting input");
	auto inp = Rect::FactoryInput(left, right, bot, top);
	auto& _path = std::get<0>(inp);
	auto& corners = std::get<1>(inp);
	//rebuild path so that corners[0] = 0;
	vector<Point> path(_path);
	std::rotate(path.begin(), path.begin() + corners[0], path.end());
	int i1 = corners[1] - corners[0]; if (i1<0) i1+=path.size();
	int i2 = corners[2] - corners[0]; if (i2<0) i2+=path.size();
	int i3 = corners[3] - corners[0]; if (i3<0) i3+=path.size();

	callback->step_after(5, "Trying approximations");
	shared_ptr<Rect> ret;
	//check if we can use rectangle approximation
	//without building conformal mapping
	if (opt.CanUseRectApprox(path, i1, i2, i3)){
		ret = Impl::RectApprox::Build(path, i1, i2, i3);
	}
	if (ret) return ret;
	
	//So we can not. 
	//Lets check if polygon is not
	//too oblong and we can use SCPACK routines
	if (opt.CanUseSCPACK(path, i1, i2, i3)){
		ret = Impl::SCPack::ToRect::Build(path, i1, i2, i3);
	}
	if (ret && std::max(ret->module(), 1.0/ret->module()) < opt.scpack_ratio_limit) return ret;
	
	//numerical solution of the Laplace problem
	auto subcaller = callback->bottom_line_subrange(90);
	return std::make_shared<Impl::ConfFem::ToRect>(
		Impl::ConfFem::ToRect::Build.UseCallback(subcaller, path, i1, i2, i3, opt));
}

shared_ptr<Impl::RectApprox>
Impl::RectApprox::Build(const vector<Point>& path, int i1, int i2, int i3){
	shared_ptr<Impl::RectApprox> ret(new RectApprox());
	//add points
	for (auto& p: path) ret->pcol.add_value(p);
	vector<Point*> allpnt = ret->pcol.pvalues();
	allpnt.push_back(allpnt[0]);
	//left contour
	ret->left = HMCont2D::Constructor::ContourFromPoints(allpnt.begin(),
			allpnt.begin()+i1+1, false);
	ret->bot = HMCont2D::Constructor::ContourFromPoints(allpnt.begin()+i1,
			allpnt.begin()+i2+1, false);
	ret->right = HMCont2D::Constructor::ContourFromPoints(allpnt.begin()+i2,
			allpnt.begin()+i3+1, false);
	ret->top = HMCont2D::Constructor::ContourFromPoints(allpnt.begin()+i3,
			allpnt.end(), false);
	ret->left.ReallyReverse();
	ret->top.ReallyReverse();

	//calc_options
	ret->_len_top = ret->top.length();
	ret->_len_bot = ret->bot.length();
	ret->_len_right = ret->right.length();
	ret->_len_left = ret->left.length();
	ret->_module = (ret->_len_top + ret->_len_bot)/(ret->_len_left+ret->_len_right);
	ret->_left_straight = ret->left.is_straight();
	ret->_bot_straight = ret->bot.is_straight();
	ret->_right_straight = ret->right.is_straight();
	ret->_top_straight = ret->top.is_straight();

	//choose a direction for approximation
	if (ret->_len_top + ret->_len_bot > ret->_len_left + ret->_len_right)
		ret->_xksi = false;
	else
		ret->_xksi = true;

	//maybe we need some checks before return?
	return ret;
}

shared_ptr<Impl::RectApprox>
Impl::RectApprox::Build(const HMCont2D::Contour& left, const HMCont2D::Contour& right,
	const HMCont2D::Contour& bot, const HMCont2D::Contour& top){
	auto prep = FactoryInput(left, right, bot, top);
	auto& vp = std::get<0>(prep);
	auto& a = std::get<1>(prep);
	return Build(vp, a[1], a[2], a[3]);
}

vector<Point> Impl::RectApprox::MapToPolygon(const vector<Point>& input) const{
	vector<Point> ret;
	//internal points
	if (_xksi) ret = MapToPolygon_xksi(input);
	else ret = MapToPolygon_yeta(input);

	//boundary points
	for (int i=0; i<input.size(); ++i){
		double ksi = input[i].x/_module;
		double eta = input[i].y;
		if (ISEQ(ksi,1)) ret[i] = HMCont2D::Contour::WeightPoint(right, eta);
		else if (ISEQ(ksi,0)) ret[i] = HMCont2D::Contour::WeightPoint(left, eta);
		else if (ISEQ(eta,0)) ret[i] = HMCont2D::Contour::WeightPoint(bot, ksi);
		else if (ISEQ(eta,1)) ret[i] = HMCont2D::Contour::WeightPoint(top, ksi);
	}
	return ret;
}

vector<Point> Impl::RectApprox::MapToPolygon_yeta(const vector<Point>& input) const{
	//gather all x and y components of input array to avoid
	//superficial computations
	// -- calculate weights
	TCoordSet xline_w;
	for (auto& p: input) xline_w.insert(p.x/_module);

	// -- calculate points from weights
	HMCont2D::PCollection bottomp =
		HMCont2D::Contour::WeightPoints(bot, vector<double>(xline_w.begin(), xline_w.end()));
	HMCont2D::PCollection topp =
		HMCont2D::Contour::WeightPoints(top, vector<double>(xline_w.begin(), xline_w.end()));
		
	// -- calculate coordinates of lines
	TCoordMap<std::pair<Point, Point>> yline; //from x coordinate
	auto topp1 = topp.begin();
	auto bottomp1 = bottomp.begin();
	for (auto& w: xline_w){
		auto& yn = yline[w*_module] = std::pair<Point, Point>(); 
		yn.first = *(*bottomp1++);
		yn.second = *(*topp1++);
	}

	// -- calculate points
	vector<Point> ret; ret.reserve(input.size());
	double ksieta[2];
	for (auto& p: input){
		auto& line2 = yline[p.x];
		double xmod = p.x/_module;
		Point pp;
		//need a special treatment for boundary points
		//to match boundaries exactly. will be done later
		if (ISEQ(xmod,1) || ISEQ(xmod,0) || ISEQ(p.y,0) || ISEQ(p.y,1)) pp = Point();
		else {
			pp = line2.first * (1-p.y) + line2.second * p.y;
		}
		ret.push_back(pp);
	}
	return ret;
}

vector<Point> Impl::RectApprox::MapToPolygon_xksi(const vector<Point>& input) const{
	//gather all x and y components of input array to avoid
	//superficial computations
	// -- calculate weights
	TCoordSet yline_w;
	for (auto& p: input) yline_w.insert(p.y);

	// -- calculate points from weights
	HMCont2D::PCollection leftp =
		HMCont2D::Contour::WeightPoints(left, vector<double>(yline_w.begin(), yline_w.end()));
	HMCont2D::PCollection rightp =
		HMCont2D::Contour::WeightPoints(right, vector<double>(yline_w.begin(), yline_w.end()));
		
	// -- calculate coordinates of lines
	TCoordMap<std::pair<Point, Point>> xline; //from y coordinate
	auto leftp1 = leftp.begin();
	auto rightp1 = rightp.begin();
	for (auto& w: yline_w){
		auto& xn = xline[w] = std::pair<Point, Point>(); 
		xn.first = *(*leftp1++);
		xn.second = *(*rightp1++);
	}

	// -- calculate points
	vector<Point> ret; ret.reserve(input.size());
	double ksieta[2];
	for (auto& p: input){
		auto& line1 = xline[p.y];
		double xmod = p.x/_module;
		Point pp;
		//need a special treatment for boundary points
		//to match boundaries exactly. will be done later
		if (ISEQ(xmod,1) || ISEQ(xmod,0) || ISEQ(p.y,0) || ISEQ(p.y,1)) pp = Point();
		else {
			pp = line1.first * (1-xmod) + line1.second * xmod;
		}
		ret.push_back(pp);
	}
	return ret;
}

Point Impl::RectApprox::MapToPolygon(const Point& input) const{
	return MapToPolygon(vector<Point> {input})[0];
}

vector<Point> Impl::RectApprox::MapToRectangle(const vector<Point>& input) const{
	vector<Point> ret; ret.reserve(input.size());
	for (auto& p: input) ret.push_back(MapToRectangle(p));
	return ret;
}

Point Impl::RectApprox::MapToRectangle(const Point& p) const{
	auto left_fnd = left.coord_at(p);
	auto right_fnd = right.coord_at(p);
	auto top_fnd = top.coord_at(p);
	auto bot_fnd = bot.coord_at(p);
	//if on side walls
	if (ISZERO(std::get<4>(left_fnd))){
		return Point(0, std::get<1>(left_fnd));
	} else if (ISZERO(std::get<4>(right_fnd))){
		return Point(_module, std::get<1>(right_fnd));
	} else if (ISZERO(std::get<4>(top_fnd))){
		return Point(_module*std::get<1>(top_fnd), 1);
	} else if (ISZERO(std::get<4>(bot_fnd))){
		return Point(_module*std::get<1>(bot_fnd), 0);
	}
	//if internal
	if (_xksi) return MapToRectangle_xksi(p);
	else return MapToRectangle_yeta(p);
}

Point Impl::RectApprox::MapToRectangle_xksi(const Point& p) const{
	//calculate eta
	auto func2 = [&](double x){
		Point start = HMCont2D::Contour::WeightPoint(left, x);
		Point end = HMCont2D::Contour::WeightPoint(right, x);
		return Point::meas_section(p, start, end);
	};
	double eta = HMMath::Compute::GoldenRatioMin(0.0, 1.0, func2, 100, 1e-16);
	assert(ISIN_EE(eta, 0, 1));

	double ksi;
	Point start, end;
	//Point coordinate by eta
	start = HMCont2D::Contour::WeightPoint(left, eta);
	end = HMCont2D::Contour::WeightPoint(right, eta);
	isOnSection(p, start, end, ksi);
	assert(ISIN_EE(ksi, 0, 1));
	return Point(ksi*_module, eta);
}

Point Impl::RectApprox::MapToRectangle_yeta(const Point& p) const{
	//calculate ksi
	auto func1 = [&](double x){
		Point start = HMCont2D::Contour::WeightPoint(bot, x);
		Point end = HMCont2D::Contour::WeightPoint(top, x);
		return Point::meas_section(p, start, end);
	};
	double ksi = HMMath::Compute::GoldenRatioMin(0.0, 1.0, func1, 100, 1e-16);
	assert(ksi>-geps && ksi<1+geps);

	double eta;
	Point start, end;
	//Point coordinate by ksi
	start = HMCont2D::Contour::WeightPoint(bot, ksi);
	end = HMCont2D::Contour::WeightPoint(top, ksi);
	isOnSection(p, start, end, eta);
	assert(eta>-geps && eta<1+geps);
	return Point(ksi*_module, eta);
}

HMCont2D::Container<HMCont2D::Contour> Rect::RectContour() const{
	auto pts = RectPoints();
	return HMCont2D::Constructor::ContourFromPoints(pts, true);
}

// ==================================== Annulus
shared_ptr<Annulus> Annulus::Factory(
		const vector<Point>& outerpnt,
		const vector<Point>& innerpnt,
		const Options& opt
){
	assert([&]()->bool{
		//anti-clockwise check
		auto oc = HMCont2D::Constructor::ContourFromPoints(outerpnt, true);
		auto ic = HMCont2D::Constructor::ContourFromPoints(innerpnt, true);
		if (HMCont2D::Area(oc)<=0) return false;
		if (HMCont2D::Area(ic)<=0) return false;
		return true;
	}());
	shared_ptr<Annulus> ret;
	//1) try DSCPACK procedure
	if (opt.CanUseDSCPACK()){
		ret = Impl::DSCPack::ToAnnulus::Build(outerpnt, innerpnt);
		if (ret) return ret;
	}

	//2) FEM procedure
	ret = Impl::ConfFem::ToAnnulus::Build(outerpnt, innerpnt, opt);
	if (ret) return ret;

	return 0;
}

HMCont2D::Container<HMCont2D::Contour> Annulus::InnerCircleContour() const{
	auto pts = InnerCirclePoints();
	return HMCont2D::Constructor::ContourFromPoints(pts, true);
}
HMCont2D::Container<HMCont2D::Contour> Annulus::OuterCircleContour() const{
	auto pts = OuterCirclePoints();
	return HMCont2D::Constructor::ContourFromPoints(pts, true);
}

