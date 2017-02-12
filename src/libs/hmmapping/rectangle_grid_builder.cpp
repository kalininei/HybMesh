#include "rectangle_grid_builder.hpp"
#include "hmmapping.hpp"
#include "hmfdm.hpp"
#include "hmconformal.hpp"
#include "buildcont.hpp"
#include "assemble2d.hpp"
#include "modcont.hpp"
#include "healgrid.hpp"
#include "buildgrid.hpp"
#include "finder2d.hpp"
#include "treverter2d.hpp"

HMCallback::FunctionWithCallback<HMMap::TOrthogonalRectGrid> HMMap::OrthogonalRectGrid;
HMCallback::FunctionWithCallback<HMMap::TLaplaceRectGrid> HMMap::LaplaceRectGrid;

namespace{

HM2D::EdgeData deepcopy(const HM2D::EdgeData& cont, bool isrev){
	HM2D::EdgeData ret;
	HM2D::DeepCopy(cont, ret);
	if (isrev) HM2D::Contour::R::ReallyRevert::Permanent(ret);
	else HM2D::Contour::R::ReallyDirect::Permanent(ret);

	return ret;
}

bool has_self_cross(const HM2D::EdgeData& cont){
	double ksieta[2];
	for (int i=0; i<(int)cont.size()-2; ++i){
		for (int j=i+2; j<cont.size(); ++j){
			Point& p1 = *cont[i]->first();
			Point& p2 = *cont[i]->last();
			Point& p3 = *cont[j]->first();
			Point& p4 = *cont[j]->last();
			SectCross(p1, p2, p3, p4, ksieta);
			if (ksieta[0]>-geps && ksieta[0]<1+geps &&
				ksieta[1]>-geps && ksieta[1]<1+geps) return true;
		}
	}
	return false;
}
bool no_cross_except_touch(const HM2D::EdgeData& c1, const HM2D::EdgeData& c2){
	//return (HM2D::Contour::Algos::Finder::CrossAll(c1, c2).size() != 1);
	int n_touches = 0;
	Point p11 = *HM2D::Contour::First(c1), p12 = *HM2D::Contour::Last(c1);
	if (p11 == *HM2D::Contour::First(c2) || p11 == *HM2D::Contour::Last(c2)){
		++n_touches;
		HM2D::Contour::First(c1)->set(c1[0]->center());
	}
	if (p12 == *HM2D::Contour::First(c2) || p12 == *HM2D::Contour::Last(c2)){
		++n_touches;
		HM2D::Contour::Last(c1)->set(c1.back()->center());
	}
	bool ret;
	if (n_touches != 1) ret = false;
	else{
		ret = !std::get<0>(HM2D::Contour::Finder::Cross(c1, c2));
	}
	HM2D::Contour::First(c1)->set(p11);
	HM2D::Contour::Last(c1)->set(p12);
	return ret;
}
bool check_for_no_cross(const HM2D::EdgeData& left, const HM2D::EdgeData& bot,
		const HM2D::EdgeData& right, const HM2D::EdgeData& top){
	if (!no_cross_except_touch(left, bot)) return false;
	if (!no_cross_except_touch(left, top)) return false;
	if (std::get<0>(HM2D::Contour::Finder::Cross(left, right))) return false;
	if (std::get<0>(HM2D::Contour::Finder::Cross(bot, top))) return false;
	if (!no_cross_except_touch(bot, right)) return false;
	if (!no_cross_except_touch(top, right)) return false;
	if (has_self_cross(right)) return false;  //only right was modified pointwisely
	return true;
}

double connect_vec_points(const HM2D::VertexData& left, const HM2D::VertexData& bot,
		const HM2D::VertexData& right, const HM2D::VertexData& top){
	Point bot_move=*left[0] - *bot[0];
	for (auto p: bot) *p += bot_move;
	Point top_move=*left.back() - *top[0];
	for (auto p: top) *p += top_move;
	Point right_move1 = *bot.back() - *right[0];
	Point right_move2 = *top.back() - *right.back();
	auto rw = HM2D::Contour::EWeights(HM2D::Contour::Assembler::Contour1(right));
	for (int i=0; i<rw.size(); ++i){
		*right[i] += (right_move1 * (1.0-rw[i]) + right_move2 * rw[i]);
	}
	//using 1.001 because if bot_move and top_move are equal than moving top is better
	return 1.001*vecLen(bot_move) + vecLen(top_move) +
		vecLen(right_move1) + vecLen(right_move2);
}

double tryopt(int opt, const HM2D::EdgeData& _left, const HM2D::EdgeData& _bot,
		const HM2D::EdgeData& _right, const HM2D::EdgeData& _top){
	bool left_rev = opt & 1;
	bool right_rev = opt & 2;
	bool top_rev = opt & 4;
	bool bot_rev = opt & 8;

	auto left = deepcopy(_left, left_rev);
	auto right = deepcopy(_right, right_rev);
	auto top = deepcopy(_top, top_rev);
	auto bot = deepcopy(_bot, bot_rev);

	double ret = connect_vec_points(
			HM2D::Contour::OrderedPoints(left),
			HM2D::Contour::OrderedPoints(bot),
			HM2D::Contour::OrderedPoints(right),
			HM2D::Contour::OrderedPoints(top));

	if (check_for_no_cross(left, bot, right, top)){
		return ret;
	} else {
		return -1;
	}
}

vector<HM2D::EdgeData> modify_contours(int opt, HM2D::EdgeData& left, HM2D::EdgeData& bot,
		HM2D::EdgeData& right, HM2D::EdgeData& top){
	ShpVector<HM2D::Contour::R::ReallyRevert> rr;
	ShpVector<HM2D::Contour::R::ReallyDirect> rd;

	if (opt & 1) rr.emplace_back(new HM2D::Contour::R::ReallyRevert(left));
	else rd.emplace_back(new HM2D::Contour::R::ReallyDirect(left));

	if (opt & 2) rr.emplace_back(new HM2D::Contour::R::ReallyRevert(right));
	else rd.emplace_back(new HM2D::Contour::R::ReallyDirect(right));

	if (opt & 4) rr.emplace_back(new HM2D::Contour::R::ReallyRevert(top));
	else rd.emplace_back(new HM2D::Contour::R::ReallyDirect(top));

	if (opt & 8) rr.emplace_back(new HM2D::Contour::R::ReallyRevert(bot));
	else rd.emplace_back(new HM2D::Contour::R::ReallyDirect(bot));

	HM2D::VertexData leftp=HM2D::Contour::OrderedPoints(left);
	HM2D::VertexData botp=HM2D::Contour::OrderedPoints(bot);
	HM2D::VertexData rightp=HM2D::Contour::OrderedPoints(right);
	HM2D::VertexData topp = HM2D::Contour::OrderedPoints(top);
	connect_vec_points(leftp, botp, rightp, topp);

	//calculate area and swap roles if needed
	HM2D::VertexData closed;
	closed.insert(closed.end(), leftp.rbegin()+1, leftp.rend());
	closed.insert(closed.end(), botp.begin()+1, botp.end());
	closed.insert(closed.end(), rightp.begin()+1, rightp.end());
	closed.insert(closed.end(), topp.rbegin()+1, topp.rend());
	double area = HM2D::Contour::Area(HM2D::Contour::Assembler::Contour1(closed, true));
	assert(fabs(area)>geps*geps);

	//assemble result
	vector<HM2D::EdgeData> ret(4);
	if (area > 0){
		HM2D::DeepCopy(left, ret[0], 0);
		HM2D::DeepCopy(bot, ret[1], 0);
		HM2D::DeepCopy(right, ret[2], 0);
		HM2D::DeepCopy(top, ret[3], 0);
	} else {
		HM2D::Contour::R::ReallyRevert r1(bot), r2(top);
		HM2D::DeepCopy(right, ret[0], 0);
		HM2D::DeepCopy(bot, ret[1], 0);
		HM2D::DeepCopy(left, ret[2], 0);
		HM2D::DeepCopy(top, ret[3], 0);
	}
	return ret;
}

vector<HM2D::EdgeData> connect_contours(HM2D::EdgeData& left, HM2D::EdgeData& bot,
		HM2D::EdgeData& right, HM2D::EdgeData& top){
	if (HM2D::Contour::IsClosed(left) || 
	    HM2D::Contour::IsClosed(bot) ||
	    HM2D::Contour::IsClosed(right) ||
	    HM2D::Contour::IsClosed(top))
		throw std::runtime_error("Closed contours are not allowed");
	if (has_self_cross(left) || has_self_cross(bot) || has_self_cross(right) || has_self_cross(top))
		throw std::runtime_error("Contour with a self cross is not allowed");
	std::map<double, int> allcases;
	int best_option = -1;
	for (int opt=0; opt<16; ++opt){
		double ans = tryopt(opt, left, bot, right, top);
		if (ISZERO(ans)) {best_option = opt; break;}
		if (ans > 0) allcases[ans] = opt;
	}
	if (best_option < 0 && allcases.size() > 0) best_option = allcases.begin()->second;
	if (best_option < 0) throw std::runtime_error("Failed to connect contours into quadrangle");

	return modify_contours(best_option, left, bot, right, top);
}

//check direction of first cell and reverse all cells if needed
void check_direction(HM2D::GridData& ret){
	if (ret.vcells.size() == 0) return;
	if (HM2D::Contour::Area(ret.vcells[0]->edges) < 0){
		for (auto e: ret.vedges) std::swap(e->left, e->right);
		for (auto c: ret.vcells) HM2D::Contour::Algos::Reverse(c->edges);
	}
	if (!HM2D::Grid::Algos::Check(ret)) throw HMMap::EInvalidGrid(std::move(ret));
}

void set_bt(const HM2D::EdgeData& from, const HM2D::EdgeData& to){
	for (int i=0; i<to.size(); ++i)
		to[i]->boundary_type = from[i]->boundary_type;
};

}

HM2D::GridData HMMap::LinearRectGrid(HM2D::EdgeData& _left, HM2D::EdgeData& _bot,
		HM2D::EdgeData& _right, HM2D::EdgeData& _top){
	if (_left.size() != _right.size() || _bot.size() != _top.size())
		throw std::runtime_error("right/top contours should have same number "
				"of nodes as left/bottom for linear rectangle grid algo");

	auto newcont = connect_contours(_left, _bot, _right, _top);
	auto& left = newcont[0]; auto pleft = HM2D::Contour::OrderedPoints(left);
	auto& bot = newcont[1]; auto pbot = HM2D::Contour::OrderedPoints(bot);
	auto& right = newcont[2]; auto pright = HM2D::Contour::OrderedPoints(right);
	auto& top = newcont[3]; auto ptop = HM2D::Contour::OrderedPoints(top);
	
	//build grid main grid connectivity
	HM2D::GridData ret = HM2D::Grid::Constructor::RectGrid01(bot.size(), left.size());

	//boundary types. use structured data given by RectGrid01 procedure.
	set_bt(bot, HM2D::Grid::Constructor::RectGridBottom(ret));
	set_bt(top, HM2D::Grid::Constructor::RectGridTop(ret));
	set_bt(left, HM2D::Grid::Constructor::RectGridLeft(ret));
	set_bt(right, HM2D::Grid::Constructor::RectGridRight(ret));

	//shift nodes
	auto allnodes = ret.vvert;
	int k=0;
	double ksieta[2];
	for (int j=0; j<left.size()+1; ++j){
		Point* lp = pleft[j].get();
		Point* rp = pright[j].get();
		for (int i=0; i<bot.size()+1; ++i){
			Point* tp = ptop[i].get();
			Point* bp = pbot[i].get();
			Point* p = allnodes[k++].get();
			if (i == 0) *p = *lp;
			else if (i == bot.size()) *p = *rp;
			else if (j == 0) *p = *bp;
			else if (j == left.size()) *p = *tp;
			else{
				bool iscr = SectCross(*lp, *rp, *bp, *tp, ksieta);
				if (!iscr) throw std::runtime_error("failed to find section cross");
				*p = Point::Weigh(*lp, *rp, ksieta[0]);
			}
		}
	}
	check_direction(ret);
	return ret;
}
namespace{
vector<Point*> topp(const HM2D::VertexData& data){
	vector<Point*> ret;
	for(auto& p: data) ret.push_back(p.get());
	return ret;
}
}

HM2D::GridData HMMap::TOrthogonalRectGrid::_run(HM2D::EdgeData& _left, HM2D::EdgeData& _bot,
		HM2D::EdgeData& _right, HM2D::EdgeData& _top){
	callback->step_after(5, "Connect contours");
	auto conres = connect_contours(_left, _bot, _right, _top);
	auto& left = conres[0];
	auto& bot = conres[1];
	auto& right = conres[2];
	auto& top = conres[3];

	//previous procedure can swap left/right contours.
	//We use input left as a source, so we need to know whether the swap occured.
	bool left_basic = (HM2D::Contour::First(left) == HM2D::Contour::First(_left) ||
			   HM2D::Contour::First(left) == HM2D::Contour::Last(_left));

	auto subcaller = callback->bottom_line_subrange(70);
	HMMap::Conformal::Options opt;
	opt.use_scpack = false;
	opt.use_rect_approx = false;
	opt.fem_nrec = std::max(std::min((size_t)10000, 2*left.size()*right.size()), (size_t)500);
	auto cmap = HMMap::Conformal::BuildRect.UseCallback(subcaller,
			left, right, bot, top, opt);
	int lsz = left_basic ? left.size() : right.size();
	HM2D::GridData ret = HM2D::Grid::Constructor::RectGrid01(bot.size(), lsz);

	//grid to conformal rectangle
	callback->step_after(15, "Map grid", 2, 1);
	vector<Point> xpts = cmap->MapToRectangle(topp(HM2D::Contour::OrderedPoints(bot)));
	vector<Point> ypts = left_basic ? cmap->MapToRectangle(topp(HM2D::Contour::OrderedPoints(left)))
	                                : cmap->MapToRectangle(topp(HM2D::Contour::OrderedPoints(right)));
	for (int i=0; i<ret.vvert.size(); ++i){
		int ni = i % (bot.size() + 1);
		int nj = i / (bot.size() + 1);
		ret.vvert[i]->set(xpts[ni].x, ypts[nj].y);
	}

	//grid to physical domain
	callback->subprocess_step_after(1);
	vector<Point> gpnt(ret.vvert.size());
	for (int i=0; i<ret.vvert.size(); ++i) gpnt[i].set(*ret.vvert[i]);
	vector<Point> physpnt = cmap->MapToPolygon(gpnt);
	for (int i=0; i<ret.vvert.size(); ++i){
		ret.vvert[i]->set(physpnt[i]);
	}

	//boundary types
	callback->step_after(5, "Assign boundary types");
	set_bt(bot, HM2D::Grid::Constructor::RectGridBottom(ret));
	auto to = HM2D::Grid::Constructor::RectGridTop(ret);
	HM2D::ECol::Algos::AssignBTypes(top, to);
	if (left_basic){
		set_bt(left, HM2D::Grid::Constructor::RectGridLeft(ret));
		to = HM2D::Grid::Constructor::RectGridRight(ret);
		HM2D::ECol::Algos::AssignBTypes(right, to);
	} else {
		set_bt(right, HM2D::Grid::Constructor::RectGridRight(ret));
		to = HM2D::Grid::Constructor::RectGridLeft(ret);
		HM2D::ECol::Algos::AssignBTypes(left, to);
	}

	callback->step_after(5, "Grid check");
	check_direction(ret);
	//move boundary nodes to their exact places
	{
		//bottom
		auto pp = HM2D::Contour::OrderedPoints(bot);
		for (int i=0; i<pp.size(); ++i){
			ret.vvert[i]->set(*pp[i]);
		}
	}
	{
		//side
		auto pp = HM2D::Contour::OrderedPoints( (left_basic) ? left : right );
		int i0 = (left_basic) ? 0 : bot.size();
		for (int i=0; i<pp.size(); ++i){
			ret.vvert[i0]->set(*pp[i]);
			i0 += bot.size()+1;
		}
	}
	return ret;
}

HM2D::GridData HMMap::FDMLaplasRectGrid(HM2D::EdgeData& _left, HM2D::EdgeData& _bot,
	HM2D::EdgeData& _right, HM2D::EdgeData& _top){
	if (_left.size() != _right.size() || _bot.size() != _top.size())
		throw std::runtime_error("right/top contours should have same number "
				"of nodes as left/bottom for laplas rectangle grid algo");

	//Cont4Connection cn = connect_rect_segments(left, bot, right, top);
	auto conres = connect_contours(_left, _bot, _right, _top);
	auto& left = conres[0];
	auto& bot = conres[1];
	auto& right = conres[2];
	auto& top = conres[3];

	//assembling fdm grid
	vector<double> x(top.size()+1, 0);
	vector<double> y(left.size()+1, 0);
	auto leftop = HM2D::Contour::OrderedPoints(left);
	auto rightop = HM2D::Contour::OrderedPoints(right);
	auto topop = HM2D::Contour::OrderedPoints(top);
	auto botop = HM2D::Contour::OrderedPoints(bot);
	for (int i=1; i<topop.size(); ++i){
		double s1 = Point::dist(*topop[i], *topop[i-1]);
		double s2 = Point::dist(*botop[i], *botop[i-1]);
		x[i] = x[i-1] + (s1 + s2)/2.0;
	}
	for (int i=1; i<leftop.size(); ++i){
		double s1 = Point::dist(*leftop[i], *leftop[i-1]);
		double s2 = Point::dist(*rightop[i], *rightop[i-1]);
		y[i] = y[i-1] + (s1 + s2)/2.0;
	}

	//fdm solver
	auto slv = HMFdm::LaplasSolver(x, y);
	//x problem
	vector<double> xcoords(x.size() * y.size(), 0);
	slv.SetBndValues(HMFdm::LaplasSolver::Bnd::Top,
			[&topop](int i, int j){ return topop[i]->x; });
	slv.SetBndValues(HMFdm::LaplasSolver::Bnd::Bottom,
			[&botop](int i, int j){ return botop[i]->x; });
	slv.SetBndValues(HMFdm::LaplasSolver::Bnd::Left,
			[&leftop](int i, int j){ return leftop[j]->x; });
	slv.SetBndValues(HMFdm::LaplasSolver::Bnd::Right,
			[&rightop](int i, int j){ return rightop[j]->x; });
	slv.Solve(xcoords);

	//y problem
	vector<double> ycoords(x.size() * y.size(), 0);
	slv.SetBndValues(HMFdm::LaplasSolver::Bnd::Top,
			[&topop](int i, int j){ return topop[i]->y; });
	slv.SetBndValues(HMFdm::LaplasSolver::Bnd::Bottom,
			[&botop](int i, int j){ return botop[i]->y; });
	slv.SetBndValues(HMFdm::LaplasSolver::Bnd::Left,
			[&leftop](int i, int j){ return leftop[j]->y; });
	slv.SetBndValues(HMFdm::LaplasSolver::Bnd::Right,
			[&rightop](int i, int j){ return rightop[j]->y; });
	slv.Solve(ycoords);

	//building resulting grid
	HM2D::GridData ret = HM2D::Grid::Constructor::RectGrid01(x.size()-1, y.size()-1);
	for (int i=0; i<ret.vvert.size(); ++i){
		ret.vvert[i]->set(xcoords[i], ycoords[i]);
	}
	set_bt(bot, HM2D::Grid::Constructor::RectGridBottom(ret));
	set_bt(top, HM2D::Grid::Constructor::RectGridTop(ret));
	set_bt(left, HM2D::Grid::Constructor::RectGridLeft(ret));
	set_bt(right, HM2D::Grid::Constructor::RectGridRight(ret));

	check_direction(ret);
	return ret;
}

HM2D::GridData HMMap::TLaplaceRectGrid::_run(HM2D::EdgeData& left, HM2D::EdgeData& bot,
		HM2D::EdgeData& right, HM2D::EdgeData& top, std::string algo){
	if (left.size() != right.size() || bot.size() != top.size())
		throw std::runtime_error("right/top contours should have same number "
				"of nodes as left/bottom for laplace rectangle grid algo");
	callback->step_after(5, "Connect contours");

	auto conres = connect_contours(left, bot, right, top);
	auto& left1 = conres[0];
	auto& bot1 = conres[1];
	auto& right1 = conres[2];
	auto& top1 = conres[3];

	//enclose
	top1[0]->vertices[0] = left1.back()->last();
	bot1[0]->vertices[0] = left1[0]->first();
	right1[0]->vertices[0] = bot1.back()->last();
	right1.back()->vertices[1] = top1.back()->last();

	//build base grid
	callback->step_after(5, "Assemble base grid");
	auto weight_contour = [](const HM2D::EdgeData& icont, bool vert, bool plus_one){
		vector<double> w = HM2D::Contour::EWeights(icont);
		vector<Point> pp; pp.reserve(w.size());
		for (auto it: w) pp.push_back(Point(it, 0));
		if (plus_one) for (auto& p: pp) p.y+=1;
		if (vert) for (auto& p: pp) std::swap(p.x, p.y);
		return HM2D::Contour::Constructor::FromPoints(pp);
	};
	auto rectleft = weight_contour(left1, true, false);
	auto rectright = weight_contour(right1, true, true);
	auto rectbot = weight_contour(bot1, false, false);
	auto recttop = weight_contour(top1, false, true);
	HM2D::GridData base = HMMap::LinearRectGrid(rectleft, rectbot, rectright, recttop);

	//boundary types for base
	set_bt(bot, HM2D::Grid::Constructor::RectGridBottom(base));
	set_bt(top, HM2D::Grid::Constructor::RectGridTop(base));
	set_bt(left, HM2D::Grid::Constructor::RectGridLeft(base));
	set_bt(right, HM2D::Grid::Constructor::RectGridRight(base));

	//map grid base and mapped points
	callback->step_after(10, "Boundary mapping");
	vector<Point> base_pnt;
	vector<Point> mapped_pnt;
	auto add_point_map = [&base_pnt, &mapped_pnt](
			const HM2D::EdgeData& bcont,
			const HM2D::EdgeData& mcont){
		auto v1 = HM2D::Contour::OrderedPoints(bcont), v2 = HM2D::Contour::OrderedPoints(mcont);
		for (int i=0; i<v1.size(); ++i){
			base_pnt.push_back(*v1[i]);
			mapped_pnt.push_back(*v2[i]);
		}
	};
	add_point_map(rectleft, left1);
	add_point_map(rectright, right1);
	add_point_map(rectbot, bot1);
	add_point_map(recttop, top1);

	//area building
	HM2D::EdgeData ecol;
	ecol.insert(ecol.end(), left1.begin(), left1.end());
	ecol.insert(ecol.end(), right1.begin(), right1.end());
	ecol.insert(ecol.end(), top1.begin(), top1.end());
	ecol.insert(ecol.end(), bot1.begin(), bot1.end());

	//options
	HMMap::Options opt(algo);
	opt.fem_nrec = 1.5*(left1.size()+1)*(bot1.size()+1);
	opt.btypes_from_contour = false;

	//calculate
	auto subcaller = callback->bottom_line_subrange(80);
	return HMMap::MapGrid.UseCallback(subcaller, base, ecol, base_pnt, mapped_pnt, false, opt);
}

HM2D::GridData HMMap::LinearTFIRectGrid(HM2D::EdgeData& left, HM2D::EdgeData& bot,
		HM2D::EdgeData& right, HM2D::EdgeData& top){
	if (left.size() != right.size() || bot.size() != top.size())
		throw std::runtime_error("right/top contours should have same number "
				"of nodes as left/bottom for tfi algo");

	auto conres = connect_contours(left, bot, right, top);
	auto& left1 = conres[0];
	auto& bot1 = conres[1];
	auto& right1 = conres[2];
	auto& top1 = conres[3];
	HM2D::VertexData leftp = HM2D::Contour::OrderedPoints(left1);
	HM2D::VertexData botp = HM2D::Contour::OrderedPoints(bot1);
	HM2D::VertexData rightp = HM2D::Contour::OrderedPoints(right1);
	HM2D::VertexData topp = HM2D::Contour::OrderedPoints(top1);
	Point p00 = *leftp[0];
	Point p10 = *rightp[0];
	Point p01 = *topp[0];
	Point p11 = *topp.back();
	vector<double> ksi = HM2D::Contour::EWeights(bot1);
	vector<double> ksi2 = HM2D::Contour::EWeights(top1);
	vector<double> eta = HM2D::Contour::EWeights(left1);
	vector<double> eta2 = HM2D::Contour::EWeights(right1);
	for (int i=0; i<ksi.size(); ++i) ksi[i] = (ksi[i] + ksi2[i])/2.0;
	for (int i=0; i<eta.size(); ++i) eta[i] = (eta[i] + eta2[i])/2.0;
	//for (int i=0; i<ksi.size(); ++i) ksi[i] = (double)i/(ksi.size()-1);
	//for (int i=0; i<eta.size(); ++i) eta[i] = (double)i/(eta.size()-1);

	HM2D::GridData U = HM2D::Grid::Constructor::RectGrid01(bot.size(), left.size());
	HM2D::GridData V = HM2D::Grid::Constructor::RectGrid01(bot.size(), left.size());
	HM2D::GridData UV = HM2D::Grid::Constructor::RectGrid01(bot.size(), left.size());

	for (int j=0; j<left.size()+1; ++j){
		//double eta = double(j)/(left.size());
		double e = eta[j];
		for (int i=0; i<bot.size()+1; ++i){
			//double ksi = double(i)/(bot.size());
			double k = ksi[i];
			int gi = i + j*(bot.size()+1);
			U.vvert[gi]->set(Point::Weigh(*leftp[j], *rightp[j], k));
			V.vvert[gi]->set(Point::Weigh(*botp[i], *topp[i], e));
			UV.vvert[gi]->set(
				p00*(1-k)*(1-e) + p11*k*e + p10*k*(1-e) + p01*(1-k)*e
			);
		}
	}
	for (int i=0; i<U.vvert.size(); ++i){
		HM2D::Vertex*  p1 = U.vvert[i].get();
		HM2D::Vertex*  p2 = V.vvert[i].get();
		HM2D::Vertex*  p3 = UV.vvert[i].get();
		p1->set(*p1 + *p2 - *p3);
	}

	//boundary types
	set_bt(bot, HM2D::Grid::Constructor::RectGridBottom(U));
	set_bt(top, HM2D::Grid::Constructor::RectGridTop(U));
	set_bt(left, HM2D::Grid::Constructor::RectGridLeft(U));
	set_bt(right, HM2D::Grid::Constructor::RectGridRight(U));

	check_direction(U);
	return U;
}

HM2D::GridData HMMap::CubicTFIRectGrid(HM2D::EdgeData& left, HM2D::EdgeData& bot,
		HM2D::EdgeData& right, HM2D::EdgeData& top, std::array<double, 4> c){
	if (left.size() != right.size() || bot.size() != top.size())
		throw std::runtime_error("right/top contours should have same number "
				"of nodes as left/bottom for tfi algo");

	auto conres = connect_contours(left, bot, right, top);
	auto& left1 = conres[0];
	auto& bot1 = conres[1];
	auto& right1 = conres[2];
	auto& top1 = conres[3];
	HM2D::VertexData leftp = HM2D::Contour::OrderedPoints(left1);
	HM2D::VertexData botp = HM2D::Contour::OrderedPoints(bot1);
	HM2D::VertexData rightp = HM2D::Contour::OrderedPoints(right1);
	HM2D::VertexData topp = HM2D::Contour::OrderedPoints(top1);
	vector<double> ksi = HM2D::Contour::EWeights(bot1);
	vector<double> ksi2 = HM2D::Contour::EWeights(top1);
	vector<double> eta = HM2D::Contour::EWeights(left1);
	vector<double> eta2 = HM2D::Contour::EWeights(right1);
	for (int i=0; i<ksi.size(); ++i) ksi[i] = (ksi[i] + ksi2[i])/2.0;
	for (int i=0; i<eta.size(); ++i) eta[i] = (eta[i] + eta2[i])/2.0;

	vector<Point> dksi_left, dksi_right, deta_bot, deta_top;
	vector<double> a01, a02, a11, a12, b01, b02, b11, b12;

	for (int i=0; i<bot.size()+1; ++i){
		//deta
		int i2, i1;
		if (i==0){
			i1 = 0; i2 = 1;
		} else if (i==bot.size()){
			i1 = botp.size()-2; i2 = botp.size()-1;
		} else {
			i1 = i-1; i2 = i+1;
		}
		Point dksi_bot = (*botp[i2] - *botp[i1])/(ksi[i2]-ksi[i1]);
		Point dksi_top = (*topp[i2] - *topp[i1])/(ksi[i2]-ksi[i1]);
		deta_bot.push_back(Point(-dksi_bot.y, dksi_bot.x)*c[1]);
		deta_top.push_back(Point(-dksi_top.y, dksi_top.x)*c[3]);
		//aij
		double k = ksi[i];
		a01.push_back(2*k*k*k - 3*k*k + 1);
		a02.push_back(-2*k*k*k + 3*k*k);
		a11.push_back(k*k*k - 2*k*k + k);
		a12.push_back(k*k*k - k*k);
	}
	for (int j=0; j<left.size()+1; ++j){
		//dksi
		int i2, i1;
		if (j==0){
			i1 = 0; i2 = 1;
		} else if (j==left.size()){
			i1 = leftp.size()-2; i2 = leftp.size()-1;
		} else {
			i1 = j-1; i2 = j+1;
		}
		Point deta_left = (*leftp[i2] - *leftp[i1])/(eta[i2]-eta[i1]);
		Point deta_right = (*rightp[i2] - *rightp[i1])/(eta[i2]-eta[i1]);
		dksi_left.push_back(Point(deta_left.y, -deta_left.x)*c[0]);
		dksi_right.push_back(Point(deta_right.y, -deta_right.x)*c[2]);
		//bij
		double e = eta[j];
		b01.push_back(2*e*e*e - 3*e*e + 1);
		b02.push_back(-2*e*e*e + 3*e*e);
		b11.push_back(e*e*e - 2*e*e + e);
		b12.push_back(e*e*e - e*e);
	}
	Point detaksi_00 = (deta_bot[1]-deta_bot[0])/(ksi[1]-ksi[0]);
	Point detaksi_01 = (deta_top[1]-deta_top[0])/(ksi[1]-ksi[0]);
	Point detaksi_10 = (deta_bot.back() - deta_bot[deta_bot.size()-2])/(ksi.back()-ksi[ksi.size()-2]);
	Point detaksi_11 = (deta_top.back() - deta_top[deta_top.size()-2])/(ksi.back()-ksi[ksi.size()-2]);

	HM2D::GridData U = HM2D::Grid::Constructor::RectGrid01(bot.size(), left.size());
	HM2D::GridData V = HM2D::Grid::Constructor::RectGrid01(bot.size(), left.size());
	HM2D::GridData UV = HM2D::Grid::Constructor::RectGrid01(bot.size(), left.size());

	for (int j=0; j<left.size()+1; ++j){
		for (int i=0; i<bot.size()+1; ++i){
			int gi = i + j*(bot.size()+1);
			U.vvert[gi]->set(
				*leftp[j]*a01[i] + *rightp[j]*a02[i]+
				dksi_left[j]*a11[i] + dksi_right[j]*a12[i]);
			V.vvert[gi]->set(
				*botp[i]*b01[j] + *topp[i]*b02[j]+
				deta_bot[i]*b11[j] + deta_top[i]*b12[j]);
		}
	}
	for (int j=0; j<left.size()+1; ++j){
		for (int i=0; i<bot.size()+1; ++i){
			int gi = i + j*(bot.size()+1);
			UV.vvert[gi]->set(
				*V.vvert[j*(bot.size()+1)]*a01[i] +
				*V.vvert[bot.size() + j*(bot.size()+1)]*a02[i]+
				(dksi_left[0]*b01[j] + dksi_left.back()*b02[j] + 
					detaksi_00*b11[j]+detaksi_01*b12[j])*a11[i]+
				(dksi_right[0]*b01[j] + dksi_right.back()*b02[j] +
					detaksi_10*b11[j]+detaksi_11*b12[j])*a12[i]
			);
		}
	}

	for (int i=0; i<U.vvert.size(); ++i){
		auto& p1 = U.vvert[i];
		auto& p2 = V.vvert[i];
		auto& p3 = UV.vvert[i];
		p1->set(*p1 + *p2 - *p3);
	}

	//boundary types
	set_bt(bot, HM2D::Grid::Constructor::RectGridBottom(U));
	set_bt(top, HM2D::Grid::Constructor::RectGridTop(U));
	set_bt(left, HM2D::Grid::Constructor::RectGridLeft(U));
	set_bt(right, HM2D::Grid::Constructor::RectGridRight(U));

	check_direction(U);
	return U;
}
