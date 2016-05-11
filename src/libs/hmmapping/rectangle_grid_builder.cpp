#include "rectangle_grid_builder.hpp"
#include "procgrid.h"
#include "hmmapping.hpp"
#include "hmfdm.hpp"
#include "hmconformal.hpp"
#include "debug_grid2d.h"

HMCallback::FunctionWithCallback<HMGMap::TOrthogonalRectGrid> HMGMap::OrthogonalRectGrid;
HMCallback::FunctionWithCallback<HMGMap::TLaplaceRectGrid> HMGMap::LaplaceRectGrid;

namespace{

HMCont2D::Container<HMCont2D::Contour> deepcopy(const HMCont2D::Contour& cont, bool isrev){
	vector<Point*> p1 = cont.ordered_points();
	vector<Point> p2(p1.size());
	for (int i=0; i<p1.size(); ++i) p2[i].set(p1[i]->x, p1[i]->y);
	if (isrev) std::reverse(p2.begin(), p2.end());
	return HMCont2D::Constructor::ContourFromPoints(p2);
}

bool has_self_cross(const HMCont2D::Contour& cont){
	double ksieta[2];
	for (int i=0; i<cont.size()-2; ++i){
		for (int j=i+2; j<cont.size(); ++j){
			Point& p1 = *cont.data[i]->pstart;
			Point& p2 = *cont.data[i]->pend;
			Point& p3 = *cont.data[j]->pstart;
			Point& p4 = *cont.data[j]->pend;
			SectCross(p1, p2, p3, p4, ksieta);
			if (ksieta[0]>-geps && ksieta[0]<1+geps &&
				ksieta[1]>-geps && ksieta[1]<1+geps) return true;
		}
	}
	return false;
}
bool check_for_no_cross(const HMCont2D::Contour& left, const HMCont2D::Contour& bot,
		const HMCont2D::Contour& right, const HMCont2D::Contour& top){
	if (HMCont2D::Algos::CrossAll(left, bot).size() != 1) return false;
	if (HMCont2D::Algos::CrossAll(left, top).size() != 1) return false;
	if (std::get<0>(HMCont2D::Algos::Cross(left, right))) return false;
	if (std::get<0>(HMCont2D::Algos::Cross(bot, top))) return false;
	if (HMCont2D::Algos::CrossAll(bot, right).size() != 1) return false;
	if (HMCont2D::Algos::CrossAll(top, right).size() != 1) return false;
	if (has_self_cross(right)) return false;  //only right was modified pointwisely
	return true;
}

double connect_vec_points(const vector<Point*>& left, const vector<Point*>& bot,
		const vector<Point*>& right, const vector<Point*>& top){
	Point bot_move=*left[0] - *bot[0];
	for (auto p: bot) *p += bot_move;
	Point top_move=*left.back() - *top[0];
	for (auto p: top) *p += top_move;
	Point right_move1 = *bot.back() - *right[0];
	Point right_move2 = *top.back() - *right.back();
	auto rw = HMCont2D::Contour::EWeights(HMCont2D::Constructor::ContourFromPoints(right));
	for (int i=0; i<rw.size(); ++i){
		*right[i] += (right_move1 * (1.0-rw[i]) + right_move2 * rw[i]);
	}
	return vecLen(bot_move) + vecLen(top_move) +
		vecLen(right_move1) + vecLen(right_move2);
}

double tryopt(int opt, const HMCont2D::Contour& _left, const HMCont2D::Contour& _bot,
		const HMCont2D::Contour& _right, const HMCont2D::Contour& _top){
	bool left_rev = opt & 1;
	bool right_rev = opt & 2;
	bool top_rev = opt & 4;
	bool bot_rev = opt & 8;

	auto left = deepcopy(_left, left_rev);
	auto right = deepcopy(_right, right_rev);
	auto top = deepcopy(_top, top_rev);
	auto bot = deepcopy(_bot, bot_rev);

	double ret = connect_vec_points(left.ordered_points(), bot.ordered_points(),
			right.ordered_points(), top.ordered_points());

	if (check_for_no_cross(left, bot, right, top)){
		return ret;
	} else {
		return -1;
	}
}

vector<HMCont2D::Contour> modify_contours(int opt, HMCont2D::Contour& left, HMCont2D::Contour& bot,
		HMCont2D::Contour& right, HMCont2D::Contour& top){
	vector<Point*> leftp=left.ordered_points();
	vector<Point*> botp=bot.ordered_points();
	vector<Point*> rightp=right.ordered_points();
	vector<Point*> topp = top.ordered_points();
	if (opt & 1) std::reverse(leftp.begin(), leftp.end());
	if (opt & 2) std::reverse(rightp.begin(), rightp.end());
	if (opt & 4) std::reverse(topp.begin(), topp.end());
	if (opt & 8) std::reverse(botp.begin(), botp.end());
	connect_vec_points(leftp, botp, rightp, topp);
	//calculate area and swap roles if needed
	vector<Point*> closed;
	closed.insert(closed.end(), leftp.rbegin()+1, leftp.rend());
	closed.insert(closed.end(), botp.begin()+1, botp.end());
	closed.insert(closed.end(), rightp.begin()+1, rightp.end());
	closed.insert(closed.end(), topp.rbegin()+1, topp.rend());
	double area = HMCont2D::Area(HMCont2D::Constructor::ContourFromPoints(closed, true));
	assert(fabs(area)>geps*geps);

	//assemble result
	vector<HMCont2D::Contour> ret;
	if (area > 0){
		ret.push_back(HMCont2D::Constructor::ContourFromPoints(leftp));
		ret.push_back(HMCont2D::Constructor::ContourFromPoints(botp));
		ret.push_back(HMCont2D::Constructor::ContourFromPoints(rightp));
		ret.push_back(HMCont2D::Constructor::ContourFromPoints(topp));
	} else {
		std::reverse(botp.begin(), botp.end());
		std::reverse(topp.begin(), topp.end());
		ret.push_back(HMCont2D::Constructor::ContourFromPoints(rightp));
		ret.push_back(HMCont2D::Constructor::ContourFromPoints(botp));
		ret.push_back(HMCont2D::Constructor::ContourFromPoints(leftp));
		ret.push_back(HMCont2D::Constructor::ContourFromPoints(topp));
	}
	return ret;
}

vector<HMCont2D::Contour> connect_contours(HMCont2D::Contour& left, HMCont2D::Contour& bot,
		HMCont2D::Contour& right, HMCont2D::Contour& top){
	if (left.is_closed() || bot.is_closed() || right.is_closed() || top.is_closed())
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


/*
struct Cont4Connection{
	Point *l1, *l2, *r1, *r2, *t1, *t2, *b1, *b2; //ordered corner points of each contour
	                                               //  l2 t1---t2 r2
	                                               //  |          |
	                                               //  |          |
	                                               //  l1 b1---b2 r1
	const HMCont2D::Contour *pleft, *pbot, *pright, *ptop;

	Cont4Connection(HMCont2D::Contour& left, HMCont2D::Contour& bot,
			HMCont2D::Contour& right, HMCont2D::Contour& top){
		pleft = &left; pbot = &bot;
		pright = &right; ptop = &top;
		//check for open contours
		if (left.is_closed() || bot.is_closed() || top.is_closed() || right.is_closed()){
			throw std::runtime_error("contours should be open");
		}
		//guess of top left, bottom left
		double m11 = Point::meas(*left.first(), *bot.first());
		double m12 = Point::meas(*left.first(), *bot.last());
		double m21 = Point::meas(*left.last(), *bot.first());
		double m22 = Point::meas(*left.last(), *bot.last());
		if (m11 <= m12 && m11 <= m21 && m11 <= m22){
			l1 = left.first(); l2 = left.last();
			b1 = bot.first(); b2 = bot.last();
		} else if (m12 <= m11 && m12 <= m21 && m12 <= m22){
			l1 = left.first(); l2 = left.last();
			b1 = bot.last(); b2 = bot.first();
		} else if (m21 <= m11 && m21 <= m12 && m21 <= m22){
			l1 = left.last(); l2 = left.first();
			b1 = bot.first(); b2 = bot.last();
		} else {
			l1 = left.last(); l2 = left.first();
			b1 = bot.last(); b2 = bot.first();
		}
		// top contour
		//if top contour is simply a copy of bottom
		if (ISZERO(Point::dist(*bot.first(), *top.first())) &&
				ISZERO(Point::dist(*bot.last(), *top.last()))){
			t1 = top.first(); t2 = top.last();
			if (bot.first() != b1) std::swap(t1, t2);
		} else {
			double m1 = Point::meas(*top.first(), *l2);
			double m2 = Point::meas(*top.last(), *l2);
			t1 = top.first(); t2 = top.last();
			if (m1 > m2) std::swap(t1, t2);
		}
		//right contour
		if (ISZERO(Point::dist(*left.first(), *right.first())) &&
				ISZERO(Point::dist(*left.last(), *right.last()))){
			r1 = right.first(); r2 = right.last();
			if (left.first() != l1) std::swap(r1, r2);
		} else {
			double m1 = Point::meas(*right.first(), *b2);
			double m2 = Point::meas(*right.last(), *b2);
			r1 = right.first(); r2 = right.last();
			if (m1 > m2) std::swap(r1, r2);
		}
	}
	HMCont2D::Contour directed_left() const{
		HMCont2D::Contour ret;
		HMCont2D::Contour::DeepCopy(*pleft, ret);
		if (ret.first() != l1){ ret.ReallyReverse(); }
		assert(ret.first() == l1);
		return ret;
	}
	HMCont2D::Contour directed_top() const{
		HMCont2D::Contour ret;
		HMCont2D::Contour::DeepCopy(*ptop, ret);
		if (ret.first() != t1){ ret.ReallyReverse(); }
		assert(ret.first() == t1);
		return ret;
	}
	HMCont2D::Contour directed_bot() const{
		HMCont2D::Contour ret;
		HMCont2D::Contour::DeepCopy(*pbot, ret);
		if (ret.first() != b1){ ret.ReallyReverse(); }
		assert(ret.first() == b1);
		return ret;
	}
	HMCont2D::Contour directed_right() const{
		HMCont2D::Contour ret;
		HMCont2D::Contour::DeepCopy(*pright, ret);
		if (ret.first() != r1){ ret.ReallyReverse(); }
		assert(ret.first() == r1);
		return ret;
	}
};

Cont4Connection connect_rect_segments(HMCont2D::Contour& left, HMCont2D::Contour& bot,
		HMCont2D::Contour& right, HMCont2D::Contour& top){
	Cont4Connection cn(left, bot, right, top);
	//move bottom contour
	Point mv1 = *cn.l1 - *cn.b1;
	for (auto p: bot.all_points()) *p += mv1;
	//move top
	Point mv2 = *cn.l2 - *cn.t1;
	for (auto p: top.all_points()) *p += mv2;
	//move right
	Point mv3 = *cn.b2 - *cn.r1;
	for (auto p: right.all_points()) *p += mv3;
	//stretch right contour to r2 = t2
	if (!ISZERO(Point::dist(*cn.t2, *cn.r2))){
		Vect mv = *cn.t2 - *cn.r2;
		bool isrev = right.first() != cn.r1;
		auto ew = HMCont2D::Contour::EWeights(right);
		if (isrev) for (auto& v: ew) v = 1 - v;
		int i=0;
		for (auto p: right.ordered_points()){
			Vect mvl = mv * ew[i++];
			*p += mvl;
		}
	}
	return cn;
}
*/

//check direction of first cell and reverse all cells if needed
void check_direction(GridGeom& ret){
	if (ret.n_cells() == 0) return;
	if (ret.get_cell(0)->area() < 0){
		GGeom::Modify::CellModify(ret, [](Cell* c){
			std::reverse(c->points.begin(), c->points.end());} );
	}
	if (!GGeom::Info::Check(ret)) throw HMGMap::MapException("Resulting grid is not valid");
}

}

GridGeom HMGMap::LinearRectGrid(HMCont2D::Contour& _left, HMCont2D::Contour& _bot,
		HMCont2D::Contour& _right, HMCont2D::Contour& _top){
	if (_left.size() != _right.size() || _bot.size() != _top.size())
		throw std::runtime_error("right/top contours should have same number "
				"of nodes as left/bottom for linear rectangle grid algo");

	//Cont4Connection cn = connect_rect_segments(left, bot, right, top);
	////assemble ordered points for each contour and their weights
	//auto pleft = left.ordered_points();
	//auto pbot = bot.ordered_points();
	//auto pright = right.ordered_points();
	//auto ptop = top.ordered_points();
	////reverse if needed
	//if (pleft[0] != cn.l1) std::reverse(pleft.begin(), pleft.end());
	//if (pbot[0] != cn.b1) std::reverse(pbot.begin(), pbot.end());
	//if (ptop[0] != cn.t1) std::reverse(ptop.begin(), ptop.end());
	//if (pright[0] != cn.r1) std::reverse(pright.begin(), pright.end());
	auto newcont = connect_contours(_left, _bot, _right, _top);
	auto& left = newcont[0]; auto pleft = left.ordered_points();
	auto& bot = newcont[1]; auto pbot = bot.ordered_points();
	auto& right = newcont[2]; auto pright = right.ordered_points();
	auto& top = newcont[3]; auto ptop = top.ordered_points();

	//build grid main grid connectivity
	GridGeom ret = GGeom::Constructor::RectGrid01(bot.size(), left.size());
	//shift nodes
	auto allnodes = GGeom::Info::SharePoints(ret);
	int k=0;
	double ksieta[2];
	for (int j=0; j<left.size()+1; ++j){
		Point* lp = pleft[j];
		Point* rp = pright[j];
		for (int i=0; i<bot.size()+1; ++i){
			Point* tp = ptop[i];
			Point* bp = pbot[i];
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

GridGeom HMGMap::TOrthogonalRectGrid::_run(HMCont2D::Contour& _left, HMCont2D::Contour& _bot,
		HMCont2D::Contour& _right, HMCont2D::Contour& _top){
	callback->step_after(5, "Connect contours");
	auto conres = connect_contours(_left, _bot, _right, _top);
	auto& left = conres[0];
	auto& bot = conres[1];
	auto& right = conres[2];
	auto& top = conres[3];

	//previous procedure can swap left/right contours.
	//We use input left as a source, so we need to know whether the swap occured.
	bool left_basic = (left.first() == _left.first() ||
			   left.first() == _left.last());

	auto subcaller = callback->bottom_line_subrange(80);
	HMMath::Conformal::Options opt;
	opt.use_scpack = false;
	opt.use_rect_approx = false;
	opt.fem_nrec = std::max(std::min(10000, left.size()*right.size()*2), 500);
	auto cmap = HMMath::Conformal::BuildRect.UseCallback(subcaller,
			left, right, bot, top, opt);
	GridGeom ret = GGeom::Constructor::RectGrid01(bot.size(), left.size());

	//grid to conformal rectangle
	callback->step_after(10, "Map grid", 2, 1);
	vector<Point> xpts = cmap->MapToRectangle(bot.ordered_points());
	vector<Point> ypts = left_basic ? cmap->MapToRectangle(left.ordered_points())
	                                : cmap->MapToRectangle(right.ordered_points());
	auto to_conf = [&xpts, &ypts, &bot](GridPoint* p){
		int i = p->get_ind() % (bot.size() + 1);
		int j = p->get_ind() / (bot.size() + 1);
		p->set(xpts[i].x, ypts[j].y);
	};
	GGeom::Modify::PointModify(ret, to_conf);

	//grid to physical domain
	callback->subprocess_step_after(1);
	vector<Point> gpnt(ret.n_points());
	for (int i=0; i<ret.n_points(); ++i) gpnt[i] = *ret.get_point(i);
	vector<Point> physpnt = cmap->MapToPolygon(gpnt);
	GGeom::Modify::PointModify(ret, [&physpnt](GridPoint* p){
			Point& pp = physpnt[p->get_ind()];
			p->set(pp.x, pp.y);
	});

	callback->step_after(5, "Grid check");
	check_direction(ret);
	return ret;
}

GridGeom HMGMap::FDMLaplasRectGrid(HMCont2D::Contour& _left, HMCont2D::Contour& _bot,
	HMCont2D::Contour& _right, HMCont2D::Contour& _top){
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
	auto leftop = left.ordered_points();
	auto rightop = right.ordered_points();
	auto topop = top.ordered_points();
	auto botop = bot.ordered_points();
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
	GridGeom ret = GGeom::Constructor::RectGrid01(x.size()-1, y.size()-1);
	GGeom::Modify::PointModify(ret, [&xcoords, &ycoords](GridPoint* p){
				p->x = xcoords[p->get_ind()];
				p->y = ycoords[p->get_ind()];
			});

	check_direction(ret);
	return ret;
}

GridGeom HMGMap::TLaplaceRectGrid::_run(HMCont2D::Contour& left, HMCont2D::Contour& bot,
		HMCont2D::Contour& right, HMCont2D::Contour& top, std::string algo){
	if (left.size() != right.size() || bot.size() != top.size())
		throw std::runtime_error("right/top contours should have same number "
				"of nodes as left/bottom for laplace rectangle grid algo");
	callback->step_after(5, "Connect contours");

	//Cont4Connection cn = connect_rect_segments(left, bot, right, top);
	//auto left1 = cn.directed_left();
	//auto top1 = cn.directed_top();
	//auto bot1 = cn.directed_bot();
	//auto right1 = cn.directed_right();

	auto conres = connect_contours(left, bot, right, top);
	auto& left1 = conres[0];
	auto& bot1 = conres[1];
	auto& right1 = conres[2];
	auto& top1 = conres[3];

	//enclose
	top1.data[0]->pstart = left1.data.back()->pend;
	bot1.data[0]->pstart = left1.data[0]->pstart;
	right1.data[0]->pstart = bot1.data.back()->pend;
	right1.data.back()->pend = top1.data.back()->pend;

	//build base grid
	callback->step_after(5, "Assemble base grid");
	auto weight_contour = [](const HMCont2D::Contour& icont, bool vert, bool plus_one){
		vector<double> w = HMCont2D::Contour::EWeights(icont);
		vector<Point> pp; pp.reserve(w.size());
		for (auto it: w) pp.push_back(Point(it, 0));
		if (plus_one) for (auto& p: pp) p.y+=1;
		if (vert) for (auto& p: pp) std::swap(p.x, p.y);
		return HMCont2D::Constructor::ContourFromPoints(pp);
	};
	auto rectleft = weight_contour(left1, true, false);
	auto rectright = weight_contour(right1, true, true);
	auto rectbot = weight_contour(bot1, false, false);
	auto recttop = weight_contour(top1, false, true);
	GridGeom base = HMGMap::LinearRectGrid(rectleft, rectbot, rectright, recttop);

	//map grid base and mapped points
	callback->step_after(10, "Boundary mapping");
	vector<Point> base_pnt;
	vector<Point> mapped_pnt;
	auto add_point_map = [&base_pnt, &mapped_pnt](
			const HMCont2D::Contour& bcont,
			const HMCont2D::Contour& mcont){
		auto v1 = bcont.ordered_points(), v2 = mcont.ordered_points();
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
	HMCont2D::ECollection ecol;
	ecol.Unite(left1);
	ecol.Unite(right1);
	ecol.Unite(top1);
	ecol.Unite(bot1);

	//options
	HMGMap::Options opt(algo);
	opt.fem_nrec = 1.5*(left1.size()+1)*(bot1.size()+1);

	//calculate
	auto subcaller = callback->bottom_line_subrange(80);
	return HMGMap::MapGrid.UseCallback(subcaller, base, ecol, base_pnt, mapped_pnt, opt);
}
