#include "rectangle_grid_builder.hpp"
#include "procgrid.h"
#include "hmmapping.hpp"
#include "hmfdm.hpp"

namespace{
struct Cont4Connection{
	Point *l1, *l2, *r1, *r2, *t1, *t2, *b1, *b2; //ordered corner points of each contour
	                                               //  l2 t1---t2 r2
	                                               //  |          |
	                                               //  |          |
	                                               //  l1 b1---b2 r1 
	Cont4Connection(HMCont2D::Contour& left, HMCont2D::Contour& bot,
			HMCont2D::Contour& right, HMCont2D::Contour& top){
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

//check direction of first cell and reverse all cells if needed
void check_direction(GridGeom& ret){
	if (ret.n_cells() == 0) return;
	if (ret.get_cell(0)->area() < 0){
		GGeom::Modify::CellModify(ret, [](Cell* c){
			std::reverse(c->points.begin(), c->points.end());} );
	}
}

}

GridGeom HMGMap::LinearRectGrid(HMCont2D::Contour& left, HMCont2D::Contour& bot,
		HMCont2D::Contour& right, HMCont2D::Contour& top){
	int nleft = left.size();
	int nbot = bot.size();
	if (nleft != right.size() || nbot != top.size())
		throw std::runtime_error("right/top contours should have same number "
				"of nodes as left/bottom for linear rectangle grid algo");

	Cont4Connection cn = connect_rect_segments(left, bot, right, top);
	//assemble ordered points for each contour and their weights
	auto pleft = left.ordered_points();
	auto pbot = bot.ordered_points();
	auto pright = right.ordered_points();
	auto ptop = top.ordered_points();
	//reverse if needed
	if (pleft[0] != cn.l1) std::reverse(pleft.begin(), pleft.end());
	if (pbot[0] != cn.b1) std::reverse(pbot.begin(), pbot.end());
	if (ptop[0] != cn.t1) std::reverse(ptop.begin(), ptop.end());
	if (pright[0] != cn.r1) std::reverse(pright.begin(), pright.end());
	//build grid main grid connectivity
	GridGeom ret = GGeom::Constructor::RectGrid01(nbot, nleft);
	//shift nodes
	auto allnodes = GGeom::Info::SharePoints(ret);
	int k=0;
	double ksieta[2];
	for (int j=0; j<nleft+1; ++j){
		Point* lp = pleft[j];
		Point* rp = pright[j];
		for (int i=0; i<nbot+1; ++i){
			Point* tp = ptop[i];
			Point* bp = pbot[i];
			Point* p = allnodes[k++].get();
			if (i == 0) *p = *lp;
			else if (i == nbot) *p = *rp;
			else if (j == 0) *p = *bp;
			else if (j == nleft) *p = *tp;
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

GridGeom HMGMap::ConformalRectGrid(HMCont2D::Contour& left, HMCont2D::Contour& bot,
		HMCont2D::Contour& right, HMCont2D::Contour& top){
	_THROW_NOT_IMP_;
}

GridGeom HMGMap::LaplasRectGrid(HMCont2D::Contour& left, HMCont2D::Contour& bot,
		HMCont2D::Contour& right, HMCont2D::Contour& top){
	Cont4Connection cn = connect_rect_segments(left, bot, right, top);

	// --- prepare for map grid invoke
	//weights and ordered points
	auto rev_if_needed = [](HMCont2D::Contour& cont, Point* p0, vector<Point*>& op, vector<double>& w){
		if (cont.first() != p0){
			std::reverse(op.begin(), op.end());
			for (auto& v: w) v = 1 - v;
			std::reverse(w.begin(), w.end());
		}
	};
	auto left_points = left.ordered_points(); auto left_w = HMCont2D::Contour::EWeights(left);
	auto right_points = right.ordered_points(); auto right_w = HMCont2D::Contour::EWeights(right);
	auto top_points = top.ordered_points(); auto top_w = HMCont2D::Contour::EWeights(top);
	auto bottom_points = bot.ordered_points(); auto bottom_w = HMCont2D::Contour::EWeights(bot);
	rev_if_needed(left, cn.l1, left_points, left_w);
	rev_if_needed(right, cn.r1, right_points, right_w);
	rev_if_needed(top, cn.t1, top_points, top_w);
	rev_if_needed(bot, cn.b1, bottom_points, bottom_w);

	// --- assemble ecollection from input contours
	HMCont2D::ECollection ecol1, ecol2, ecol3, ecol4;
	HMCont2D::ECollection::DeepCopy(left, ecol1);
	HMCont2D::ECollection::DeepCopy(bot, ecol2);
	HMCont2D::ECollection::DeepCopy(right, ecol3);
	HMCont2D::ECollection::DeepCopy(top, ecol4);

	// --- change b1->l1; t1->l2; r1->b2; r2->t2
	auto change_end_points = [](Point* from, Point* to, HMCont2D::ECollection& ecol){
		if (ecol.data[0]->pstart == from) ecol.data[0]->pstart = to;
		else if (ecol.data[0]->pend == from) ecol.data[0]->pend = to;
		else if (ecol.data.back()->pstart == from) ecol.data.back()->pend = to;
		else if (ecol.data.back()->pend == from) ecol.data.back()->pend = to;
		else assert(false);
	};
	change_end_points(cn.b1, cn.l1, ecol2);
	change_end_points(cn.t1, cn.l2, ecol4);
	change_end_points(cn.r1, cn.b2, ecol3);
	change_end_points(cn.r2, cn.t2, ecol3);

	// --- gather all collection in one
	HMCont2D::ECollection ecol;
	ecol.Unite(ecol1);
	ecol.Unite(ecol2);
	ecol.Unite(ecol3);
	ecol.Unite(ecol4);

	//basic rectangle grid grid
	int nleft = left.size();
	int nbot = bot.size();
	if (nleft != right.size() || nbot != top.size())
		throw std::runtime_error("right/top contours should have same number "
				"of nodes as left/bottom for laplas rectangle grid algo");
	//create rectangle of four lines with each lines partition equals input data
	vector<Point> left_aux, bot_aux, right_aux, top_aux;
	for (int i=0; i < nleft+1; ++i){
		left_aux.push_back(Point(0, left_w[i]));
		right_aux.push_back(Point(1, right_w[i]));
	}
	for (int i=0; i < nbot+1; ++i){
		bot_aux.push_back(Point(bottom_w[i], 0));
		top_aux.push_back(Point(top_w[i], 1));
	}
	auto cleft_aux = HMCont2D::Constructor::ContourFromPoints(left_aux);
	auto cbot_aux = HMCont2D::Constructor::ContourFromPoints(bot_aux);
	auto cright_aux = HMCont2D::Constructor::ContourFromPoints(right_aux);
	auto ctop_aux = HMCont2D::Constructor::ContourFromPoints(top_aux);
	GridGeom rg = LinearRectGrid(cleft_aux, cbot_aux, cright_aux, ctop_aux);
	
	//points mappints
	vector<Point> base_points, mapped_points;
	for (int i=0; i<left_points.size(); ++i){
		mapped_points.push_back(*left_points[i]);
		base_points.push_back(*rg.get_point((nbot+1) * i));
		mapped_points.push_back(*right_points[i]);
		base_points.push_back(*rg.get_point((nbot+1) * i + nbot));
	}
	for (int i=0; i<bottom_points.size(); ++i){
		mapped_points.push_back(*bottom_points[i]);
		base_points.push_back(*rg.get_point(i));
		mapped_points.push_back(*top_points[i]);
		base_points.push_back(*rg.get_point(nleft*(nbot+1) + i));
	}

	//calculate
	return HMGMap::MapGrid(rg, ecol, base_points, mapped_points);
}

GridGeom HMGMap::FDMLaplasRectGrid(HMCont2D::Contour& left, HMCont2D::Contour& bot,
	HMCont2D::Contour& right, HMCont2D::Contour& top){
	if (left.size() != right.size() || bot.size() != top.size())
		throw std::runtime_error("right/top contours should have same number "
				"of nodes as left/bottom for laplas rectangle grid algo");

	Cont4Connection cn = connect_rect_segments(left, bot, right, top);

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
