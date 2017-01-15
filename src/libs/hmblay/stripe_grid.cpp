#include "contour.hpp"
#include "hmblay.hpp"
#include "constructor.hpp"
#include "buildgrid.hpp"
#include "healgrid.hpp"
#include "modgrid.hpp"

HM2D::GridData HalfCirc(HM2D::EdgeData& b1, HM2D::EdgeData& b2, double sz){
	Point b1f = *HM2D::Contour::First(b1);
	Point b1l = *HM2D::Contour::Last(b1);
	Point b2f = *HM2D::Contour::First(b2);
	Point b2l = *HM2D::Contour::Last(b2);
	//calculate size
	bool has_center = b1f == b2f;
	double rad = Point::dist(b1l, b2l)/2.0;
	int narc = 2*M_PI*rad/sz;
	if (narc<8) narc = 8;
	if (narc % 2 != 0) narc+=1;
	sz = 2*M_PI*rad/narc; 
	
	int nrad = b1.size();
	//if (!has_center) ++nrad;
	Point pc;
	if (has_center) pc = b1f;
	else pc = (b1f + b2f)/2.0;

	//build circular grid
	HM2D::GridData ret = HM2D::Grid::Constructor::Circle(Point(0, 0), rad, narc, nrad, true);

	//change radii
	std::vector<double> rads;
	auto pts1 = HM2D::Contour::OrderedPoints(b1);
	auto pts2 = HM2D::Contour::OrderedPoints(b2);
	if (!has_center) pts1[0]->set(pc);
	if (!has_center) pts2[0]->set(pc);
	for (int i=1; i<pts1.size(); ++i){
		rads.push_back(Point::dist(pc, *pts1[i]));
	}
	for (int j=0; j<nrad; ++j){
		for (int i=0; i<narc; ++i){
			int k = (nrad-1-j)*narc + i;
			vecSetLen(*ret.vvert[k], rads[j]);
		}
	}

	//cut
	for (int j=0; j<nrad; ++j){
		for (int i=narc/2; i<narc; ++i){
			int k = (nrad-1-j)*narc+i;
			ret.vcells[k] = nullptr;
		}
	}
	HM2D::Grid::Algos::RestoreFromCells(ret);
	
	//translate and rotate
	double angle = atan2(b1l.y - b1f.y, b1l.x - b1f.x);
	for (auto& v: ret.vvert){
		v->set(vecRotate(*v, angle));
		*v += pc;
	}
	//snap to input boundaries
	for (int j=0; j<nrad; ++j){
		int k1 = (nrad-1-j)*(narc/2+1) + 0;
		int k2 = (nrad-1-j)*(narc/2+1) + narc/2;
		ret.vvert[k1]->set(*pts1[j+1]);
		ret.vvert[k2]->set(*pts2[j+1]);
	}

	//connect to parent
	if (!has_center){
		//move
		Vect mv = b1l - b1f;
		mv = vecRotate(mv, M_PI/2.0);
		vecSetLen(mv, sz);
		for (auto v: ret.vvert) *v += mv;

		//additional section
		ShpVector<HM2D::Vertex> apoints, ap1, ap2;
		for (auto p: pts1) aa::add_shared(ap1, HM2D::Vertex(*p));
		for (auto p: pts2) aa::add_shared(ap2, HM2D::Vertex(*p));
		ap1[0]->set(b1f); ap2[0]->set(b2f);
		DeepCopy(ret.vvert, apoints);
		std::copy(ap1.begin(), ap1.end(), std::back_inserter(apoints));
		std::copy(ap2.begin(), ap2.end(), std::back_inserter(apoints));
		aa::enumerate_ids_pvec(apoints);
		vector<vector<int>> cell_vert;
		for (int j=0; j<nrad; ++j){
			int k1 = (nrad-j)*(narc/2+1) + 0;
			int k2 = (nrad-j-1)*(narc/2+1) + 0;
			cell_vert.emplace_back();
			cell_vert.back().push_back(ap1[j]->id);
			cell_vert.back().push_back(ap1[j+1]->id);
			cell_vert.back().push_back(apoints[k2]->id);
			cell_vert.back().push_back(apoints[k1]->id);
		}
		for (int j=0; j<nrad; ++j){
			int k1 = (nrad-j)*(narc/2+1);
			int k2 = (nrad-j-1)*(narc/2+1) + narc/2;
			if (j!=0) k1 += narc/2;
			cell_vert.emplace_back();
			cell_vert.back().push_back(ap2[j+1]->id);
			cell_vert.back().push_back(ap2[j]->id);
			cell_vert.back().push_back(apoints[k1]->id);
			cell_vert.back().push_back(apoints[k2]->id);
		}
		//center triangle
		cell_vert.emplace_back();
		cell_vert.back().push_back(apoints[nrad*(narc/2+1)]->id);
		cell_vert.back().push_back(ap2[0]->id);
		cell_vert.back().push_back(ap1[0]->id);
		HM2D::GridData ag = HM2D::Grid::Constructor::FromTab(apoints, cell_vert);
		HM2D::Grid::Algos::MergeTo(ag, ret);
	}
	return ret;
}

HMCallback::FunctionWithCallback<HMBlay::TBuildStripeGrid> HMBlay::BuildStripeGrid;

HM2D::GridData HMBlay::TBuildStripeGrid::_run(const HM2D::EdgeData& cont,
		const std::vector<double>& partition,
		int tip_algo, Point& bl, Point& br, Point& tr, Point& tl){
	HM2D::EdgeData cp(cont);
	HMBlay::Input opt;
	opt.partition = partition;
	if (opt.partition[0] > 0) opt.partition.insert(opt.partition.begin(), 0);
	opt.bnd_step_method = HMBlay::BndStepMethod::NO_BND_STEPPING;
	opt.edges = &cp;
	opt.start = *HM2D::Contour::First(cont);
	opt.end = *HM2D::Contour::Last(cont);

	callback->step_after(45, "Upper grid");
	opt.direction = HMBlay::Direction::INNER;
	HM2D::GridData g1 = HMBlay::BuildBLayerGrid({opt});
	HM2D::EdgeData cleft1, cright1;
	{
		vector<Point> _pv1, _pv2;
		if (HM2D::Contour::IsOpen(cont)){
			int Nplen = cont.size()+1;
			for (int i=0; i<opt.partition.size(); ++i){
				_pv1.push_back(*g1.vvert[i*Nplen]);
				_pv2.push_back(*g1.vvert[i*Nplen+Nplen-1]);
			}
			cleft1 = HM2D::Contour::Constructor::FromPoints(_pv1);
			cright1 = HM2D::Contour::Constructor::FromPoints(_pv2);
			tl = *HM2D::Contour::Last(cleft1);
			tr = *HM2D::Contour::Last(cright1);
		} else { tl = *g1.vvert.back(); tr = tl; }
	}

	callback->step_after(45, "Lower grid");
	opt.direction = HMBlay::Direction::OUTER;
	HM2D::GridData g2 = HMBlay::BuildBLayerGrid({opt});
	HM2D::EdgeData cleft2, cright2;
	{
		vector<Point> _pv1, _pv2;
		if (HM2D::Contour::IsOpen(cont)){
			int Nplen = cont.size()+1;
			for (int i=0; i<opt.partition.size(); ++i){
				_pv2.push_back(*g2.vvert[i*Nplen]);
				_pv1.push_back(*g2.vvert[i*Nplen+Nplen-1]);
			}
			cleft2 = HM2D::Contour::Constructor::FromPoints(_pv1);
			cright2 = HM2D::Contour::Constructor::FromPoints(_pv2);
			bl = *HM2D::Contour::Last(cleft2);
			br = *HM2D::Contour::Last(cright2);
		} else { bl = *g2.vvert.back(); br=bl; }
	}

	//if no center line
	if (partition[0] > 0){
		int N = cont.size(); if (HM2D::Contour::IsOpen(cont)) ++N;

		vector<int> lowpts1, lowpts2, lowcls2;
		for (int i=0; i<N; ++i) lowpts1.push_back(i);
		for (int i=0; i<N; ++i) lowpts2.push_back(i+N);
		if (HM2D::Contour::IsOpen(cont)) std::reverse(lowpts2.begin(), lowpts2.end());
		else std::reverse(lowpts2.begin()+1, lowpts2.end());
		for (int i=0; i<cont.size(); ++i) lowcls2.push_back(i);

		for (int i=0; i<N; ++i){
			g1.vvert[lowpts1[i]]->set(*g2.vvert[lowpts2[i]]);
		}

		for (auto i: lowcls2) g2.vcells[i] = nullptr;
		HM2D::Grid::Algos::RestoreFromCells(g2);

		if (HM2D::Contour::IsOpen(cont)){
			cleft1.erase(cleft1.begin());
			cleft2.erase(cleft2.begin());
			cright1.erase(cright1.begin());
			cright2.erase(cright2.begin());
		}
		
	}

	callback->step_after(5, "Tip grids");
	shared_ptr<HM2D::GridData> tip1, tip2;

	if (HM2D::Contour::IsClosed(cont)){ 
	} else if (tip_algo == 0){
	} else if (tip_algo == 1){
		if (cleft1.size()>0){
			tip1.reset(new HM2D::GridData(HalfCirc(cleft1, cleft2, cont[0]->length())));
			tip2.reset(new HM2D::GridData(HalfCirc(cright2, cright1, cont.back()->length())));
		} else {
			HM2D::VertexData pts;
			Vect av1 = tl - *HM2D::Contour::First(cont);
			vecSetLen(av1, (vecLen(av1) + cont[0]->length())/2.0);
			av1 = vecRotate(av1, M_PI/2.);
			aa::add_shared(pts, HM2D::Vertex(tl));
			aa::add_shared(pts, HM2D::Vertex( *HM2D::Contour::First(cont) + av1  ));
			aa::add_shared(pts, HM2D::Vertex(bl));
			tip1 = std::make_shared<HM2D::GridData>(
				HM2D::Grid::Constructor::FromTab(pts, {{0, 1, 2}}));

			pts.clear();
			Vect av2 = br - *HM2D::Contour::Last(cont);
			vecSetLen(av2, (vecLen(av2) + cont.back()->length())/2.0);
			av2 = vecRotate(av2, M_PI/2.);
			aa::add_shared(pts, HM2D::Vertex(br));
			aa::add_shared(pts, HM2D::Vertex( *HM2D::Contour::Last(cont) + av2  ));
			aa::add_shared(pts, HM2D::Vertex(tr));
			tip2 = std::make_shared<HM2D::GridData>(
				HM2D::Grid::Constructor::FromTab(pts, {{0, 1, 2}}));
		}
	}

	callback->step_after(5, "Merge");
	HM2D::GridData ret = g1;
	HM2D::Grid::Algos::MergeTo(g2, ret);
	if (tip1) HM2D::Grid::Algos::MergeTo(*tip1, ret);
	if (tip2) HM2D::Grid::Algos::MergeTo(*tip2, ret);

	return ret;
}
