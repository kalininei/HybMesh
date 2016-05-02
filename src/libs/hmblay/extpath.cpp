#include "extpath.hpp"

using namespace HMBlay::Impl;


void PathPntData::fill(Point* p1, Point* p2, Point* p3){
	p = p2;
	set_exact_angle(0, true);
};

//smoothing takes place with length = HCOEF * (boundary depth)
const double HCOEF = 2.1;
void PathPntData::set_smooth_angle(int dir, bool redefine_type){
	double h = HCOEF*opt->partition.back();
	if (dir == 0){
		Vect direct1 = HMCont2D::Algos::SmoothedDirection2(*opt->get_full_source(), p,  1, h, 0);
		Vect direct2 = HMCont2D::Algos::SmoothedDirection2(*opt->get_full_source(), p, -1, 0, h);
		if (direct1.x == 0 && direct1.y == 0) direct1 = direct2 * -1;
		if (direct2.x == 0 && direct2.y == 0) direct2 = direct1 * -1;
		angle = Angle(direct2, Point(0,0), direct1);
		if (direct1.x != 0 || direct1.y != 0){
			normal = vecRotate(direct1, angle/2.0);
		} else {
			normal = vecRotate(direct2, -angle/2.0);
		}
	} else if (dir == 1){
		Vect direct = HMCont2D::Algos::SmoothedDirection2(*opt->get_full_source(), p,  1, h, 0);
		angle = M_PI;
		normal = vecRotate(direct, angle/2.0);
	} else if (dir == -1){
		Vect direct = HMCont2D::Algos::SmoothedDirection2(*opt->get_full_source(), p,  1, 0, h);
		angle = M_PI;
		normal = vecRotate(direct, angle/2.0);
	}
	vecNormalize(normal);

	if (redefine_type){
		//angles which were non-straight but become straight
		//are neglectable
		bool negl = tp != CornerTp::STRAIGHT;
		tp = opt->CornerType(angle);
		if (tp == CornerTp::STRAIGHT && negl) tp = CornerTp::NEGLECTABLE;
	}
}

void PathPntData::set_exact_angle(int dir, bool redefine_type){
	Point *p1, *p2, *p3;
	p2 = p;
	p1 = (dir != 1) ? opt->get_full_source()->prev_point(p2) : 0;
	p3 = (dir != -1) ? opt->get_full_source()->next_point(p2) : 0;

	Vect v;
	if (p1 == 0 || p3 == 0){
		assert( p1 != 0 || p3 != 0);
		angle = M_PI;
		if (p1 == 0) v = *p3 - *p2;
		else v = *p2 - *p1;
	} else {
		angle = Angle(*p1, *p2, *p3);
		v = *p3 - *p2;
	}
	normal = vecRotate(v, angle/2.0);
	vecNormalize(normal);

	if (redefine_type){
		tp = opt->CornerType(angle);
	}
}

double ExtPath::largest_depth() const{
	vector<double> h; h.reserve(ext_data.size());
	for (auto& e: ext_data) h.push_back(e.opt->partition.back());
	if (!is_closed()) h.resize(h.size() - 1);
	return *std::max_element(h.begin(), h.end());
}

int ExtPath::largest_vpart_size() const{
	vector<int> h(ext_data.size());
	for (int i=0; i<ext_data.size(); ++i) h[i] = ext_data[i].opt->partition.size();
	if (!is_closed()) h.resize(h.size() - 1);
	return *std::max_element(h.begin(), h.end());
}

ExtPath ExtPath::Assemble(const vector<Options*>& data){
	//data[i]->path are in strict sequence.
	//all data may have different options
	ExtPath ret;
	vector<Options*> edgeopt;
	for (auto d: data){
		auto p = d->get_path();
		//check for ordering
		assert(ret.size() == 0 || ret.last() == p->first());
		//add edges
		ret.Unite(*p);
		for (int i=0; i<p->size(); ++i) edgeopt.push_back(d);
	}
	//add empty options
	for (int i=0; i<edgeopt.size(); ++i){
		ret.ext_data.push_back(PathPntData(edgeopt[i]));
	}
	//fill options
	auto p = ret.ordered_points();
	for (int i=0; i<edgeopt.size(); ++i){
		Point* pprev = (i==0) ? 0 : p[i-1];
		Point* pnext = p[i+1];
		ret.ext_data[i].fill(pprev, p[i], pnext);
	}
	//for closed contours take into account first-last connection
	if (ret.is_closed()){
		Point* p0 = p[p.size()-2];
		Point* p1 = p[p.size()-1];
		Point* p2 = p[0];
		Point* p3 = p[1];
		ret.ext_data[0].fill(p0,p2,p3);
		//ret.ext_data.back().fill(p0,p1,p);
	}
	//add one entry to match ordered_points() length
	if (ret.is_closed()){
		ret.ext_data.push_back(ret.ext_data[0]);
	} else {
		ret.ext_data.push_back(ret.ext_data.back());
		Point *pprev = ret.ext_data[0].opt->get_full_source()->point_siblings(p[0])[0];
		Point *pnext = ret.ext_data.back().opt->get_full_source()->point_siblings(p.back())[2];
		ret.ext_data[0].fill(pprev, p[0], p[1]);
		ret.ext_data.back().fill(p[p.size()-2], p.back(), pnext);
	}
		
	ret.FillEndConditions();
	return ret;
}

void ExtPath::FillEndConditions(){
	assert(size()>0);
	full_source = ext_data[0].opt->get_full_source();
	leftbc.clear();
	rightbc.clear();
	// ====== start point angle adjustment
	switch (ext_data[0].tp){
		case CornerTp::ACUTE:
			ext_data[0].set_smooth_angle(1, false);
			break;
		case CornerTp::REENTRANT: case CornerTp::ROUND:
			//Calculating exact angle since smoothing
			//gives bad result on a circle like curves.
			//This is a temporary solution and may give
			//awful results in some cases. Better algo should be done here.
			ext_data[0].set_exact_angle(1, false);
			break;
		default:
			break;  // no need to adjust
	}
	// ===== end point angle adjustment
	switch (ext_data.back().tp){
		case CornerTp::ACUTE:
			ext_data.back().set_smooth_angle(-1, false);
			break;
		case CornerTp::REENTRANT: case CornerTp::ROUND:
			ext_data.back().set_exact_angle(-1, false);
			break;
		default:
			break;  // no need to adjust
	}

	//building start point boundary
	if (ext_data[0].tp == CornerTp::RIGHT){
		double h = 2*ext_data[0].opt->partition.back();
		leftbc = HMCont2D::Constructor::CutContour(*full_source, *first(), -1, h);
		//if leftbc contains only single edge its direction is not defined
		//we have to guarantee direction of bnd.
		if (*first() != *leftbc.first()) leftbc.ReallyReverse();
	} else {
		PerpendicularStart();
	}
	//end point boundary
	if (ext_data.back().tp == CornerTp::RIGHT){
		double h = 2*ext_data.back().opt->partition.back();
		rightbc = HMCont2D::Constructor::CutContour(*full_source, *last(), 1, h);
		if (*rightbc.first() != *last()) rightbc.ReallyReverse();
	} else {
		PerpendicularEnd();
	}
}

void ExtPath::PerpendicularStart(){
	leftbc.clear();
	Vect v = ext_data[0].normal * 100;
	auto c = HMCont2D::Constructor::ContourFromPoints({*first(), *first()+v});
	leftbc.Unite(c);
}

void ExtPath::PerpendicularEnd(){
	rightbc.clear();
	Vect v = ext_data.back().normal * 100;
	auto c = HMCont2D::Constructor::ContourFromPoints({*last(), *last()+v});
	rightbc.Unite(c);
}


vector<ExtPath> ExtPath::DivideByAngle(const ExtPath& pth, CornerTp tp){
	return DivideByAngle(pth, vector<CornerTp>({tp}));
}

vector<ExtPath> ExtPath::DivideByAngle(const ExtPath& pth, vector<CornerTp> tps){
	vector<ExtPath> ret;
	if (pth.size()==0) return ret;
	if (pth.size()==1) return {pth};
	auto info = pth.ordered_info();
	ret.push_back(ExtPath());
	auto istps = [&](CornerTp t){
		return std::find(tps.begin(), tps.end(), t)
				!= tps.end();
	};
	for (int i=0; i<info.size()-1; ++i){
		auto& it=info[i];
		if (i>0 && istps(pth.ext_data[it.index].tp)){
			ret.back().ext_data.push_back(pth.ext_data[it.index]);
			ret.push_back(ExtPath());
		}
		ret.back().add_value(it.enext);
		ret.back().ext_data.push_back(pth.ext_data[it.index]);
	}
	ret.back().ext_data.push_back(pth.ext_data[info.back().index]);
	//if this is a closed contour get rid of division at the first point
	if (pth.is_closed() && ret.size()>1 && !istps(ret[0].ext_data[0].tp)){
		ExtPath np(ret.back());
		np.Unite(ret[0]);
		for (int i=1; i<ret[0].ext_data.size(); ++i){
			np.ext_data.push_back(ret[0].ext_data[i]);
		}
		vector<ExtPath> nret { np };
		for (int i=1; i<ret.size()-1; ++i) nret.push_back(ret[i]);
		std::swap(ret, nret);
	}

	for (auto& r: ret){
		r.FillEndConditions();
		assert(r.ext_data.size() == r.ordered_points().size());
	}
	return ret;
}

vector<ExtPath> ExtPath::DivideByHalf(const ExtPath& pth){
	// === find half point
	assert(pth.size() > 2);
	Point* hp = HMCont2D::ECollection::FindClosestNode(
			pth, HMCont2D::Contour::WeightPoint(pth, 0.5));
	vector<Point*> dp = pth.ordered_points();
	int hind = std::find(dp.begin(), dp.end(), hp) - dp.begin();
	if (hind == 0) hind = 1;
	if (hind == dp.size()-1) hind = dp.size()-2;
	// === write data
	vector<ExtPath> ret(2);
	//first edges
	ret[0].add_value(HMCont2D::Edge(dp[0], dp[1]));
	for (int i=2; i<hind+1; ++i) ret[0].AddLastPoint(dp[i]);
	//second edges
	ret[1].add_value(HMCont2D::Edge(dp[hind], dp[hind+1]));
	for (int i = hind+2; i<dp.size(); ++i) ret[1].AddLastPoint(dp[i]);
	//first ext_data
	for (int i=0; i<hind+1; ++i) ret[0].ext_data.push_back(pth.ext_data[i]);
	//second ext_data
	for (int i=hind; i<dp.size(); ++i) ret[1].ext_data.push_back(pth.ext_data[i]);
	// === fill ends and return
	ret[0].FillEndConditions();
	ret[1].FillEndConditions();
	return ret;
}

void ExtPath::ReinterpretCornerTp(ExtPath& pth){
	vector<Point*> op = pth.ordered_points();
	//calculate smoothed angles and place neglactable feature
	for (int i=0; i<op.size(); ++i){
		auto ct = pth.ext_data[i].tp;
		if ((!pth.is_closed() && (i==0 || i==op.size()-1)) ||
				ct == CornerTp::REENTRANT ||
				ct == CornerTp::NO ||
				ct == CornerTp::ROUND ||
				ct == CornerTp::ACUTE  ||
				ct == CornerTp::RIGHT){
			pth.ext_data[i].set_smooth_angle(0, true);
		}
	}
}

vector<double> ExtPath::PathPartition(double len1, double len2) const{
	HMCont2D::PCollection apoints;

	//1) create a sub-contour [len1->len2]
	ExtPath subpath = SubPath(len1, len2, apoints);

	//2) create sub-sub-contours whith same partition options
	vector<ExtPath> subs = DivideByBndPart(subpath);

	//3) make a partition of each sub-sub-contour
	vector<Contour> conts2; conts2.reserve(subs.size());
	for (auto& p: subs){
		conts2.push_back(p.Partition(apoints));
	}

	//4) get lengths of each subcontour segments and gather result
	vector<double> ret {len1};
	for (auto& c: conts2){
		bool is0 = true;
		for (auto& p: c.ordered_points()){
			if (is0) {is0=false; continue; }
			auto coord = coord_at(*p);
			ret.push_back(std::get<0>(coord));
		}
	}
	return ret;
}

HMCont2D::Contour ExtPath::Partition(HMCont2D::PCollection& apoints) const{
	auto& opt = *ext_data[0].opt;
	switch (opt.bnd_step_method){
		case BndStepMethod::NO_BND_STEPPING:
			return *this;
		case BndStepMethod::CONST_BND_STEP:
			return HMCont2D::Algos::Partition(opt.bnd_step, *this,
					apoints, HMCont2D::PartitionTp::IGNORE_ALL);
		case BndStepMethod::CONST_BND_STEP_KEEP_SHAPE:
			return HMCont2D::Algos::Partition(opt.bnd_step, *this,
					apoints, HMCont2D::PartitionTp::KEEP_SHAPE);
		case BndStepMethod::CONST_BND_STEP_KEEP_ALL:
			return HMCont2D::Algos::Partition(opt.bnd_step, *this,
					apoints, HMCont2D::PartitionTp::KEEP_ALL);
		case BndStepMethod::INCREMENTAL:{
			double wstart = std::get<1>(opt.get_full_source()->coord_at(*first()));
			double wend = std::get<1>(opt.get_full_source()->coord_at(*last()));
			double w1 = std::get<1>(opt.get_full_source()->coord_at(opt.bnd_step_basis[0].first));
			double w2 = std::get<1>(opt.get_full_source()->coord_at(opt.bnd_step_basis[1].first));
			assert(!ISZERO(wstart-wend) && !ISZERO(w1-w2));
			if (opt.get_full_source()->is_closed()){
				if (wstart>wend) wend+=1.0;
				while (w1-wstart>geps) w1-=1.0;   //to get: w1<=wstart
				while (wend-w2>geps) w2+=1.0;         //to get: w2>=wend
			}
			double b1 = (wstart-w1)*opt.bnd_step_basis[1].second+(w2-wstart)*opt.bnd_step_basis[0].second;
			double b2 = (wend-w1)*opt.bnd_step_basis[1].second+(w2-wend)*opt.bnd_step_basis[0].second;
			b1/=(w2-w1); b2/=(w2-w1);
			std::map<double, double> wbas;
			wbas[0] = b1; wbas[1] = b2;
			return HMCont2D::Algos::WeightedPartition(wbas, *this,
					apoints, HMCont2D::PartitionTp::IGNORE_ALL);
		}
		default:
			throw std::runtime_error("Unsupported partition algorithm");
	}
}

vector<ExtPath> ExtPath::DivideByBndPart(const ExtPath& pth){
	vector<Point*> ap = pth.ordered_points();
	vector<Point*> rpoints {ap[0]};
	auto last_method = pth.ext_data[0].opt->bnd_step_method;
	auto last_step = pth.ext_data[0].opt->bnd_step;
	for (int i=0; i<ap.size(); ++i){
		auto meth = pth.ext_data[i].opt->bnd_step_method;
		auto step = pth.ext_data[i].opt->bnd_step;
		if (( meth != last_method) ||
		    (last_method != BndStepMethod::NO_BND_STEPPING && !ISEQ(step, last_step))){
			rpoints.push_back(ap[i]);
			last_method = meth;
			last_step = step;
		}
	}
	if (rpoints.size() == 1) { return vector<ExtPath> {pth}; }
	if (rpoints.back() != ap.back()) rpoints.push_back(ap.back());
	vector<ExtPath> subs;
	for (int i=0; i<rpoints.size()-1; ++i){
		subs.push_back(pth.SubPath(rpoints[i], rpoints[i+1]));
	}
	return subs;
}

ExtPath ExtPath::SubPath(const Point* p1, const Point* p2) const{
	HMCont2D::Contour c = HMCont2D::Assembler::Contour1(*this, p1, p2);
	ExtPath ret;
	ret.Unite(c);
	auto orig_pts = ordered_points();
	for (auto p: c.ordered_points()){
		int i = std::find(orig_pts.begin(), orig_pts.end(), p) - orig_pts.begin();
		ret.ext_data.push_back(ext_data[i]);
	}
	ret.FillEndConditions();
	return ret;
}

ExtPath ExtPath::SubPath(double len1, double len2, HMCont2D::PCollection& apoints) const{
	ExtPath ret(*this);
	Point* p1 = ret.AddPoint(len1, apoints);
	Point* p2 = ret.AddPoint(len2, apoints);
	return ret.SubPath(p1, p2);
}

Point* ExtPath::AddPoint(double len, HMCont2D::PCollection& apoints){
	auto w = WeightPointsByLen(*this, {len});
	auto g = GuaranteePoint(*w.point(0), apoints);
	Point* ret = std::get<1>(g);
	//if point already in path -> return it
	if (std::get<0>(g) == false) return ret;
	//if not -> add new entry to ext_data vector
	vector<Point*> op = this->ordered_points();
	auto ifnd = std::find(op.begin(), op.end(), ret) - op.begin();
	if (ifnd!=0) ext_data.insert(ext_data.begin() + ifnd, ext_data[ifnd-1]);
	else ext_data.insert(ext_data.begin(), ext_data[0]);
	return ret;
}

vector<double> ExtPath::VerticalPartition(double len) const{
	int i=0;
	double used_len = 0;
	while (ISGREATER(len, used_len) && i<size()-1) {used_len += edge(i)->length(); ++i;}
	if (ISGREATER(used_len, len)) --i;
	Options* ret;
	//if point lies exactly on edge division point choose longest partition
	if (ISEQ(used_len, len)){
		Options *cand1 = 0, *cand2 = 0;
		cand1 = ext_data[i].opt;
		if (i == 0){
			if (is_closed()){
				cand2 = ext_data[size()-1].opt;
			} else {
				cand2 = 0;
			}
		} else{
			cand2 = ext_data[i-1].opt;
		}
		if (cand2 == 0) ret = cand1;
		else{
			double L1 = cand1->partition.back();
			double L2 = cand2->partition.back();
			if (L1>L2) ret = cand1;
			else ret = cand2;
		}
	} else ret = ext_data[i].opt;

	return ret->partition;
}




