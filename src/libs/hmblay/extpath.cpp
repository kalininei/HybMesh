#include "extpath.hpp"

using namespace HMBlay::Impl;

//smoothing takes place with length = HCOEF * (boundary depth)
const double HCOEF = 2.1;

void PathPntData::fill(Point* p1, Point* p2, Point* p3){
	//building of object for 'ext_data' field
	if (p1 == 0 || p3 == 0){
		tp = CornerTp::NO;
	} else {
		double angle = Angle(*p1, *p2, *p3);
		tp = opt->CornerType(angle);
	}
	p = p2;
};

double ExtPath::largest_depth(){
	double k = 0;
	for (auto& e: ext_data){
		if (e.opt->partition.back() > k)
			k = e.opt->partition.back();
	}
	return k;
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
		Point* p0 = p[p.size()-3];
		Point* p1 = p[p.size()-2];
		Point* p2 = p[0];
		Point* p3 = p[1];
		ret.ext_data[0].fill(p1,p2,p3);
		ret.ext_data.back().fill(p0,p1,p2);
	}
	//add one fictive entry to match ordered_points() length
	if (ret.is_closed()){
		ret.ext_data.push_back(ret.ext_data[0]);
	} else {
		ret.ext_data.push_back(ret.ext_data.back());
		ret.ext_data.back().fill(p[p.size()-2], p.back(), 0);
	}
		
	ret.FillEndConditions();
	return ret;
}

void ExtPath::FillEndConditions(){
	assert(size()>0);
	full_source = ext_data[0].opt->get_full_source();
	leftbc.clear();
	rightbc.clear();
	//if contour is closed draw an infinite perpendicular as end conditions
	if (is_closed()){ PerpendicularStart(); PerpendicularEnd(); return;}

	// ============ start point
	double h1 = ext_data[0].opt->partition.back();
	Point* p1 = first();
	if (!full_source->is_closed() && (p1 == full_source->first() || p1 == full_source->last()))
		PerpendicularStart();
	else{
		Vect direct1 = Contour::SmoothedDirection(*full_source, p1,  1, HCOEF*h1);
		Vect direct2 = Contour::SmoothedDirection(*full_source, p1, -1, HCOEF*h1);
		double angle = Angle(direct2, Point(0,0), direct1);
		CornerTp tp = ext_data[0].opt->CornerType(angle);
		switch (tp){
			case CornerTp::ZERO: case CornerTp::NO: case CornerTp::SHARP:
			case CornerTp::OBTUSE: case CornerTp::NEGLECTABLE:
				PerpendicularStart();
				break;
			case CornerTp::REGULAR:
				PerpendicularStart(angle/2.0);
				break;
			case CornerTp::CORNER:
				leftbc = HMCont2D::Constructor::CutContour(*full_source, *p1, -1, HCOEF*h1);
				//if leftbc contains only single edge its direction is not defined
				//we have to guarantee direction of bnd.
				if (*p1 != *leftbc.first()) leftbc.ReallyReverse();
				break;
		}
	}

	// ============= end point
	double h2 = ext_data.back().opt->partition.back();
	Point* p2 = last();
	if (!full_source->is_closed() && (p2 == full_source->first() || p2 == full_source->last()))
		PerpendicularEnd();
	else{
		Vect direct1 = Contour::SmoothedDirection(*full_source, p2,  1, HCOEF*h2);
		Vect direct2 = Contour::SmoothedDirection(*full_source, p2, -1, HCOEF*h2);
		double angle = Angle(direct2, Point(0,0), direct1);
		CornerTp tp = ext_data.back().opt->CornerType(angle);
		switch (tp){
			case CornerTp::ZERO: case CornerTp::NO: case CornerTp::SHARP:
			case CornerTp::OBTUSE: case CornerTp::NEGLECTABLE:
				PerpendicularEnd();
				break;
			case CornerTp::REGULAR:
				PerpendicularEnd(angle/2.0);
				break;
			case CornerTp::CORNER:
				rightbc = HMCont2D::Constructor::CutContour(*full_source, *p2, 1, HCOEF*h2);
				if (*rightbc.first() != *p2) rightbc.ReallyReverse();
				break;
		}
	};
}

void ExtPath::PerpendicularStart(double a){
	double h1 = ext_data[0].opt->partition.back();
	Vect v = HMCont2D::Contour::SmoothedDirection(*full_source, first(), 1, HCOEF*h1);
	//this is infinite since we use normalized geometry
	v = vecRotate(v*100.0, a);
	auto c = HMCont2D::Constructor::ContourFromPoints({*first(), *first()+v});
	leftbc.clear();
	leftbc.Unite(c);
}

void ExtPath::PerpendicularEnd(double a){
	double h1 = ext_data.back().opt->partition.back();
	Vect v = HMCont2D::Contour::SmoothedDirection(*full_source, last(), -1, HCOEF*h1);
	//this is infinite since we use normalized geometry
	v = vecRotate(v*100.0, 2*M_PI-a);
	auto c = HMCont2D::Constructor::ContourFromPoints({*last(), *last()+v});
	rightbc.clear();
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
	// === wright data
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
		if (pth.ext_data[i].tp == CornerTp::OBTUSE ||
				pth.ext_data[i].tp == CornerTp::SHARP  ||
				pth.ext_data[i].tp == CornerTp::CORNER){
			
			double h = pth.ext_data[i].opt->partition.back();
			Vect direct1 = SmoothedDirection(pth, op[i],  1, HCOEF*h);
			Vect direct2 = SmoothedDirection(pth, op[i], -1, HCOEF*h);
			double angle = Angle(direct2, Point(0, 0), direct1);
			CornerTp realtp = pth.ext_data[i].opt->CornerType(angle);
			if (realtp == CornerTp::ZERO || realtp == CornerTp::REGULAR)
				pth.ext_data[i].tp = CornerTp::NEGLECTABLE;
			else if (realtp != pth.ext_data[i].tp)
				//using sharp here as the most straightforward method
				pth.ext_data[i].tp = CornerTp::SHARP;
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
	for (auto& p: conts2){
		vector<double> lens = HMCont2D::ECollection::ELengths(p);
		std::copy(lens.begin(), lens.end(), std::back_inserter(ret));
	}
	std::partial_sum(ret.begin(), ret.end(), ret.begin());

	assert(ISEQ(ret.back(), len2));
	return ret;
}

HMCont2D::Contour ExtPath::Partition(HMCont2D::PCollection& apoints) const{
	auto& opt = *ext_data[0].opt;
	switch (opt.bnd_step_method){
		case BndStepMethod::NO_BND_STEPPING:
			return *this;
		case BndStepMethod::CONST_BND_STEP:
			return HMCont2D::Contour::Partition(opt.bnd_step, *this,
					apoints, HMCont2D::PartitionTp::IGNORE_ALL);
		case BndStepMethod::CONST_BND_STEP_KEEP_SHAPE:
			return HMCont2D::Contour::Partition(opt.bnd_step, *this,
					apoints, HMCont2D::PartitionTp::KEEP_SHAPE);
		case BndStepMethod::CONST_BND_STEP_KEEP_ALL:
			return HMCont2D::Contour::Partition(opt.bnd_step, *this,
					apoints, HMCont2D::PartitionTp::KEEP_ALL);
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
	HMCont2D::Contour c = HMCont2D::Contour::Assemble(*this, p1, p2);
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




