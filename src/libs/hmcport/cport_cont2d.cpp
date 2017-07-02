#include "cport_cont2d.h"
#include "c2cpp_helper.hpp"
#include "tscaler.hpp"
#include "primitives2d.hpp"
#include "contour_tree.hpp"
#include "buildcont.hpp"
#include "modcont.hpp"
#include "clipdomain.hpp"
#include "treverter2d.hpp"
#include "partcont.hpp"
#include "assemble2d.hpp"
#include "finder2d.hpp"
#include "export2d_hm.hpp"
#include "partition01.hpp"


int c2_dims(void* obj, int* ret){
	try{
		auto c = static_cast<HM2D::EdgeData*>(obj);
		ret[1] = c->size();
		ret[0] = HM2D::AllVertices(*c).size();
		return HMSUCCESS;
	} catch (std::exception& e){
		add_error_message(e.what());
		return HMERROR;
	}
}

int c2_deepcopy(void* obj, void** ret){
	try{
		HM2D::EdgeData* ret_ = new HM2D::EdgeData();
		HM2D::DeepCopy(*static_cast<HM2D::EdgeData*>(obj), *ret_);
		*ret = ret_;
		return HMSUCCESS;
	} catch (std::exception& e){
		add_error_message(e.what());
		return HMERROR;
	}
}

int c2_free(void* obj){
	try{
		if (obj != nullptr) delete static_cast<HM2D::EdgeData*>(obj);
		return HMSUCCESS;
	} catch (std::exception& e){
		add_error_message(e.what());
		return HMERROR;
	}
}

int c2_area(void* obj, double* ret){
	try{
		auto c = static_cast<HM2D::EdgeData*>(obj);
		auto tree = HM2D::Contour::Tree::Assemble(*c);
		*ret = tree.area();
		return HMSUCCESS;
	} catch (std::exception& e){
		add_error_message(e.what());
		return HMERROR;
	}
}

int c2_length(void* obj, double* ret){
	try{
		*ret = HM2D::Contour::Length(*static_cast<HM2D::EdgeData*>(obj));
		return HMSUCCESS;
	} catch (std::exception& e){
		add_error_message(e.what());
		return HMERROR;
	}
}

int c2_move(void* obj, double* dp){
	try{
		auto cont = static_cast<HM2D::EdgeData*>(obj);
		for (auto p: HM2D::AllVertices(*cont)){
			p->x += dp[0];
			p->y += dp[1];
		}
		return HMSUCCESS;
	} catch (std::exception& e){
		add_error_message(e.what());
		return HMERROR;
	}
}

int c2_scale(void* obj, double* pc, double* p0){
	try{
		auto cont = static_cast<HM2D::EdgeData*>(obj);
		double xc = pc[0] / 100., yc = pc[1] / 100.;
		for (auto& p: HM2D::AllVertices(*cont)){
			p->x -= p0[0]; p->x *= xc;
			p->y -= p0[1]; p->y *= yc;
		}
		return HMSUCCESS;
	} catch (std::exception& e){
		add_error_message(e.what());
		return HMERROR;
	}
}

int c2_reflect(void* obj, double* v0, double* v1){
	try{
		auto cont = static_cast<HM2D::EdgeData*>(obj);
		double lx = v1[0] - v0[0], ly = v1[1] - v0[1];
		double r2 = sqrt(lx*lx + ly*ly);
		lx /= r2; ly /= r2;
		double A[4] = {lx*lx - ly*ly, 2*lx*ly,
		               2*lx*ly, ly*ly-lx*lx};
		for (auto& p: HM2D::AllVertices(*cont)){
			double a = p->x - v0[0];
			double b = p->y - v0[1];
			p->x = A[0]*a + A[1]*b + v0[0];
			p->y = A[2]*a + A[3]*b + v0[1];
		}
		for (auto& e: *cont){ std::swap(e->left, e->right); }
		return HMSUCCESS;
	} catch (std::exception& e){
		add_error_message(e.what());
		return HMERROR;
	}
}

int c2_rotate(void* obj, double* p0, double a){
	try{
		auto cont = static_cast<HM2D::EdgeData*>(obj);
		double cs = cos(a / 180. * M_PI), sn = sin(a / 180. * M_PI);
		for (auto& p: HM2D::AllVertices(*cont)){
			double a = p->x - p0[0];
			double b = p->y - p0[1];
			p->x = a*cs - b*sn + p0[0];
			p->y = a*sn + b*cs + p0[1];
		}
		return HMSUCCESS;
	} catch (std::exception& e){
		add_error_message(e.what());
		return HMERROR;
	}
}

int c2_frompoints(int npts, double* pts, int* bnds, int force_closed, void** ret){
	try{
		vector<double> pp(pts, pts+npts*2);
		HM2D::EdgeData ret_ = HM2D::Contour::Constructor::FromPoints(pp, force_closed);
		for (int i=0; i<ret_.size(); ++i){
			ret_[i]->boundary_type = bnds[i];
		}
		c2cpp::to_pp(ret_, ret);
		return HMSUCCESS;
	} catch (std::exception& e){
		add_error_message(e.what());
		return HMERROR;
	}
}

int c2_simplify_self(void* obj, double angle, int keep_btypes){
	try{
		auto cont = static_cast<HM2D::EdgeData*>(obj);
		Autoscale::D2 sc(cont);
		auto ret1 = HM2D::ECol::Algos::Simplified(*cont, angle, keep_btypes);
		cont->clear();
		cont->insert(cont->end(), ret1.begin(), ret1.end());
		return HMSUCCESS;
	} catch (std::exception& e){
		add_error_message(e.what());
		return HMERROR;
	}
}

int c2_quick_separate(void* obj, int* nret, void*** ret){
	try{
		auto cont = static_cast<HM2D::EdgeData*>(obj);
		HM2D::EdgeData cp;
		HM2D::DeepCopy(*cont, cp);
		vector<HM2D::EdgeData> sep = HM2D::SplitData(cp);
		c2cpp::to_ppp(sep, nret, ret);
		return HMSUCCESS;
	} catch (std::exception& e){
		add_error_message(e.what());
		return HMERROR;
	}
}


int c2_unite(int nobjs, void** objs, void** ret){
	try{
		auto conts = c2cpp::to_pvec<HM2D::EdgeData>(nobjs, objs);
		Autoscale::D2 sc(conts);
		HM2D::EdgeData* ret_ = new HM2D::EdgeData();
		for (auto& c: conts){ HM2D::DeepCopy(*c, *ret_); }
		HM2D::ECol::Algos::MergePoints(*ret_);
		sc.unscale(ret_);
		*ret = ret_;
		return HMSUCCESS;
	} catch (std::exception& e){
		add_error_message(e.what());
		return HMERROR;
	}
}


int c2_decompose(void* obj, int* nret, void*** ret){
	try{
		auto cont = static_cast<HM2D::EdgeData*>(obj);
		Autoscale::D2 sc(cont);
		vector<HM2D::EdgeData> sep = HM2D::Contour::Constructor::ExtendedSeparate(*cont);
		vector<HM2D::EdgeData> ret_(sep.size());
		for (int i=0; i<sep.size(); ++i){
			HM2D::DeepCopy(sep[i], ret_[i]);
			sc.unscale(&ret_[i]);
		}
		c2cpp::to_ppp(ret_, nret, ret);
		return HMSUCCESS;
	} catch (std::exception& e){
		add_error_message(e.what());
		return HMERROR;
	}
}


int c2_spline(int np, double* _pts, int ned, int* bnds, void** ret){
	try{
		vector<Point> pts = c2cpp::to_points2(np, _pts);
		Autoscale::D2 sc(pts);
		HM2D::EdgeData spline = HM2D::Contour::Constructor::Spline(pts, ned);

		int icur = 0;
		for (int i=0; i<spline.size(); ++i){
			if (*spline[i]->pfirst() == pts[icur+1]){
				++icur;
			}
			spline[i]->boundary_type = bnds[icur];
		}

		sc.unscale(&spline);
		c2cpp::to_pp(spline, ret);
		return HMSUCCESS;
	} catch (std::exception& e){
		add_error_message(e.what());
		return HMERROR;
	}
}

int c2_clip_domain(void* obj1, void* obj2, const char* op, int simplify, void** ret){
	try{
		auto cont1 = static_cast<HM2D::EdgeData*>(obj1);
		auto cont2 = static_cast<HM2D::EdgeData*>(obj2);
		Autoscale::D2 sc(vector<HM2D::EdgeData*> {cont1, cont2});

		auto tree1 = HM2D::Contour::Tree::Assemble(*cont1);
		auto tree2 = HM2D::Contour::Tree::Assemble(*cont2);

		if (tree1.bound_contours().size() == 0) throw std::runtime_error("not a closed contour");
		if (tree2.bound_contours().size() == 0) throw std::runtime_error("not a closed contour");

		HM2D::Contour::Tree t;
		if (c2cpp::eqstring(op, "union")){
			t = HM2D::Contour::Clip::Union(tree1, tree2);
		} else if (c2cpp::eqstring(op, "difference")){
			t = HM2D::Contour::Clip::Difference(tree1, tree2);
		} else if (c2cpp::eqstring(op, "intersection")){
			t = HM2D::Contour::Clip::Intersection(tree1, tree2);
		} else if (c2cpp::eqstring(op, "xor")){
			t = HM2D::Contour::Clip::XOR(tree1, tree2);
		} else {
			throw std::runtime_error("unknown operation");
		}
	
		//heal to remove 0 angled sections
		HM2D::Contour::Clip::Heal(t);
		auto ae = t.alledges();
		//restore non-significant points
		if (!simplify && t.nodes.size() > 0){
			vector<Point*> allpnt;
			for (auto& v: HM2D::AllVertices(*cont1)) allpnt.push_back(v.get());
			for (auto& v: HM2D::AllVertices(*cont2)) allpnt.push_back(v.get());
			for (auto p: allpnt){
				auto fnd = HM2D::Finder::ClosestEdge(ae, *p);
				if (std::get<0>(fnd)<0) continue;
				auto ed = ae[std::get<0>(fnd)];
				//if point lies on edge and does not equal edge end points
				if (std::get<1>(fnd)<geps && ISIN_NN(std::get<2>(fnd), geps, 1-geps)){
					auto cont = t.find_node(ed.get());
					HM2D::Contour::Algos::GuaranteePoint(cont->contour, *p);
					ae = t.alledges();
				}
			}
		}
		//unscale and return
		sc.unscale(&ae);
		c2cpp::to_pp(ae, ret);
		return HMSUCCESS;
	} catch (std::exception& e){
		add_error_message(e.what());
		return HMERROR;
	}
}

namespace _c2part{
// c2_partion helper functions

vector<int> required_nedges(vector<HM2D::EdgeData>& part_contours, int nedges,
		const HM2D::VertexData& vd){
//splits int nedges by vector where each entry corresponds to i-th subcontour
	if (nedges<0) return vector<int>(part_contours.size(), -1);
	if (part_contours.size() == 1) return {nedges};

	vector<int> mined(part_contours.size(), 1);
	//each non-first and non-last vd increases number of needed edges
	for (int i=0; i<part_contours.size(); ++i){
		auto av = HM2D::AllVertices(part_contours[i]);
		for (auto p: vd){
			if (!HM2D::Finder::Contains(av, p.get()))
				continue;
			if (HM2D::Contour::IsClosed(part_contours[i])){
				++mined[i];
				continue;
			}
			if (HM2D::Contour::First(part_contours[i]) == p)
				continue;
			if (HM2D::Contour::Last(part_contours[i]) == p)
				continue;
			++mined[i];
		}
	}
	//for closed contours no less the 3 edges
	for (int i=0; i<part_contours.size(); ++i){
		if (mined[i]<3 && HM2D::Contour::IsClosed(part_contours[i])){
			mined[i] = 3;
		}
	}

	//nums double
	vector<double> lens;
	for (int i=0; i<part_contours.size(); ++i){
		lens.push_back(HM2D::Contour::Length(part_contours[i]));
	}
	vector<double> nums_double;
	double full_len = std::accumulate(lens.begin(), lens.end(), 0.0);
	for (int i=0; i<lens.size(); ++i){
		nums_double.push_back(lens[i]/full_len*nedges);
	}
	//nums int
	return HMMath::RoundVector(nums_double, mined);
}

void place_keep(HM2D::EdgeData& data, Point p, HM2D::VertexData& keep){
	if (HM2D::Contour::IsContour(data)){
		auto gp = HM2D::Contour::Algos::GuaranteePoint(data, p);
		keep.push_back(std::get<1>(gp));
		return;
	}
	auto fnd = HM2D::Finder::ClosestEdge(data, p);
	double locx = std::get<2>(fnd);
	if (ISEQ(locx, 0) || ISEQ(locx, 1)) return;
	shared_ptr<HM2D::Edge> e1 = data[std::get<0>(fnd)];
	shared_ptr<HM2D::Edge> e2 = std::make_shared<HM2D::Edge>(*e1);
	shared_ptr<HM2D::Vertex> v1=std::make_shared<HM2D::Vertex>(
			Point::Weigh(*e1->pfirst(), *e1->plast(), locx));
	e1->vertices[1] = v1;
	e2->vertices[0] = v1;
	data.push_back(e2);
	keep.push_back(v1);
}
void place_cross(HM2D::EdgeData& data, const HM2D::EdgeData& cc,
		HM2D::VertexData& keep){
	assert(HM2D::Contour::IsContour(data));
	assert(HM2D::Contour::IsContour(cc));
	auto cres = HM2D::Contour::Finder::CrossAll(data, cc);
	for (auto cross: cres){
		auto gp = HM2D::Contour::Algos::GuaranteePoint(data, std::get<1>(cross));
		keep.push_back(std::get<1>(gp));
	}
}
void place_angle(HM2D::EdgeData& data, double a0, bool keepbnd, HM2D::VertexData& keep){
	if (a0 >= 0){
		auto sc = HM2D::ECol::Algos::Simplified(data, a0, keepbnd, keep);
		for (auto v: AllVertices(sc)) keep.push_back(v);
	} else {
		for (auto v: AllVertices(data)) keep.push_back(v);
	}
	aa::no_duplicates(keep);
}
std::map<double, double> ref_weights_to_weights(const HM2D::EdgeData& data, const vector<double>& step){
	std::map<double, double> ret;
	for (int i=0; i<step.size(); i+=2){
		ret[step[i+1]] = step[i];
	}
	return ret;
}
std::map<double, double> ref_points_to_weights(const HM2D::EdgeData& data, const vector<double>& step){
	vector<double> step1;
	for (int i=0; i<step.size(); i+=3){
		Point p(step[i+1], step[i+2]);
		auto cr = HM2D::Contour::CoordAt(data, p);
		step1.push_back(step[i]);
		step1.push_back(std::get<1>(cr));
	}
	return ref_weights_to_weights(data, step1);
}
std::map<double, double> ref_lengths_to_weights(const HM2D::EdgeData& data, vector<double> step){
	double len = HM2D::Contour::Length(data);
	for (int i=0; i<step.size(); i+=2){
		step[i+1]/=len;
		if (step[i+1] < 0) step[i+1] = 1 + step[i+1];
	}
	return ref_weights_to_weights(data, step);
}
void force_start(HM2D::EdgeData& cont, Point start, HM2D::EdgeData& nonprocessed){
	auto av = HM2D::AllVertices(cont);
	int i1 = std::get<0>(HM2D::Finder::ClosestPoint(av, start));
	auto ret = HM2D::Contour::Assembler::Contour1(cont, av[i1].get(), 0);
	aa::constant_ids_pvec(cont, 0);
	aa::constant_ids_pvec(ret, 1);
	for (auto e: cont) if (e->id == 0){
		nonprocessed.push_back(e);
	}
	cont = ret;
}
void cut_by_startend(HM2D::EdgeData& cont, Point start, Point end, HM2D::EdgeData& nonprocessed){
	auto av = HM2D::AllVertices(cont);
	int i1 = std::get<0>(HM2D::Finder::ClosestPoint(av, start));
	int i2 = std::get<0>(HM2D::Finder::ClosestPoint(av, end));
	if (i1 == i2) force_start(cont, start, nonprocessed);
	auto ret = HM2D::Contour::Assembler::Contour1(cont, av[i1].get(), av[i2].get());
	aa::constant_ids_pvec(cont, 0);
	aa::constant_ids_pvec(ret, 1);
	for (auto e: cont) if (e->id == 0){
		nonprocessed.push_back(e);
	}
	std::swap(cont, ret);
}
HM2D::EdgeData const_partition(HM2D::EdgeData& data, double step, double a0,
		bool keepbnd, int nedges,
		const vector<HM2D::EdgeData>& ccont, const vector<Point>& keeppts){
	//place keeppts 
	HM2D::VertexData keep;
	for (auto& p: keeppts) place_keep(data, p, keep);
	//assemble simple contours
	vector<HM2D::EdgeData> simpc = HM2D::Contour::Assembler::SimpleContours(data);
	//place crosses
	for (auto& sc: simpc)
	for (auto& c: ccont) place_cross(sc, c, keep);
	//place angle points and boundary significant points
	for (auto& d: simpc) place_angle(d, a0, keepbnd, keep);
	//divide edges number with respect to simple contours lengths
	vector<int> ned(simpc.size(), -1);
	if (nedges > 0){
		ned = required_nedges(simpc, nedges, keep);
	}

	//building
	std::map<double, double> m; m[0]=step;
	HM2D::EdgeData ret;
	for (int i=0; i<simpc.size(); ++i){
		auto r = HM2D::Contour::Algos::WeightedPartition(m, simpc[i], ned[i], keep);
		ret.insert(ret.end(), r.begin(), r.end());
	}
	return ret;
}
HM2D::EdgeData ref_partition(HM2D::EdgeData& data, const char* algo,
		vector<double> step, double a0,
		bool keepbnd, int nedges,
		const vector<HM2D::EdgeData>& ccont,
		const vector<Point>& keeppts){
	//check for single contour
	if (!HM2D::Contour::IsContour(data)){
		throw std::runtime_error("Only singly connected contours "
				"are valid for ref_* partitions");
	}
	//place keeppts, crosses
	HM2D::VertexData keep;
	for (auto& p: keeppts) place_keep(data, p, keep);
	for (auto& c: ccont) place_cross(data, c, keep);
	//place angle points and boundary significant points
	place_angle(data, a0, keepbnd, keep);
	//modify step depending on algorithm
	std::map<double, double> weights;

	switch (std::map<std::string, int>({
			{"const", 1},
			{"ref_points", 2},
			{"ref_lengths", 3},
			{"ref_weights", 4}})[algo]){
	case 2:
		weights = ref_points_to_weights(data, step);
		break;
	case 4:
		weights = ref_weights_to_weights(data, step);
		break;
	case 3:
		weights = ref_lengths_to_weights(data, step);
		break;
	default:
		throw std::runtime_error("Unknown partition algorithm");
	}
	//do partition
	return (nedges<=0) ? HM2D::Contour::Algos::WeightedPartition(weights, data, keep)
	                   : HM2D::Contour::Algos::WeightedPartition(weights, data, nedges, keep);
}
} // c2_partion helper functions

int c2_partition(void* obj, const char* algo, int nstep, double* step, double a0,
		int keepbnd, int nedges, int ncrosses, void** crosses, int nkeeppts, double* _keeppts,
		double* _start, double* _end, void** ret){
	try{
		//read data
		auto _cont = static_cast<HM2D::EdgeData*>(obj);
		auto _ccont = c2cpp::to_pvec<HM2D::EdgeData>(ncrosses, crosses);
		vector<Point> keeppts = c2cpp::to_points2(nkeeppts, _keeppts);
		vector<double> basis = c2cpp::to_vec(nstep, step);
		Point start(0, 0), end(0, 0);
		if (_start) start.set(_start[0], _start[1]);
		if (_end) end.set(_end[0], _end[1]);

		//deepcopy
		HM2D::EdgeData cont;
		HM2D::DeepCopy(*_cont, cont);
		HM2D::EdgeData _tmp_ccont;
		for (int i=0; i<_ccont.size(); ++i){
			HM2D::DeepCopy(*_ccont[i], _tmp_ccont);
		}
		vector<HM2D::EdgeData> ccont = HM2D::Contour::Assembler::SimpleContours(_tmp_ccont);

		//scale
		ScaleBase sc = HM2D::Scale01(cont);
		for (int i=0; i<ccont.size(); ++i)
			HM2D::Scale(ccont[i], sc);
		switch (std::map<std::string, int>({
				{"const", 1},
				{"ref_points", 2},
				{"ref_lengths", 3},
				{"ref_weights", 4}})[algo]){
		case 2:
			for (int i=0; i<basis.size(); i+=3){
				basis[i] /= sc.L;
				basis[i+1] = (basis[i+1] - sc.p0.x)/sc.L;
				basis[i+2] = (basis[i+2] - sc.p0.y)/sc.L;
			}
			break;
		case 1: case 3:
			for (int i=0; i<basis.size(); ++i) basis[i] /= sc.L; 
			break;
		case 4:
			for (int i=0; i<basis.size(); i+=2) basis[i] /= sc.L; 
			break;
		default:
			throw std::runtime_error("unknown parition algorithm");
		}
		sc.scale(keeppts.begin(), keeppts.end());
		//cut by start/end
		HM2D::EdgeData nonprocessed;
		if (_start && (!_end || start == end)){
			_c2part::force_start(cont, start, nonprocessed);
		} else if (_start && _end){
			_c2part::cut_by_startend(cont, start, end, nonprocessed);
		}

		//build algorithm-depending result
		HM2D::EdgeData ret_;
		switch (std::map<std::string, int>({
				{"const", 1},
				{"ref_points", 2},
				{"ref_lengths", 3},
				{"ref_weights", 4}})[algo]){
		case 1:
			ret_ = _c2part::const_partition(cont, basis[0], a0,
					keepbnd, nedges, ccont, keeppts);
			break;
		case 2: case 3: case 4:
			ret_ = _c2part::ref_partition(cont, algo, basis, a0,
					keepbnd, nedges, ccont, keeppts);
			break;
		default:
			throw std::runtime_error("unknown parition algorithm");
		}
		//add non-processed
		ret_.insert(ret_.end(), nonprocessed.begin(), nonprocessed.end());
		HM2D::ECol::Algos::MergePoints(ret_);

		//unscale and return
		HM2D::Unscale(ret_, sc);
		c2cpp::to_pp(ret_, ret);
		return HMSUCCESS;
	} catch (std::exception& e){
		add_error_message(e.what());
		return HMERROR;
	}
}

int c2_matched_partition(void* obj, int nconds, void** conds, int npts, double* pts, double step,
                         double infdist, double power, double a0, void** ret){
	try{
		//scale
		Autoscale::D2 sc(static_cast<HM2D::EdgeData*>(obj));
		sc.add_data(c2cpp::to_pvec<HM2D::EdgeData>(nconds, conds));
		sc.scale(infdist);
		sc.scale(step);
	
		//read point conditions
		vector<std::pair<Point, double>> pconditions;
		for (int i=0; i<npts; ++i){
			Point p1(pts[3*i+1], pts[3*i+2]);
			pconditions.emplace_back(p1, pts[3*i]);
			sc.scale(pconditions.back().first);
			sc.scale(pconditions.back().second);
		}

		//assemble input data to contours
		vector<HM2D::EdgeData> input;
		{
			auto ac = HM2D::Contour::Assembler::SimpleContours(
					*static_cast<HM2D::EdgeData*>(obj));
			for (auto& c: ac) input.push_back(c);
		}
		vector<HM2D::EdgeData> conditions;
		for (int i=0; i<nconds; ++i){
			auto ac = HM2D::Contour::Assembler::SimpleContours(
					*static_cast<HM2D::EdgeData*>(conds[i]));
			for (auto& c: ac) conditions.push_back(c);
		}

		//process each input contour
		HM2D::EdgeData ret_;
		for (auto& icont: input){
			//set angle points
			HM2D::VertexData fixpoints;
			_c2part::place_angle(icont, a0, true, fixpoints);
			//set cross points
			for (auto& cond: conditions) _c2part::place_cross(icont, cond, fixpoints);
			//build partition
			auto r = HM2D::Contour::Algos::ConditionalPartition(icont, step, infdist,
					conditions, pconditions, power, fixpoints);
			//add to answer
			HM2D::DeepCopy(r, ret_);
		}
		sc.unscale(&ret_);
		c2cpp::to_pp(ret_, ret);
		return HMSUCCESS;
	} catch (std::exception& e){
		add_error_message(e.what());
		return HMERROR;
	}
}

int c2_segment_partition(double start, double end, double hstart, double hend,
		int ninternal, double* hinternal, int* nret, double** ret){
	try{
		double len = end - start;
		std::map<double, double> w;
		w[0] = hstart / len;
		w[1.] = hend / len;
		for (int i=0; i<ninternal; ++i){
			double c = hinternal[2*i];
			double s = hinternal[2*i+1];
			c = (c - start) / len;
			s/=len;
			w[c] = s;
		}
		vector<double> ret_ = HMMath::Partition01(w);
		for (auto& v: ret_) v = (v*len) + start;
		*ret = new double[ret_.size()];
		std::copy(ret_.begin(), ret_.end(), *ret);
		*nret = ret_.size();
		return HMSUCCESS;
	} catch (std::exception& e){
		add_error_message(e.what());
		return HMERROR;
	}
}

int c2_extract_subcontours(void* obj, int nplist, double* plist, void** ret){
	try{
		vector<Point> p0 = c2cpp::to_points2(nplist, plist);
		HM2D::EdgeData* ss = static_cast<HM2D::EdgeData*>(obj);
		Autoscale::D2 sc(ss);
		sc.scale(p0);

		vector<Point*> p1(p0.size(), nullptr);
		if (p0.size()<2) throw std::runtime_error("insufficient number of base points");
		vector<HM2D::EdgeData> et = HM2D::Contour::Assembler::SimpleContours(*ss);
		//1) find source contour
		HM2D::EdgeData* src=nullptr;
		double mind = 1e32;
		for (auto& c: et){
			auto ce1 = HM2D::Finder::ClosestEdge(c, p0[0]);
			if (std::get<0>(ce1) >= 0)
			if (std::get<1>(ce1) < mind){
				mind = std::get<1>(ce1);
				src = &c;
			}
		}
		if (src == nullptr) throw std::runtime_error("source contour was not found");
		//2) project base points
		for (int i=0; i<p0.size(); ++i){
			auto gp0 = HM2D::Contour::Algos::GuaranteePoint(*src, p0[i]);
			p0[i].set(*std::get<1>(gp0));
			p1[i] = std::get<1>(gp0).get();
		}
		//3) Reverse p0 if needed
		bool reversed_order = false;
		HM2D::Contour::R::Clockwise trev(*src, false);
		if (HM2D::Contour::IsClosed(*src)){
			if (p0.size() > 2){
				double w1 = std::get<1>(HM2D::Contour::CoordAt(*src, p0[0]));
				double w2 = std::get<1>(HM2D::Contour::CoordAt(*src, p0[1]));
				double w3 = std::get<1>(HM2D::Contour::CoordAt(*src, p0[2]));
				if (w2 <= w1) w2 += 1.0;
				if (w3 <= w1) w3 += 1.0;
				if (w3 < w2) reversed_order = true;
			}
		}
		//4) assemble
		for (int i=0; i<nplist-1; ++i){
			int i1 = i, i2 = i+1;
			if (reversed_order) std::swap(i1, i2);
			auto c = HM2D::Contour::Assembler::ShrinkContour(*src, p1[i1], p1[i2]);
			auto econt = new HM2D::EdgeData(); 
			HM2D::DeepCopy(c, *econt);
			sc.unscale(econt);
			ret[i] = econt;
		}
		return HMSUCCESS;
	} catch (std::exception& e){
		add_error_message(e.what());
		return HMERROR;
	}
}

int c2_connect_subcontours(int nobjs, void** objs, int nfix, int* fix, int shift,
		const char* close, void** ret){
	try{
		vector<HM2D::EdgeData*> ecols = c2cpp::to_pvec<HM2D::EdgeData>(nobjs, objs);
		Autoscale::D2 sc(ecols);

		std::string close_method(close);
		//assemble contours
		vector<HM2D::EdgeData> vconts;
		for (int i=0; i<nobjs; ++i) if (ecols[i]->size()>0){
			vconts.push_back(HM2D::Contour::Assembler::Contour1(*ecols[i], (*ecols[i])[0]->first().get()));
		}
		//construct new contour
		std::set<int> fixset(fix, fix+nfix);
		HM2D::EdgeData c = HM2D::Contour::Constructor::FromContours(
			vconts, close_method=="yes", shift==1, fixset);

		if (close_method == "force" && HM2D::Contour::IsOpen(c)){
			HM2D::Contour::Algos::AddLastPoint(c, HM2D::Contour::First(c));
		}
		//deepcopy
		HM2D::EdgeData ret_;
		HM2D::DeepCopy(c, ret_);
		sc.unscale(&ret_);
		c2cpp::to_pp(ret_, ret);
		return HMSUCCESS;
	} catch (std::exception& e){
		add_error_message(e.what());
		return HMERROR;
	}
}

//set boundary types to contour
int c2_assign_boundary_types(void* obj, int* bnd, int** revdif){
	try{
		auto cont = static_cast<HM2D::EdgeData*>(obj);
		auto whole_assign = [&cont](int bt, std::map<int, int>& mp)->void{
			for (int k=0; k<cont->size(); ++k){
				mp[k] = bt;
			}
		};
		auto bt_by_index = [&cont](int ind)->int&{
			return (*cont)[ind]->boundary_type;
		};
		c2cpp::assign_boundary_types(bnd, revdif, whole_assign, bt_by_index);
		return HMSUCCESS;
	} catch (std::exception& e){
		add_error_message(e.what());
		return HMERROR;
	}
}

int c2_contour_type(void* obj, int* ret){
	try{
		auto cont = static_cast<HM2D::EdgeData*>(obj);
		if (HM2D::Contour::IsContour(*cont)){
			if (HM2D::Contour::IsClosed(*cont)) *ret = 1;
			else *ret = 0;
			return HMSUCCESS;
		}

		auto scont = HM2D::Contour::Assembler::SimpleContours(*cont);
		if (scont.size() == 1){
			if (HM2D::Contour::IsClosed(scont[0])) *ret = 1;
			else *ret = 0;
			return HMSUCCESS;
		}
		*ret = 3;
		for (auto& s: scont){
			if (HM2D::Contour::IsOpen(s)) *ret = 4;
		}
		return HMSUCCESS;
	} catch (std::exception& e){
		add_error_message(e.what());
		return HMERROR;
	}
}

int c2_concatenate(int nobjs, void** objs, void** ret){
	try{
		HM2D::EdgeData ret_;
		for (auto& it: c2cpp::to_pvec<HM2D::EdgeData>(nobjs, objs)){
			ret_.insert(ret_.end(), it->begin(), it->end());
		}
		c2cpp::to_pp(ret_, ret);
		return HMSUCCESS;
	} catch (std::exception& e){
		add_error_message(e.what());
		return HMERROR;
	}
}

int c2_closest_points(void* obj, int npts, double* pts, const char* proj, double* ret){
	try{
		auto cont = static_cast<HM2D::EdgeData*>(obj);
		vector<Point> points = c2cpp::to_points2(npts, pts);
		vector<Point> ret_;
		if (c2cpp::eqstring(proj, "vertex")){
			auto av = HM2D::AllVertices(*cont);
			for (auto& p: points){
				int fnd = std::get<0>(HM2D::Finder::ClosestPoint(av, p));
				ret_.push_back(*av[fnd]);
			}
		} else if (c2cpp::eqstring(proj, "edge")){
			for (auto& p: points){
				auto fnd = HM2D::Finder::ClosestEdge(*cont, p);
				HM2D::Edge& e = *(*cont)[std::get<0>(fnd)];
				double w = std::get<2>(fnd);
				ret_.push_back(Point::Weigh(*e.pfirst(), *e.plast(), w));
			}
		}
		else throw std::runtime_error("unknown projection option");

		for (int i=0; i<ret_.size(); ++i){
			ret[2*i] = ret_[i].x;
			ret[2*i+1] = ret_[i].y;
		}
		return HMSUCCESS;
	} catch (std::exception& e){
		add_error_message(e.what());
		return HMERROR;
	}
}

int c2_point_at(void* obj, int index, double* ret){
	try{
		auto cont = static_cast<HM2D::EdgeData*>(obj);
		auto av = HM2D::AllVertices(*cont);
		if (index >= av.size()) throw std::runtime_error("index is out of range");
		ret[0] = av[index]->x;
		ret[1] = av[index]->y;
		return HMSUCCESS;
	} catch (std::exception& e){
		add_error_message(e.what());
		return HMERROR;
	}
}

int c2_end_points(void* obj, double* start, double* end){
	try{
		auto sc = HM2D::Contour::Assembler::Contour1(
			*static_cast<HM2D::EdgeData*>(obj));
		if (HM2D::Contour::IsClosed(sc)){
			throw std::runtime_error("closed contour has no end points");
		}
		start[0] = HM2D::Contour::First(sc)->x;
		start[1] = HM2D::Contour::First(sc)->y;
		end[0] = HM2D::Contour::Last(sc)->x;
		end[1] = HM2D::Contour::Last(sc)->y;
		return HMSUCCESS;
	} catch (std::exception& e){
		add_error_message(e.what());
		return HMERROR;
	}
}

int c2_tab_edgevert(void* obj, int* ret){
	try{
		auto cont = static_cast<HM2D::EdgeData*>(obj);
		auto av = HM2D::AllVertices(*cont);
		aa::enumerate_ids_pvec(av);
		for (auto& e: *cont){
			*ret++ = e->pfirst()->id;
			*ret++ = e->plast()->id;
		}
		return HMSUCCESS;
	} catch (std::exception& e){
		add_error_message(e.what());
		return HMERROR;
	}
}
int c2_tab_vertices(void* obj, double* ret){
	try{
		auto cont = static_cast<HM2D::EdgeData*>(obj);
		auto av = HM2D::AllVertices(*cont);
		for (int i=0; i<av.size(); ++i){
			*ret++ = av[i]->x;
			*ret++ = av[i]->y;
		}
		return HMSUCCESS;
	} catch (std::exception& e){
		add_error_message(e.what());
		return HMERROR;
	}
}
int c2_tab_btypes(void* obj, int* ret){
	try{
		auto cont = static_cast<HM2D::EdgeData*>(obj);
		for (int i=0; i<cont->size(); ++i){
			ret[i] = (*cont)[i]->boundary_type;
		}
		return HMSUCCESS;
	} catch (std::exception& e){
		add_error_message(e.what());
		return HMERROR;
	}
}

int c2_to_hm(void* doc, void* node, void* obj, const char* name, const char* fmt){
	try{
		HMXML::ReaderA* wr = static_cast<HMXML::ReaderA*>(doc);
		HMXML::Reader* sn = static_cast<HMXML::Reader*>(node);
		HM2D::EdgeData* c = static_cast<HM2D::EdgeData*>(obj);
		HM2D::Export::EColWriter(*c, wr, sn, name, fmt);
		return HMSUCCESS;
	} catch (std::exception& e){
		add_error_message(e.what());
		return HMERROR;
	}
}
