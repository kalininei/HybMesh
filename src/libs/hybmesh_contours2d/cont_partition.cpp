#include "cont_partition.hpp"
#include "piecewise.hpp"
#include "partition01.hpp"
#include "cont_assembler.hpp"
#include "treverter2d.hpp"

using namespace HM2D;
using namespace HM2D::Contour;
namespace cns = Algos;

namespace{
typedef std::list<shared_ptr<Vertex>> Vlist;
typedef std::set<shared_ptr<Vertex>> Vset;

vector<double> partition_new_points_w(double step, const EdgeData& contour){
	double len = Length(contour);
	if (len<1.5*step){
		if (IsOpen(contour)) return {0.0, 1.0};
		else step = len/3.1;
	}
	//calculates new points
	int n = (int)round(len/step);
	vector<double> w;
	for (int i=0; i<n+1; ++i) w.push_back((double)i/n);
	return w;
}
vector<double> partition_new_points_w(std::map<double, double> basis, const EdgeData& contour){
	assert(basis.size()>0);
	bool isconst = true;
	for (auto& v: basis) if (!ISEQ(v.second, basis.begin()->second)){
		isconst = false; break;
	}
	if (isconst) return partition_new_points_w(basis.begin()->second, contour);
	double len = Length(contour);
	for (auto& m: basis) m.second/=len;
	//for closed contour values at 0 and 1 should match
	if (IsClosed(contour)){
		double x0 = basis.begin()->first;
		double x1 = basis.rbegin()->first;
		double v0 = basis.begin()->second;
		double v1 = basis.rbegin()->second;
		if (ISEQ(x0, 0) && !ISEQ(x1, 1.0)) basis[1] = v0;
		else if (ISEQ(x1, 1.0) && !ISEQ(x0, 0)) basis[0] = v1;
		else if (ISEQ(x0, 0.0) && ISEQ(x1, 1.0)){
			basis.begin()->second = (v0 + v1)/2.0;
			basis.rbegin()->second = (v0 + v1)/2.0;
		} else {
			//no values at 0 and 1 and more then 1 entries in basis
			x1 -= 1.0;
			double t = -x0/(x1-x0);
			double v = (1-t)*v0+t*v1;
			basis[0] = v;
			basis[1] = v;
		}
	}
	auto ret = HMMath::Partition01(basis);
	if (IsClosed(contour) && ret.size()<4) return {0, 1.0/3.0, 2.0/3.0, 1.0};
	else return ret;
}

double build_substep(double step, const EdgeData&, shared_ptr<Vertex>, shared_ptr<Vertex>){ return step; }

struct Conditions2D{
	static const int LINEARIZE = 50;
	vector<EdgeData> contcond;
	vector<std::pair<Point, double>> pointcond;
	double default_step;
	double influence_dist;
	int Ncond() const {return contcond.size() + pointcond.size(); }
	double pw;

	std::vector<BoundingBox> boxes;
	
	bool is_out_of_box(int i, const Point& p){
		if (boxes.size() == 0){
			for (auto& cnt: contcond){
				boxes.push_back(BBox(cnt, influence_dist));
			}
			for (auto& p: pointcond){
				boxes.push_back(BoundingBox(p.first, influence_dist));
			}
		};
		return boxes[i].whereis(p) == OUTSIDE;
	}

	//returns step size and weight
	std::pair<double, double> cond_for_point(int i, const Point& p){
		double dist = Point::dist(pointcond[i].first, p);
		if (dist > influence_dist) return std::make_pair(1.0, 0.0);
		double h = pointcond[i].second;
		double w = (influence_dist - dist)/influence_dist;
		return std::make_pair(h, w);
	}

	//returns step size and weight
	std::pair<double, double> cond_for_contour(int i, const Point& p){
		auto cle = FindClosestEdge(contcond[i], p);
		double dist = std::get<1>(cle);
		if (dist > influence_dist) return std::make_pair(1.0, 0.0);
		double h = contcond[i][std::get<0>(cle)]->length();
		double w = (influence_dist - dist)/influence_dist;
		return std::make_pair(h, w);
	}

	double aver(double w1, double x1, double x2){
		double w = pow(w1, 1.0/pw);
		return w*x1 + (1-w)*x2;
	}

	//ws - weight-step pair
	double aver(std::vector<std::pair<double, double>>& ws){
		vector<double> w;
		for (auto& it: ws) w.push_back(pow(it.first, 1.0/pw));
		double ret = 0;
		for (int i=0; i<w.size(); ++i) ret += ws[i].second*w[i];
		return ret/std::accumulate(w.begin(), w.end(), 0.0);
	}

	double stepfrom(int i, const Point& p){
		if (is_out_of_box(i, p)) return default_step;

		std::pair<double, double> step_delta;
		if (i < contcond.size()){
			step_delta = cond_for_contour(i, p);
		} else {
			step_delta = cond_for_point(i-contcond.size(), p);
		}
		return aver(step_delta.second, step_delta.first, default_step);
	}

	double Value(const Point& p){
		//vector<double> hs;
		//for (int i=0; i<Ncond(); ++i) hs.push_back(stepfrom(i, p));
		//return *min_element(hs.begin(), hs.end());
		std::vector<std::pair<double, double>> ws;  //weight-step
		for (int i=0; i<Ncond(); ++i){
			std::pair<double, double> step_delta;
			if (i < contcond.size()){
				step_delta = cond_for_contour(i, p);
			} else {
				step_delta = cond_for_point(i-contcond.size(), p);
			}
			if (step_delta.second > geps){
				ws.emplace_back(
					step_delta.second,
					aver(step_delta.second, step_delta.first, default_step));
			}
		}

		if (ws.size()==0) return default_step;
		if (ws.size()==1) return ws.begin()->second;
		return aver(ws);
	}

	std::map<double, double> linear_conditions(const EdgeData& cont){
		std::map<double, double> ret;
		std::vector<double> ws(LINEARIZE);
		for (int i=0; i<ws.size(); ++i) ws[i] = (double)i/(ws.size()-1);

		vector<Point> wp = WeightPoints(cont, ws);

		for (int i=0; i<wp.size(); ++i){
			double v = Value(wp[i]);
			if (ret.size()>1){
				auto mit = std::prev(ret.end());
				auto mmit = std::prev(std::prev(ret.end()));
				if (ISEQ(mit->second, mmit->second) && ISEQ(mit->second, v)){
					ret.erase(mit);
				}
			}
			ret.emplace(ws[i], v);
		}
		return ret;
	}
};
Conditions2D build_substep(Conditions2D step, const EdgeData&,
		shared_ptr<Vertex>, shared_ptr<Vertex>){
	return step;
}

vector<double> partition_new_points_w(Conditions2D step, const EdgeData& cont){
	auto basis = step.linear_conditions(cont);
	return partition_new_points_w(basis, cont);
}
std::pair<
	std::map<double, double>::iterator,
	bool
> insert_into_basismap(std::map<double, double>& m, double v){
	auto ins = m.emplace(v, -1);
	if (ins.second){
		if (ins.first == m.begin()) ins.first->second = (++m.begin())->second;
		else if (ins.first == --m.end()) ins.first->second = (++m.rbegin())->second;
		else{
			auto it0 = std::prev(ins.first);
			auto it1 = ins.first;
			auto it2 = std::next(ins.first);
			double t = (it1->first - it0->first)/(it2->first - it0->first);
			it1->second = (1-t)*it0->second + t*it2->second;
		}
	}
	return ins;
}
std::map<double, double> shrink01_basismap(std::map<double, double>& m,
		std::map<double, double>::iterator it0, 
		std::map<double, double>::iterator it1){
	assert(it0->first < it1->first && it1!=m.end());
	m.erase(m.begin(), it0);
	m.erase(++it1, m.end());
	double start = m.begin()->first;
	double len = m.rbegin()->first - m.begin()->first;
	std::map<double, double> ret;
	for (auto it: m) ret.emplace((it.first-start)/len, it.second);
	return ret;
}
std::map<double, double> build_substep(std::map<double, double> step,
		const EdgeData& cont, shared_ptr<Vertex> p0, shared_ptr<Vertex> p1){
	int i0 = PInfo(cont, p0.get()).index, i1 = PInfo(cont, p1.get()).index;
	if (IsOpen(cont)){
		assert(i1 > i0);
		if (i0 == 0 && i1 == cont.size()){
			insert_into_basismap(step, 0);
			insert_into_basismap(step, 1);
			return step;
		}
	} else {
		if (i0 == 0 && i1 == 0) return step;
	}
	auto ew = EWeights(cont);
	double w0 = ew[i0], w1 = ew[i1];
	if (IsClosed(cont)){
		//enlarge step by adding fictive entries to its end
		for (auto it = step.begin(); it!=step.end(); ++it){
			if (it->first >= 1) break;
			step[it->first + 1] = it->second;
			step[it->first - 1] = it->second;
		}
		if (ISEQLOWER(w1, w0)) w1 += 1;
	}
	auto ins0 = insert_into_basismap(step, w0);
	auto ins1 = insert_into_basismap(step, w1);
	auto ret = shrink01_basismap(step, ins0.first, ins1.first);
	return ret;
}

//build a new contour based on input contour begin/end points
template<class A>
EdgeData partition_core(A& step, const EdgeData& contour){
	vector<double> w = partition_new_points_w(step, contour);
	if (w.size() == 2){
		EdgeData ret;
		ret.emplace_back(new Edge{First(contour), Last(contour)});
		return ret;
	}
	if (w.size()<2){
		partition_core(step, contour);
	}
	vector<double> w2(w.begin()+1, w.end()-1);
	vector<Point> wp = WeightPoints(contour, w2);
	//Construct new contour
	VertexData cpoints;
	cpoints.push_back(First(contour));
	std::transform(wp.begin(), wp.end(), std::back_inserter(cpoints),
			[](const Point& x){ return std::make_shared<Vertex>(x); });
	cpoints.push_back(Last(contour));
	return Assembler::Contour1(cpoints);
}

template<class A>
EdgeData partition_section(A& step, const EdgeData& cont, shared_ptr<Vertex> pstart, shared_ptr<Vertex> pend){
	EdgeData sub = Assembler::Contour1(cont, pstart.get(), pend.get());
	A substep = build_substep(step, cont, pstart, pend);
	//temporary reverse contour with single edge if needed
	std::unique_ptr<ReallyRevert> trev;
	if (First(sub) != pstart) trev.reset(new ReallyRevert(sub));
	assert(First(sub) == pstart);
	return partition_core(substep, sub);
}

//returns repartitioned copy of a contour
template<class A>
EdgeData partition_core(A& step, const EdgeData& contour, const Vlist& keep){
	auto it0 = keep.begin(), it1 = std::next(it0);
	EdgeData ret;
	while (it1 != keep.end()){
		Connect(ret, partition_section(step, contour, *it0, *it1));
		//if sub.size() == 1 then its direction is not defined
		//so we need to check the resulting direction.
		//We did it after second union when direction matters
		if (it0 == std::next(keep.begin())){
			if (Last(ret) != *it1) Reverse(ret);
		}

		++it0; ++it1;
	}
	return ret;
}

Vlist build_sorted_pnt(const EdgeData& contour, const VertexData& keepit){
	//sort points in keepit and add start, end point there
	VertexData orig_pnt = Contour::OrderedPoints(contour);
	Vset keepset(keepit.begin(), keepit.end());

	Vlist keep_sorted;
	std::copy_if(orig_pnt.begin(), orig_pnt.end(), std::back_inserter(keep_sorted),
		[&keepset](shared_ptr<Vertex>& p){ return keepset.find(p) != keepset.end(); }
	);

	//add first, last points
	if (!IsClosed(contour)){
		if (keep_sorted.size() == 0 || *keep_sorted.begin() != orig_pnt[0])
			keep_sorted.push_front(orig_pnt[0]);
		if (*keep_sorted.rbegin() != orig_pnt.back())
			keep_sorted.push_back(orig_pnt.back());
	} else if (keep_sorted.size() == 1){
		keep_sorted.push_back(*keep_sorted.begin());
	} else if (keep_sorted.size() > 0 && *keep_sorted.begin() != *keep_sorted.rbegin()){
		keep_sorted.push_back(*keep_sorted.begin());
	} else if (keep_sorted.size() == 0){
		keep_sorted.push_back(orig_pnt[0]);
		keep_sorted.push_back(orig_pnt[0]);
	}

	return keep_sorted;
}
template<class A>
EdgeData partition_with_keepit(A& step, const EdgeData& contour,
		const VertexData& keepit){
	Vlist keep_sorted = build_sorted_pnt(contour, keepit);
	//call core procedure
	return partition_core(step, contour, keep_sorted);
}
}//namespace

EdgeData cns::Partition(double step, const EdgeData& contour, PartitionTp tp){
	switch (tp){
		case PartitionTp::IGNORE_ALL:
			return Partition(step, contour);
		case PartitionTp::KEEP_ALL:
			return Partition(step, contour, AllVertices(contour));
		case PartitionTp::KEEP_SHAPE:
			return Partition(step, contour, CornerPoints(contour));
	};
}
EdgeData cns::WeightedPartition(const std::map<double, double>& basis,
		const EdgeData& contour, PartitionTp tp){
	switch (tp){
		case PartitionTp::IGNORE_ALL:
			return WeightedPartition(basis, contour);
		case PartitionTp::KEEP_ALL:
			return WeightedPartition(basis, contour, AllVertices(contour));
		case PartitionTp::KEEP_SHAPE:
			return WeightedPartition(basis, contour, CornerPoints(contour));
	};
}
EdgeData cns::Partition(double step, const EdgeData& contour,
		const VertexData& keepit){
	return partition_with_keepit(step, contour, keepit);
}
EdgeData cns::WeightedPartition(const std::map<double, double>& basis,
		const EdgeData& contour,
		const VertexData& keepit){
	return partition_with_keepit(basis, contour, keepit);
}
EdgeData cns::WeightedPartition(std::map<double, double> basis,
		const EdgeData& contour, int nedges,
		const VertexData& keepit){
	if (nedges<=0) return WeightedPartition(basis, contour, keepit);
	Vlist keep_sorted = build_sorted_pnt(contour, keepit);
	if (IsClosed(contour) && (keep_sorted.size()>nedges || nedges<3))
		throw std::runtime_error("Failed to satisfy forced number of edges property");
	if (IsOpen(contour) && keep_sorted.size()>nedges+1)
		throw std::runtime_error("Failed to satisfy forced number of edges property");
	if (keep_sorted.size() == nedges+1){
		VertexData p(keep_sorted.begin(), keep_sorted.end());
		return Assembler::Contour1(p, false);
	} else {
		//lengths of keep points
		std::vector<double> lengths;
		for (auto it = keep_sorted.begin(); it!=--keep_sorted.end(); ++it){
			auto c = CoordAt(contour, **it);
			lengths.push_back(std::get<0>(c));
		}
		lengths.push_back(Length(contour));
		//h(len) function
		HMMath::LinearPiecewise pw;
		for (auto& v: basis){
			double x = v.first*lengths.back();
			double y = v.second;
			pw.add_point(x, y);
		}
		if (!ISEQ(basis.begin()->first, 0)){
			if (IsClosed(contour)){
				pw.add_point((basis.rbegin()->first - 1)*lengths.back(), basis.rbegin()->second);
			} else {
				pw.add_point(0, basis.begin()->second);
			}
		}
		if (!ISEQ(basis.rbegin()->first, 1.0)){
			if (IsClosed(contour)){
				pw.add_point((basis.begin()->first + 1)*lengths.back(), basis.begin()->second);
			} else {
				pw.add_point(lengths.back(), basis.rbegin()->second);
			}
		}
		//modify data to fit nedges
		double haver = pw.Average(0, lengths.back());
		double naver = lengths.back()/haver;
		for (auto& v: basis) v.second*=(naver/nedges);
		pw *= (naver/nedges);
		//build approximations for each subcontour partition number
		std::vector<double> nums_double;
		for (int i=0; i<lengths.size()-1; ++i){
			double len = lengths[i+1] - lengths[i];
			double av = pw.Average(lengths[i], lengths[i+1]);
			nums_double.push_back(len/av);
		}
		double sum_double = std::accumulate(nums_double.begin(), nums_double.end(), 0.0);
		for (auto& v: nums_double) v*=(double(nedges)/sum_double);
		vector<int> nums = HMMath::RoundVector(nums_double, vector<int>(nums_double.size(), 1));
		//build subcontours
		auto it0 = keep_sorted.begin(), it1 = std::next(it0);
		auto itn = nums.begin();
		EdgeData ret;
		while (it1 != keep_sorted.end()){
			EdgeData psub;
			double leftx=-std::numeric_limits<double>::max(), rightx=std::numeric_limits<double>::max();
			int lefty=0, righty=0;
			int tries=0;
			double coef = 1;
			if (*itn == 1){
				psub = Assembler::Contour1(VertexData{*it0, *it1});
			} else for (tries=0; tries<100; ++tries){
				//try
				auto b = basis;
				for (auto& it: b) it.second*=coef;
				psub = partition_section(b, contour, *it0, *it1);
				if (psub.size() == *itn) break;
				//approximate
				if (psub.size() < *itn){
					leftx = coef; lefty = psub.size();
				} else{
					rightx = coef; righty = psub.size();
				}
				if (leftx == -std::numeric_limits<double>::max()){
					coef = rightx * righty / (*itn);
				} else if (rightx == std::numeric_limits<double>::max()){
					coef = leftx * lefty / (*itn);
				} else{
					double t = double(*itn - lefty)/(righty - lefty);
					coef = (1.0-t)*leftx + t*rightx;
				}
			}
			if (tries >= 100) throw std::runtime_error(
				"Failed to satisfy forced number of edges property");
			Connect(ret, psub);
			if (it0 == std::next(keep_sorted.begin())){
				if (Last(ret) != *it1) Reverse(ret);
			}
			++it0; ++it1; ++itn;
		}
		assert(ret.size() == nedges);
		return ret;
	}
}
EdgeData cns::ConditionalPartition(const EdgeData& input, double step, double influence,
		const vector<EdgeData>& condconts,
		const vector<std::pair<Point, double>>& condpoints,
		double pw, const VertexData& keepit){
	Conditions2D cond;
	cond.contcond = condconts;
	cond.pointcond = condpoints;
	cond.default_step = step;
	cond.influence_dist = influence;
	cond.pw = pw;
	return partition_with_keepit(cond, input, keepit);
}

