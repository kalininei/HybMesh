#include "cont_partition.hpp"
#include "constructor.hpp"

using namespace HMCont2D;
namespace cns = Algos;

// ===================================== Partition implementation
namespace {

class Partition01{
// builds partition of [0, 1] segment using weights associated with
// section lengths
	const double Epsilon = 1e-5;
	std::map<double, double> data;    //weight -> value data map
	double w0() const{ return data.begin()->first; }
	double w1() const{ return data.rbegin()->first; }
	double val0() const{ return data.begin()->second; }
	double val1() const{ return data.rbegin()->second; }
	int nseg() const{ return data.size() - 1; }

	Partition01(){}

	Partition01(const std::map<double, double>& inp): data(inp){
		//only positive data
		auto it = data.begin();
		while(it!=data.end()){
			if (ISEQLOWER(it->second, 0.0)) it = data.erase(it);
			else ++it;
		}
		if (data.size() == 0) throw GeomError("Invalid segment partition input");
		//if single value add dummy entry
		if (inp.size() == 1) data[inp.begin()->first + 1]  = inp.begin()->second;
		//leave only [0->1] segments
		Crop(0, 1);
	}
	
	void Crop(double v0, double v1){
		if (data.find(v0) == data.end()){
			data.emplace(v0, Value(v0));
		}
		if (data.find(v1) == data.end()){
			data.emplace(v1, Value(v1));
		}
		auto fnd0 = data.find(v0);
		auto fnd1 = data.find(v1);
		++fnd1;
		data.erase(data.begin(), fnd0);
		data.erase(fnd1, data.end());
	}

	std::pair<
		std::map<double, double>::const_iterator, 
		double
	> _Geti(double v) const{
		if (ISLOWER(v, w0())) return std::make_pair(data.end(), w0() - v);   //negative weight for <v0
		if (ISGREATER(v, w1())) return std::make_pair(data.end(), v - w1()); //positive w for >v1
		if (ISEQ(v, w0())) return std::make_pair(data.begin(), 0);
		if (ISEQ(v, w1())) return std::make_pair(----data.end(), 1);
		auto fndu = data.upper_bound(v);
		auto fndl = fndu; --fndl;
		double t = (v - fndl->first)/(fndu->first - fndl->first);
		return std::make_pair(fndl, t);
	}
	double _Integr(std::map<double, double>::const_iterator it, double t0, double t1) const{
		auto itnext = it; ++itnext;
		double len = (itnext->first - it->first) * (t1 - t0);
		double p1 = (1-t0)*it->second + t0*itnext->second;
		double p2 = (1-t1)*it->second + t1*itnext->second;
		return (p1 + p2)/2*len;
	}

	double Integrate(double x0, double x1) const{
		if (x1<x0) return -Integrate(x1, x0);
		if (ISEQ(x1, x0)) return 0;
		auto f1 = _Geti(x0);
		auto f2 = _Geti(x1);
		bool before1 = f1.first == data.end() && f1.second < 0;
		bool before2 = f2.first == data.end() && f2.second < 0;
		bool after1 = f1.first == data.end() && f1.second > 0;
		bool after2 = f2.first == data.end() && f2.second > 0;
		if (before1 && before2) return val0() * (x1-x0);
		if (after1 && after2) return val1() * (x1-x0);
		double ret = 0;
		if (before1){
			ret += val0() * (w0() - x0);
			x0 = w0();
			f1.first = data.begin();
			f1.second = 0.0;
		}
		if (after2){
			ret += val1() * (x1 - w1());
			x1 = w1();
			f2.first = ----data.end();
			f2.second = 1.0;
		}
		if (f1.first == f2.first){
			ret += _Integr(f1.first, f1.second, f2.second);
		} else {
			ret +=_Integr(f1.first, f1.second, 1.0);
			ret +=_Integr(f2.first, 0.0, f2.second);
			for (auto it = ++f1.first; it!=f2.first; ++it){
				ret += _Integr(it, 0.0, 1.0);
			}
		}
		return ret;
	}

	double Value(double x) const{
		auto f = _Geti(x);
		if (f.first == data.end()){
			return f.second>0 ? val1() : val0();
		}
		double v0 = f.first->second;
		double v1 = (++f.first)->second;
		double t = f.second;
		return (1-t)*v0 + t*v1;
	}

	vector<double> _RunForward() const{
		vector<double> ret {w0()};
		while (1){
			double s0 = Value(ret.back());
			ret.push_back(ret.back() + s0);
			while (1){
				double& m1 = ret[ret.size() - 1];
				double& m2 = ret[ret.size() - 2];
				double len = m1 - m2;
				double sr = Integrate(m2, m1)/len;
				double diff = fabs(len - sr);
				if (diff < Partition01::Epsilon) break;
				else m1 = m2 + sr;
			}
			if (ISEQGREATER(ret.back(), w1())) break;
		}
		if (ret.size() > 2 && (ret.back() - w1()) > 2 * (w1()-ret[ret.size()-2])){
			ret.resize(ret.size()-1);
		}
		return ret;
	}

	vector<double> _RunBackward() const{
		auto revcoord = [&](double x){ return w0() + w1() - x; };
		//revert data
		Partition01 rdata;
		for (auto it=data.rbegin(); it!=data.rend(); ++it)
			rdata.data.emplace(revcoord(it->first), it->second);
		//calculate using reverted data
		vector<double> ret = rdata._RunForward();
		//revert result back
		for (auto& r: ret) r = revcoord(r);
		std::reverse(ret.begin(), ret.end());
		//return
		return ret;
	}

	static vector<double> _AverageFB(const vector<double>& f, const vector<double>& b){
		assert(f.size() == b.size());
		vector<double> ret;
		for (int i=0; i<f.size(); ++i){
			double x = f[i], y = b[i];
			double t = (float)i/(f.size() - 1);
			ret.push_back((1-t)*x + t*y);
		}
		return ret;
	}
	vector<double> _ConnectFB(const vector<double>& f, const vector<double>& b, int fi, int bi) const{
		double cp = (f[fi] + b[bi])/2.0;
		vector<double> ret;
		double cf = (cp - w0())/(f[fi] - w0());
		for (int i=0; i<fi; ++i){
			ret.push_back(cf*(f[i]-w0()) + w0());
		}
		double cb = (cp - w1())/(b[bi] - w1());
		for (int i=bi; i<b.size(); ++i){
			ret.push_back(cb*(b[i]-w1()) + w1());
		}
		return ret;
	}
	vector<double> _ConnectFB(const vector<double>& f, const vector<double>& b) const{
		//find closest points
		double dist0=gbig;
		int ic=-1, jc=-1;
		for (int i=1; i<f.size(); ++i){
			for (int j=0; j<b.size()-1; ++j){
				double d = fabs(f[i] - b[j]);
				if (ISLOWER(d, dist0)){
					dist0 = d;
					ic = i;
					jc = j;
				}
			}
		}
		assert(ic>=0 && jc>=0);
		auto ret = _ConnectFB(f, b, ic, jc);
		//if input data is symmetric then force output to be symmetric
		bool is_sym = (f.size() == b.size() &&
			ISEQ(f[ic] - w0(), w1() - b[b.size() - 1 - ic]) &&
			ISEQ(w1() - b[jc], f[f.size() - 1 - jc] - w0()));
		if (is_sym){
			vector<double> ret2;
			for (auto it = ret.rbegin(); it!=ret.rend(); ++it){
				ret2.push_back(w0() + w1() - *it);
			}
			ret = _AverageFB(ret, ret2);
		}
		//return
		return ret;
	}

	double Quality(const vector<double>& ans) const{
		double ret = 0;
		for (int i=0; i<(int)ans.size() - 1; ++i){
			double d = ans[i+1] - ans[i];
			if (d<=0) return gbig;
			double s = Integrate(ans[i], ans[i+1]) / d;
			double err = fabs(d-s)/s;
			ret += err;
		}
		return ret;
	}

	vector<double> Run() const{
		auto ansf = _RunForward();
		auto ansb = _RunBackward();
		if (ansf.size() == 2 && ansb.size() == 2) return {w0(), w1()};
		auto ret = _ConnectFB(ansf, ansb);
		double q = Quality(ret);
		//try one more algorithm
		if (ansf.size() == ansb.size()){
			auto ret2 = _AverageFB(ansf, ansb);
			double q2 = Quality(ret2);
			if (q2 < q) { ret = ret2; q = q2; }
		}
		return ret;
	}
public:
	static vector<double> Build(const std::map<double, double>& wts){
		Partition01 part(wts);
		return part.Run();
	}
};

vector<double> partition_new_points_w(double step, const Contour& contour){
	double len = contour.length();
	if (len<1.5*step){
		if (!contour.is_closed()) return {0.0, 1.0};
		else step = len/3.1;
	}
	//calculates new points
	int n = (int)round(len/step);
	vector<double> w;
	for (int i=0; i<n+1; ++i) w.push_back((double)i/n);
	return w;
}
vector<double> partition_new_points_w(std::map<double, double> basis, const Contour& contour){
	assert(basis.size()>0);
	if (basis.size() == 1) return partition_new_points_w(basis.begin()->second, contour);
	double len = contour.length();
	for (auto& m: basis) m.second/=len;
	//for closed contour values at 0 and 1 should match
	if (contour.is_closed()){
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
	auto ret = Partition01::Build(basis);
	if (contour.is_closed() && ret.size()<4) return {0, 1.0/3.0, 2.0/3.0, 1.0};
	else return ret;

	//auto funcstep = [&](double x)->double{
		//if (x<=basis.begin()->first) return basis.begin()->second;
		//if (x>=basis.rbegin()->first) return basis.rbegin()->second;
		//auto it0 = basis.begin();
		//while (x>it0->first) ++it0;
		//auto it1 = it0--;
		//double ret = (it0->second*(it1->first-x) + it1->second*(x-it0->first));
		//return x + ret/(it1->first - it0->first);
	//};
	//vector<double> ret;
	//ret.push_back(0);
	//while (ret.back()<1.0) ret.push_back(funcstep(ret.back()));
	////analyze last entry
	//if ((ret.back() - 1.0) > 1.0 - ret[ret.size() - 1]) ret.resize(ret.size()-1);
	//if (ret.size()<3) return {0, 1};
	//double stretch_coef = 1.0/ret.back();
	//for (auto& x: ret) x*=stretch_coef;
	//return ret;
}

double build_substep(double step, const Contour&, Point*, Point*){ return step; }

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
			it1->second = (1-t)*it0->second + t*it1->second;
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
		const Contour& cont, Point* p0, Point* p1){
	int i0 = cont.pinfo(p0).index, i1 = cont.pinfo(p1).index;
	if (!cont.is_closed()){
		assert(i1 > i0);
		if (i0 == 0 && i1 == cont.size()) return step;
	} else {
		if (i0 == 0 && i1 == 0) return step;
	}
	auto ew = HMCont2D::Contour::EWeights(cont);
	double w0 = ew[i0], w1 = ew[i1];
	if (cont.is_closed()){
		//enlarge step by adding fictive entries to its end
		for (auto it = step.begin(); it!=step.end(); ++it){
			if (it->first >= 1) break;
			step[it->first + 1] = it->second;
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
Contour partition_core(A& step, const Contour& contour, PCollection& pstore){
	vector<double> w = partition_new_points_w(step, contour);
	if (w.size() == 2){
		Contour ret;
		ret.add_value(Edge{contour.first(), contour.last()});
		return ret;
	}
	if (w.size()<2){
		partition_core(step, contour, pstore);
	}
	vector<double> w2(w.begin()+1, w.end()-1);
	PCollection wp = Contour::WeightPoints(contour, w2);
	//add new points to pstore
	pstore.Unite(wp);
	//Construct new contour
	vector<Point*> cpoints;
	cpoints.push_back(contour.first());
	std::transform(wp.begin(), wp.end(), std::back_inserter(cpoints),
			[](shared_ptr<Point> x){ return x.get(); });
	cpoints.push_back(contour.last());
	return Constructor::ContourFromPoints(cpoints);
}
//returns repartitioned copy of a contour
template<class A>
Contour partition_core(A& step, const Contour& contour, PCollection& pstore, const std::list<Point*>& keep){
	auto it0 = keep.begin(), it1 = std::next(it0);
	Contour ret;
	while (it1 != keep.end()){
		Contour sub = Assembler::Contour1(contour, *it0, *it1);
		A substep = build_substep(step, contour, *it0, *it1);
		Contour psub = partition_core(substep, sub, pstore);
		ret.Unite(psub);
		//if sub.size() == 1 then its direction is not defined
		//so we need to check the resulting direction.
		//We did it after second Unition when direction matters
		if (it0 == std::next(keep.begin())){
			if (ret.last() != *it1) ret.Reverse();
		}

		++it0; ++it1;
	}
	return ret;
}
template<class A>
Contour partition_with_keepit(A& step, const Contour& contour, PCollection& pstore,
		const std::vector<Point*>& keepit){
	//sort points in keepit and add start, end point there
	vector<Point*> orig_pnt = contour.ordered_points();
	std::set<Point*> keepset(keepit.begin(), keepit.end());

	std::list<Point*> keep_sorted;
	std::copy_if(orig_pnt.begin(), orig_pnt.end(), std::back_inserter(keep_sorted),
		[&keepset](Point* p){ return keepset.find(p) != keepset.end(); }
	);

	//add first, last points
	if (!contour.is_closed()){
		if (*keep_sorted.begin() != orig_pnt[0])
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

	//call core procedure
	return partition_core(step, contour, pstore, keep_sorted);
}
}//namespace

Contour cns::Partition(double step, const Contour& contour, PCollection& pstore, PartitionTp tp){
	switch (tp){
		case PartitionTp::IGNORE_ALL:
			return Partition(step, contour, pstore);
		case PartitionTp::KEEP_ALL:
			return Partition(step, contour, pstore, contour.all_points());
		case PartitionTp::KEEP_SHAPE:
			return Partition(step, contour, pstore, contour.corner_points());
		default:
			assert(false);
	};
}
Contour cns::WeightedPartition(const std::map<double, double>& basis,
		const Contour& contour, PCollection& pstore, PartitionTp tp){
	switch (tp){
		case PartitionTp::IGNORE_ALL:
			return WeightedPartition(basis, contour, pstore);
		case PartitionTp::KEEP_ALL:
			return WeightedPartition(basis, contour, pstore, contour.all_points());
		case PartitionTp::KEEP_SHAPE:
			return WeightedPartition(basis, contour, pstore, contour.corner_points());
		default:
			assert(false);
	};
}

Contour cns::Partition(double step, const Contour& contour, PCollection& pstore,
		const std::vector<Point*>& keepit){
	return partition_with_keepit(step, contour, pstore, keepit);
}
Contour cns::WeightedPartition(const std::map<double, double>& basis,
		const Contour& contour, PCollection& pstore,
		const std::vector<Point*>& keepit){
	return partition_with_keepit(basis, contour, pstore, keepit);
}
