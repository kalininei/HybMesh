#include "partition01.hpp"

namespace{
// builds partition of [0, 1] segment using weights associated with
// section lengths
class partition01{
	const double Epsilon = 1e-5;
	std::map<double, double> data;    //weight -> value data map
	double w0() const{ return data.begin()->first; }
	double w1() const{ return data.rbegin()->first; }
	double val0() const{ return data.begin()->second; }
	double val1() const{ return data.rbegin()->second; }
	int nseg() const{ return data.size() - 1; }

	partition01(){}

	partition01(const std::map<double, double>& inp): data(inp){
		//only positive data
		auto it = data.begin();
		while(it!=data.end()){
			if (ISEQLOWER(it->second, 0.0)) it = data.erase(it);
			else ++it;
		}
		if (data.size() == 0) throw std::runtime_error("Invalid segment partition input");
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

	std::vector<double> _RunForward() const{
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
				if (diff < partition01::Epsilon) break;
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
		partition01 rdata;
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
		partition01 part(wts);
		return part.Run();
	}
};

}

vector<double> HMMath::Partition01(const std::map<double, double>& wts){
	return partition01::Build(wts);
}


vector<int> HMMath::RoundVector(const vector<double>& vect, const vector<int>& minsize){
	int sum = (int)std::round(std::accumulate(vect.begin(), vect.end(), 0.0));
	int minnum = std::accumulate(minsize.begin(), minsize.end(), 0);
	if (sum < minnum) throw std::runtime_error("Failed to satisfy forced number of edges property");
	if (sum == minnum) return minsize;
	std::vector<int> nums;
	std::vector<double> positive_errors, negative_errors;
	for (int i=0; i<vect.size(); ++i){
		double n = vect[i];
		int nr = std::max(minsize[i], (int)std::round(n));
		nums.push_back(nr);
		positive_errors.push_back((nr - n + 1)/n);
		negative_errors.push_back(nr > minsize[i] ? (n - nr + 1)/n :
			       std::numeric_limits<double>::max());
	}
	//fit approximation to equal nedges
	int sumn = std::accumulate(nums.begin(), nums.end(), 0);
	while (sumn > sum){
		int index = std::min_element(negative_errors.begin(), negative_errors.end())-
			negative_errors.begin();
		assert(negative_errors[index] != std::numeric_limits<double>::max());
		--nums[index]; --sumn;
		if (nums[index] > minsize[index])
			negative_errors[index] = (vect[index]-nums[index]+1)/vect[index];
		else negative_errors[index] = std::numeric_limits<double>::max();
	}
	while (sumn < sum){
		int index = std::min_element(positive_errors.begin(), positive_errors.end())-
			positive_errors.begin();
		++nums[index]; ++sumn;
		positive_errors[index] = (nums[index] - vect[index] + 1)/vect[index];
	}
	return nums;
}
