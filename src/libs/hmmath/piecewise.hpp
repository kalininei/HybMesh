#ifndef HMMATH_PIECEWISE_HPP
#define HMMATH_PIECEWISE_HPP
#include "hmproject.h"
#include "addalgo.hpp"

namespace HMMath{


template<class A>
class Piecewise{
protected:
	aa::DoubleMap<A> _data;
	virtual double _section_integral(typename std::map<double, A>::const_iterator, double x1, double x2) const = 0;
public:
	struct Error: public std::runtime_error{
		Error(): std::runtime_error("Piecewise operation error"){};
	};

	static constexpr double inf() { return std::numeric_limits<double>::max()/10.0; }
	static constexpr double eps() { return geps; }

	//information
	double xstart() const { return _data.begin()->first; }
	double xend() const { return _data.rbegin()->first; }
	int n() const { return (_data.size() == 0) ? 0 : _data.size() - 1; }

	//calculations
	A& value(double x);
	const A& value(double x) const;
	double Integral(double x1, double x2) const;
	double Average(double x1, double x2) const { return Integral(x1, x2)/fabs(x2 - x1); }
	
	//data manipulations
	void add_section(double x1, double x2, const A& v);
	void clear(){ _data.clear(); }

	Piecewise():_data(eps()){}
};

template<class A>
A& Piecewise<A>::value(double x){
	if (n() == 0 || x<xstart()-eps() || x>xend()+eps()) throw Error();
	auto fnd = _data.upper_bound(x);
	return std::prev(fnd)->second;
}
template<class A>
const A& Piecewise<A>::value(double x) const{
	if (n() == 0 || x<xstart()-eps() || x>xend()+eps()) throw Error();
	auto fnd = _data.upper_bound(x);
	return std::prev(fnd)->second;
}

template<class A>
void Piecewise<A>::add_section(double x1, double x2, const A& v){
	if (fabs(x2 - x1) < eps()) return;
	else if (x2<x1) std::swap(x1, x2);
	if (n() == 0){
		_data.emplace(x1, v);
		_data.emplace(x2, v);
	} else {
		A v2 = value(x2);
		auto em1 = _data.emplace(x1, v);
		auto em2 = _data.emplace(x2, v);
		if (!em1.second) em1.first->second = v;
		if (em2.second){
			auto f2 = std::prev(em2.first);
			em2.first->second = v2;
		}
		auto it = std::next(em1.first);
		while (it != em2.first && it != _data.end()){
			it = _data.erase(it);
		}
	}
}

template<class A>
double Piecewise<A>::Integral(double x1, double x2) const{
	if (fabs(x2 - x1) < eps()) return 0;
	else if (x2<x1) std::swap(x1, x2);
	if (n() == 0 || x1<xstart()-eps() || x2>xend()+eps()) throw Error();

	auto fnd1 = _data.upper_bound(x1);
	auto fnd2 = _data.upper_bound(x2);
	--fnd1; --fnd2;
	if (fnd1 == fnd2) return _section_integral(fnd1, x1, x2);
	else {
		double ret = 0;
		ret += _section_integral(fnd1++, x1, inf());
		ret += _section_integral(fnd2, -inf(), x2);
		while (fnd1 != fnd2){
			ret += _section_integral(fnd1++, -inf(), inf());
		}
		return ret;
	}
}

// ============================================ Step Function
template<class A>
class StepFunction: public Piecewise<A>{
	double _section_integral(typename std::map<double, A>::const_iterator it, double x1, double x2) const override{
		double a = it->first;
		double b = std::next(it)->first;
		x1 = std::max(x1, a);
		x2 = std::min(x2, b);
		return (x2 - x1) * it->second;
	}
public:
	using Piecewise<A>::inf;
	using Piecewise<A>::value;
	using Piecewise<A>::add_section;
	StepFunction(A init): Piecewise<A>(){
		add_section(-inf(), inf(), init);
	}
	A operator()(double x) const { return value(x); }

	template<class B>
	friend std::ostream& operator<<(std::ostream& str, const StepFunction<B>& f);
};

template<class A>
std::ostream& operator<<(std::ostream& str, const StepFunction<A>& f){
	if (f.n() == 0) return str;
	str<<"Step Function:";
	auto it = f._data.begin();
	for (int i=0; i<f.n(); ++i){
		str<<std::endl;
		auto it2 = std::next(it);
		double x1 = it->first;
		double x2 = it2->first;
		A v = it->second;
		std::string s1 = (x1 == -StepFunction<A>::inf()) ? "-inf" : std::to_string(x1);
		std::string s2 = (x2 ==  StepFunction<A>::inf()) ? "inf" : std::to_string(x2);
		str<<"["<<s1<<", "<<s2<<"] --> "<<v;
		++it;
	}
	return str;
}

struct LinFunc{
	double a, b;
	double operator()(double x) const {return a*x + b; }
	double integral(double x1, double x2) const{
		return a/2.0*(x2*x2 - x1*x1) + b*(x2 - x1);
	}
	LinFunc(double x1, double v1, double x2, double v2){
		a = (v2 - v1)/(x2 - x1);
		b = v1 - (v2 - v1)*x1/(x2 - x1);
	}
};

class LinearPiecewise: public Piecewise<LinFunc>{
	double _section_integral(typename std::map<double, LinFunc>::const_iterator it, double x1, double x2) const override{
		double a = it->first;
		double b = std::next(it)->first;
		x1 = std::max(x1, a);
		x2 = std::min(x2, b);
		return it->second.integral(x1, x2);
	}
public:
	using Piecewise<LinFunc>::inf;
	using Piecewise<LinFunc>::value;
	void add_point(double x, double v){
		if (n() == 0){
			add_section(-inf(), inf(), LinFunc(-1, v, 1, v));
			add_section(-inf(), x, LinFunc(-1, v, 1, v));
			add_section(x, inf(), LinFunc(-1, v, 1, v));
		} else {
			auto f2 = _data.upper_bound(x);
			auto f1 = std::prev(f2);
			assert(!ISEQ(f1->first, x) && !ISEQ(f2->first, x));
			//before
			if (f1 == _data.begin()){
				auto fn = std::next(f1);
				f1->second = LinFunc(x, v, fn->first, fn->second(fn->first));
			} else {
				LinFunc f(f1->first, f1->second(f1->first), x, v);
				f1->second = f;
				if (f1 == ++_data.begin()) _data.begin()->second = f;
			}
			//after
			if (f2 == std::prev(_data.end())){
				auto fp = std::prev(f2);
				add_section(x, inf(), fp->second);
			} else{
				LinFunc f(x, v, f2->first, f2->second(f2->first));
				add_section(x, f2->first, f);
				if (f2 == std::prev(_data.end(), 2)) f2->second = f;
			}
		}
	}
	std::map<double, double> to_map() const{
		std::map<double, double> ret;
		for (auto& it: _data){
			if (fabs(it.first) != inf()){
				double val = it.second(it.first);
				ret.emplace(it.first, val);
			}
		}
		return ret;
	}
	double operator()(double x) const { return value(x)(x); }
	void operator*=(double x){ for (auto& it: _data){ it.second.a*=x; it.second.b*=x; }; }

	friend std::ostream& operator<<(std::ostream& str, const LinearPiecewise& f);
};

inline std::ostream& operator<<(std::ostream& str, const LinearPiecewise& f){
	if (f.n() == 0) return str;
	str<<"Linear Piecewise Function:";
	auto it = std::next(f._data.begin());
	if (f.n() == 2){
		str<<std::endl<<it->first<<" --> "<<it->second(it->first);
	} else for (int i=1; i<f.n(); ++i){
		str<<std::endl<<it->first<<" --> "<<it->second(it->first);
		++it;
	}
	return str;
}

}
#endif
