#include "hmtimer.hpp"
#include <map>
#include <string>
#include <vector>
#include <iostream>

using namespace HMTimer;

namespace{

struct _tclass{
	std::map<std::string, TicToc> data;
	~_tclass(){ FinReport(); }
	bool has(std::string s){ return data.find(s) != data.end(); }
	TicToc& get(std::string s){
		auto fnd = data.find(s);
		if (fnd == data.end()){
			data[s] = TicToc(s.c_str(), 0);
		}
		return data[s];
	}
	std::vector<std::string> keys(){
		std::vector<std::string> ret;
		for (auto& v: data) ret.push_back(v.first);
		return ret;
	}
	void erase(std::string s){
		auto fnd = data.find(s);
		if (fnd != data.end()) data.erase(fnd);
	}
};
_tclass _alltimers;
}

void HMTimer::Tic(std::string s){
	if (s.size() == 0){
		for (int i=0; i<99999; ++i){
			std::string nm = "Timer" + std::to_string(i);
			if (!_alltimers.has(nm)) return Toc(nm);
		}
	} else {
		auto& tm = _alltimers.get(s);
		tm.tic();
	}
}
void HMTimer::Tic1(std::string s){
	HMTimer::Toc();
	HMTimer::Tic(s);
}

void HMTimer::Toc(std::string s){
	if (s.size() == 0){
		for (auto k: _alltimers.keys()) Toc(k);
	} else {
		if (_alltimers.has(s)){
			_alltimers.get(s).toc();
		}
	}
}

void HMTimer::Report(std::string s){
	if (s.size() == 0){
		for (auto k: _alltimers.keys()) Report(k);
	} else {
		if (_alltimers.has(s)){
			_alltimers.get(s).report();
		}
	}
}
void HMTimer::FinReport(std::string s){
	if (s.size() == 0){
		for (auto k: _alltimers.keys()) FinReport(k);
	} else {
		if (_alltimers.has(s)){
			_alltimers.get(s).fintoc();
		}
		_alltimers.erase(s);
	}
}

TicToc::TicToc(const char* _name, bool start):name(_name), is_working(start), dur(TDuration::zero()){
	if (start) tic();
}
void TicToc::init(){
	dur = TDuration::zero();
	tp = TClock::now();
}
void TicToc::tic(){
	if (!is_working){
		is_working = true;
		tp = TClock::now();
	}
}
void TicToc::toc(){
	if (is_working){
		is_working = false;
		dur += std::chrono::duration_cast<TDuration>(TClock::now() - tp);
	}
}
void TicToc::report() const{
	std::cout<<name<<":  "<<elapsed()<<" seconds"<<std::endl;
}
double TicToc::elapsed() const{
	if (!is_working) return dur.count();
	else return  (dur + std::chrono::duration_cast<TDuration>(TClock::now() - tp)).count();
}


