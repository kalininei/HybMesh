#ifndef HMPROJECT_TIMER_HPP
#define HMPROJECT_TIMER_HPP

#include <chrono>
#include <string>

namespace HMTimer{

//Timer
struct TicToc{
	//usage:
	//TicToc timer();
	//......
	//timer.fintoc();
	TicToc(const char* name = "Time duration", bool start=true);
	~TicToc(){}
	void init();
	void tic();
	void toc();
	void report() const;
	void fintoc(){ toc(); report(); }
	double elapsed() const;
private:
	typedef std::chrono::high_resolution_clock TClock;
	typedef TClock::time_point TTimePoint;
	typedef std::chrono::duration<double> TDuration;
	const char* name;
	bool is_working;
	TTimePoint tp;
	TDuration dur;
};

void Tic(std::string s="");        //starts a timer with string id. default is "timer1,2,3 etc"
void Tic1(std::string s="");       //starts a timer and stops all others
void Toc(std::string s="");        //stops it or all if ""
void Report(std::string s="");     //reports timer with string id or all
void FinReport(std::string s="");  //reports and deletes timer with id or all

}

#endif
