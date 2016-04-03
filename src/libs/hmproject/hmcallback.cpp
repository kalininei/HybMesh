#include "hmcallback.hpp"
#include <iostream>

using namespace HMCallback;

namespace{

//callbacks
//default callback function
int default_callback2(const char* proc, const char* subproc, double percent, double subpercent){
	std::cout<<proc;
	if (percent>=0) std::cout<<" - "<<int(100*percent)<<"%";
	std::cout<<";";
	if (subproc!=0){
		std::cout<<" ("<<subproc;
		if (subpercent>=0) std::cout<<" - "<<int(100*subpercent)<<"%";
		std::cout<<")";
	}
	std::cout<<std::endl;
	return OK;
}

int silent_callback2(const char* proc, const char* subproc, double percent, double subpercent){
	return OK;
}

}

Fun2 HMCallback::to_cout2 = default_callback2;
Fun2 HMCallback::silent2 = silent_callback2;

Caller2::Caller2(std::string proc_name, double proc_duration, Fun2 func){
	call = func;
	reset(proc_name, proc_duration);
}

void Caller2::reset(std::string proc_name, double proc_duration){
	before1(false); before2(false);
	name1 = proc_name;
	dur1 = proc_duration;
	prog1 = 0;
	new_sub_process("", -1);
}

void Caller2::new_sub_process(std::string subproc_name, double subproc_duration){
	name2 = subproc_name;
	dur2 = subproc_duration;
	prog2 = 0;
	step_before2 = 0;
}

void Caller2::flush(){
	double w1 = (dur1>0) ? prog1/dur1 : -1;
	double w2 = (dur2>0) ? prog2/dur2 : -1;
	int res = call(name1.c_str(), name2.c_str(), w1, w2);
	if (res == CANCEL) throw Cancelled(name1);
}

void Caller2::before1(bool use=true){
	if (use) prog1 += step_before1;
	step_before1 = 0;
}

void Caller2::before2(bool use=true){
	if (use) prog2 += step_before2;
	step_before2 = 0;
}

void Caller2::silent_move_now(double progress, std::string subproc_name, double subproc_duration){
	before1(false);
	prog1 = progress;
	new_sub_process(subproc_name, subproc_duration);
}
void Caller2::silent_step_now(double progress, std::string subproc_name, double subproc_duration){
	before1(true);
	prog1 += progress;
	new_sub_process(subproc_name, subproc_duration);
}

void Caller2::silent_move_after(double progress, std::string subproc_name, double subproc_duration, double subproc_progress){
	before1(true);
	new_sub_process(subproc_name, subproc_duration);
	step_before1 = progress-prog1;
	step_before2 = subproc_progress;
}
void Caller2::silent_step_after(double progress, std::string subproc_name, double subproc_duration, double subproc_progress){
	before1(true);
	new_sub_process(subproc_name, subproc_duration);
	step_before1 = progress;
	step_before2 = subproc_progress;
}
void Caller2::silent_subprocess_move_now(double progress){
	before2(false);
	prog2 = progress;
}
void Caller2::silent_subprocess_step_now(double progress){
	before2(true);
	prog2 += progress;
}
void Caller2::silent_subprocess_move_after(double progress){
	before2(true);
	step_before2 = progress - prog2;
}
void Caller2::silent_subprocess_step_after(double progress){
	before2(true);
	step_before2 = progress;
}

void Caller2::move_now(double progress, std::string subproc_name, double subproc_duration){
	silent_move_now(progress, subproc_name, subproc_duration);
	flush();
}
void Caller2::step_now(double progress, std::string subproc_name, double subproc_duration){
	silent_step_now(progress, subproc_name, subproc_duration);
	flush();
}
void Caller2::move_after(double progress, std::string subproc_name, double subproc_duration, double subproc_move_after){
	silent_move_after(progress, subproc_name, subproc_duration, subproc_move_after);
	flush();
}
void Caller2::step_after(double progress, std::string subproc_name, double subproc_duration, double subproc_step_after){
	silent_step_after(progress, subproc_name, subproc_duration, subproc_step_after);
	flush();
}

void Caller2::subprocess_move_now(double progress){
	silent_subprocess_move_now(progress);
	flush();
}
void Caller2::subprocess_step_now(double progress){
	silent_subprocess_step_now(progress);
	flush();
}
void Caller2::subprocess_move_after(double progress){
	silent_subprocess_move_after(progress);
	flush();
}
void Caller2::subprocess_step_after(double progress){
	silent_subprocess_step_after(progress);
	flush();
}


void Caller2::fin(){
	before1(false); before2(false);
	prog1 = dur1;
	name2 = "Done";
	prog2 = 1;
	dur2 = 1;
	flush();
}
void Caller2::subprocess_fin(){
	before2(false);
	prog2 = dur2;
	flush();
}

Caller2 Caller2::subrange(double parent_duration, double child_duration){
	before1(true);
	double parent_start = prog1;
	prog1 += parent_duration;
	Fun2& f=call;
	auto sf = [&f, parent_start, parent_duration](const char* s1, const char* s2, double p1, double p2)->int{
		return f(s1, s2, parent_start + parent_duration*p1, p2);
	};
	return Caller2(name1, child_duration, sf);
}

LoopCaller2 Caller2::looper(double parent_duration, double nloops, std::string caption){
	before1(true);
	double parent_start = prog1;
	prog1 += parent_duration;
	Fun2& f=call;
	auto sf = [&f, parent_start, parent_duration](const char* s1, const char* s2, double p1, double p2)->int{
		return f(s1, s2, parent_start + parent_duration*p1, p2);
	};
	return LoopCaller2(name1, caption, nloops, sf);
}

LoopCaller2::LoopCaller2(std::string proc_name, std::string loop_caption, double nloops, Fun2 call):
	Caller2(proc_name, nloops, call), caption(loop_caption), iter(0){}


void LoopCaller2::new_iteration(double subproc_duration){
	std::string nm = caption + ": " + std::to_string(iter+1) + "/" + std::to_string((int)dur1);
	move_now(iter, nm, subproc_duration);
	++iter;
}
