#ifndef HYBMESH_CALLBACK_HPP
#define HYBMESH_CALLBACK_HPP
#include <functional>

namespace HMCallback{
struct Caller1;
struct Caller2;
struct LoopCaller2;

const int OK = 0;
const int CANCEL = 1;

class Cancelled: public std::runtime_error{
	static std::string errmsg(std::string proc){
		return "User interrupt at " + proc;
	}
public:
	Cancelled(std::string proc) noexcept: std::runtime_error(errmsg(proc)){}
};

//Function of the following arguments:
//	(global procedure name, local procedure name, global percentage, local percentage)
//percentages are doubles in [0, 1] range. If it is less then, then percentages of procedures is not tracking.
//normally returns OK. Should return CANCEL for cancellation require
typedef std::function<int(const char*, const char*, double, double)> Fun2;

//same with single progress bar
typedef std::function<int(const char*, double)> Fun1;

extern Fun1 to_cout1;  //callback to std::cout
extern Fun2 to_cout2;
extern Fun1 silent1;   //no callback at all
extern Fun2 silent2;

struct Caller1{
	void move(double progress);
	void step(double progress);
	void fin();
};

struct Caller2{
	Caller2(std::string proc_name, double proc_duration, Fun2 func);
	void reset(std::string proc_name, double proc_duration);

	//Moving/steppin progress
	void silent_move_now(double progress, std::string subproc_name, double subproc_duration=-1);
	void silent_step_now(double progress, std::string subproc_name, double subproc_duration=-1);
	void silent_move_after(double progress, std::string subproc_name, double subproc_duration=-1, double subproc_move_after=0);
	void silent_step_after(double progress, std::string subproc_name, double subproc_duration=-1, double subproc_step_after=0);

	void silent_subprocess_move_now(double progress);
	void silent_subprocess_step_now(double progress);
	void silent_subprocess_move_after(double progress);
	void silent_subprocess_step_after(double progress);

	void move_now(double progress, std::string subproc_name, double subproc_duration=-1);
	void step_now(double progress, std::string subproc_name, double subproc_duration=-1);
	void move_after(double progress, std::string subproc_name, double subproc_duration=-1, double subproc_move_after=0);
	void step_after(double progress, std::string subproc_name, double subproc_duration=-1, double subproc_step_after=0);

	void subprocess_move_now(double progress);
	void subprocess_step_now(double progress);
	void subprocess_move_after(double progress);
	void subprocess_step_after(double progress);

	Caller2 subrange(double parent_duration, double child_duration);
	LoopCaller2 looper(double parent_duration, double nloops, std::string caption);

	//output 100%
	void fin();
	void subprocess_fin();
protected:
	//currant state data
	std::string name1, name2;
	double prog1, prog2;
	double dur1, dur2;
	double step_before1, step_before2;

	//construct data
	Fun2 call; 

	void flush();
	void new_sub_process(std::string subproc_name, double subproc_duration);
	void before1(bool);
	void before2(bool);
};

struct LoopCaller2: public Caller2{
	LoopCaller2(std::string proc_name, std::string loop_caption, double nloops, Fun2 call);
	void new_iteration(double subproc_duration=-1);
private:
	std::string caption;
	int iter;
};

}



#endif
