#ifndef HYBMESH_CALLBACK_HPP
#define HYBMESH_CALLBACK_HPP
#include <functional>

namespace HMCallback{
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

extern Fun2 to_cout2;         //callback to std::cout
extern Fun2 to_cout2_timer;  //callback to std::cout with timer
extern Fun2 silent2;          //no callback at all

struct Caller2{
	Caller2(std::string proc_name="", double proc_duration=0, Fun2 func=HMCallback::silent2);
	void reset(std::string proc_name, double proc_duration);
	void setfun(Fun2 func);

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


class Singleton2: public HMCallback::Caller2{
	Singleton2();
public:
	static Singleton2& init(std::string s, double duration);
	static Singleton2& get();

	struct Beholder{
		Beholder(HMCallback::Fun2& cb){ Singleton2::get().call = cb; }
		~Beholder(){ Singleton2::get().call = HMCallback::silent2; }
	};
	static Beholder enable(HMCallback::Fun2& cb){ return Beholder(cb); }
};

template<class Fun, class... Args>
void WithCallback(HMCallback::Fun2 cb, Fun&& f, Args &&... args){
	auto _tmp = Singleton2::enable(cb);
	f(std::forward<Args>(args)...);
	Singleton2::get().fin();
}

template<class TExecutor> class FunctionWithCallback;

template<int N>
struct TDuration{ static constexpr int value = N; };

template<typename T, class... Args>
struct GetDuration{
	static const int value = decltype( T::duration(0, std::declval<Args>()...) )::value;
};

//base class for functors with callback
class ExecutorBase{
	//used/assigned by FunctionWithCallback<>
	void init(const char* nm, double dur){ callback.reset(nm, dur); }
	void fin(){ callback.fin(); }
protected:
	//used by _run(...) procedures
	HMCallback::Caller2 callback;

	ExecutorBase(): callback("", 0, silent2){}

	void set_callback(Fun2 func){ callback.setfun(func); }
	void set_callback(HMCallback::Caller2 func){ callback = func; }
	void swap_callback(Caller2& f){ std::swap(callback, f); }

	template<class X>
	friend class FunctionWithCallback;
};

// ===== Macros which should present in a functor definition
#define HMCB_SET_DEFAULT_DURATION(X) \
	static constexpr HMCallback::TDuration<X> duration(...);
#define HMCB_SET_PROCNAME(X) \
	static constexpr const char* procname() { return X; }
//set custom duration for arguments lists.
//Should be called multiple times for functions with default arguments
//
//HMCB_DURATION(10, int, double);
//HMCB_DURATION(10, int);
//HMCB_DURATION(10);
//_run(int x=0, double y=0);
#define HMCB_SET_DURATION(X, ...) \
	static constexpr HMCallback::TDuration<X> duration(bool, ##__VA_ARGS__);

// ===== Get duration by functor type and arguments list
#define HMCB_DURATION(FUNCTOR, ...) \
	HMCallback::GetDuration<FUNCTOR,  ##__VA_ARGS__>::value

template<class TExecutor>
class FunctionWithCallback{
	template<class... Args>
	using TRet = decltype( std::declval<TExecutor>()._run(std::declval<Args>()...) );
	TExecutor exe;

	//execution with callback function reset
	template<class... Args>
	struct Beholder{
		TExecutor* e;
		static constexpr const char* nm = TExecutor::procname();
		static constexpr double dr = HMCB_DURATION(TExecutor, Args...);
		Beholder(TExecutor* _e):e(_e){ e->init(nm, dr); }
		~Beholder() { e->fin(); }
	};

	template<class... Args>
	TRet<Args...> invoke(Args&&... arg){
		Beholder<Args...> b(&exe);  //to call exe->fin() before return;
		return exe._run(std::forward<Args>(arg)...);
	}

	//execution with exeisted callback function without reset
	template<class... Args>
	struct Beholder1{
		TExecutor* e;
		Caller2* cb;
		Beholder1(Caller2& _cb, TExecutor* _e):e(_e), cb(&_cb){ e->swap_callback(_cb); }
		~Beholder1() { e->swap_callback(*cb); }
	};
	template<class... Args>
	TRet<Args...> invoke1(Caller2& cb, Args&&... arg){
		Beholder1<Args...> b(cb, &exe);  //to call exe->fin() before return;
		return exe._run(std::forward<Args>(arg)...);
	}
public:
	FunctionWithCallback(): exe(){}

	TExecutor& functor(){ return exe; }
	const TExecutor& functor() const{ return exe; }

	//call with inactive callback
	template<class... Args>
	TRet<Args...> operator()(Args&&... arg){
		exe.set_callback(silent2);
		return invoke(std::forward<Args>(arg)...);
	}

	//call with inactive callback
	template<class... Args>
	TRet<Args...> Silent(Args&&... arg){
		exe.set_callback(silent2);
		return invoke(std::forward<Args>(arg)...);
	}

	//call with cout callback
	template<class... Args>
	TRet<Args...> ToCout(Args&&... arg){
		exe.set_callback(to_cout2);
		return invoke(std::forward<Args>(arg)...);
	}

	//call with callback with timer
	template<class... Args>
	TRet<Args...> WTimer(Args&&... arg){
		exe.set_callback(to_cout2_timer);
		return invoke(std::forward<Args>(arg)...);
	}

	//call with defined callback with its initializing
	template<class TCallback, class... Args>
	TRet<Args...> WithCallback(TCallback&& cb, Args&&... arg){
		exe.set_callback(std::forward<TCallback>(cb));
		return invoke(std::forward<Args>(arg)...);
	}

	//run with existing callback without it reinitializing
	template<class... Args>
	TRet<Args...> MoveCallback(Caller2& cb, Args&&... arg){
		return invoke1(cb, std::forward<Args>(arg)...);
	}

};



}

#endif
