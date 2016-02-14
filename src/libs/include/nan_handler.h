//http://www.yolinux.com/TUTORIALS/C++Signals.html
#ifndef NAN_HANDLER_H
#define NAN_HANDLER_H

//class for catching nan signals in unix systems.
#ifndef WIN32

#include <stdexcept>
#include <fenv.h>
#include <signal.h>
#include <errno.h>

class NanSignalHandler{
	//true if signal was catched
	bool mbGotExitSignal;
	//initialization
	void setupSignalHandlers(){
		if (signal((int) SIGFPE, NanSignalHandler::exitSignalHandler) == SIG_ERR){
			throw; 
		}
	}

	//exception which will be generated on floating point error
	struct E_SignalException : public std::runtime_error{
		E_SignalException() noexcept
			:std::runtime_error("NanSignalHandler: Floating point error"){}
	};

	/**
	* Default Contructor.
	*/
	NanSignalHandler(){
		try{
			mbGotExitSignal = false;
			feenableexcept(FE_DIVBYZERO|FE_INVALID|FE_OVERFLOW);
			setupSignalHandlers();
		} catch (...){
			throw std::runtime_error("Error during NanSignalHandler initialisation");
		}
	};

	/**
	* Returns the bool flag indicating whether we received an exit signal
	* @return Flag indicating shutdown of program
	*/
	static bool gotExitSignal(){
		return Instance().mbGotExitSignal;
	}
	
	/**
	* Sets the bool flag indicating whether we received an exit signal
	*/
	static void setExitSignal(bool _bExitSignal){
		Instance().mbGotExitSignal = _bExitSignal;
	}
	
	/**
	* Sets exit signal to true.
	* @param[in] _ignored Not used but required by function prototype
	*                     to match required handler.
	*/
	static void exitSignalHandler(int _ignored){
		Instance().mbGotExitSignal = true;
		throw E_SignalException();
	}
	
	//get singleton instance
	static NanSignalHandler& Instance(){
		static NanSignalHandler hand;
		return hand;
	}
public:
	NanSignalHandler(const NanSignalHandler&) = delete;
	NanSignalHandler& operator=(const NanSignalHandler& ) = delete;

	static void StartCheck(){
		Instance();
	}
};

#else
class NanSignalHandler{
public:
	void StartCheck(){};
};  //empty class for windows platform
#endif


#endif
          


















