#ifndef NDEBUG

#include "hmdebug.hpp"
#include "stdarg.h"
#include "nan_handler.h"
#include <stdio.h>
#include <execinfo.h>
#include <signal.h>
#include <stdlib.h>
#include <unistd.h>


namespace{
struct _D{ _D(){HMDebug::get();} };
_D _d;
}


int HMDebug::tabs = 0;

void HMDebug::pre(){
	NanSignalHandler::StartCheck();
}

void HMDebug::Print(const char* fmt, ...){
	std::string s;
	for (int i=0; i<tabs; ++i) s+="        ";
	s+=fmt;
	va_list va;
	const char* format = s.c_str();
	va_start(va, fmt);
	vprintf(format, va);
	va_end(va);
}

std::ostream& HMDebug::Cout(){
	for (int i=0; i<tabs; ++i) std::cout<<"        ";
	return std::cout;
}

void HMDebug::print_trace(int d){
	void *array[d];
	size_t size;
	//get void*'s for all entries on the stack
	size = backtrace(array, d);
	backtrace_symbols_fd(array, size, STDOUT_FILENO);
}
#endif
