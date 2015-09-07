#ifndef CROSSGRID_CROSSGRID_H
#define CROSSGRID_CROSSGRID_H

namespace CrossGridCallback{

// === set default callbacks
//crossgrid callback pointer is the function of the following arguments:
//	(global procedure name, local procedure name, global percentage, local percentage)
//percentages are doubles in [0, 1] range. If it is less then, then percentages of procedures is not tracking.
//normally returns OK. Should return CANCEL for cancellation require
const int OK = 0;
const int CANCEL = 1;
typedef int (*func)(const char*, const char*, double, double);

//set global callbacks procedures
func to_cout();  //callback to std::cout
func silent();   //no callback at all

}

//base grid class
class Grid{
public:
	virtual ~Grid(){}
};

//base contour class
class Cont{
public:
	virtual ~Cont(){}
};

extern "C"{

// === testing
void crossgrid_internal_tests();

}; //extern C


#endif
