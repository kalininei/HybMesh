#ifndef CROSSGRID_CROSSGRID_H
#define CROSSGRID_CROSSGRID_H

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
