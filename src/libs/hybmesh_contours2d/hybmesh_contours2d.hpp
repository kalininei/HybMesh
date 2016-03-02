#ifndef HYBMESH_CONTOURS2D_HPP
#define HYBMESH_CONTOURS2D_HPP
#include "hmproject.h"
#include "bgeom2d.h"
#include "edges.hpp"
#include "collections.hpp"
#include "containers.hpp"
#include "contour.hpp"
#include "tree.hpp"
#include "constructor.hpp"
#include "cont_assembler.hpp"
#include "contclipping.hpp"
#include "algos.hpp"
#include "debug_cont2d.hpp"

namespace HMCont2D{


// =================== Algorithms functions
template<class C>
void SaveVtk(const C& obj, const char* fn){
	C::SaveVtk(obj, fn);
}

template<class C>
double Area(const C& obj){
	return C::Area(obj);
}

template<class C>
ScaleBase Scale01(C& obj){
	return C::Scale01(obj);
}
template<class C>
void Scale(C& obj, const ScaleBase& sc){
	C::Scale(obj, sc);
}
template<class C>
void Unscale(C& obj, const ScaleBase& sc){
	C::Unscale(obj, sc);
}


}//HMCont2D

#endif
