#ifndef NDEBUG

#include "debug_cont2d.hpp"
#include "stdarg.h"

using namespace HMCont2D;

Debug _dbg;
int Debug::tabs = 0;

void Debug::Print(const char* fmt, ...){
	std::string s;
	for (int i=0; i<tabs; ++i) s+="        ";
	s+=fmt;
	va_list va;
	const char* format = s.c_str();
	va_start(va, fmt);
	vprintf(format, va);
	va_end(va);
}

std::ostream& Debug::Cout(){
	for (int i=0; i<tabs; ++i) std::cout<<"        ";
	return std::cout;
}

void Debug::info_contour(const Contour& c){
	Cout()<<"+++ Contour at "<<&c<<". "<<c.size()<<" edges. ";
	if (c.is_closed()) std::cout<<"Closed. "<<"Area = "<<Contour::Area(c)<<std::endl;
	else std::cout<<"Open."<<std::endl;
	int i = 0;
	for (auto p: c.ordered_points()){
		Print("Point %p:  (%10.6f, %10.6f)", p, p->x, p->y);
		if (i>=0) printf("   --> Edge %p\n", c.edge(i));
		else printf("\n");
		++i;
		if (i==c.size()) i = c.is_closed()?0:-1;
	}
	Cout()<<"+++++++++++++++++++++++++++++++++++++"<<std::endl;
}

void Debug::info_ecollection(const ECollection& c){
	Cout()<<"+++ ECollection at "<<&c<<". "<<c.size()<<" edges."<<std::endl;
	for (auto e: c.data){
		Print("Points %p -- %p, (%10.6f, %10.6f) -- (%10.6f, %10.6f) --> Edge %p\n",
			e->pstart, e->pend, e->pstart->x, e->pstart->y,
			e->pend->x, e->pend->y, e.get());
	}
	Cout()<<"+++++++++++++++++++++++++++++++++++++"<<std::endl;
}

void Debug::info_tree(const ContourTree& c){
	Print("+++ ContourTree at %p. %i Contours. Area = %f.\n", &c, c.nodes.size(),
			ContourTree::Area(c));
	tabs += 1;
	for (auto r: c.roots()){
		info_tree_node(*r);
	}
	tabs -= 1;
	Cout()<<"+++++++++++++++++++++++++++++++++++++"<<std::endl;
}


void Debug::info_tree_node(const ContourTree::TreeNode& c){
	Print("+++ TreeNode at %p. Parent = %p. %i Children\n", &c, c.parent, c.children.size());
	info_contour(c);
	tabs += 1;
	for (auto c2: c.children){
		info_tree_node(*c2);
	}
	tabs -= 1;
	Cout()<<"+++++++++++++++++++++++++++++++++++++"<<std::endl;
}




#endif
