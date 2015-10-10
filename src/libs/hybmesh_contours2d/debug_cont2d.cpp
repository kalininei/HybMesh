#ifndef NDEBUG

#include "debug_cont2d.hpp"
#include "stdarg.h"
#include "constructor.hpp"
#include "nan_handler.h"

using namespace HMCont2D;

Debug HMCont2D::_dbg;
int Debug::tabs = 0;

void Debug::pre(){
	NanSignalHandler::StartCheck();
}

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

#include <time.h>
std::string random_id(int len){
	static bool k = true;
	if (k){
		srand (time(NULL));
		k = false;
	}
	std::string ret;
	for (int i=0; i<len; ++i) ret += char(97 + rand()%26);
	return ret;
}

void Debug::geogebra_contour(const Contour& c){
	vector<Point*> op = c.ordered_points();
	std::string id = random_id(3);
	for (int i=0; i<op.size(); ++i){
		if (i == op.size() - 1 && c.is_closed()) break;
		printf("P%s_{%i} = (%10.8f, %10.8f)\n", id.c_str(), i+1, op[i]->x, op[i]->y);
	}
	if (op.size() == 2){
		std::cout<<"L"+id<<" = Segment[";
		std::cout<<"P"+id+"_{1}, ";
		std::cout<<"P"+id+"_{2}]"<<std::endl;
	} else if (!c.is_closed()){
		std::cout<<"L"+id<<" = Polyline[";
		for (int i=0; i<op.size(); ++i){
			if (i>0) std::cout<<", ";
			std::cout<<"P"+id+"_{"<<i+1<<"}";
		}
		std::cout<<"]"<<std::endl;
	} else {
		std::cout<<"L"+id<<" = Polygon[";
		for (int i=0; i<op.size()-1; ++i){
			if (i>0) std::cout<<", ";
			std::cout<<"P"+id+"_{"<<i+1<<"}";
		}
		std::cout<<"]"<<std::endl;
	}
}

void Debug::geogebra_tree(const ContourTree& c){
	for (auto n: c.nodes) geogebra_contour(*n);
}

void Debug::geogebra_box(const BoundingBox& c){
	return geogebra_contour(HMCont2D::Constructor::ContourFromPoints(c.FourPoints(), true));
}

void Debug::geogebra_etree(const ExtendedTree& c){
	geogebra_tree(c);
	for (auto n: c.open_contours) geogebra_contour(*n);
}



#endif
