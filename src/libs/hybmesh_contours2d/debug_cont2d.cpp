#ifndef NDEBUG

#include "debug_cont2d.hpp"
#include "constructor.hpp"
#include <fstream>

using namespace HMCont2D;

void Debug::info_contour(const Contour& c){
	Cout()<<"+++ Contour at "<<&c<<". "<<c.size()<<" edges. ";
	if (c.is_closed()) Cout()<<"Closed. "<<"Area = "<<Contour::Area(c)<<std::endl;
	else Cout()<<"Open."<<std::endl;
	Cout()<<"+++ Length = "<<c.length()<<std::endl;
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

void Debug::geogebra_ecollection(const ECollection& ecol){
	auto etree = HMCont2D::Assembler::ETree(ecol);
	geogebra_etree(etree);
}

void Debug::vtk_contours(const ShpVector<HMCont2D::Contour>& c, const char* fn){
	HMCont2D::ECollection ecol;
	for (auto& x: c) ecol.Unite(*x);
	HMCont2D::ECollection::SaveVtk(ecol, fn);
}

void Debug::vtk_contour(const ECollection& col, std::map<Point*, double>& pdata, double defval, const char* fn){
	HMCont2D::ECollection::SaveVtk(col, fn);
	std::ofstream fs(fn, std::ios_base::app);
	auto ap = col.all_points();
	fs<<"POINT_DATA "<<ap.size()<<std::endl;
	fs<<"SCALARS debug_data float 1"<<std::endl;
	fs<<"LOOKUP_TABLE default"<<std::endl;
	for (auto p: ap){
		double val;
		if (pdata.find(p) != pdata.end()) val = pdata[p];
		else val = defval;
		fs<<val<<std::endl;
	}
}

void Debug::vtk_contour(const ECollection& col, std::map<Edge*, double>& edata, double defval, const char* fn){
	HMCont2D::ECollection::SaveVtk(col, fn);
	std::ofstream fs(fn, std::ios_base::app);
	fs<<"CELL_DATA "<<col.size()<<std::endl;
	fs<<"SCALARS debug_data float 1"<<std::endl;
	fs<<"LOOKUP_TABLE default"<<std::endl;
	for (auto e: col.data){
		double val;
		if (edata.find(e.get()) != edata.end()) val = edata[e.get()];
		else val = defval;
		fs<<val<<std::endl;
	}
}


#endif
