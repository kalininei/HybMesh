#ifndef NDEBUG

#include "debug2d.hpp"
#include "constructor.hpp"
#include <fstream>
#include "vtk_export2d.hpp"

using namespace HM2D;

void Debug::info_contour(const EdgeData& c){
	Cout()<<"+++ EdgeData at "<<&c<<". "<<c.size()<<" edges. ";
	if (Contour::IsContour(c)){
		if (Contour::IsClosed(c)) Cout()<<"Closed contour. "<<"Area = "<<Contour::Area(c)<<std::endl;
		else Cout()<<"Open contour"<<std::endl;
	} else {
		Cout()<<"Shattered edges."<<std::endl;
	}
	Cout()<<"+++ Length = "<<Length(c)<<std::endl;
	int i = 0;
	for (auto e: c){
		Print("%i: Points %p -- %p, (%10.6f, %10.6f) -- (%10.6f, %10.6f) --> Edge %p\n",
			i,
			e->first().get(), e->last().get(), e->first()->x, e->first()->y,
			e->last()->x, e->last()->y, e.get());
		++i;
	}
	Cout()<<"+++++++++++++++++++++++++++++++++++++"<<std::endl;
}

void Debug::info_tree(const Contour::Tree& c){
	Print("+++ ContourTree at %p. %i Contours. Area = %f.\n", &c, c.nodes.size(),
			c.area());
	tabs += 1;
	for (auto r: c.roots()){
		info_tree_node(*r);
	}
	tabs -= 1;
	Cout()<<"+++++++++++++++++++++++++++++++++++++"<<std::endl;
}


void Debug::info_tree_node(const Contour::Tree::TNode& c){
	Print("+++ TreeNode at %p. Parent = %p. %i Children\n", &c, c.parent.lock().get(), c.children.size());
	info_contour(c.contour);
	tabs += 1;
	for (auto c2: c.children){
		info_tree_node(*c2.lock());
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

void Debug::geogebra_contour(const EdgeData& c){
	auto op = Contour::OrderedPoints(c);
	std::string id = random_id(3);
	for (int i=0; i<op.size(); ++i){
		if (i == op.size() - 1 && Contour::IsClosed(c)) break;
		printf("P%s_{%i} = (%10.8f, %10.8f)\n", id.c_str(), i+1, op[i]->x, op[i]->y);
	}
	if (op.size() == 2){
		std::cout<<"L"+id<<" = Segment[";
		std::cout<<"P"+id+"_{1}, ";
		std::cout<<"P"+id+"_{2}]"<<std::endl;
	} else if (!Contour::IsClosed(c)){
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


void Debug::geogebra_tree(const Contour::Tree& c){
	for (auto n: c.nodes) geogebra_contour(n->contour);
}
void Debug::geogebra_ecollection(const EdgeData& ecol){
	auto etree = HM2D::Contour::Tree::Assemble(ecol);
	geogebra_tree(etree);
}


double Debug::hash(const EdgeData& ecol){
	double sum;
	int i=0;
	for (auto e: ecol){
		auto vert1 = e->first();
		auto vert2 = e->last();
		sum += 0.333*sin(i*vert1->x+3);
		sum -= 0.333*cos(i*vert1->y+4);
		sum += 0.333*sin(i*vert2->x+5);
		sum -= 0.333*cos(i*vert2->y+6);
		++i;
	}
	return sum;
}

VertexData Debug::allvertices(const EdgeData& c){
	VertexData ret;
	for (auto e: c){
		ret.push_back(e->first());
		ret.push_back(e->last());
	}
	return aa::no_duplicates(ret);
}

void Debug::save_edges_vtk(const EdgeData& c){
	auto av = allvertices(c);
	aa::RestoreIds<EdgeData> r1(c);
	aa::RestoreIds<VertexData> r2(av);

	HM2D::Export::ContourVTK(c, "_dbgout.vtk");
}

void Debug::save_edges_vtk(const EdgeData& c, const std::map<Point*, double>& data){
	save_edges_vtk(c);
	auto av = allvertices(c);
	vector<float> values(av.size(), -1);
	for (int i=0; i<values.size(); ++i){
		auto fnd = data.find(av[i].get());
		if (fnd == data.end()) continue;
		values[i] = fnd->second;
	}
	std::ofstream fs("_dbgout.vtk", std::ios::app);
	fs<<"POINT_DATA "<<av.size()<<std::endl;
	fs<<"SCALARS mapped_values float 1"<<std::endl;
	fs<<"LOOKUP_TABLE default"<<std::endl;
	for (auto x: values) fs<<x<<std::endl;
}

void Debug::save_vertices_vtk(const VertexData& c){
	aa::RestoreIds<VertexData> r1(c);

	HM2D::Export::VerticesVTK(c, "_dbgout.vtk");
}

#endif
