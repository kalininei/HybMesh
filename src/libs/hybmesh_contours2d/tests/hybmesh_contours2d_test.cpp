#include "c_hybmesh_contours2d.h"
#include "hybmesh_contours2d.hpp"
#include <type_traits>


using namespace HMCont2D;

int FAILED_CHECKS = 0;

void add_check(bool ex, std::string info){
	if (info.size()==0){
		std::cout<<"\tunknown check: ";
	} else{
		std::cout<<"\t"<<info;
	}
	if (ex){
		std::cout<<": True"<<std::endl;
	} else {
		++FAILED_CHECKS;
		std::cout<<": False <<<<<<<<<<<<<<<<<<<"<<std::endl;
	}
};

void test1(){
	std::cout<<"Offset closed polygon"<<std::endl;
	auto c6 = Constructor::Circle(6, 1.0, Point(0,0));
	auto oret = Offset(c6, 0.2, OffsetTp::CLOSED_POLY);
	add_check(fabs(Area(oret) - 3.923) < 0.01, "6-sided closed polygon offsetting");
}

void test2(){
	std::cout<<"Deep copy procedures"<<std::endl;
	//edges and points pool
	auto c = Constructor::Circle(3, 1, Point(5,4));
	//copy only edges
	ECollection ne;
	ECollection::DeepCopy(c, ne);
	try{
		for (auto e=ne.begin(); e!=ne.end(); ++e){
			//new edges are not in c edges
			if (c.contains((*e).get())) throw 1;
			//points in new edges are all in e.pdata
			if (!c.pdata.contains((*e)->pstart)) throw 1;
			if (!c.pdata.contains((*e)->pend)) throw 1;
		}
		add_check(true, "Deep copy of edges");
	} catch (int) {add_check(false, "Deep copy of edges");}

	//copy only points
	PCollection pc;
	PCollection::DeepCopy(c.pdata, pc);
	try{
		for (auto p=pc.begin(); p!=pc.end(); ++p){
			//new points are not in c.pdata
			if (c.pdata.contains((*p).get())) throw 1;
		}
		add_check(true, "Deep copy of points");
	} catch (int) { add_check(false, "Deep copy of points"); }

	//copy whole container
	Container<Contour> nc;
	Container<Contour>::DeepCopy(c, nc);
	try{
		for (auto e=nc.begin(); e!=nc.end(); ++e){
			//edges of nc are not in c edges
			if (c.contains((*e).get())) throw 1;
			//points of nc are not in in e.pdata
			if (c.pdata.contains((*e)->pstart)) throw 1;
			if (c.pdata.contains((*e)->pend)) throw 1;
			//points of nc are in nc
			if (!nc.pdata.contains((*e)->pstart)) throw 1;
			if (!nc.pdata.contains((*e)->pend)) throw 1;
		}
		add_check(true, "Deep copy of a whole container");
	} catch (int) { add_check(false, "Deep copy of a whole container"); }
	
	//copy of a whole container from edge collection
	Container<ECollection> nc2;
	Container<ECollection>::DeepCopy(ne, nc2);
	try{
		for (auto e=nc2.begin(); e!=nc2.end(); ++e){
			//edges of nc2 are not in ne edges
			if (ne.contains((*e).get())) throw 1;
		}
		for (auto e=ne.begin(); e!=ne.end(); ++e){
			//points of ne edges are not in nc2.pdata
			if (nc2.pdata.contains((*e)->pstart)) throw 1;
			if (nc2.pdata.contains((*e)->pend)) throw 1;
		}
		add_check(true, "Deep copy of a whole container from edge collection");
	} catch (int) { add_check(false, "Deep copy of a whole container from edge collection"); }
}


void test3(){
	std::cout<<"Partition algorithm"<<std::endl;
	auto c12 = Constructor::Circle(12, 10, Point(5,4));
	auto c6 = Constructor::ContourFromPoints({0.0,0.0, 0.2,0, 0.33,0, 0.7,0.2, 0.9,0.6, 1.0,0.0 });
	double len12 = c12.length();
	double len6 = c6.length();
	//add_check(fabs(len12 - 62) - 0.5, "12-sided contour length");
	//1) ignore all algorithm
	auto r1 = Contour::Partition(len12/5.0, c12, PartitionTp::IGNORE_ALL);
	add_check(fabs(r1.length()-57.5067)<1e-3, "ignore all");
	
	//2) keep all
	auto r2 = Contour::Partition(len12/5.0, c12, PartitionTp::KEEP_ALL);
	add_check(fabs(r2.length()-len12)<1e-3 && r2.size() == c12.size(), "keep all for closed contour");

	PCollection pc3;
	auto r3 = Contour::Partition(0.11, c6, pc3, PartitionTp::KEEP_ALL);
	add_check(pc3.size() == 12 && fabs(r3.length() - len6)<1e-12, "keep all for a polyline");
	
	//3) keep shape
	PCollection pc4;
	auto r4 = Contour::Partition(1.0, c6, pc4, PartitionTp::KEEP_SHAPE);
	add_check(r4.size() == 4 && fabs(r4.length() -len6)<1e-12 && pc4.size() == 0, 
				"keep shape coarse");
	auto r5 = Contour::Partition(0.1, c6, pc4, PartitionTp::KEEP_SHAPE);
	add_check(fabs(r5.length() - len6)<1e-12 && pc4.size() == r5.size() + 1 - 5, "keep shape fine");

}

void test4(){
	std::cout<<"Partition direction"<<std::endl;
	//ignore all
	//auto c10 = Constructor::Circle(10, 1, Point(-5,4));
	//auto r1 = Contour::Partition(100, c10, PartitionTp::IGNORE_ALL);
	//add_check(r1.size() == 3 && Area(r1) > 0, "positive closed path, ignore all");
	//c10.Reverse();
	//auto r2 = Contour::Partition(0.1, c10, PartitionTp::IGNORE_ALL);
	//add_check(r2.size() == std::lround(c10.length()/0.1) && Area(r2) < 0, "negative closed path, ignore all");
	
	//keep shape
	auto in4 = Constructor::Circle(10, 0.2, Point(2,3));
	in4.Reverse();
	auto out4 = Contour::Partition(0.1, in4, PartitionTp::KEEP_SHAPE);
	add_check(fabs(Area(in4) - Area(out4))<1e-8, "negative closed path, keep shape");
	std::cout<<Area(in4)<<",  "<<Area(out4)<<std::endl;
	SaveVtk(in4, "1.vtk");
	SaveVtk(out4, "2.vtk");

}

void test5(){
	std::cout<<"Tree structure"<<std::endl;
	ECollection collection;

	//1) closed contour
	auto c1 = Constructor::Circle(13, 5, Point(0, 0));
	c1.Reverse();
	ECollection::DeepCopy(c1, collection);
	auto etree1 = ExtendedTree::Assemble(collection);
	add_check(fabs(Area(etree1) + Area(c1))<1e-12, "1 closed contour");

	//2) + inside contour
	auto c2 = Constructor::Circle(20, 0.3, Point(1,1));
	ECollection::DeepCopy(c2, collection);
	auto etree2 = ExtendedTree::Assemble(collection);
	add_check(fabs(Area(etree2) - (-Area(c1)-Area(c2)))<1e-8, "2 closed");

	
	//3) + sibling contour
	auto c3 = Constructor::Circle(3, 0.2, Point(-2.3, -3));
	ECollection::DeepCopy(c3, collection);
	auto etree3 = ExtendedTree::Assemble(collection);
	add_check(fabs(Area(etree3) - (-Area(c1)-Area(c2)-Area(c3)))<1e-8, "3 closed");

	//4) + open path
	auto c4 = Constructor::ContourFromPoints({2,2, 3,3, 2.5,3.5});
	ECollection::DeepCopy(c4, collection);
	auto etree4 = ExtendedTree::Assemble(collection);
	double area4=(-Area(c1)-Area(c2)-Area(c3));
	double len4 = c1.length() + c2.length() + c3.length() + c4.length();
	add_check(fabs(Area(etree4) - area4)<1e-12 && fabs(etree4.length() - len4)<1e-12, "3 closed, 1 open");

	//5) + big contour
	auto c5 = Constructor::Circle(100, 30, Point(1,1));
	ECollection::DeepCopy(c5, collection);
	auto etree5 = ExtendedTree::Assemble(collection);
	double area5=Area(c5) - (-Area(c1)-Area(c2)-Area(c3));
	double len5 = c1.length() + c2.length() + c3.length() + c4.length() + c5.length();
	add_check(fabs(Area(etree5) - area5)<1e-8 && fabs(etree5.length() - len5)<1e-12, "4 closed, 1 open");

	//6) + another contour
	auto c6 = Constructor::Circle(13, 4.5, Point(0, 0));
	ECollection::DeepCopy(c6, collection);
	auto etree6 = ExtendedTree::Assemble(collection);
	double area6=Area(c5) - (-Area(c1)) + Area(c6) - Area(c2) - Area(c3);
	double len6 = c6.length() + c1.length() + c2.length() + c3.length() + c4.length() + c5.length();
	add_check(fabs(Area(etree6) - area6)<1e-8 && fabs(etree6.length() - len6)<1e-12, "5 closed, 1 open");
}

void test6(){
	std::cout<<"Contour unite"<<std::endl;
	
	Container<Contour> c1 = Constructor::ContourFromPoints({0,0, 2,1, 4,3, 5,1, 6,2, 7,0});
	Contour a1 = Contour::ShallowCopy(c1, 0, 1);
	Contour a2 = Contour::ShallowCopy(c1, 2, 2);
	a1.Unite(a2);
	add_check(a1.ordered_points().size() == 4, "unite open lines");

	try{
		Contour a3 = Contour::ShallowCopy(c1, 1, 1);
		a1.Unite(a3);
		add_check(false, "try to add not connected line");
	} catch (GeomError& g){
		add_check(true, "try to add not connected line");
	}
	
	a1.add_value(Edge{c1.pdata.point(0), c1.pdata.point(3)});
	add_check(a1.is_closed(), "closing contour");

	try{
		Contour a4 = Contour::ShallowCopy(c1);
		a1.Unite(a4);
		add_check(false, "try to add to a closed contour");
	} catch (GeomError& g){
		add_check(true, "try to add to a closed contour");
	}


}

int main(){
	std::cout<<"hybmesh_contours2d testing"<<std::endl;
	if (hybmesh_contours2d_ping(1) == 2) 
		std::cout<<"Ping OK"<<std::endl;

	test1();
	test2();
	test3();
	test4();
	test5();
	test6();
	

	if (FAILED_CHECKS == 1){
		std::cout<<FAILED_CHECKS<<" test failed <<<<<<<<<<<<<<<<<<<"<<std::endl;
	} else if (FAILED_CHECKS > 1) {
		std::cout<<FAILED_CHECKS<<" tests failed <<<<<<<<<<<<<<<<<<<"<<std::endl;
	} else {
		std::cout<<"All tests passed"<<std::endl;
	}
	std::cout<<"DONE"<<std::endl;
}
