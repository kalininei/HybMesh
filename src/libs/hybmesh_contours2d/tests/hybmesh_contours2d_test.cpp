#include <type_traits>
#include "hmtesting.hpp"
#include "primitives2d.hpp"
#include "constructor.hpp"
#include "algos.hpp"
#include "cont_partition.hpp"
#include "contclipping.hpp"
#include "cont_assembler.hpp"
#include "vtk_export2d.hpp"
#include "cont_repart.hpp"

using HMTesting::add_check;

using namespace HM2D;

void test1(){
	std::cout<<"1. Offset closed polygon"<<std::endl;
	auto c6 = Contour::Constructor::Circle(6, 1.0, Point(0,0));
	auto oret = Contour::Algos::Offset(c6, 0.2, Contour::Algos::OffsetTp::RC_CLOSED_POLY);
	add_check(fabs(oret.area() - 3.923) < 0.01, "6-sided closed polygon offsetting");
}

void test2(){
	std::cout<<"2. Deep copy procedures"<<std::endl;
	//edges and points pool
	auto c = Contour::Constructor::Circle(3, 1, Point(5,4));
	//copy only edges
	EdgeData ne;
	DeepCopy(c, ne, 0);
	VertexData cvert = AllVertices(c);
	try{
		for (auto e=ne.begin(); e!=ne.end(); ++e){
			//new edges are not in c edges
			if (std::find(c.begin(), c.end(), *e) != c.end()) throw 1;
			//points in new edges are all in e.pdata
			if (std::find(cvert.begin(), cvert.end(), (*e)->first()) == cvert.end()) throw 1;
			if (std::find(cvert.begin(), cvert.end(), (*e)->last()) == cvert.end()) throw 1;
		}
		add_check(true, "Deep copy of edges");
	} catch (int) {add_check(false, "Deep copy of edges");}

	//copy only points
	VertexData pc;
	DeepCopy(cvert, pc);
	try{
		for (auto p=pc.begin(); p!=pc.end(); ++p){
			//new points are not in c.pdata
			if (std::find(cvert.begin(), cvert.end(), *p) != cvert.end()) throw 1;
		}
		add_check(true, "Deep copy of points");
	} catch (int) { add_check(false, "Deep copy of points"); }

	//copy whole container
	EdgeData nc;
	VertexData ncvert = AllVertices(nc);
	DeepCopy(c, nc, 1);
	try{
		for (auto e=nc.begin(); e!=nc.end(); ++e){
			//edges of nc are not in c edges
			if (std::find(c.begin(), c.end(), *e) != c.end()) throw 1;
			//points of nc are not in in e.pdata
			if (std::find(cvert.begin(), cvert.end(), (*e)->first()) != cvert.end()) throw 1;
			if (std::find(cvert.begin(), cvert.end(), (*e)->last()) != cvert.end()) throw 1;
			//points of nc are in nc
			if (std::find(ncvert.begin(), ncvert.end(), (*e)->first()) == ncvert.end()) throw 1;
			if (std::find(ncvert.begin(), ncvert.end(), (*e)->last()) == ncvert.end()) throw 1;
		}
		add_check(true, "Deep copy of a whole container");
	} catch (int) { add_check(false, "Deep copy of a whole container"); }
}

void test3(){
	std::cout<<"3. Partition algorithm"<<std::endl;
	auto c12 = Contour::Constructor::Circle(12, 10, Point(5,4));
	auto c6 = Contour::Constructor::FromPoints({0.0,0.0, 0.2,0, 0.33,0, 0.7,0.2, 0.9,0.6, 1.0,0.0 });
	double len12 = Length(c12);
	double len6 = Length(c6);
	//add_check(fabs(len12 - 62) - 0.5, "12-sided contour length");
	//1) ignore all algorithm
	auto r1 = Contour::Algos::Partition(len12/5.0, c12, Contour::Algos::PartitionTp::IGNORE_ALL);
	add_check(fabs(Length(r1)-57.5067)<1e-3, "ignore all");

	//2) keep all
	auto r2 = Contour::Algos::Partition(len12/5.0, c12, Contour::Algos::PartitionTp::KEEP_ALL);
	add_check(fabs(Length(r2)-len12)<1e-3 && r2.size() == c12.size(), "keep all for closed contour");

	auto r3 = Contour::Algos::Partition(0.11, c6, Contour::Algos::PartitionTp::KEEP_ALL);
	add_check(fabs(Length(r3)-len6)<1e-12, "keep all for a polyline");

	//3) keep shape
	auto r4 = Contour::Algos::Partition(1.0, c6, Contour::Algos::PartitionTp::KEEP_SHAPE);
	add_check(r4.size() == 4 && fabs(Length(r4)-len6)<1e-12,
				"keep shape coarse");
	auto r5 = Contour::Algos::Partition(0.1, c6, Contour::Algos::PartitionTp::KEEP_SHAPE);
	add_check(fabs(Length(r5)-len6)<1e-12, "keep shape fine");
}

void test4(){
	std::cout<<"4. Partition direction"<<std::endl;
	//ignore all
	auto c10 = Contour::Constructor::Circle(10, 1, Point(-5,4));
	auto r1 = Contour::Algos::Partition(100, c10, Contour::Algos::PartitionTp::IGNORE_ALL);
	add_check(r1.size() == 3 && Contour::Area(r1) > 0, "positive closed path, ignore all");
	Contour::Reverse(c10);
	auto r2 = Contour::Algos::Partition(0.1, c10, Contour::Algos::PartitionTp::IGNORE_ALL);
	add_check(r2.size() == std::lround(Length(c10)/0.1) && Contour::Area(r2) < 0, "negative closed path, ignore all");

	//keep shape
	auto in4 = Contour::Constructor::Circle(10, 0.2, Point(2,3));
	Contour::Reverse(in4);
	auto out4 = Contour::Algos::Partition(0.1, in4, Contour::Algos::PartitionTp::KEEP_SHAPE);
	add_check(fabs(Contour::Area(in4) - Contour::Area(out4))<1e-8, "negative closed path, keep shape");
}

void test5(){
	std::cout<<"5. Tree structure"<<std::endl;
	EdgeData collection;

	//1) closed contour
	auto c1 = Contour::Constructor::Circle(13, 5, Point(0, 0));
	Contour::Reverse(c1);
	DeepCopy(c1, collection);
	auto etree1 = Contour::Tree::Assemble(collection);
	add_check(fabs(etree1.area() + Contour::Area(c1))<1e-12, "1 closed contour");

	//2) + inside contour
	auto c2 = Contour::Constructor::Circle(20, 0.3, Point(1,1));
	DeepCopy(c2, collection);
	auto etree2 = Contour::Tree::Assemble(collection);
	add_check(fabs(etree2.area() - (-Contour::Area(c1)-Contour::Area(c2)))<1e-8, "2 closed");

	//3) + sibling contour
	auto c3 = Contour::Constructor::Circle(3, 0.2, Point(-2.3, -3));
	DeepCopy(c3, collection);
	auto etree3 = Contour::Tree::Assemble(collection);
	add_check(fabs(etree3.area() - (-Contour::Area(c1)-Contour::Area(c2)-Contour::Area(c3)))<1e-8, "3 closed");

	//4) + open path
	auto c4 = Contour::Constructor::FromPoints({2,2, 3,3, 2.5,3.5});
	DeepCopy(c4, collection);
	auto etree4 = Contour::Tree::Assemble(collection);
	double area4=(-Contour::Area(c1)-Contour::Area(c2)-Contour::Area(c3));
	double len4 = Length(c1) + Length(c2) + Length(c3) + Length(c4);
	add_check(fabs(etree4.area() - area4)<1e-12 && fabs(Length(etree4.alledges()) - len4)<1e-12, "3 closed, 1 open");

	//5) + big contour
	auto c5 = Contour::Constructor::Circle(100, 30, Point(1,1));
	DeepCopy(c5, collection);
	auto etree5 = Contour::Tree::Assemble(collection);
	double area5=Contour::Area(c5) - (-Contour::Area(c1)-Contour::Area(c2)-Contour::Area(c3));
	double len5 = Length(c1) + Length(c2) + Length(c3) + Length(c4) + Length(c5);
	add_check(fabs(etree5.area() - area5)<1e-8 && fabs(Length(etree5.alledges()) - len5)<1e-12, "4 closed, 1 open");

	//6) + another contour
	auto c6 = Contour::Constructor::Circle(13, 4.5, Point(0, 0));
	DeepCopy(c6, collection);
	auto etree6 = Contour::Tree::Assemble(collection);
	double area6=Contour::Area(c5) - (-Contour::Area(c1)) + Contour::Area(c6) - Contour::Area(c2) - Contour::Area(c3);
	double len6 = Length(c6) + Length(c1) + Length(c2) + Length(c3) + Length(c4) + Length(c5);
	add_check(fabs(etree6.area() - area6)<1e-8 && fabs(Length(etree6.alledges()) - len6)<1e-12, "5 closed, 1 open");
}

void test6(){
	std::cout<<"6. Contour unite"<<std::endl;

	EdgeData c1 = Contour::Constructor::FromPoints({0,0, 2,1, 4,3, 5,1, 6,2, 7,0});
	EdgeData a1 {c1[0], c1[1]};
	EdgeData a2 {c1[2]};
	Contour::Connect(a1, a2);
	add_check(Contour::OrderedPoints(a1).size() == 4, "unite open lines");

	try{
		EdgeData a3 {c1[1]};
		Contour::Connect(a1, a3);
		add_check(false, "try to add not connected line");
	} catch (std::runtime_error& g){
		add_check(true, "try to add not connected line");
	}

	a1.emplace_back(new Edge(c1[0]->first(), c1.back()->last()));
	add_check(Contour::IsClosed(a1), "closing contour");

	try{
		Contour::Connect(a1, c1);
		add_check(false, "try to add to a closed contour");
	} catch (std::runtime_error& g){
		add_check(true, "try to add to a closed contour");
	}


}

void test7(){
	std::cout<<"7. Contour point coordinate by weight"<<std::endl;
	auto top = Contour::Constructor::FromPoints({
		Point{0.43701993, 0.07659573},
		Point{0.44410482, 0.07840499},
		Point{0.45201250, 0.07960025},
		Point{0.46000000, 0.08000000},
		Point{1.00000000, 0.08000000}});

	Point p(0.519522452, 0.08);
	auto coord = Contour::CoordAt(top, p);
	double w = std::get<1>(coord);
	Point p1 = Contour::WeightPoint(top, w);
	Point diff = p-p1;
	add_check(ISZERO(diff.x) && ISZERO(diff.y),
			"weight by point -> point by weight");
}

void test8(){
	std::cout<<"8. Contours Cross"<<std::endl;
	auto c1 = Contour::Constructor::FromPoints(
		{0,0, 1,0, 1,1, 0,1}, true);
	auto c2 = Contour::Constructor::FromPoints(
		{0.5,0.5, 1,3, 0.5,2}, true);
	auto res1 = Contour::Clip::Intersection(c1, c2);
	add_check(ISEQ(res1.area(), 0.025), "square cross triangle");

	auto c3 = Contour::Constructor::FromPoints(
		{0.5,1, 1,3, 0.5,2}, true);
	auto res2 = Contour::Clip::Intersection(c1, c3);
	add_check(ISEQ(res2.area(), 0.0), "Only one common point");

	auto c4 = Contour::Constructor::FromPoints(
		{0.5, 0.5, 0.6,0.6, 0.5,0.6}, true);
	auto res3 = Contour::Clip::Intersection(c4, c1);
	add_check(ISEQ(res3.area(), Contour::Area(c4)), "One within the other");

	auto c5 = Contour::Constructor::FromPoints(
		{1.1, 0.5, 0.5, 1.1, 2,2, 2,-1, 0,-0.1}, true);
	auto res4 = Contour::Clip::Intersection(c5, c1);
	add_check(res4.nodes.size()==2 && fabs(res4.area() - 0.2618939)<1e-7,
			"Multiple contours in result");
}

void test9(){
	std::cout<<"9. Contours union"<<std::endl;
	auto c1 = Contour::Constructor::FromPoints(
		{0,0, 1,0, 1,1, 0,1}, true);
	auto c2 = Contour::Constructor::FromPoints(
		{0.5,0.5, 1,3, 0.5,2}, true);
	auto res1 = Contour::Clip::Union(c1, c2);
	add_check(ISEQ(res1.area(), 1.35), "square union triangle");

	auto c3 = Contour::Constructor::FromPoints(
		{1,0, 2,0, 2,1, 1,1}, true);
	auto res2 = Contour::Clip::Union({c1, c2, c3});
	add_check(fabs(res2.area()-2.35)<1e-7, "squares with common edge");

	auto c4 = Contour::Constructor::FromPoints(
		{1.5,1, 2,3, 1.5,2}, true);
	auto res3 = Contour::Clip::Union({c1, c2, c3, c4});
	add_check(fabs(res3.area() - 2.6)<1e-7, "common point");
	vector<EdgeData> sep3 = Contour::Assembler::SimpleContours(res3.alledges());
	add_check([&](){
		if (sep3.size() != 2) return false;
		if (Contour::IsOpen(sep3[0])) return false;
		if (Contour::IsOpen(sep3[1])) return false;
		double v1 = fabs(Contour::Area(sep3[0]));
		double v2 = fabs(Contour::Area(sep3[1]));
		return ISZERO(res3.area() - v1 - v2);
	}(), "common point result separate");

	auto c5 = Contour::Constructor::FromPoints(
		{0.2,0.2, 0.7,0.2, 0.7,0.7, 0.2,0.7}, true);

	Contour::Tree tree;
	tree.AddContour(c1);
	tree.AddContour(c5);

	auto res4 = Contour::Clip::Union(tree, c2);
	add_check(fabs(res4.area() - 1.104)<1e-7, "union with a tree");
}

void test10(){
	std::cout<<"10. Contours substruct"<<std::endl;
	auto c1 = Contour::Constructor::FromPoints(
		{0,0, 1,0, 1,1, 0,1}, true);
	auto c2 = Contour::Constructor::FromPoints(
		{0.5,0.5, 1,3, 0.5,2}, true);
	auto res1 = Contour::Clip::Difference(c1, c2);
	add_check(ISEQ(res1.area(), 0.975), "square minus triangle");


	auto c5 = Contour::Constructor::FromPoints(
		{0.2,0.2, 0.7,0.2, 0.7,0.7, 0.2,0.7}, true);
	Contour::Tree tree;
	tree.AddContour(c1);
	tree.AddContour(c5);

	auto res2 = Contour::Clip::Difference(tree, c2);
	add_check(fabs(res2.area() - 0.729)<1e-7, "substruction from a tree");
}

void test11(){
	std::cout<<"11. Extended tree assembling"<<std::endl;
	//common point
	VertexData pcol;
	pcol.emplace_back(new Vertex(0, 0));
	pcol.emplace_back(new Vertex(1, 0));
	pcol.emplace_back(new Vertex(2, 0));
	pcol.emplace_back(new Vertex(0, 1));
	pcol.emplace_back(new Vertex(1, 1));
	pcol.emplace_back(new Vertex(2, 1));
	pcol.emplace_back(new Vertex(0, 2));
	pcol.emplace_back(new Vertex(1, 2));
	pcol.emplace_back(new Vertex(2, 2));

	Edge ed[13];
	//basic
	ed[0] = Edge(pcol[0], pcol[1]);
	ed[1] = Edge(pcol[1], pcol[2]);
	ed[2] = Edge(pcol[3], pcol[4]);
	ed[3] = Edge(pcol[4], pcol[5]);
	ed[4] = Edge(pcol[6], pcol[7]);
	ed[5] = Edge(pcol[7], pcol[8]);
	ed[6] = Edge(pcol[0], pcol[3]);
	ed[7] = Edge(pcol[1], pcol[4]);
	ed[8] = Edge(pcol[2], pcol[5]);
	ed[9] = Edge(pcol[3], pcol[6]);
	ed[10] = Edge(pcol[4], pcol[7]);
	ed[11] = Edge(pcol[5], pcol[8]);
	//deep copies
	ed[12] = Edge(pcol[0], pcol[1]);

	//1. 1 side
	EdgeData ecol1;
	ecol1.emplace_back(new Edge(ed[0]));
	auto tree1 = Contour::Tree::Assemble(ecol1);
	add_check(tree1.nodes.size() == 1 && Length(tree1.open_contours()[0]->contour) == 1.0, "one side");

	//2. 2 divergent edges
	EdgeData ecol2;
	ecol2.emplace_back(new Edge(ed[0]));
	ecol2.emplace_back(new Edge(ed[11]));
	auto tree2 = Contour::Tree::Assemble(ecol2);
	add_check(tree2.nodes.size() == 2 && Length(tree2.open_contours()[0]->contour) == 1.0 &&
		Length(tree2.open_contours()[1]->contour) == 1.0, "two divergent sides");

	//3. 2 connected edges
	EdgeData ecol3;
	ecol3.emplace_back(new Edge(ed[0]));
	ecol3.emplace_back(new Edge(ed[1]));
	auto tree3 = Contour::Tree::Assemble(ecol3);
	add_check(tree3.nodes.size() == 1 &&
			Length(tree3.open_contours()[0]->contour) == 2.0, "two connected edges");

	//4. 2 connected and 1 disconnected
	EdgeData ecol4;
	ecol4.emplace_back(new Edge(ed[3]));
	ecol4.emplace_back(new Edge(ed[0]));
	ecol4.emplace_back(new Edge(ed[2]));
	auto tree4 = Contour::Tree::Assemble(ecol4);
	add_check(tree4.nodes.size() == 2 && [&]()->bool{
		auto oc1 = tree4.open_contours()[0]->contour;
		auto oc2 = tree4.open_contours()[1]->contour;
		if (oc1.size() == 2) std::swap(oc1, oc2);
		if (oc1.size() != 1 || Length(oc1) != 1.0) return false;
		if (oc2.size() != 2 || Length(oc2) != 2.0) return false;
		return true;
	}(), "2 connected and 1 disconnected edges");

	//5. 3 connected edges
	EdgeData ecol5;
	ecol5.emplace_back(new Edge(ed[3]));
	ecol5.emplace_back(new Edge(ed[11]));
	ecol5.emplace_back(new Edge(ed[2]));
	auto tree5 = Contour::Tree::Assemble(ecol5);
	add_check(tree5.nodes.size() == 1 && Length(tree5.open_contours()[0]->contour) == 3.0,
			"3 connected edges");

	//6. Closed contour
	EdgeData ecol6;
	ecol6.emplace_back(new Edge(ed[0]));
	ecol6.emplace_back(new Edge(ed[1]));
	ecol6.emplace_back(new Edge(ed[8]));
	ecol6.emplace_back(new Edge(ed[11]));
	ecol6.emplace_back(new Edge(ed[5]));
	ecol6.emplace_back(new Edge(ed[4]));
	ecol6.emplace_back(new Edge(ed[9]));
	ecol6.emplace_back(new Edge(ed[6]));
	auto tree6_1 = Contour::Tree::Assemble(ecol6);
	std::reverse(ecol6.begin(), ecol6.end());
	auto tree6_2 = Contour::Tree::Assemble(ecol6);
	add_check(tree6_1.nodes.size() == 1 && tree6_1.area() == 4.0 &&
		  tree6_2.nodes.size() == 1 && tree6_2.area() == 4.0,
		  "closed contour");

	//7. 2 closed contours
	EdgeData ecol7;
	ecol7.emplace_back(new Edge(ed[0]));
	ecol7.emplace_back(new Edge(ed[12]));
	ecol7.emplace_back(new Edge(ed[5]));
	ecol7.emplace_back(new Edge(ed[11]));
	ecol7.emplace_back(new Edge(ed[3]));
	ecol7.emplace_back(new Edge(ed[10]));
	auto tree7 = Contour::Tree::Assemble(ecol7);
	add_check(tree7.nodes.size() == 2 && tree7.area() == 1.0,
		  "2 closed contours");
	
	//8. 2 closed contours with common point
	EdgeData ecol8;
	ecol8.emplace_back(new Edge(ed[0]));
	ecol8.emplace_back(new Edge(ed[2]));
	ecol8.emplace_back(new Edge(ed[3]));
	ecol8.emplace_back(new Edge(ed[11]));
	ecol8.emplace_back(new Edge(ed[5]));
	ecol8.emplace_back(new Edge(ed[6]));
	ecol8.emplace_back(new Edge(ed[10]));
	ecol8.emplace_back(new Edge(ed[7]));
	auto tree8 = Contour::Tree::Assemble(ecol8);
	add_check(tree8.nodes.size() == 2 && tree8.area() == 2.0,
		  "2 closed contours with a  common point");

	//9. Closed contour with unclosed contour
	EdgeData ecol9;
	ecol9.emplace_back(new Edge(ed[4]));
	ecol9.emplace_back(new Edge(ed[10]));
	ecol9.emplace_back(new Edge(ed[2]));
	ecol9.emplace_back(new Edge(ed[9]));
	ecol9.emplace_back(new Edge(ed[5]));
	ecol9.emplace_back(new Edge(ed[11]));
	ecol9.emplace_back(new Edge(ed[8]));
	auto tree9 = Contour::Tree::Assemble(ecol9);
	add_check(tree9.nodes.size() == 2 && tree9.area() == 1.0 &&
		Length(tree9.open_contours()[0]->contour) == 3.0,
		"closed contour + open contour with common points");

	//10. Swastika sign 
	EdgeData ecol10;
	ecol10.emplace_back(new Edge(ed[6]));
	ecol10.emplace_back(new Edge(ed[2]));
	ecol10.emplace_back(new Edge(ed[3]));
	ecol10.emplace_back(new Edge(ed[11]));
	ecol10.emplace_back(new Edge(ed[1]));
	ecol10.emplace_back(new Edge(ed[7]));
	ecol10.emplace_back(new Edge(ed[4]));
	ecol10.emplace_back(new Edge(ed[10]));

	auto tree10 = Contour::Tree::Assemble(ecol10);
	add_check(tree10.nodes.size() == 2 &&
		Length(tree10.open_contours()[0]->contour) == 4.0 &&
		Length(tree10.open_contours()[1]->contour) == 4.0,
		"two crossed open contours");

	//11. All in one. Result here is ambiguous.
	//    in current version it tends to built smaller open contours.
	//    but this should not be used in further algorithms
	EdgeData ecol11;
	for (int i=0; i<12; ++i) ecol11.emplace_back(new Edge(ed[i]));
	auto tree11 = Contour::Tree::Assemble(ecol11);
	add_check(tree11.nodes.size() == 4 && tree11.area() == 2.0, "all in one");
}

void test12(){
	std::cout<<"12. Contour smoothed direction"<<std::endl;
	//auto cn1 = Contour::Constructor::Circle(16, 1, Point(0, 0));
	//auto vec1 = Contour::Algos::SmoothedDirection(cn1, cn1.point(0), 1, 0.1);
	//auto an1 = atan2(vec1.y, vec1.x)/M_PI*180;
	//std::cout<<an1<<std::endl;

}

void test13(){
	std::cout<<"13. Sorting points out"<<std::endl;
	auto cn1 = Contour::Constructor::FromPoints({0,0, 1,0, 1,1, 0,1}, true);
	vector<Point> p1 = {Point(-0.5, 0.5), Point(0, 0.5), Point(0.5, 0.5), Point(1, 0.5), Point(1, 0)};
	auto ans1 = Contour::Algos::SortOutPoints(cn1, p1);

	add_check([&]{
		if (ans1[0] != OUTSIDE) return false;
		if (ans1[1] != BOUND) return false;
		if (ans1[2] != INSIDE) return false;
		if (ans1[3] != BOUND) return false;
		if (ans1[4] != BOUND) return false;
		return true;
	}(), "square contour");

	auto cn2 = Contour::Constructor::FromPoints({0.1,0.1, 0.9,0.1, 0.9,0.9, 0.1,0.9}, true);
	Contour::Tree ctree;
	ctree.AddContour(cn1);
	ctree.AddContour(cn2);

	vector<Point> p2 = {Point(-0.5, 0.5), Point(0, 0.5), Point(0.05, 0.5), Point(0.1, 0.5), Point(0.5, 0.5)};
	auto ans2 = Contour::Algos::SortOutPoints(ctree, p2);
	add_check([&]{
		if (ans2[0] != OUTSIDE) return false;
		if (ans2[1] != BOUND) return false;
		if (ans2[2] != INSIDE) return false;
		if (ans2[3] != BOUND) return false;
		if (ans2[4] != OUTSIDE) return false;
		return true;
	}(), "square contour with a hole");

	auto rotate = [](double angle, Point& pnt){
		pnt.set(vecRotate(pnt, angle/180.0*M_PI));
	};
	for (auto& p: AllVertices(cn1)) rotate(37, *p);
	for (auto& p: AllVertices(cn2)) rotate(37, *p);
	for (auto& p: p1) rotate(37, p);
	for (auto& p: p2) rotate(37, p);

	auto ans3 = Contour::Algos::SortOutPoints(cn1, p1);
	auto ans4 = Contour::Algos::SortOutPoints(ctree, p2);
	add_check([&]{
		if (ans3[0] != OUTSIDE) return false;
		if (ans3[1] != BOUND) return false;
		if (ans3[2] != INSIDE) return false;
		if (ans3[3] != BOUND) return false;
		if (ans3[4] != BOUND) return false;
		if (ans4[0] != OUTSIDE) return false;
		if (ans4[1] != BOUND) return false;
		if (ans4[2] != INSIDE) return false;
		if (ans4[3] != BOUND) return false;
		if (ans4[4] != OUTSIDE) return false;
		return true;
	}(), "after rotation");

}

void test14(){
	std::cout<<"14. Weighted Partition"<<std::endl;
	auto line01 = Contour::Constructor::FromPoints({-0.3,0, 1.7,0});
	std::map<double, double> w1;
	w1[0] = 0.2;
	w1[0.5] = 0.6;
	w1[0.75] = 0.4;
	w1[1.0] = 0.2;
	auto ans1 = Contour::Algos::WeightedPartition(w1, line01);
	add_check(ans1.size() == 5 && fabs(Contour::OrderedPoints(ans1)[1]->x - 0.154386 * 2.0 + 0.3)<1e-4,
			"symmetrical");
	std::map<double, double> w2;
	w2[0] = 0.4; w2[1] = 0.2;
	w2[0.2] = 0.01; w2[0.25] = 0.05;
	w2[0.6] = 0.03;
	auto ans2 = Contour::Algos::WeightedPartition(w2, line01);
	add_check(ans2.size() == 34 && fabs(Contour::OrderedPoints(ans2)[6]->x - 0.224817*2.0 + 0.3)<1e-4,
			"fine grid");
	std::map<double, double> w3;
	w3[0] = 1.4; w3[1] = 0.2;
	auto ans3 = Contour::Algos::WeightedPartition(w3, line01);
	add_check(ans3.size() == 3 && fabs(Contour::OrderedPoints(ans3)[2]->x - 0.842779*2.0 + 0.3)<1e-4,
			"coarse grid");
	std::map<double, double> w4;
	w4[0] = 166; w4[1] = 182;
	auto ans4 = Contour::Algos::WeightedPartition(w4, line01);
	add_check(ans4.size() == 1, "very coarse grid");

	auto circ1 = Contour::Constructor::Circle(256, 1, Point(0, 0));
	std::map<double, double> w5;
	w5[0] = 0.4; w5[0.5] = 0.2;
	auto ans5 = Contour::Algos::WeightedPartition(w5, circ1);
	add_check(ans5.size() == 22 && fabs(Contour::OrderedPoints(ans5)[11]->x + 1.0)<1e-12,
			"closed contour");
	std::map<double, double> w6;
	w6[0.25] = 0.4; w6[0.75] = 0.2;
	auto ans6 = Contour::Algos::WeightedPartition(w6, circ1);
	add_check(ans6.size() == 22 && fabs(Contour::OrderedPoints(ans6)[5]->x + 0.180351)<1e-4,
			"closed contour without 0, 1 definition");

	std::map<double, double> w7;
	w7[0.25] = 18; w7[0.75] = 19;
	auto ans7 = Contour::Algos::WeightedPartition(w7, circ1);
	add_check(ans7.size() == 3,
			"coarse closed contour");

}

void test15(){
	std::cout<<"15. Contour repartition"<<std::endl;
	auto watch_res = [](const vector<VertexData>& inp, bool isclosed){
		VertexData allpnt;
		for (auto& v1: inp) for (auto& v2: v1) allpnt.push_back(v2);
		VertexData vertpnt;
		for (auto& v1: inp) vertpnt.push_back(v1[0]);
		auto c1 = Contour::Assembler::Contour1(allpnt, isclosed);
		auto c2 = Contour::Assembler::Contour1(vertpnt, isclosed);
		Export::ContourVTK(c1, "res1.vtk");
		VertexData pcol = AllVertices(c2);
		Export::VerticesVTK(pcol, "res2.vtk");
	};
	auto totpnts = [](const vector<VertexData>& inp){
		int ret = 0;
		for (auto& v: inp) ret+=v.size();
		return ret;
	};
	auto vpnts = [](const vector<VertexData>& inp){ return inp.size(); };
	{
		auto c1 = Contour::Constructor::FromPoints({0,0, 0.99,0, 1,0, 1,1, 0,0.973}, true);
		shared_ptr<Vertex> mp = Contour::OrderedPoints(c1)[1];
		auto res1 = Contour::Algos::Coarsening(c1, {}, 0.312678, M_PI/4, 0.4);
		add_check(totpnts(res1) == 13 && vpnts(res1) == 12, "coarse square");
		auto res2 = Contour::Algos::Coarsening(c1, {mp}, 0.312678, M_PI/4, 0.4);
		add_check(totpnts(res2) == 21 && vpnts(res2) == 21, "coarse square with a mandatory point");
	}
	{
		auto c1 = Contour::Constructor::FromPoints({0,0, 1,0, 1,1});
		std::map<double, double> m;
		m[0] = 0.2; m[0.5] = 0.1; m[0.75] = 0.01; m[1]=0.2;
		auto c2 = Contour::Algos::WeightedPartition(m, c1, {Contour::OrderedPoints(c1)[1]});
		auto res = Contour::Algos::Coarsening(c2, {}, 0.2, M_PI/4, 0.2);
		add_check(totpnts(res) == 30 && vpnts(res) == 11 , "open contour");
	}
	{
		auto c1 = Contour::Constructor::Circle(133, 1, Point(0, 0));
		for (auto p: AllVertices(c1)) p->y/=4;
		auto res = Contour::Algos::Coarsening(c1, {}, 0.1, M_PI/4, 0.2);
		add_check(vpnts(res) == 43 && totpnts(res) == 135, "ellipse");
	}
	{
		auto c1 = Contour::Constructor::Circle(112, 1, Point(0, 0));
		for (auto p: AllVertices(c1)) p->y/=5;
		auto c2 = Contour::Constructor::Circle(256, 0.07, Point(-0.95, 0));
		auto c3 = Contour::Clip::Difference(c1, c2);
		auto res = Contour::Algos::Coarsening(c3.nodes[0]->contour, {}, 0.2, M_PI/4, 0.2);
		add_check(vpnts(res) == 29 && totpnts(res) == 216, "ellipse minus circ");
	}

};

void test16(){
	std::cout<<"16. Contour partition with fixed nedges"<<std::endl;
	{
		auto inp1 = Contour::Constructor::FromPoints({0,0, 1,0});
		std::map<double, double> m1;
		m1.clear(); m1[0]=0.1; m1[1.0]=0.4;
		auto ans1 = Contour::Algos::WeightedPartition(m1, inp1, 10, {});
		m1.clear(); m1[0.1]=1; m1[0.4]=5; m1[0.8]=2;
		auto ans2 = Contour::Algos::WeightedPartition(m1, inp1, 50, {});
		m1.clear(); m1[0] = 1;
		auto ans3 = Contour::Algos::WeightedPartition(m1, inp1, 4, {});
		add_check(ans1.size() == 10 && ans2.size() == 50 && ans3.size()==4, "single segment");
	}
	{
		auto inp1 = Contour::Constructor::FromPoints({0,0, 2,0, 2,1, 0,1}, true);
		std::map<double, double> m;
		VertexData ap = Contour::OrderedPoints(inp1);
		m.clear(); m[1.0/6.0]=5; m[4.0/6.0]=1;
		auto ans1 = Contour::Algos::WeightedPartition(m, inp1, 50, ap);
		add_check(ans1.size() == 50 &&
				ISEQ(ELengths(ans1)[0], ELengths(ans1)[9]), "closed contour");
	}
	{
		auto inp1 = Contour::Constructor::FromPoints({0,0, 0.001,0, 1.99,0, 2,0, 2,1, 0,1}, true);
		std::map<double, double> m;
		VertexData ap = Contour::OrderedPoints(inp1);
		ap.erase(ap.begin()+1);
		m.clear(); m[1.0/6.0]=5; m[4.0/6.0]=1;
		auto ans1 = Contour::Algos::WeightedPartition(m, inp1, 20, ap);
		Export::ContourVTK(ans1, "ans1.vtk");
		add_check(ans1.size() == 20 &&
				ISEQ(ELengths(ans1)[5], ELengths(ans1)[19]) &&
				ISEQ(ELengths(ans1)[4], 0.01), "closed contour with restriction");
	}
}


int main(){
	std::cout<<"hybmesh_contours2d testing"<<std::endl;
	test1();
	test2();
	test3();
	test4();
	test5();
	test6();
	test7();
	test8();
	test9();
	test10();
	test11();
	test12();
	test13();
	test14();
	test15();
	test16();


	HMTesting::check_final_report();
	std::cout<<"DONE"<<std::endl;
}
