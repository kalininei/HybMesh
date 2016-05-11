#include "c_hybmesh_contours2d.h"
#include "hybmesh_contours2d.hpp"
#include <type_traits>
#include "hmtesting.hpp"

using HMTesting::add_check;

using namespace HMCont2D;

void test1(){
	std::cout<<"Offset closed polygon"<<std::endl;
	auto c6 = Constructor::Circle(6, 1.0, Point(0,0));
	auto oret = Algos::Offset(c6, 0.2, OffsetTp::RC_CLOSED_POLY);
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
	auto r1 = Algos::Partition(len12/5.0, c12, PartitionTp::IGNORE_ALL);
	add_check(fabs(r1.length()-57.5067)<1e-3, "ignore all");

	//2) keep all
	auto r2 = Algos::Partition(len12/5.0, c12, PartitionTp::KEEP_ALL);
	add_check(fabs(r2.length()-len12)<1e-3 && r2.size() == c12.size(), "keep all for closed contour");

	PCollection pc3;
	auto r3 = Algos::Partition(0.11, c6, pc3, PartitionTp::KEEP_ALL);
	add_check(pc3.size() == 12 && fabs(r3.length() - len6)<1e-12, "keep all for a polyline");

	//3) keep shape
	PCollection pc4;
	auto r4 = Algos::Partition(1.0, c6, pc4, PartitionTp::KEEP_SHAPE);
	add_check(r4.size() == 4 && fabs(r4.length() -len6)<1e-12 && pc4.size() == 0,
				"keep shape coarse");
	auto r5 = Algos::Partition(0.1, c6, pc4, PartitionTp::KEEP_SHAPE);
	add_check(fabs(r5.length() - len6)<1e-12 && pc4.size() == r5.size() + 1 - 5, "keep shape fine");

}

void test4(){
	std::cout<<"Partition direction"<<std::endl;
	//ignore all
	auto c10 = Constructor::Circle(10, 1, Point(-5,4));
	auto r1 = Algos::Partition(100, c10, PartitionTp::IGNORE_ALL);
	add_check(r1.size() == 3 && Area(r1) > 0, "positive closed path, ignore all");
	c10.Reverse();
	auto r2 = Algos::Partition(0.1, c10, PartitionTp::IGNORE_ALL);
	add_check(r2.size() == std::lround(c10.length()/0.1) && Area(r2) < 0, "negative closed path, ignore all");

	//keep shape
	auto in4 = Constructor::Circle(10, 0.2, Point(2,3));
	in4.Reverse();
	auto out4 = Algos::Partition(0.1, in4, PartitionTp::KEEP_SHAPE);
	add_check(fabs(Area(in4) - Area(out4))<1e-8, "negative closed path, keep shape");
}

void test5(){
	std::cout<<"Tree structure"<<std::endl;
	ECollection collection;

	//1) closed contour
	auto c1 = Constructor::Circle(13, 5, Point(0, 0));
	c1.Reverse();
	ECollection::DeepCopy(c1, collection);
	auto etree1 = Assembler::ETree(collection);
	add_check(fabs(Area(etree1) + Area(c1))<1e-12, "1 closed contour");

	//2) + inside contour
	auto c2 = Constructor::Circle(20, 0.3, Point(1,1));
	ECollection::DeepCopy(c2, collection);
	auto etree2 = Assembler::ETree(collection);
	add_check(fabs(Area(etree2) - (-Area(c1)-Area(c2)))<1e-8, "2 closed");


	//3) + sibling contour
	auto c3 = Constructor::Circle(3, 0.2, Point(-2.3, -3));
	ECollection::DeepCopy(c3, collection);
	auto etree3 = Assembler::ETree(collection);
	add_check(fabs(Area(etree3) - (-Area(c1)-Area(c2)-Area(c3)))<1e-8, "3 closed");

	//4) + open path
	auto c4 = Constructor::ContourFromPoints({2,2, 3,3, 2.5,3.5});
	ECollection::DeepCopy(c4, collection);
	auto etree4 = Assembler::ETree(collection);
	double area4=(-Area(c1)-Area(c2)-Area(c3));
	double len4 = c1.length() + c2.length() + c3.length() + c4.length();
	add_check(fabs(Area(etree4) - area4)<1e-12 && fabs(etree4.length() - len4)<1e-12, "3 closed, 1 open");

	//5) + big contour
	auto c5 = Constructor::Circle(100, 30, Point(1,1));
	ECollection::DeepCopy(c5, collection);
	auto etree5 = Assembler::ETree(collection);
	double area5=Area(c5) - (-Area(c1)-Area(c2)-Area(c3));
	double len5 = c1.length() + c2.length() + c3.length() + c4.length() + c5.length();
	add_check(fabs(Area(etree5) - area5)<1e-8 && fabs(etree5.length() - len5)<1e-12, "4 closed, 1 open");

	//6) + another contour
	auto c6 = Constructor::Circle(13, 4.5, Point(0, 0));
	ECollection::DeepCopy(c6, collection);
	auto etree6 = Assembler::ETree(collection);
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

void test7(){
	std::cout<<"Contour point coordinate by weight"<<std::endl;
	auto top = HMCont2D::Constructor::ContourFromPoints({
		Point{0.43701993, 0.07659573},
		Point{0.44410482, 0.07840499},
		Point{0.45201250, 0.07960025},
		Point{0.46000000, 0.08000000},
		Point{1.00000000, 0.08000000}});

	Point p(0.519522452, 0.08);
	auto coord = top.coord_at(p);
	double w = std::get<1>(coord);
	Point p1 = Contour::WeightPoint(top, w);
	Point diff = p-p1;
	add_check(ISZERO(diff.x) && ISZERO(diff.y),
			"weight by point -> point by weight");
}

void test8(){
	std::cout<<"Contours Cross"<<std::endl;
	auto c1 = HMCont2D::Constructor::ContourFromPoints(
		{0,0, 1,0, 1,1, 0,1}, true);
	auto c2 = HMCont2D::Constructor::ContourFromPoints(
		{0.5,0.5, 1,3, 0.5,2}, true);
	auto res1 = HMCont2D::Clip::Intersection(c1, c2);
	add_check(ISEQ(HMCont2D::Area(res1), 0.025), "square cross triangle");

	auto c3 = HMCont2D::Constructor::ContourFromPoints(
		{0.5,1, 1,3, 0.5,2}, true);
	auto res2 = HMCont2D::Clip::Intersection(c1, c3);
	add_check(ISEQ(HMCont2D::Area(res2), 0.0), "Only one common point");

	auto c4 = HMCont2D::Constructor::ContourFromPoints(
		{0.5, 0.5, 0.6,0.6, 0.5,0.6}, true);
	auto res3 = HMCont2D::Clip::Intersection(c4, c1);
	add_check(ISEQ(HMCont2D::Area(res3), HMCont2D::Area(c4)), "One within the other");

	auto c5 = HMCont2D::Constructor::ContourFromPoints(
		{1.1, 0.5, 0.5, 1.1, 2,2, 2,-1, 0,-0.1}, true);
	auto res4 = HMCont2D::Clip::Intersection(c5, c1);
	add_check(res4.cont_count()==2 && fabs(HMCont2D::Area(res4) - 0.2618939)<1e-7,
			"Multiple contours in result");
}

void test9(){
	std::cout<<"Contours union"<<std::endl;
	auto c1 = HMCont2D::Constructor::ContourFromPoints(
		{0,0, 1,0, 1,1, 0,1}, true);
	auto c2 = HMCont2D::Constructor::ContourFromPoints(
		{0.5,0.5, 1,3, 0.5,2}, true);
	auto res1 = HMCont2D::Clip::Union(c1, c2);
	add_check(ISEQ(HMCont2D::Area(res1), 1.35), "square union triangle");

	auto c3 = HMCont2D::Constructor::ContourFromPoints(
		{1,0, 2,0, 2,1, 1,1}, true);
	auto res2 = HMCont2D::Clip::Union({c1, c2, c3});
	add_check(fabs(HMCont2D::Area(res2)-2.35)<1e-7, "squares with common edge");

	auto c4 = HMCont2D::Constructor::ContourFromPoints(
		{1.5,1, 2,3, 1.5,2}, true);
	auto res3 = HMCont2D::Clip::Union({c1, c2, c3, c4});
	add_check(fabs(HMCont2D::Area(res3) - 2.6)<1e-7, "common point");

	auto c5 = HMCont2D::Constructor::ContourFromPoints(
		{0.2,0.2, 0.7,0.2, 0.7,0.7, 0.2,0.7}, true);

	HMCont2D::ContourTree tree;
	auto p1 = shared_ptr<HMCont2D::Contour>(new HMCont2D::Contour(c1));
	auto p2 = shared_ptr<HMCont2D::Contour>(new HMCont2D::Contour(c5));
	tree.AddContour(p1);
	tree.AddContour(p2);

	auto res4 = HMCont2D::Clip::Union(tree, c2);
	add_check(fabs(HMCont2D::Area(res4) - 1.104)<1e-7, "union with a tree");
}

void test10(){
	std::cout<<"Contours substruct"<<std::endl;
	auto c1 = HMCont2D::Constructor::ContourFromPoints(
		{0,0, 1,0, 1,1, 0,1}, true);
	auto c2 = HMCont2D::Constructor::ContourFromPoints(
		{0.5,0.5, 1,3, 0.5,2}, true);
	auto res1 = HMCont2D::Clip::Difference(c1, c2);
	add_check(ISEQ(HMCont2D::Area(res1), 0.975), "square minus triangle");


	auto c5 = HMCont2D::Constructor::ContourFromPoints(
		{0.2,0.2, 0.7,0.2, 0.7,0.7, 0.2,0.7}, true);
	HMCont2D::ContourTree tree;
	auto p1 = shared_ptr<HMCont2D::Contour>(new HMCont2D::Contour(c1));
	auto p2 = shared_ptr<HMCont2D::Contour>(new HMCont2D::Contour(c5));
	tree.AddContour(p1);
	tree.AddContour(p2);

	auto res2 = HMCont2D::Clip::Difference(tree, c2);
	add_check(fabs(HMCont2D::Area(res2) - 0.729)<1e-7, "substruction from a tree");
}

void test11(){
	std::cout<<"Extended tree assembling"<<std::endl;
	//common point
	HMCont2D::PCollection pcol;
	pcol.add_value(Point(0, 0));
	pcol.add_value(Point(1, 0));
	pcol.add_value(Point(2, 0));
	pcol.add_value(Point(0, 1));
	pcol.add_value(Point(1, 1));
	pcol.add_value(Point(2, 1));
	pcol.add_value(Point(0, 2));
	pcol.add_value(Point(1, 2));
	pcol.add_value(Point(2, 2));

	HMCont2D::Edge ed[13];
	//basic
	ed[0] = HMCont2D::Edge(pcol.pvalue(0), pcol.pvalue(1));
	ed[1] = HMCont2D::Edge(pcol.pvalue(1), pcol.pvalue(2));
	ed[2] = HMCont2D::Edge(pcol.pvalue(3), pcol.pvalue(4));
	ed[3] = HMCont2D::Edge(pcol.pvalue(4), pcol.pvalue(5));
	ed[4] = HMCont2D::Edge(pcol.pvalue(6), pcol.pvalue(7));
	ed[5] = HMCont2D::Edge(pcol.pvalue(7), pcol.pvalue(8));
	ed[6] = HMCont2D::Edge(pcol.pvalue(0), pcol.pvalue(3));
	ed[7] = HMCont2D::Edge(pcol.pvalue(1), pcol.pvalue(4));
	ed[8] = HMCont2D::Edge(pcol.pvalue(2), pcol.pvalue(5));
	ed[9] = HMCont2D::Edge(pcol.pvalue(3), pcol.pvalue(6));
	ed[10] = HMCont2D::Edge(pcol.pvalue(4), pcol.pvalue(7));
	ed[11] = HMCont2D::Edge(pcol.pvalue(5), pcol.pvalue(8));
	//deep copies
	ed[12] = HMCont2D::Edge(pcol.pvalue(0), pcol.pvalue(1));

	//1. 1 side
	HMCont2D::ECollection ecol1;
	ecol1.add_value(ed[0]);
	auto tree1 = HMCont2D::Assembler::ETree(ecol1);
	add_check(tree1.cont_count() == 1 && tree1.open_contours[0]->length() == 1.0, "one side");

	//2. 2 divergent edges
	HMCont2D::ECollection ecol2;
	ecol2.add_value(ed[0]);
	ecol2.add_value(ed[11]);
	auto tree2 = HMCont2D::Assembler::ETree(ecol2);
	add_check(tree2.cont_count() == 2 && tree2.open_contours[0]->length() == 1.0 &&
		tree2.open_contours[1]->length() == 1.0, "two divergent sides");

	//3. 2 connected edges
	HMCont2D::ECollection ecol3;
	ecol3.add_value(ed[0]);
	ecol3.add_value(ed[1]);
	auto tree3 = HMCont2D::Assembler::ETree(ecol3);
	add_check(tree3.cont_count() == 1 && tree3.open_contours[0]->length() == 2.0, "two connected edges");

	//4. 2 connected and 1 disconnected
	HMCont2D::ECollection ecol4;
	ecol4.add_value(ed[3]);
	ecol4.add_value(ed[0]);
	ecol4.add_value(ed[2]);
	auto tree4 = HMCont2D::Assembler::ETree(ecol4);
	add_check(tree4.cont_count() == 2 && [&]()->bool{
		auto oc1 = tree4.open_contours[0].get();
		auto oc2 = tree4.open_contours[1].get();
		if (oc1->size() == 2) std::swap(oc1, oc2);
		if (oc1->size() != 1 || oc1->length() != 1.0) return false;
		if (oc2->size() != 2 || oc2->length() != 2.0) return false;
		return true;
	}(), "2 connected and 1 disconnected edges");

	//5. 3 connected edges
	HMCont2D::ECollection ecol5;
	ecol5.add_value(ed[3]);
	ecol5.add_value(ed[11]);
	ecol5.add_value(ed[2]);
	auto tree5 = HMCont2D::Assembler::ETree(ecol5);
	add_check(tree5.cont_count() == 1 && tree5.open_contours[0]->length() == 3.0, "3 connected edges");

	//6. Closed contour
	HMCont2D::ECollection ecol6;
	ecol6.add_value(ed[0]);
	ecol6.add_value(ed[1]);
	ecol6.add_value(ed[8]);
	ecol6.add_value(ed[11]);
	ecol6.add_value(ed[5]);
	ecol6.add_value(ed[4]);
	ecol6.add_value(ed[9]);
	ecol6.add_value(ed[6]);
	auto tree6_1 = HMCont2D::Assembler::ETree(ecol6);
	std::reverse(ecol6.data.begin(), ecol6.data.end());
	auto tree6_2 = HMCont2D::Assembler::ETree(ecol6);
	add_check(tree6_1.cont_count() == 1 && HMCont2D::Area(tree6_1) == 4.0 &&
		  tree6_2.cont_count() == 1 && HMCont2D::Area(tree6_2) == 4.0,
		  "closed contour");

	//7. 2 closed contours
	HMCont2D::ECollection ecol7;
	ecol7.add_value(ed[0]);
	ecol7.add_value(ed[12]);
	ecol7.add_value(ed[5]);
	ecol7.add_value(ed[11]);
	ecol7.add_value(ed[3]);
	ecol7.add_value(ed[10]);
	auto tree7 = HMCont2D::Assembler::ETree(ecol7);
	add_check(tree7.cont_count() == 2 && HMCont2D::Area(tree7) == 1.0,
		  "2 closed contours");
	
	//8. 2 closed contours with common point
	HMCont2D::ECollection ecol8;
	ecol8.add_value(ed[0]);
	ecol8.add_value(ed[2]);
	ecol8.add_value(ed[3]);
	ecol8.add_value(ed[11]);
	ecol8.add_value(ed[5]);
	ecol8.add_value(ed[6]);
	ecol8.add_value(ed[10]);
	ecol8.add_value(ed[7]);
	auto tree8 = HMCont2D::Assembler::ETree(ecol8);
	add_check(tree8.cont_count() == 2 && HMCont2D::Area(tree8) == 2.0,
		  "2 closed contours with a  common point");

	//9. Closed contour with unclosed contour
	HMCont2D::ECollection ecol9;
	ecol9.add_value(ed[4]);
	ecol9.add_value(ed[10]);
	ecol9.add_value(ed[2]);
	ecol9.add_value(ed[9]);
	ecol9.add_value(ed[5]);
	ecol9.add_value(ed[11]);
	ecol9.add_value(ed[8]);
	auto tree9 = HMCont2D::Assembler::ETree(ecol9);
	add_check(tree9.cont_count() == 2 && HMCont2D::Area(tree9) == 1.0 &&
		tree9.open_contours[0]->length() == 3.0,
		"closed contour + open contour with common points");

	//10. Swastika sign 
	HMCont2D::ECollection ecol10;
	ecol10.add_value(ed[6]);
	ecol10.add_value(ed[2]);
	ecol10.add_value(ed[3]);
	ecol10.add_value(ed[11]);
	ecol10.add_value(ed[1]);
	ecol10.add_value(ed[7]);
	ecol10.add_value(ed[4]);
	ecol10.add_value(ed[10]);

	auto tree10 = HMCont2D::Assembler::ETree(ecol10);
	add_check(tree10.cont_count() == 2 &&
		tree10.open_contours[0]->length() == 4.0 &&
		tree10.open_contours[1]->length() == 4.0,
		"two crossed open contours");

	//11. All in one. Result here is ambiguous.
	//    in current version it tends to built smaller open contours.
	//    but this should not be used in further algorithms
	HMCont2D::ECollection ecol11;
	for (int i=0; i<12; ++i) ecol11.add_value(ed[i]);
	auto tree11 = HMCont2D::Assembler::ETree(ecol11);
	add_check(tree11.cont_count() == 4 && HMCont2D::Area(tree11) == 2.0, "all in one");
	HMCont2D::SaveVtk(tree11, "dbgout.vtk");
	
}

void test12(){
	std::cout<<"Contour smoothed direction"<<std::endl;
	//auto cn1 = HMCont2D::Constructor::Circle(16, 1, Point(0, 0));
	//auto vec1 = HMCont2D::Algos::SmoothedDirection(cn1, cn1.point(0), 1, 0.1);
	//auto an1 = atan2(vec1.y, vec1.x)/M_PI*180;
	//std::cout<<an1<<std::endl;

}

void test13(){
	std::cout<<"Sorting points out"<<std::endl;
	auto cn1 = HMCont2D::Constructor::ContourFromPoints({0,0, 1,0, 1,1, 0,1}, true);
	vector<Point> p1 = {Point(-0.5, 0.5), Point(0, 0.5), Point(0.5, 0.5), Point(1, 0.5), Point(1, 0)};
	auto ans1 = HMCont2D::Algos::SortOutPoints(cn1, p1);

	add_check([&]{
		if (ans1[0] != OUTSIDE) return false;
		if (ans1[1] != BOUND) return false;
		if (ans1[2] != INSIDE) return false;
		if (ans1[3] != BOUND) return false;
		if (ans1[4] != BOUND) return false;
		return true;
	}(), "square contour");

	auto cn2 = HMCont2D::Constructor::ContourFromPoints({0.1,0.1, 0.9,0.1, 0.9,0.9, 0.1,0.9}, true);
	HMCont2D::ContourTree ctree;
	ctree.AddContour(cn1);
	ctree.AddContour(cn2);

	vector<Point> p2 = {Point(-0.5, 0.5), Point(0, 0.5), Point(0.05, 0.5), Point(0.1, 0.5), Point(0.5, 0.5)};
	auto ans2 = HMCont2D::Algos::SortOutPoints(ctree, p2);
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
	for (auto& p: cn1.all_points()) rotate(37, *p);
	for (auto& p: cn2.all_points()) rotate(37, *p);
	for (auto& p: p1) rotate(37, p);
	for (auto& p: p2) rotate(37, p);

	auto ans3 = HMCont2D::Algos::SortOutPoints(cn1, p1);
	auto ans4 = HMCont2D::Algos::SortOutPoints(ctree, p2);
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
	std::cout<<"Weighted Partition"<<std::endl;
	HMCont2D::PCollection pstore;
	auto line01 = HMCont2D::Constructor::ContourFromPoints({-0.3,0, 1.7,0});
	std::map<double, double> w1;
	w1[0] = 0.2;
	w1[0.5] = 0.6;
	w1[0.75] = 0.4;
	w1[1.0] = 0.2;
	auto ans1 = HMCont2D::Algos::WeightedPartition(w1, line01, pstore);
	add_check(ans1.size() == 5 && fabs(ans1.ordered_points()[1]->x - 0.154386 * 2.0 + 0.3)<1e-4,
			"symmetrical");
	std::map<double, double> w2;
	w2[0] = 0.4; w2[1] = 0.2;
	w2[0.2] = 0.01; w2[0.25] = 0.05;
	w2[0.6] = 0.03;
	auto ans2 = HMCont2D::Algos::WeightedPartition(w2, line01, pstore);
	add_check(ans2.size() == 34 && fabs(ans2.ordered_points()[6]->x - 0.224817*2.0 + 0.3)<1e-4,
			"fine grid");
	std::map<double, double> w3;
	w3[0] = 1.4; w3[1] = 0.2;
	auto ans3 = HMCont2D::Algos::WeightedPartition(w3, line01, pstore);
	add_check(ans3.size() == 3 && fabs(ans3.ordered_points()[2]->x - 0.842779*2.0 + 0.3)<1e-4,
			"coarse grid");
	std::map<double, double> w4;
	w4[0] = 166; w4[1] = 182;
	auto ans4 = HMCont2D::Algos::WeightedPartition(w4, line01, pstore);
	add_check(ans4.size() == 1,
			"very coarse grid");

	auto circ1 = HMCont2D::Constructor::Circle(256, 1, Point(0, 0));
	std::map<double, double> w5;
	w5[0] = 0.4; w5[0.5] = 0.2;
	auto ans5 = HMCont2D::Algos::WeightedPartition(w5, circ1, pstore);
	add_check(ans5.size() == 22 && fabs(ans5.ordered_points()[11]->x + 1.0)<1e-12,
			"closed contour");
	std::map<double, double> w6;
	w6[0.25] = 0.4; w6[0.75] = 0.2;
	auto ans6 = HMCont2D::Algos::WeightedPartition(w6, circ1, pstore);
	add_check(ans6.size() == 22 && fabs(ans6.ordered_points()[5]->x + 0.180351)<1e-4,
			"closed contour without 0, 1 definition");

	std::map<double, double> w7;
	w7[0.25] = 18; w7[0.75] = 19;
	auto ans7 = HMCont2D::Algos::WeightedPartition(w7, circ1, pstore);
	add_check(ans7.size() == 3,
			"coarse closed contour");

}

void test15(){
	std::cout<<"15. Contour repartition"<<std::endl;
	auto watch_res = [](const vector<vector<Point*>>& inp, bool isclosed){
		vector<Point*> allpnt;
		for (auto& v1: inp) for (auto& v2: v1) allpnt.push_back(v2);
		vector<Point*> vertpnt;
		for (auto& v1: inp) vertpnt.push_back(v1[0]);
		auto c1 = HMCont2D::Constructor::ContourFromPoints(allpnt, isclosed);
		auto c2 = HMCont2D::Constructor::ContourFromPoints(vertpnt, isclosed);
		HMCont2D::SaveVtk(c1, "res1.vtk");
		PCollection pcol;
		for (auto p: c2.all_points()) pcol.add_value(*p);
		HMCont2D::SaveVtk(pcol, "res2.vtk");
	};
	auto totpnts = [](const vector<vector<Point*>>& inp){
		int ret = 0;
		for (auto& v: inp) ret+=v.size();
		return ret;
	};
	auto vpnts = [](const vector<vector<Point*>>& inp){ return inp.size(); };
	{
		PCollection pcol;
		auto c1 = HMCont2D::Constructor::ContourFromPoints({0,0, 0.99,0, 1,0, 1,1, 0,0.973}, true);
		Point* mp = c1.ordered_points()[1];
		auto res1 = HMCont2D::Algos::Coarsening(c1, {}, pcol, 0.312678, M_PI/4, 0.4);
		add_check(totpnts(res1) == 13 && vpnts(res1) == 12 && pcol.size() == 8, "coarse square");
		pcol.clear();
		auto res2 = HMCont2D::Algos::Coarsening(c1, {mp}, pcol, 0.312678, M_PI/4, 0.4);
		add_check(totpnts(res2) == 21 && vpnts(res2) == 21 && pcol.size() == 16, "coarse square with a mandatory point");
	}
	{
		PCollection pcol;
		auto c1 = HMCont2D::Constructor::ContourFromPoints({0,0, 1,0, 1,1});
		std::map<double, double> m;
		m[0] = 0.2; m[0.5] = 0.1; m[0.75] = 0.01; m[1]=0.2;
		auto c2 = HMCont2D::Algos::WeightedPartition(m, c1, c1.pdata, {c1.ordered_points()[1]});
		auto res = HMCont2D::Algos::Coarsening(c2, {}, pcol, 0.2, M_PI/4, 0.2);
		add_check(pcol.size()==1 && totpnts(res) == 30 && vpnts(res) == 11 , "open contour");
	}
	{
		PCollection pcol;
		auto c1 = HMCont2D::Constructor::Circle(133, 1, Point(0, 0));
		for (auto p: c1.pdata) p->y/=4;
		auto res = HMCont2D::Algos::Coarsening(c1, {}, pcol, 0.1, M_PI/4, 0.2);
		add_check(pcol.size()==2 && vpnts(res) == 43 && totpnts(res) == 135, "ellipse");
	}
	{
		PCollection pcol;
		auto c1 = HMCont2D::Constructor::Circle(112, 1, Point(0, 0));
		for (auto p: c1.pdata) p->y/=5;
		auto c2 = HMCont2D::Constructor::Circle(256, 0.07, Point(-0.95, 0));
		auto c3 = HMCont2D::Clip::Difference(c1, c2);
		auto res = HMCont2D::Algos::Coarsening(*c3.nodes[0], {}, pcol, 0.2, M_PI/4, 0.2);
		add_check(pcol.size()==0 && vpnts(res) == 29 && totpnts(res) == 216, "ellipse minus circ");
	}

};


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
	test7();
	test8();
	test9();
	test10();
	test11();
	test12();
	test13();
	test14();
	test15();


	HMTesting::check_final_report();
	std::cout<<"DONE"<<std::endl;
}
