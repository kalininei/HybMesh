#ifndef HMCONT2D_ALGOS_HPP
#define HMCONT2D_ALGOS_HPP
#include "contour.hpp"
#include "tree.hpp"

namespace HMCont2D{

//Offset contour
enum class OffsetTp{
	CLOSED_POLY,
	OPEN_ROUND,
	OPEN_BUTT,
};

enum class PartitionTp{
	IGNORE_ALL,   //ignore all source vertices
	KEEP_ALL,     //keep all source vertices
	KEEP_SHAPE    //keep only shape significant points
};

namespace Algos{
// ================================== Offset
//takes into account direction of source and sign of delta:
//all positives -> offsets to the left etc.
Container<ContourTree> Offset(const Contour& source, double delta, OffsetTp tp);
//forces singly connected output contour. tp = CLOSED_POLY or OPEN_ROUND
Container<Contour> Offset1(const Contour& source, double delta);

// ================================== Contour Partition
//step - step of partitioning
//pstore - point collection where to put new generated points with defined ShpGenerator
Contour Partition(double step, const Contour& contour, PCollection& pstore, PartitionTp tp);
//make a partition keeping points defined in keepit
Contour Partition(double step, const Contour& contour, PCollection& pstore,
		const std::vector<Point*>& keepit = {});
//this can be called with Container<Contour> object.
template<class TContainer,
	class = Tpp::IsBase<Contour, typename TContainer::TParent>>
static TContainer Partition(double step, const TContainer& container, PartitionTp tp){
	PCollection pc;
	auto x = Partition(step, container, pc, tp);
	TContainer ret;
	TContainer::DeepCopy(x, ret);
	return ret;
}
template<class TContainer,
	class = Tpp::IsBase<Contour, typename TContainer::TParent>>
static TContainer Partition(double step, const TContainer& container,
		const vector<Point*>& keepit = {}){
	return Partition(step, container, container.pdata, keepit);
}
//here the partition step will be calculated according to basis:
//  dictionary that maps  (contour weight in [0,1]) -> (required partition size)
Contour WeightedPartition(const std::map<double, double>& basis,
		const Contour& contour, PCollection& pstore, PartitionTp tp);
Contour WeightedPartition(const std::map<double, double>& basis,
		const Contour& contour, PCollection& pstore,
		const std::vector<Point*>& keepit = {});

// ============================ Crosses and intersections
//finds first cross (with respect to length of c1) of contours c1, c2.
//returns <0>: if cross was found
//        <1>: cross point
//        <2,3>: normalized length coordinate of intersection
std::tuple<bool, Point, double, double> 
Cross(const Contour& c1, const Contour& c2);

vector<std::tuple<bool, Point, double, double>>
CrossAll(const Contour& c1, const Contour& c2);

//returns true if c1 and c2 have common area (not a point, but may be an edge)
bool DoIntersect(const Contour& c1, const Contour& c2);
bool DoIntersect(const ContourTree& t1, const Contour& c2);

//returns true if c1 and c2 have common area (not a point, not an edge)
//is not realiable if  Area(c1) >> Area(c2) 
bool DoReallyIntersect(const Contour& c1, const Contour& c2);



// =========================== Simplifications
//remove points which lie on the same edge
//removes zero length edges
Contour Simplified(const Contour& cont);
ContourTree Simplified(const ContourTree& t1);
ExtendedTree Simplified(const ExtendedTree& t1);

// =========================== Smoothing
//calculates vector representing direction of contour at point p smoother by lengh len
//if p doesn't lie on c -> project it to c and calculate
Vect SmoothedDirection(const Contour& c, Point* p, int direction, double len);


}
};


#endif
