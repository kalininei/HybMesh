#ifndef HMCONT2D_CONT_PARTITION_HPP
#define HMCONT2D_CONT_PARTITION_HPP

#include "contour.hpp"
#include "tree.hpp"

namespace HMCont2D{

enum class PartitionTp{
	IGNORE_ALL,   //ignore all source vertices
	KEEP_ALL,     //keep all source vertices
	KEEP_SHAPE    //keep only shape significant points
};


// ================================== Contour Partition
namespace Algos{

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
Contour WeightedPartition(std::map<double, double> basis,
		const Contour& contour, PCollection& pstore,
		int nedges, const std::vector<Point*>& keepit = {});

//rounds vector keeping constant vector sum
vector<int> RoundVector(const vector<double>& vect, const vector<int>& minsize);

}}
#endif
