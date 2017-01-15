#ifndef HMCONT2D_CONT_PARTITION_HPP
#define HMCONT2D_CONT_PARTITION_HPP

#include "contour.hpp"

namespace HM2D{ namespace Contour{ namespace Algos{

enum class PartitionTp{
	IGNORE_ALL,   //ignore all source vertices
	KEEP_ALL,     //keep all source vertices
	KEEP_SHAPE    //keep only shape significant points
};

//step - step of partitioning
//builds new edges, but all points which should be kept will present in resulting edges
EdgeData Partition(double step, const EdgeData& contour, PartitionTp tp);
EdgeData Partition(double step, const EdgeData& contour, const VertexData& keepit = {});

//here the partition step will be calculated according to basis:
//  dictionary that maps  (contour weight in [0,1]) -> (required partition size)
EdgeData WeightedPartition(const std::map<double, double>& basis,
		const EdgeData& contour,
		PartitionTp tp);
EdgeData WeightedPartition(const std::map<double, double>& basis,
		const EdgeData& contour,
		const VertexData& keepit = {});
EdgeData WeightedPartition(std::map<double, double> basis,
		const EdgeData& contour, int nedges,
		const VertexData& keepit = {});
//partition with respect to other contour partitions
EdgeData ConditionalPartition(const EdgeData& input, double step, double influence,
		const vector<EdgeData>& condconts,
		const vector<std::pair<Point, double>>& condpoints,
		double pw,
		const VertexData& keepit = {});

}}}

#endif
