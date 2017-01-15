#ifndef HYBMESH_GRAPH_ALGOS_HPP
#define HYBMESH_GRAPH_ALGOS_HPP
#include <vector>

namespace HMMath{ namespace Graph{ 

//split bidirectional graph into clusters.
std::vector<std::vector<int>> SplitGraph(const std::vector<std::vector<int>>& graph);


}}

#endif
