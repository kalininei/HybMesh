#include "hmgraph.hpp"
#include <algorithm>
using namespace HMMath;

namespace{

std::vector<int> subgraph(std::vector<bool>& used,
		const std::vector<std::vector<int>>& graph, int i1){
	std::vector<int> ret(1, i1);
	used[i1] = true;
	int uu = 0;
	while (uu<ret.size()){
		for (auto a: graph[ret[uu]]) if (!used[a]){
			ret.push_back(a);
			used[a] = true;
		}
		++uu;
	}
	return ret;
}

};

std::vector<std::vector<int>>
Graph::SplitGraph(const std::vector<std::vector<int>>& graph){
	std::vector<std::vector<int>> ret;
	std::vector<bool> used(graph.size(), false);
	while (1){
		int i1 = std::find(used.begin(), used.end(), false) - used.begin();
		if (i1 >= graph.size()) break;
		ret.push_back(subgraph(used, graph, i1));
	};

	return ret;
}
