#include "hmblay.hpp"
#include <memory>
#include "bgrid.hpp"

int hmblay_ping(int a){
	return 2*a;
}

using namespace HMBlay::Impl;

// ========================= Main Algo
GridGeom HMBlay::BuildBLayerGrid(const vector<HMBlay::Input>& orig_opt){
	if (orig_opt.size() == 0) return GridGeom(0, 0, 0, 0);
	//0) check input
	for (auto& o: orig_opt){
		if (o.edges == 0 || o.edges->size() == 0){
			throw EBuildError("source edges were not set");
		}
		if (o.partition.size() < 2 || !ISZERO(o.partition[0])){
			throw EBuildError("invalid partition");
		}
	}

	//1,2) internal, non-dimensional copy of orig_options
	std::vector<Options> opt = Options::CreateFromParent(orig_opt);

	//4) Sequences of PathOptions
	//Organize options in a sequences depending on start/end points
	auto seqvec = Options::BuildSequence(opt);

	//5) Build Grid for each sequence
	ShpVector<BGrid> gg;
	for (auto& seq: seqvec){
		auto a = BGrid::MeshSequence(seq);
		gg.push_back(a);
	}

	//6) impose boundary grids
	auto impres = BGrid::ImposeBGrids(gg);

	//7) to original size
	impres->undo_scale(*opt[0].get_scaling());

	return GridGeom::sum({impres.get()});
}

