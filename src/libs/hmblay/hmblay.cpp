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
		for (int i=1; i<(int)o.partition.size(); ++i)
			if (o.partition[i]<=o.partition[i-1]) throw EBuildError("invalid partition");

		if (o.acute_angle>180) throw EBuildError("invalid acute angle");
		if (o.right_angle>180) throw EBuildError("invalid right angle");
		if (o.reentrant_angle<185) throw EBuildError("invalid reentrant angle");
		if (o.acute_angle>o.right_angle ||
				o.right_angle>o.straight_angle ||
				o.straight_angle>o.reentrant_angle)
			throw EBuildError("invalid angles set");
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

