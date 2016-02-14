#include "bgrid.hpp"
#include "canonic_bgrid.hpp"
#include "trigrid.h"
#include "bgrid_impose.hpp"
#include "connectors.hpp"

using namespace HMBlay::Impl;

shared_ptr<BGrid> BGrid::MeshFullPath(const ExtPath& epath){
	//1. divide by angles
	vector<ExtPath> pths = ExtPath::DivideByAngle(epath,
			{CornerTp::RIGHT,
			CornerTp::REENTRANT, 
			CornerTp::ACUTE,
			CornerTp::ROUND});
	//if closed path with one special corner divide it into two
	if (epath.is_closed() && pths.size() == 1 &&
			(pths[0].ext_data.back().tp == CornerTp::RIGHT ||
			 pths[0].ext_data.back().tp == CornerTp::REENTRANT ||
			 pths[0].ext_data.back().tp == CornerTp::ROUND ||
			 pths[0].ext_data.back().tp == CornerTp::ACUTE)){
		pths = ExtPath::DivideByHalf(pths[0]);
	}

	//2. build conform mapping for each subpath
	ShpVector<MappedRect> mps;
	bool use_rect_approx = !epath.ext_data[0].opt->force_conformal;
	for (auto& p: pths){
		double h = p.largest_depth();
		mps.push_back(MappedRect::Factory(p.leftbc, p.rightbc, p, h,
					use_rect_approx));
	}

	//3. build rectangular meshers
	ShpVector<MappedMesher> mesher4;
	for (auto& m: mps){
		mesher4.push_back(shared_ptr<MappedMesher>());
		mesher4.back().reset(new MappedMesher(m.get(), 0, 1));
	}

	//4. Connection rules
	ShpVector<MConnector> connectors;
	for (int i=0; i<mps.size()-1; ++i){
		CornerTp t = pths[i].ext_data.back().tp;
		connectors.push_back(MConnector::Build(t, mesher4[i].get(), mesher4[i+1].get()));
	}
	//add first->last connection
	if (epath.is_closed()){
		CornerTp t = pths[0].ext_data[0].tp;
		if (mps.size() > 1){
			connectors.push_back(MConnector::Build(t, mesher4.back().get(), mesher4[0].get()));
		}
	}


	//5. build rectangular meshes
	for (int i=0; i<mesher4.size(); ++i){
		ExtPath& ipth = pths[i];
		double ilen = ipth.length();
		double depth = ipth.largest_depth();
		auto bpart = [&](double w1, double w2)->vector<double>{
			vector<double> reallen = ipth.PathPartition(w1*ilen, w2*ilen);
			for (auto& x: reallen) x/=ilen;
			return reallen;
		};
		auto wpart = [&](double w1)->vector<double>{
			vector<double> realdep = ipth.VerticalPartition(w1*ilen);
			for (auto& x: realdep) x/=depth;
			return realdep;
		};
		mesher4[i]->Fill(bpart, wpart, 2);
	}

	//6. Connection procedures
	for (auto& c: connectors) c->Apply();
	
	//7. Gather all resulting meshes
	shared_ptr<BGrid> g(new BGrid());
	if (connectors.size() > 0){
		for (auto& c: connectors){
			//add_previous with grid is empty
			bool with_prev = (g->n_cells() == 0);
			//add next always except if this
			//is an end connector for closed path
			bool with_next = (!epath.is_closed() || &c != &connectors.back());
			c->Add(g, with_prev, with_next);
		}
		g->merge_congruent_points();
	} else {
		assert(mesher4.size() == 1);
		g.reset(new BGrid(mesher4[0]->result));
	}

	return g;
}

shared_ptr<BGrid> BGrid::MeshSequence(vector<Options*>& data){ 
	auto ret = shared_ptr<BGrid>();

	//1) assemble extended path = path + boundaries + angles
	ExtPath fullpath = ExtPath::Assemble(data);

	//2) if corner angle section is very short then
	//   make it sharp or regular depending on adjacent corner types
	ExtPath::ReinterpretCornerTp(fullpath);

	//3) build grid for a path
	ret = BGrid::MeshFullPath(fullpath);

	//4) guarantee no self-intersections.
	//   Includes acute angle postprocessing.
	ret = NoSelfIntersections(ret, fullpath);

	return ret;
}

shared_ptr<BGrid> BGrid::NoSelfIntersections(shared_ptr<BGrid> g, const HMCont2D::Contour& source){
	//does grid have intersections
	if (!GGeom::Repair::HasSelfIntersections(*g)) return g;
	//else do imposition on the basis of cells weights
	auto wfun = [&](const Cell* c){
		int w = g->get_weight(c);
		if (w == 0) return 1e3;
		else return 1.0/g->get_weight(c);
	};
	return BGridImpose(g, wfun, source);
}

shared_ptr<BGrid> BGrid::ImposeBGrids(ShpVector<BGrid>& gg){
	if (gg.size() == 0) return shared_ptr<BGrid>();
	if (gg.size() == 1) return gg[0];
	_THROW_NOT_IMP_;
}

void BGrid::AddWeights(const std::map<const Cell*, int>& w){
	for (auto it: w){
		auto fnd = std::find_if(cells.begin(), cells.end(),
				[&it](shared_ptr<Cell> c){ return c.get() == it.first; });
		if (fnd != cells.end()) weight[it.first] = it.second;
	}
}
void BGrid::AddSourceFeat(const std::map<const Cell*, shared_ptr<int>>& f){
	for (auto it: f){
		auto fnd = std::find_if(cells.begin(), cells.end(),
				[&it](shared_ptr<Cell> c){ return c.get() == it.first; });
		if (fnd != cells.end()) source_feat[it.first] = it.second;
	}
}

void BGrid::ShallowAdd(const BGrid& g){
	GGeom::Modify::ShallowAdd(&g, this);
	AddWeights(g.weight);
	AddSourceFeat(g.source_feat);
}

void BGrid::ShallowAddCell(shared_ptr<Cell> c, const Cell* same_feat_cell){
	cells.push_back(c);
	if (same_feat_cell != 0){
		auto fnd1 = weight.find(same_feat_cell);
		auto fnd2 = source_feat.find(same_feat_cell);
		if (fnd1 != weight.end()) weight[c.get()] = fnd1->second;
		if (fnd2 != source_feat.end()) source_feat[c.get()] = fnd2->second;
	}
	set_indicies();
}

void BGrid::RemoveFeatures(const Cell* c){
	auto fnd1 = weight.find(c);
	auto fnd2 = source_feat.find(c);
	if (fnd1 != weight.end()) weight.erase(fnd1);
	if (fnd2 != source_feat.end()) source_feat.erase(fnd2);
}