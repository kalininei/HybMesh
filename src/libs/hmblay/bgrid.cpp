#include "bgrid.hpp"
#include "canonic_bgrid.hpp"
#include "bgrid_impose.hpp"
#include "connectors.hpp"
#include "hmtimer.hpp"
#include "contclipping.hpp"
#include "debug2d.hpp"
#include "healgrid.hpp"
#include "modgrid.hpp"

using namespace HMBlay::Impl;

BGrid BGrid::MeshFullPath(const ExtPath& epath){
	//1. divide by angles
	vector<ExtPath> pths = ExtPath::DivideByAngle(epath,
			{CornerTp::RIGHT,
			CornerTp::REENTRANT, 
			CornerTp::ACUTE,
			CornerTp::ROUND});
	//if closed path with one special corner divide it into two
	if (HM2D::Contour::IsClosed(epath) && pths.size() == 1 &&
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
		//estimate vertical and horizontal partition
		int isz = p.PathPartition(0, HM2D::Length(p)).size();
		int jsz = p.largest_vpart_size();

		double h = p.largest_depth();
		mps.push_back(MappedRect::Factory(p.leftbc, p.rightbc, p, h,
					isz*jsz*1.5, use_rect_approx));
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
	if (HM2D::Contour::IsClosed(epath)){
		CornerTp t = pths[0].ext_data[0].tp;
		if (mps.size() > 1){
			connectors.push_back(MConnector::Build(t, mesher4.back().get(), mesher4[0].get()));
		}
	}

	//5. build rectangular meshes
	for (int i=0; i<mesher4.size(); ++i){
		ExtPath& ipth = pths[i];
		double ilen = HM2D::Length(ipth);
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
	for (auto& c: connectors){
		c->Apply();
	}
	
	//7. Gather all resulting meshes
	BGrid g;
	if (connectors.size() > 0){
		for (auto& c: connectors){
			//add_previous with grid is empty
			bool with_prev = (g.vcells.size() == 0);
			//add next always except if this
			//is an end connector for closed path
			bool with_next = (!HM2D::Contour::IsClosed(epath) || &c != &connectors.back());
			c->Add(g, with_prev, with_next);
		}
		HM2D::Grid::Algos::Heal(g);
	} else {
		assert(mesher4.size() == 1);
		g = std::move(mesher4[0]->result);
	}

	return g;
}

BGrid BGrid::MeshSequence(vector<Options*>& data){ 
	auto ret = BGrid();

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

BGrid BGrid::NoSelfIntersections(BGrid& g, const HM2D::EdgeData& source){
	//does grid have intersections
	if (HM2D::Grid::Algos::Check(g)) return g;
	//else do imposition on the basis of cells weights
	auto wfun = [&](const HM2D::Cell* c){
		int w = g.get_weight(c);
		if (w == 0) return 1e3;
		else return 1.0/g.get_weight(c);
	};
	BGridImpose(g, wfun, source);
	return g;
}

BGrid BGrid::ImposeBGrids(ShpVector<BGrid>& gg){
	if (gg.size() == 0) return BGrid();
	//we can omit self intersections checks for single grid
	//because it should be done in ::NoSelfIntersections procedure before
	if (gg.size() == 1) return *gg[0];
	//straight grid addition. This may cause new self intersections which
	//should be treated correctly
	bool has_intersections = false;
	for (int i=1; i<(int)gg.size(); ++i){
		if (!has_intersections){
			auto c1 = HM2D::Contour::Tree::GridBoundary(*gg[0]);
			auto c2 = HM2D::Contour::Tree::GridBoundary(*gg[i]);
			auto inters = HM2D::Contour::Clip::Intersection(c1, c2);
			if (inters.nodes.size() > 0) has_intersections = true;
		}
		gg[0]->add_grid(*gg[i]);
	}
	if (!has_intersections) return std::move(*gg[0]);
	//TODO: here should be a procedure which
	//      make impositions if grid has self intersections
	std::cout<<"Self intersected grids treatment is not done yet"<<std::endl;
	_THROW_NOT_IMP_;
}

void BGrid::add_weights(const std::map<const HM2D::Cell*, int>& w){
	for (auto& it: w) it.first->id = -1;
	aa::constant_ids_pvec(vcells, 0);
	for (auto& it: w) if (it.first->id == 0){
		weight.insert(it);
	}
}
void BGrid::add_source_feat(const std::map<const HM2D::Cell*, shared_ptr<int>>& f){
	for (auto& it: f) it.first->id = -1;
	aa::constant_ids_pvec(vcells, 0);
	for (auto& it: f) if (it.first->id == 0){
		source_feat.insert(it);
	}
}

void BGrid::remove_cells(const vector<int>& bad_cells){
	for (int ic: bad_cells){
		weight.erase(vcells[ic].get());
		source_feat.erase(vcells[ic].get());
	}
	HM2D::Grid::Algos::RemoveCells(*this, bad_cells);
}

void BGrid::add_grid(const BGrid& g){
	vvert.insert(vvert.end(), g.vvert.begin(), g.vvert.end());
	vedges.insert(vedges.end(), g.vedges.begin(), g.vedges.end());
	vcells.insert(vcells.end(), g.vcells.begin(), g.vcells.end());

	add_weights(g.weight);
	add_source_feat(g.source_feat);
}

void BGrid::add_cell(shared_ptr<HM2D::Cell> c, const HM2D::Cell* same_feat_cell){
	vcells.push_back(c);
	if (same_feat_cell != 0){
		auto fnd1 = weight.find(same_feat_cell);
		auto fnd2 = source_feat.find(same_feat_cell);
		if (fnd1 != weight.end()) weight[c.get()] = fnd1->second;
		if (fnd2 != source_feat.end()) source_feat[c.get()] = fnd2->second;
	}
}

void BGrid::remove_features(const HM2D::Cell* c){
	auto fnd1 = weight.find(c);
	auto fnd2 = source_feat.find(c);
	if (fnd1 != weight.end()) weight.erase(fnd1);
	if (fnd2 != source_feat.end()) source_feat.erase(fnd2);
}

shared_ptr<BGrid> BGrid::MoveFrom1(HM2D::GridData&& gg){
	shared_ptr<BGrid> ret = std::make_shared<BGrid>();
	ret->vcells = std::move(gg.vcells);
	ret->vedges = std::move(gg.vedges);
	ret->vvert = std::move(gg.vvert);
	return ret;
}

BGrid BGrid::MoveFrom2(HM2D::GridData&& gg){
	BGrid ret;
	ret.vcells = std::move(gg.vcells);
	ret.vedges = std::move(gg.vedges);
	ret.vvert = std::move(gg.vvert);
	return ret;
}
