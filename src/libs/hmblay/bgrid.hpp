#ifndef HYBMESH_BGRID_HPP
#define HYBMESH_BGRID_HPP

#include "options.hpp"
#include "extpath.hpp"
#include "primitives2d.hpp"

namespace HMBlay{
namespace Impl{

class BGrid: public HM2D::GridData{
	static ExtPath AssembleExtendedPath(vector<Options*>& data);
	static BGrid NoSelfIntersections(BGrid& g, const HM2D::EdgeData& source);
	static BGrid MeshFullPath(const ExtPath& p);

public:
	//includes deletion of features, edges and vertices
	void remove_cells(const vector<int>& bad_cells);
	//adds cell to the end of cell list.
	//gets features (weight, source_feat) from same_feat_cell
	//!!! does not add edges and vertices.
	void add_cell(shared_ptr<HM2D::Cell> c, const HM2D::Cell* same_feat_cell=0);
	void add_grid(const BGrid& g);

	//Removes all source features
	void remove_features(const HM2D::Cell* c);

	//index of layer starting from source
	//used for bgrid imposition to calculate cell priority
	std::map<const HM2D::Cell*, int> weight;
	int get_weight(const HM2D::Cell* c){
		auto fnd = weight.find(c);
		return (fnd == weight.end()) ? 1e3 : fnd->second;
	}
	int get_weight(int i){ return get_weight(vcells[i].get()); }
	void add_weights(const std::map<const HM2D::Cell*, int>& w);

	//all cells which were created from the same source
	//has unique address at source_feat. Value by itself doesn't matter.
	//HM2D::Cells are considered to belong to same source if
	//source_feat[c1].get() == source_feat[c2].get()
	std::map<const HM2D::Cell*, shared_ptr<int>> source_feat;
	void add_source_feat(const std::map<const HM2D::Cell*, shared_ptr<int>>& f);

	bool is_from_source(const HM2D::Cell* c1) const{
		return source_feat.find(c1) != source_feat.end();
	}
	bool is_from_same_source(const HM2D::Cell* c1, const HM2D::Cell* c2) const{
		auto fnd1 = source_feat.find(c1);
		auto fnd2 = source_feat.find(c2);
		if (fnd1 == source_feat.end() || 
			fnd2 == source_feat.end()) return false;
		else return (fnd1->second.get() == fnd2->second.get());
	}
	bool is_from_same_source (int i1, int i2) const{
		return is_from_same_source(vcells[i1].get(), vcells[i2].get());
	}
	
	
	static BGrid MeshSequence(vector<Options*>& data);
	static BGrid ImposeBGrids(ShpVector<BGrid>& gg);
	static shared_ptr<BGrid> MoveFrom1(HM2D::GridData&& gg);
	static BGrid MoveFrom2(HM2D::GridData&& gg);
};



}//Impl

}//HMBlay



#endif
