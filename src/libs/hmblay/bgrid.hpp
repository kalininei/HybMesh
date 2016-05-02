#ifndef HYBMESH_BGRID_HPP
#define HYBMESH_BGRID_HPP

#include "grid.h"
#include "hybmesh_contours2d.hpp"
#include "options.hpp"
#include "extpath.hpp"

namespace HMBlay{
namespace Impl{

class BGrid: public GridGeom{
protected:
	typedef HMCont2D::Contour Pth;
	typedef HMCont2D::Edge Ed;

	static ExtPath AssembleExtendedPath(vector<Options*>& data);
	static shared_ptr<BGrid> NoSelfIntersections(shared_ptr<BGrid> g, const HMCont2D::Contour& source);
	static shared_ptr<BGrid> MeshFullPath(const ExtPath& p);

	//remove cells overriden to include deletion of features
	void remove_cells(const vector<int>& bad_cells) override;
public:
	//adds cell to the end of cell list.
	//gets features (weight, source_feat) from same_feat_cell
	void ShallowAddCell(shared_ptr<Cell> c, const Cell* same_feat_cell=0);
	void ShallowAddNode(shared_ptr<GridPoint> p) { points.push_back(p); set_indicies(); }
	void ShallowAddNodes(const ShpVector<GridPoint>& p) {
		std::copy(p.begin(), p.end(), std::back_inserter(points));
		set_indicies();
	}

	//Removes all source features
	void RemoveFeatures(const Cell* c);

	//index of layer starting from source
	//used for bgrid imposition to calculate cell priority
	std::map<const Cell*, int> weight;
	int get_weight(const Cell* c){
		auto fnd = weight.find(c);
		return (fnd == weight.end()) ? 1e3 : fnd->second;
	}
	int get_weight(int i){ return get_weight(get_cell(i)); }
	void AddWeights(const std::map<const Cell*, int>& w);

	//all cells which were created from the same source
	//has unique address at source_feat. Value by itself doesn't matter.
	//Cells are considered to belong to same source if
	//source_feat[c1].get() == source_feat[c2].get()
	std::map<const Cell*, shared_ptr<int>> source_feat;
	void AddSourceFeat(const std::map<const Cell*, shared_ptr<int>>& f);

	bool is_from_source(const Cell* c1) const{
		return source_feat.find(c1) != source_feat.end();
	}
	bool is_from_same_source(const Cell* c1, const Cell* c2) const{
		auto fnd1 = source_feat.find(c1);
		auto fnd2 = source_feat.find(c2);
		if (fnd1 == source_feat.end() || 
			fnd2 == source_feat.end()) return false;
		else return (fnd1->second.get() == fnd2->second.get());
	}
	bool is_from_same_source (int i1, int i2) const{
		return is_from_same_source(get_cell(i1), get_cell(i2));
	}
	
	void ShallowAdd(const BGrid& g);
	
	static shared_ptr<BGrid> MeshSequence(vector<Options*>& data);
	static shared_ptr<BGrid> ImposeBGrids(ShpVector<BGrid>& gg);

};



}//Impl

}//HMBlay



#endif
