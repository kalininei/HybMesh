#include "infogrid.hpp"
#include "tree.hpp"

using namespace HM2D;
namespace hg=HM2D::Grid;

double hg::Area(const GridData& grid){
	auto vt = Contour::Tree::GridBoundary01(grid);
	double ret = 0;
	for (auto it: vt) ret += it.area();
	return ret;
}

vector<double> hg::CellAreas(const GridData& grid){
	vector<double> ret(grid.vcells.size());
	for (int i=0; i<ret.size(); ++i){
		ret[i] = HM2D::Contour::Area(grid.vcells[i]->edges);
	}
	return ret;
}

//calculate skewness
vector<double> hg::Skewness(const GridData& grid){
	vector<double> ret(grid.vcells.size());
	for (int i=0; i<grid.vcells.size(); ++i){
		const Cell& c = *grid.vcells[i];
		int dim = c.edges.size();
		if (dim < 3){
			ret[i] = 1.0;  //very bad cell anyway
			continue;
		}
		vector<double> angles(dim);
		auto op = Contour::OrderedPoints(c.edges);
		for (int j=1; j<dim; ++j){
			const Point& p0 = *op[j-1];
			const Point& p1 = *op[j];
			const Point& p2 = *op[j+1];
			angles[j] = Angle(p0, p1, p2);
		}
		angles[0] = M_PI *(dim - 2) - 
			std::accumulate(angles.begin() + 1, angles.begin() + dim, 0.0);
		auto minmax = std::minmax_element(angles.begin(), angles.end());
		double minv = *minmax.first;
		double maxv = *minmax.second;
		double refan = M_PI * (dim - 2) / dim;
		ret[i] = std::max( (maxv-refan)/(M_PI-refan), (refan-minv)/refan );
		if (ret[i] > 1.0) ret[i] = 1.0;   //for non-convex cells
	}
	return ret;
}
