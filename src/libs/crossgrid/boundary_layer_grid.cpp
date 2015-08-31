#include "boundary_layer_grid.h"
namespace{

GridGeom BuildSquareGrid(const vector<double>& part_x,vector<double>& part_y){
	int Npts = part_x.size()*part_y.size();
	int Ncls = (part_x.size()-1)*(part_y.size()-1);
	vector<double> points;
	for (int j=0; j<part_y.size(); ++j){
		for (int i=0; i<part_x.size(); ++i){
			points.push_back(part_x[i]);
			points.push_back(part_y[j]);
		}
	}
	vector<int> cells;
	for (int j=0; j<part_y.size() - 1; ++j){
		for (int i=0; i<part_x.size() - 1; ++i){
			int i0 = j*part_x.size() + i;
			int i1 = j*part_x.size() + i + 1;
			int i2 = (j + 1) * part_x.size() + i + 1;
			int i3 = (j + 1) * part_x.size() + i;
			cells.push_back(4);
			cells.push_back(i0);
			cells.push_back(i1);
			cells.push_back(i2);
			cells.push_back(i3);
		}
	}
	return GridGeom(Npts, Ncls, &points[0], &cells[0]);
}

}

BLayerGrid::BLayerGrid(const BLayerGridInput& opt): GridGeom(){
	std::cout<<"DUMMY BLayerGrid"<<std::endl;
	//0) check if options are valid
	if (opt.round_off || opt.bnd_step_method != BLayerGridInput::CONST_BND_STEP)
		throw std::runtime_error("Not implemented option");
	//1) get square sizes
	BoundingBox b = opt.tree->BuildBoundingBox();
	Point a0 = b.BottomLeft();
	Point a1 = b.TopRight();
	//2) bnd steps: from zero
	vector<double> h = opt.partition;

	//3) build a grid around a square
	CGBoundingBox innerb(0, 0, 0, 0);
	vector<double> partx, party;
	if (opt.direction == INSIDE){
		double lenx = a1.x - a0.x;
		double leny = a1.y - a0.y;
		double mainx = lenx - 2*h.back();
		double mainy = leny - 2*h.back();
		if (mainx<=geps || mainy<=geps)
			throw std::runtime_error("Delta is too big");
		Point r = a1 - a0;
		//x
		for (auto x: h) partx.push_back(x);
		double hx = mainx / int(mainx/opt.bnd_step);
		for (int i=0; i<mainx/hx-1;++i) partx.push_back(partx.back() + hx);
		for (int i=h.size()-1; i>=0; --i) partx.push_back(r.x - h[i]);
		for (auto& x: partx) x+=a0.x;
		//y
		for (auto x: h) party.push_back(x);
		double hy = mainy / int(mainy/opt.bnd_step);
		for (int i=0; i<mainy/hy-1;++i) party.push_back(party.back() + hy);
		for (int i=h.size()-1; i>=0; --i) party.push_back(r.y - h[i]);
		for (auto& x: party) x+=a0.y;
		//bounding box
		innerb = CGBoundingBox(partx[0] + h.back(), party[0] + h.back(),
				partx.back() - h.back(), party.back() - h.back());
	} else {
		double lenx = a1.x - a0.x;
		double leny = a1.y - a0.y;
		double mainx = lenx + 2*h.back();
		double mainy = leny + 2*h.back();
		Point r = a1 - a0;
		//x
		partx = h;
		std::reverse(partx.begin(), partx.end());
		for (auto& a: partx) a = -(a-h.back());
		double hx = lenx / int(mainx/opt.bnd_step);
		for (int i=0; i<lenx/hx-1;++i) partx.push_back(partx.back() + hx);
		for (int i=0; i<h.size(); ++i) partx.push_back(lenx + h.back() + h[i]);
		for (auto& x: partx) x+=a0.x - h.back();
		//y
		party = h;
		std::reverse(party.begin(), party.end());
		for (auto& a: party) a = -(a-h.back());
		double hy = leny / int(mainy/opt.bnd_step);
		for (int i=0; i<leny/hy-1;++i) party.push_back(party.back() + hy);
		for (int i=0; i<h.size(); ++i) party.push_back(leny + h.back() + h[i]);
		for (auto& x: party) x+=a0.y - h.back();
		//bounding box
		innerb = CGBoundingBox(partx[0] + h.back(), party[0] + h.back(),
				partx.back() - h.back(), party.back() - h.back());
	}
	//build filled grids
	GridGeom g = BuildSquareGrid(partx, party);
	swap_data(*this, g);
	//remove all superfluous cells
	vector<int> badcells;
	for (auto c: cells){
		Point p = (*c->get_point(0) + *c->get_point(1) + 
			*c->get_point(2) + *c->get_point(3) ) / 4; 
		if (innerb.whereis(p) == INSIDE)
			badcells.push_back(c->get_ind());
	}
	remove_cells(badcells);
}





