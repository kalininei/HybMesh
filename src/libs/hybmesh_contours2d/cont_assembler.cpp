#include "cont_assembler.hpp"

using namespace HMCont2D;
namespace cns = Assembler;

std::vector<Contour> cns::AllContours(const HMCont2D::ECollection& input){
	auto ap = input.all_points();
	std::set<Point*> unusedpnts(ap.begin(), ap.end());
	//1) assemble all possible contours
	vector<Contour> conts;
	while (unusedpnts.size() > 0){
		Point* p1 = *unusedpnts.begin();
		conts.push_back(Contour::Assemble(input, p1));
		for (Point* p: conts.back().all_points()) unusedpnts.erase(p);
	}
	return conts;
}
