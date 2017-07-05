#include "connectors.hpp"
#include "bgrid_impose.hpp"
#include "clipdomain.hpp"
#include "modcont.hpp"
#include "buildcont.hpp"
#include "treverter2d.hpp"
#include "buildgrid.hpp"
#include "assemble2d.hpp"
#include "healgrid.hpp"
#include "finder2d.hpp"

using namespace HMBlay::Impl;
// =============================== Factory
shared_ptr<MConnector> MConnector::Build(CornerTp tp, MappedMesher* prev, MappedMesher* next){
	switch (tp){
		case CornerTp::ACUTE:  return std::make_shared<AcuteConnector>(prev, next);
		case CornerTp::REENTRANT:  return std::make_shared<ReentrantConnector>(prev, next);
		case CornerTp::RIGHT: return std::make_shared<RightConnector>(prev, next);
		case CornerTp::ROUND: return std::make_shared<RoundConnector>(prev, next);
		case CornerTp::ZERO: case CornerTp::STRAIGHT:
			return std::make_shared<PlainConnector>(prev, next);
		default: assert(false);
	}
}

void MConnector::Add(BGrid& g, bool with_prev, bool with_next){
	//Simply adds cells. Imposition will be done in the postprocessing
	//routine
	if (with_prev) g.add_grid(prev->result);
	if (with_next) g.add_grid(next->result);
	//add connection grid if it exists
	BGrid* cg = ConnectionGrid();
	if (cg) g.add_grid(*cg);
}

// =============================== Plain
//TODO: set same source feature to next and prev mesh cells
void PlainConnector::ModifyAdjacents() {
	//make left points of next equal to right points of prev
	int imin = std::min(prev->left_points.size(), next->right_points.size());
	for (int i=0; i<imin; ++i){
		*next->left_points[i] = *prev->right_points[i];
	}
}

// ================================= Acute
void AcuteConnector::BuildInternals(){
	//small triangle is calculated at intersection of first level cells
	//it is used to fill acute corner.
	//1) find two cells which contain corner point
	Point pc = *HM2D::Contour::Last(prev->rect->BottomContour());
	auto get_cell = [](BGrid& g, Point& pc)->const HM2D::Cell*{
		for (auto& c: g.vcells)
		for (auto& e: c->edges){
			if (*e->first() == pc) return c.get();
			if (*e->last() == pc) return c.get();
		}
		return 0;
	};
	const HM2D::Cell* c1 = get_cell(prev->result, pc);
	const HM2D::Cell* c2 = get_cell(next->result, pc);
	if (c1 == 0 || c2 == 0) return;
	//2) intersect
	auto icont = HM2D::Contour::Clip::Intersection(c1->edges, c2->edges);
	HM2D::Contour::Clip::Heal(icont);
	HM2D::Contour::R::RevertTree::Permanent(icont);
	auto _av = HM2D::AllVertices(icont.alledges());
	auto _fn = HM2D::Finder::ClosestPoint(_av, pc);
	Point* pc2=_av[std::get<0>(_fn)].get();
	//return if intersection was not found.
	//Anyway we can safely proceed without this triangle
	if (icont.nodes.size() != 1 || *pc2 != pc) return;
	//3) get 3 points for a triangle
	auto pc2info = HM2D::Contour::PInfo(icont.nodes[0]->contour, pc2);
	HM2D::VertexData cl {pc2info.pprev, pc2info.p, pc2info.pnext};
	//4) build grid with one triangle
	filler = BGrid::MoveFrom1(HM2D::Grid::Constructor::FromTab(cl, {{0, 1, 2}}));
	filler->vedges[0]->boundary_type = pc2info.eprev->boundary_type;
	filler->vedges[1]->boundary_type = pc2info.enext->boundary_type;
	//5) set highest priority
	filler->weight[filler->vcells[0].get()] = 0;
};

void AcuteConnector::ModifyAdjacents(){
	//point of intersection to grid point
	auto crosses = HM2D::Contour::Finder::CrossAll(prev->rect->TopContour(),
			next->rect->TopContour());
	if (crosses.size() == 0) return;
	auto prevcont = HM2D::ECol::Assembler::GridBoundary(prev->result);
	auto nextcont = HM2D::ECol::Assembler::GridBoundary(next->result);
	auto prevcontv = HM2D::AllVertices(prevcont);
	auto nextcontv = HM2D::AllVertices(nextcont);
	for (auto c: crosses){
		Point p = std::get<1>(c);
		auto fp1 = HM2D::Finder::ClosestPoint(prevcontv, p);
		auto fp2 = HM2D::Finder::ClosestPoint(nextcontv, p);
		Point* p1 = prevcontv[std::get<0>(fp1)].get();
		Point* p2 = nextcontv[std::get<0>(fp2)].get();
		p = Point::Weigh(*p1, *p2, 0.5);
		p1->set(p.x, p.y);
		p2->set(p.x, p.y);
	}
}

// =========================== Right connector
RightConnector::RightConnector(MappedMesher* _prev, MappedMesher* _next): MConnector(_prev, _next){
	//find cross point between top lines of previous and
	//next mesh_areas
	HM2D::EdgeData prev_top = prev->rect->TopContour();
	HM2D::EdgeData next_top = next->rect->TopContour();
	//first cross from start of next_top contour (essential)
	std::tuple<bool, Point, double, double> cross =
		HM2D::Contour::Finder::Cross(next_top, prev_top);
	assert(std::get<0>(cross));
	//find crosspoint in conformal rectangle
	Point pcross_prev = prev->rect->MapToSquare(std::get<1>(cross));
	Point pcross_next = next->rect->MapToSquare(std::get<1>(cross));
	//cut meshed area by found weights
	prev->wend = pcross_prev.x;
	next->wstart = pcross_next.x;
}

void RightConnector::BuildInternals() {
	Point realpoint1 = prev->rect->MapToReal(Point(prev->wend, 0));
	Point realpoint2 = next->rect->MapToReal(Point(next->wstart, 0));
	//1) assemble borders
	HM2D::EdgeData _left = prev->rect->BottomContour();
	left = HM2D::Contour::Constructor::CutContour(_left, realpoint1, *HM2D::Contour::Last(_left)); 
	HM2D::Contour::R::ReallyRevert::Permanent(left);
	HM2D::EdgeData _bot = next->rect->BottomContour();
	bot = HM2D::Contour::Constructor::CutContour(_bot, *HM2D::Contour::First(_bot), realpoint2); 
	HM2D::EdgeData right = next->LeftContour();
	HM2D::EdgeData top = prev->RightContour();
	//get rid of numerical errors
	Point pc = (*HM2D::Contour::Last(right) + *HM2D::Contour::Last(top))/2.0;
	HM2D::Contour::Last(right)->set(pc);
	HM2D::Contour::Last(top)->set(pc);


	//2) assemble connector
	connection_area = MappedRect::Factory(left, right, bot, top,
			top.size()*right.size()*1.5, prev->rect->use_rect_approx());
	connection_grid.reset(new MappedMesher(connection_area.get(), 0.0, 1.0));

	//3) build a grid there
	auto f1out = HM2D::Contour::EWeights(top);
	for (auto& x: f1out) x = connection_area->top2bot(x); 
	auto bot_part=[&f1out](double, double)->vector<double>{ return f1out; };

	//using conform weights instead of real weighs to match adjacent partition
	vector<double> f2out = HM2D::Contour::EWeights(right);
	for(auto& x: f2out) x = connection_area->right2conf(x);
	auto vert_part=[&f2out](double)->vector<double>{ return f2out; };

	connection_grid->Fill(bot_part, vert_part, 1);

	//snap to top and right contours
	//1. top
	auto toppnt = HM2D::Contour::OrderedPoints(top);
	for (int i=0; i<toppnt.size(); ++i){
		int gridind = right.size()*toppnt.size() + i;
		connection_grid->result.vvert[gridind]->set(*toppnt[i]);
	}
	//2. right
	auto rpnt = HM2D::Contour::OrderedPoints(right);
	for (int i=0; i<rpnt.size(); ++i){
		int gridind = i*(top.size()+1) + top.size();
		connection_grid->result.vvert[gridind]->set(*rpnt[i]);
	}
}

// ============================ Reentrant connector
void ReentrantConnector::BuildInternals(){
	//1) assemble left/bot:
	left = prev->RightContour();
	bot = next->LeftContour();
	// left/bot are straight lines. Otherwise there was an error
	// during ExtPath division
	assert(ISEQ( HM2D::Contour::Length(left), Point::dist(*HM2D::Contour::First(left), *HM2D::Contour::Last(left)) ));
	assert(ISEQ( HM2D::Contour::Length(bot),  Point::dist(*HM2D::Contour::First(bot), *HM2D::Contour::Last(bot))  ));
	//int sz = std::min(left.size(), bot.size());
	//if (left.size()>sz) left = HMCont2D::ECollection::ShallowCopy(left, 0, sz);
	//if (bot.size()>sz) bot = HMCont2D::ECollection::ShallowCopy(bot, 0, sz);

	//2) assembling top/right
	// Find a corner point and draw lines from left/bot to it
	auto leftlast = HM2D::Contour::Last(left);
	auto leftfirst = HM2D::Contour::First(left);
	auto botlast = HM2D::Contour::Last(bot);
	auto botfirst = HM2D::Contour::First(bot);
	Vect v1 = *leftlast - *leftfirst;
	v1 = vecRotate(v1, 3*M_PI/2); 
	Vect v2 = *botlast - *botfirst;
	v2 = vecRotate(v2, M_PI/2);
	double ksieta[2];
	SectCross(*leftlast, *leftlast + v1, *botlast, *botlast + v2, ksieta);
	top_right_pnt = std::make_shared<HM2D::Vertex>(Point::Weigh(*leftlast, *leftlast + v1, ksieta[0]));
	top.clear(); right.clear();
	top.emplace_back(new HM2D::Edge ( leftlast, top_right_pnt));
	right.emplace_back(new HM2D::Edge ( botlast, top_right_pnt));

	//3) assemble connector
	bool force_rect = (
			prev->rect->use_rect_approx() &&
			HM2D::Contour::CornerPoints(left).size() == 2 &&
			HM2D::Contour::CornerPoints(right).size() == 2);
	connection_area.reset(
			new RectForOpenArea(left, right, bot, top, true, force_rect)
			);
	connection_grid.reset(new MappedMesher(connection_area.get(), 0.0, 1.0));

	//3) build a grid there
	auto bot_part=[&](double, double)->vector<double>{
		return HM2D::Contour::EWeights(bot);
	};

	//calculate before because this function is called multiple times
	//using conform weights instead of real weighs to match adjacent partition
	vector<double> f2out = HM2D::Contour::EWeights(left);
	for(auto& x: f2out) x = connection_area->left2conf(x);
	auto vert_part=[&f2out](double)->vector<double>{ return f2out; };

	connection_grid->Fill(bot_part, vert_part, 1);
}

// =========================== RoundConnector
void RoundConnector::BuildInternals(){
	assert(next->left_points.size()>1);
	assert(prev->right_points.size()>1);
	assert(*next->left_points[0] == *prev->right_points[0]);
	Point cp = *next->left_points[0];
	int Nlay = std::max(next->left_points.size(), prev->right_points.size()) - 1;
	assert(Nlay>0);

	//start and end points of circles
	vector<Point> pstart, pend;
	for (int i=0; i<Nlay; ++i) pstart.push_back(*prev->right_points[1+i]);
	for (int i=0; i<Nlay; ++i) pend.push_back(*next->left_points[1+i]);

	vector<HM2D::EdgeData> circ_start, circ_end;
	//build circles
	for (int i=0; i<Nlay; ++i){
		circ_start.push_back(
			HM2D::Contour::Constructor::Circle(36, cp, pstart[i])
		);
		circ_end.push_back(
			HM2D::Contour::Constructor::Circle(36, cp, pend[i])
		);
		HM2D::Contour::Algos::Reverse(circ_start.back());
		HM2D::Contour::Algos::Reverse(circ_end.back());
	}

	//contours from lowest to highest
	vector<HM2D::EdgeData> conts;
	for (int i=0; i<Nlay; ++i){
		conts.push_back(
			ContoursWeight(circ_start[i], pstart[i], circ_end[i], pend[i])
		);
	}

	//weights
	vector<double> w;
	double lennext = Point::dist(*next->left_points[Nlay],
			*next->left_points[Nlay-1]);
	double lenprev = Point::dist(*prev->right_points[Nlay],
			*prev->right_points[Nlay-1]);
	double lenav = (lennext+lenprev)/2.0;
	int N = std::ceil(HM2D::Contour::Length(conts.back()) / lenav);
	for (int i=0; i<=N; ++i) w.push_back((double)i/N);

	//get points
	vector<vector<Point>> allpts;
	for (auto& c: conts) allpts.push_back(HM2D::Contour::WeightPoints(c, w));

	//build grid points from geometric points
	vector<ShpVector<HM2D::Vertex>> pallpts(allpts.size());
	for (int i=0; i<allpts.size(); ++i){
		pallpts[i].resize(allpts[i].size());
		for (int j=0; j<allpts[i].size(); ++j)
			pallpts[i][j].reset(new HM2D::Vertex(allpts[i][j]));
		
	}
	shared_ptr<HM2D::Vertex> cpoint(new HM2D::Vertex(cp));

	//assemble procedure
	AssembleGrid(pallpts, cpoint);
}

void RoundConnector::AssembleGrid(vector<HM2D::VertexData>& pallpts, shared_ptr<HM2D::Vertex> cpoint){
	//assemble grid
	shared_ptr<int> feat(new int);
	HM2D::VertexData ap;
	for (auto& c: pallpts) ap.insert(ap.end(), c.begin(), c.end());
	ap.push_back(cpoint);
	aa::enumerate_ids_pvec(ap);
	vector<vector<int>> cell_vert;
	vector<int> cell_weights;
	vector<shared_ptr<int>> cell_feat;
	//rectangular cells
	for (int i=0; i<pallpts.size()-1; ++i){
		auto& pts1 = pallpts[i];
		auto& pts2 = pallpts[i+1];
		for (int j=0; j<pts1.size()-1; ++j){
			cell_vert.emplace_back();
			cell_vert.back().push_back(pts1[j]->id);
			cell_vert.back().push_back(pts1[j+1]->id);
			cell_vert.back().push_back(pts2[j+1]->id);
			cell_vert.back().push_back(pts2[j]->id);
			cell_weights.push_back(i+2);
			cell_feat.push_back(feat);
		}
	}
	//triangle cells
	for (int j=0; j<pallpts[0].size()-1; ++j){
		cell_vert.emplace_back();
		cell_vert.back().push_back(ap.back()->id);
		cell_vert.back().push_back(pallpts[0][j]->id);
		cell_vert.back().push_back(pallpts[0][j+1]->id);
		cell_weights.push_back(1);
		cell_feat.push_back(feat);
	}
	HM2D::Grid::Constructor::FixCellVert(ap, cell_vert);
	filler = BGrid::MoveFrom1(HM2D::Grid::Constructor::FromTab(ap, cell_vert));
	for (int i=0; i<filler->vcells.size(); ++i){
		filler->weight.emplace(filler->vcells[i].get(), cell_weights[i]);
		filler->source_feat.emplace(filler->vcells[i].get(), cell_feat[i]);
	}
	HM2D::Grid::Algos::Heal(*filler);
}
