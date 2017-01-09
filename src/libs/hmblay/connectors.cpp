#include "connectors.hpp"
#include "bgrid_impose.hpp"
#include "contclipping.hpp"
#include "algos.hpp"
#include "constructor.hpp"
#include "treverter2d.hpp"

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

void MConnector::Add(shared_ptr<BGrid>& g, bool with_prev, bool with_next){
	//Simply adds cells. Imposition will be done in the postprocessing
	//routine
	if (with_prev) g->ShallowAdd(prev->result);
	if (with_next) g->ShallowAdd(next->result);
	//add connection grid if it exists
	BGrid* cg = ConnectionGrid();
	if (cg) g->ShallowAdd(*cg);
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
	//small triangle is found at intersection of first level cells
	//it is used to fill acute corner.
	//1) find two cells which contain corner point
	Point pc = *HM2D::Contour::Last(prev->rect->BottomContour());
	auto get_cell = [](BGrid& g, Point& pc)->const Cell*{
		for (int i=0; i<g.n_cells(); ++i){
			const Cell* c = g.get_cell(i);
			for (int i=0; i<c->dim(); ++i){
				if (*c->points[i] == pc) return c;
			}
		}
		return 0;
	};
	const Cell* c1 = get_cell(prev->result, pc);
	const Cell* c2 = get_cell(next->result, pc);
	if (c1 == 0 || c2 == 0) return;
	//2) intersect
	auto cc1 = Cell2Cont(c1), cc2 = Cell2Cont(c2);
	auto icont = HM2D::Contour::Clip::Intersection(cc1, cc2);
	HM2D::Contour::Clip::Heal(icont);
	HM2D::Contour::R::RevertTree::Permanent(icont);
	auto _av = HM2D::AllVertices(icont.alledges());
	auto _fn = HM2D::FindClosestNode(_av, pc);
	Point* pc2=_av[std::get<0>(_fn)].get();
	//return if intersection was not found.
	//Anyway we can safely proceed without this triangle
	if (icont.nodes.size() != 1 || *pc2 != pc) return;
	//3) get 3 points for a triangle
	auto pc2info = HM2D::Contour::PInfo(icont.nodes[0]->contour, pc2);
	std::array<Point*, 3> cl {pc2info.pprev.get(), pc2info.p.get(), pc2info.pnext.get()};
	//4) build grid with one triangle
	filler.reset(new BGrid());
	GGeom::Modify::AddCell(*filler, {*cl[0], *cl[1], *cl[2]});
	//5) set highest priority
	filler->weight[filler->get_cell(0)] = 0;
};

namespace{
bool pntless(const Point* p1, const Point* p2){ return *p1 < *p2; }
}
void AcuteConnector::ModifyAdjacents(){
	//point of intersection to grid point
	auto crosses = HM2D::Contour::Algos::CrossAll(prev->rect->TopContour(),
			next->rect->TopContour());
	if (crosses.size() == 0) return;
	auto prevcont = GGeom::Info::Contour1(prev->result);
	auto nextcont = GGeom::Info::Contour1(next->result);
	auto prevset1 = prev->result.get_bnd_points();
	auto nextset1 = next->result.get_bnd_points();
	typedef std::set<const Point*, bool(*)(const Point* p1, const Point* p2)> Tpset;
	Tpset prevset2(prevset1.begin(), prevset1.end(), pntless);
	Tpset nextset2(nextset1.begin(), nextset1.end(), pntless);
	auto _pv = HM2D::AllVertices(prevcont);
	auto _pn = HM2D::AllVertices(nextcont);
	for (auto c: crosses){
		Point p = std::get<1>(c);
		auto _f1 = HM2D::FindClosestNode(_pv, p);
		auto _f2 = HM2D::FindClosestNode(_pn, p);
		Point* p1 = _pv[std::get<0>(_f1)].get();
		Point* p2 = _pn[std::get<0>(_f2)].get();
		p = Point::Weigh(*p1, *p2, 0.5);
		auto fnd1 = prevset2.find(p1);
		auto fnd2 = nextset2.find(p2);
		assert(fnd1 != prevset2.end() && fnd2 != nextset2.end());
		GridPoint* gp1 = const_cast<GridPoint*>(static_cast<const GridPoint*>(*fnd1));
		GridPoint* gp2 = const_cast<GridPoint*>(static_cast<const GridPoint*>(*fnd2));
		prevset2.erase(fnd1);
		nextset2.erase(fnd2);
		gp1->set(p.x, p.y);
		gp2->set(p.x, p.y);
		p1->set(p.x, p.y);
		p2->set(p.x, p.y);
		prevset2.insert(gp1);
		nextset2.insert(gp2);
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
		HM2D::Contour::Algos::Cross(next_top, prev_top);
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
		connection_grid->result.get_point(gridind)->set(*toppnt[i]);
	}
	//2. right
	auto rpnt = HM2D::Contour::OrderedPoints(right);
	for (int i=0; i<rpnt.size(); ++i){
		int gridind = i*(top.size()+1) + top.size();
		connection_grid->result.get_point(gridind)->set(*rpnt[i]);
	}
}

// ============================ Reentrant connector
void ReentrantConnector::BuildInternals(){
	//1) assemble left/bot:
	left = prev->RightContour();
	bot = next->LeftContour();
	// left/bot are straight lines. Otherwise there was an error
	// during ExtPath division
	assert(ISEQ( HM2D::Length(left), Point::dist(*HM2D::Contour::First(left), *HM2D::Contour::Last(left)) ));
	assert(ISEQ( HM2D::Length(bot),  Point::dist(*HM2D::Contour::First(bot), *HM2D::Contour::Last(bot))  ));
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
		HM2D::Contour::Reverse(circ_start.back());
		HM2D::Contour::Reverse(circ_end.back());
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
	int N = std::ceil(HM2D::Length(conts.back()) / lenav);
	for (int i=0; i<=N; ++i) w.push_back((double)i/N);

	//get points
	vector<vector<Point>> allpts;
	for (auto& c: conts) allpts.push_back(HM2D::Contour::WeightPoints(c, w));

	//build grid points from geometric points
	vector<ShpVector<GridPoint>> pallpts(allpts.size());
	for (int i=0; i<allpts.size(); ++i){
		pallpts[i].resize(allpts[i].size());
		for (int j=0; j<allpts[i].size(); ++j)
			pallpts[i][j].reset(new GridPoint(allpts[i][j]));
		
	}
	shared_ptr<GridPoint> cpoint(new GridPoint(cp));

	//assemble procedure
	AssembleGrid(pallpts, cpoint);
}

void RoundConnector::AssembleGrid(vector<ShpVector<GridPoint>>& pallpts, shared_ptr<GridPoint> cpoint){
	//assemble grid
	filler.reset(new BGrid());
	for (auto& c: pallpts) filler->ShallowAddNodes(c);
	shared_ptr<int> feat(new int);
	//rectangular cells
	for (int i=0; i<pallpts.size()-1; ++i){
		auto& pts1 = pallpts[i];
		auto& pts2 = pallpts[i+1];
		for (int j=0; j<pts1.size()-1; ++j){
			shared_ptr<Cell> cell(new Cell);
			cell->points.push_back(pts1[j].get());
			cell->points.push_back(pts1[j+1].get());
			cell->points.push_back(pts2[j+1].get());
			cell->points.push_back(pts2[j].get());
			filler->ShallowAddCell(cell);
			filler->weight[cell.get()] = i+2;
			filler->source_feat[cell.get()] = feat;
		}
	}
	//triangle cells
	filler->ShallowAddNode(cpoint);
	for (int j=0; j<pallpts[0].size()-1; ++j){
		shared_ptr<Cell> cell(new Cell);
		cell->points.push_back(cpoint.get());
		cell->points.push_back(pallpts[0][j].get());
		cell->points.push_back(pallpts[0][j+1].get());
		filler->ShallowAddCell(cell);
		filler->weight[cell.get()] = 1;
		filler->source_feat[cell.get()] = feat;
	}
	GGeom::Repair::Heal(*filler);
}
