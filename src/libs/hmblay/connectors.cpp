#include "connectors.hpp"
#include "bgrid_impose.hpp"

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
		*next->left_points.point(i) =
			*prev->right_points.point(i);
	}
}

// ================================= Acute
void AcuteConnector::BuildInternals(){
	//small triangle is found at intersection of first level cells
	//it is used to fill acute corner.
	//1) find two cells which contain corner point
	Point pc = *prev->rect->BottomContour().last();
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
	auto icont = HMCont2D::Clip::Intersection(cc1, cc2);
	Point* pc2=HMCont2D::ECollection::FindClosestNode(icont, pc);
	//return if intersection was not found.
	//Anyway we can safely proceed without this triangle
	if (icont.cont_count() != 1 || *pc2 != pc) return;
	//3) get 3 points for a triangle
	std::array<Point*, 3> cl = icont.nodes[0]->point_siblings(pc2);
	//4) build grid with one triangle
	filler.reset(new BGrid());
	GGeom::Modify::AddCell(*filler, {*cl[0], *cl[1], *cl[2]});
	//5) set highest priority
	filler->weight[filler->get_cell(0)] = 0;
};

void AcuteConnector::ModifyAdjacents(){
	//point of intersection to grid point
	auto crosses = HMCont2D::Contour::CrossAll(prev->rect->TopContour(),
			next->rect->TopContour());
	auto prevcont = GGeom::Info::Contour1(prev->result);
	auto nextcont = GGeom::Info::Contour1(next->result);
	for (auto c: crosses){
		Point p = std::get<1>(c);
		Point* p1 = HMCont2D::ECollection::FindClosestNode(prevcont, p);
		Point* p2 = HMCont2D::ECollection::FindClosestNode(nextcont, p);
		p = Point::Weigh(*p1, *p2, 0.5);
		p1->set(p.x, p.y);
		p2->set(p.x, p.y);
	}
}

// =========================== Right connector
RightConnector::RightConnector(MappedMesher* _prev, MappedMesher* _next): MConnector(_prev, _next){
	//find cross point between top lines of previous and
	//next mesh_areas
	HMCont2D::Contour prev_top = prev->rect->TopContour();
	HMCont2D::Contour next_top = next->rect->TopContour();
	//first cross from start of next_top contour (essential)
	std::tuple<bool, Point, double, double> cross =
		HMCont2D::Contour::Cross(next_top, prev_top);
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
	HMCont2D::Contour _left = prev->rect->BottomContour();
	left = HMCont2D::Constructor::CutContour(_left, realpoint1, *_left.last()); 
	left.ReallyReverse();
	HMCont2D::Contour _bot = next->rect->BottomContour();
	bot = HMCont2D::Constructor::CutContour(_bot, *_bot.first(), realpoint2); 
	HMCont2D::Contour right = next->LeftContour();
	HMCont2D::Contour top = prev->RightContour();

	//2) assemble connector
	connection_area = MappedRect::Factory(left, right, bot, top,
			prev->rect->use_rect_approx());
	connection_grid.reset(new MappedMesher(connection_area.get(), 0.0, 1.0));

	//3) build a grid there
	auto f1out = HMCont2D::Contour::EWeights(top);
	for (auto& x: f1out) x = connection_area->top2bot(x); 
	auto bot_part=[&f1out](double, double)->vector<double>{ return f1out; };

	//using conform weights instead of real weighs to match adjacent partition
	vector<double> f2out = HMCont2D::Contour::EWeights(right);
	for(auto& x: f2out) x = connection_area->right2conf(x);
	auto vert_part=[&f2out](double)->vector<double>{ return f2out; };

	connection_grid->Fill(bot_part, vert_part, 1);
}

// ============================ Reentrant connector
void ReentrantConnector::BuildInternals(){
	//1) assemble left/bot:
	left = prev->RightContour();
	bot = next->LeftContour();
	// left/bot are straight lines. Otherwise there was an error
	// during ExtPath division
	assert(ISEQ( left.length(), Point::dist(*left.first(), *left.last()) ));
	assert(ISEQ( bot.length(),  Point::dist(*bot.first(), *bot.last())  ));
	//why is it here?
	//int sz = std::min(left.size(), bot.size());
	//if (left.size()>sz) left = HMCont2D::ECollection::ShallowCopy(left, 0, sz);
	//if (bot.size()>sz) bot = HMCont2D::ECollection::ShallowCopy(bot, 0, sz);

	//2) assembling top/right
	// Find a corner point and draw lines from left/bot to it
	Vect v1 = *left.last() - *left.first();
	v1 = vecRotate(v1, 3*M_PI/2); 
	Vect v2 = *bot.last() - *bot.first();
	v2 = vecRotate(v2, M_PI/2);
	double ksieta[2];
	SectCross(*left.last(), *left.last() + v1, *bot.last(), *bot.last() + v2, ksieta);
	top_right_pnt = Point::Weigh(*left.last(), *left.last() + v1, ksieta[0]);
	top.clear(); right.clear();
	top.add_value(HMCont2D::Edge { left.last(), &top_right_pnt});
	right.add_value(HMCont2D::Edge { bot.last(), &top_right_pnt});

	//3) assemble connector
	bool force_rect = (left.is_straight() && right.is_straight() &&
			prev->rect->use_rect_approx());
	connection_area.reset(
			new RectForOpenArea(left, right, bot, top, true, force_rect)
			);
	connection_grid.reset(new MappedMesher(connection_area.get(), 0.0, 1.0));

	//3) build a grid there
	auto bot_part=[&](double, double)->vector<double>{
		return HMCont2D::Contour::EWeights(bot);
	};

	//calculate before because this function is called multiple times
	//using conform weights instead of real weighs to match adjacent partition
	vector<double> f2out = HMCont2D::Contour::EWeights(left);
	for(auto& x: f2out) x = connection_area->left2conf(x);
	auto vert_part=[&f2out](double)->vector<double>{ return f2out; };

	connection_grid->Fill(bot_part, vert_part, 1);
}

// =========================== RoundConnector
void RoundConnector::BuildInternals(){
	assert(next->left_points.size()>1);
	assert(prev->right_points.size()>1);
	assert(*next->left_points.point(0) == *prev->right_points.point(0));
	Point cp = *next->left_points.point(0);
	int Nlay = std::max(next->left_points.size(), prev->right_points.size()) - 1;
	assert(Nlay>0);

	//start and end points of circles
	vector<Point> pstart, pend;
	for (int i=0; i<Nlay; ++i) pstart.push_back(*prev->right_points.point(1+i));
	for (int i=0; i<Nlay; ++i) pend.push_back(*next->left_points.point(1+i));

	vector<HMCont2D::Container<HMCont2D::Contour>> circ_start, circ_end;
	//build circles
	for (int i=0; i<Nlay; ++i){
		circ_start.push_back(
			HMCont2D::Constructor::Circle(36, cp, pstart[i])
		);
		circ_end.push_back(
			HMCont2D::Constructor::Circle(36, cp, pend[i])
		);
		circ_start.back().Reverse();
		circ_end.back().Reverse();
	}

	//contours from lowest to highest
	vector<HMCont2D::Container<HMCont2D::Contour>> conts;
	for (int i=0; i<Nlay; ++i){
		conts.push_back(
			ContoursWeight(circ_start[i], pstart[i], circ_end[i], pend[i])
		);
	}

	//weights
	vector<double> w;
	double lennext = Point::dist(next->left_points.value(Nlay),
			next->left_points.value(Nlay-1));
	double lenprev = Point::dist(prev->right_points.value(Nlay),
			prev->right_points.value(Nlay-1));
	double lenav = (lennext+lenprev)/2.0;
	int N = std::ceil(conts.back().length() / lenav);
	for (int i=0; i<=N; ++i) w.push_back((double)i/N);

	//get points
	vector<HMCont2D::PCollection> allpts;
	for (auto& c: conts) allpts.push_back(HMCont2D::Contour::WeightPoints(c, w));

	//build grid points from geometric points
	vector<ShpVector<GridPoint>> pallpts(allpts.size());
	for (int i=0; i<allpts.size(); ++i){
		pallpts[i].resize(allpts[i].size());
		for (int j=0; j<allpts[i].size(); ++j)
			pallpts[i][j].reset(new GridPoint(*allpts[i].point(j)));
		
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




