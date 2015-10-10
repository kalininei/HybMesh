#include "bgrid.hpp"
#include "canonic_bgrid.hpp"
#include "trigrid.h"
#include "bgrid_impose.hpp"

using namespace HMBlay::Impl;

namespace{
struct MConnector{
private:
	virtual void ModifyAdjacents(){};
	virtual void BuildInternals(){};
protected:
	MappedMesher *prev, *next;
	MConnector(MappedMesher* _prev, MappedMesher* _next): prev(_prev), next(_next){}
public:
	static shared_ptr<MConnector>
	Build(CornerTp tp, MappedMesher* prev, MappedMesher* next);
	
	//modifies prev, next meshes,
	//build connection mesh if nessessary.
	void Apply(){ ModifyAdjacents(); BuildInternals();}

	//adds next grid along with connection section to g
	virtual void Add(shared_ptr<BGrid>& g, bool with_prev, bool with_next) = 0;
};

struct PlainConnector: public MConnector{
//connects grids which have congruent left|right nodes
private:
	//TODO: set same source feature to next and prev mesh cells
	void ModifyAdjacents() override{
		//make left points of next equal to right points of prev
		int imin = std::min(prev->left_points.size(), next->right_points.size());
		for (int i=0; i<imin; ++i){
			*next->left_points.point(i) =
				*prev->right_points.point(i);
		}
	}
public:
	PlainConnector(MappedMesher* prev, MappedMesher* next): MConnector(prev, next){};
	void Add(shared_ptr<BGrid>& g, bool with_prev, bool with_next){
		if (with_prev) g->ShallowAdd(prev->result);
		if (with_next) g->ShallowAdd(next->result);
	}
};

struct SharpConnector: public MConnector{
private:
	shared_ptr<BGrid> filler;
	void BuildInternals(){
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

	void ModifyAdjacents() override{
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
public:
	SharpConnector(MappedMesher* _prev, MappedMesher* _next): MConnector(_prev, _next){}

	void Add(shared_ptr<BGrid>& g, bool with_prev, bool with_next) override{
		//Simply adds cells. Imposition will be done in the postprocessing
		//routine
		if (with_prev) g->ShallowAdd(prev->result);
		if (with_next) g->ShallowAdd(next->result);
		//triangle cell to fill the acute angle with highest priority
		if (filler) g->ShallowAdd(*filler);
	}
};


struct FillerConnector: public MConnector{
protected:
	shared_ptr<MappedMesher> connection_grid;
	shared_ptr<MappedRect> connection_area;
	bool use_rect_approx;
	FillerConnector(MappedMesher* _prev, MappedMesher* _next): MConnector(_prev, _next),
		use_rect_approx(_prev->rect->use_rect_approx()){}
public:
	void Add(shared_ptr<BGrid>& g, bool with_prev, bool with_next) override{
		//1) add previous grid if it has not been add yet
		if (with_prev) g->ShallowAdd(prev->result);
		//2) add connection grid
		g->ShallowAdd(connection_grid->result);
		//3) add next grid
		if (with_next) g->ShallowAdd(next->result);
	}
};

struct RightConnector: public FillerConnector{
private:
	HMCont2D::Container<HMCont2D::Contour> left, bot;
	HMCont2D::Contour right, top;
	void BuildInternals() override {
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
				use_rect_approx);
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
public:
	RightConnector(MappedMesher* _prev, MappedMesher* _next): FillerConnector(_prev, _next){
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
};

struct ObtuseConnector: public FillerConnector{
private:
	Point top_right_pnt;
	HMCont2D::Contour left, bot, right, top;
	void BuildInternals() override {
		//1) assemble left/bot:
		// left/bot are straight lines. Otherwise there was an error
		// during ExtPath division
		left = prev->RightContour();
		bot = next->LeftContour();
		assert(ISEQ( left.length(), Point::dist(*left.first(), *left.last()) ));
		assert(ISEQ( bot.length(),  Point::dist(*bot.first(), *bot.last())  ));
		int sz = std::min(left.size(), bot.size());
		if (left.size()>sz) left = HMCont2D::ECollection::ShallowCopy(left, 0, sz);
		if (bot.size()>sz) bot = HMCont2D::ECollection::ShallowCopy(bot, 0, sz);

		//2) assembling top/right
		// Find a corner point and draw lines from left/bot to it
		Vect v1 = *left.last() - *left.first();
		v1 = vecRotate(v1, 3*M_PI/2); 
		Vect v2 = *bot.last() - *bot.first();
		v2 = vecRotate(v2, M_PI/2);
		double ksieta[2];
		SectCross(*left.last(), *left.last() + v1, *bot.last(), *bot.last() + v2, ksieta);
		if (ksieta[0]<2 && ksieta[1]<2){
			top_right_pnt = Point::Weigh(*left.last(), *left.last() + v1, ksieta[0]);
		} else {
			//Circle grid is needed
			_THROW_NOT_IMP_;
		}
		top.clear(); right.clear();
		top.add_value(HMCont2D::Edge { left.last(), &top_right_pnt});
		right.add_value(HMCont2D::Edge { bot.last(), &top_right_pnt});

		//3) assemble connector
		connection_area = MappedRect::Factory(left, right, bot, top, use_rect_approx);
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
public:
	ObtuseConnector(MappedMesher* _prev, MappedMesher* _next): FillerConnector(_prev, _next){
	}
};

shared_ptr<MConnector> MConnector::Build(CornerTp tp, MappedMesher* prev, MappedMesher* next){
	switch (tp){
		case CornerTp::SHARP:  return std::make_shared<SharpConnector>(prev, next);
		case CornerTp::OBTUSE:  return std::make_shared<ObtuseConnector>(prev, next);
		case CornerTp::CORNER: return std::make_shared<RightConnector>(prev, next);
		case CornerTp::ZERO: case CornerTp::REGULAR:
			return std::make_shared<PlainConnector>(prev, next);
		default: assert(false);
	}
}

}

shared_ptr<BGrid> BGrid::MeshFullPath(const ExtPath& epath){
	//1. divide by angles
	vector<ExtPath> pths = ExtPath::DivideByAngle(epath,
			{CornerTp::CORNER, CornerTp::OBTUSE, CornerTp::SHARP});
	//if closed path with one special corner divide it into two
	if (epath.is_closed() && pths.size() == 1 &&
			(pths[0].ext_data.back().tp == CornerTp::CORNER ||
			 pths[0].ext_data.back().tp == CornerTp::OBTUSE ||
			 pths[0].ext_data.back().tp == CornerTp::SHARP)){
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
