#include "bgrid.hpp"
#include "fileproc.h"
#include "canonic_bgrid.hpp"

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
		if (with_prev) GGeom::Modify::ShallowAdd(&prev->result, g.get());
		if (with_next) GGeom::Modify::ShallowAdd(&next->result, g.get());
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
		if (with_prev) GGeom::Modify::ShallowAdd(&prev->result, g.get());
		//2) add connection grid
		GGeom::Modify::ShallowAdd(&connection_grid->result, g.get());
		//3) add next grid
		if (with_next) GGeom::Modify::ShallowAdd(&next->result, g.get());
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

		connection_grid->Fill(bot_part, vert_part);
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

		connection_grid->Fill(bot_part, vert_part);
	}
public:
	ObtuseConnector(MappedMesher* _prev, MappedMesher* _next): FillerConnector(_prev, _next){
	}
};

struct SharpConnector: public MConnector{
	SharpConnector(MappedMesher* _prev, MappedMesher* _next): MConnector(_prev, _next){
		_THROW_NOT_IMP_;
	}
	void Add(shared_ptr<BGrid>& g, bool with_prev, bool with_next) override{ _THROW_NOT_IMP_; }
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
		mesher4[i]->Fill(bpart, wpart);
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

	//4) guarantee no self-intersections
	ret = NoSelfIntersections(ret);


	return ret;
}

shared_ptr<BGrid> BGrid::NoSelfIntersections(shared_ptr<BGrid> g){
	return g;
	auto cont = GGeom::Info::Contour(*g);
	//check for no self-contacts in resulting tree
	assert( cont.cont_count() > 0 && HMCont2D::ContourTree::CheckNoContact(cont) );
	//area as sum of cells areas
	double area1 = g->area();
	//area as bnd area
	double area2 = HMCont2D::ContourTree::Area(cont);
	if (ISEQ(area1, area2)) return g;
	else{
		_THROW_NOT_IMP_;
	}
}

shared_ptr<BGrid> BGrid::ImposeBGrids(ShpVector<BGrid>& gg){
	if (gg.size() == 0) return shared_ptr<BGrid>();
	if (gg.size() == 1) return gg[0];
	_THROW_NOT_IMP_;
}
