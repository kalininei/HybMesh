#include "bgrid.hpp"
#include "simple_bgrid.hpp"
#include "canonic_bgrid.hpp"

using namespace HMBlay::Impl;

namespace{
struct MConnector{
private:
	virtual void ModifyAdjecents(){};
	virtual void BuildInternals(){};
protected:
	MappedMesher *prev, *next;
	MConnector(MappedMesher* _prev, MappedMesher* _next): prev(_prev), next(_next){}
public:
	static shared_ptr<MConnector>
	Build(CornerTp tp, MappedMesher* prev, MappedMesher* next);
	
	//modifies prev, next meshes,
	//build connection mesh if nessessary.
	void Apply(){ ModifyAdjecents(); BuildInternals();}

	//adds next grid along with connection section to g
	virtual void Add(shared_ptr<BGrid>& g) = 0;
};

struct FillerConnector: public MConnector{
protected:
	shared_ptr<MappedMesher> connection_grid;
	shared_ptr<MappedRect> connection_area;
	FillerConnector(MappedMesher* _prev, MappedMesher* _next): MConnector(_prev, _next){}
public:
	void Add(shared_ptr<BGrid>& g) override{
		//1) add previous grid if it has not been add yet
		if (!g) g.reset(new BGrid());
		if (g->n_cells() == 0) GGeom::Modify::ShallowAdd(&prev->result, g.get());
		//2) add connection grid
		GGeom::Modify::ShallowAdd(&connection_grid->result, g.get());
		//3) add next grid
		GGeom::Modify::ShallowAdd(&next->result, g.get());
	}
};

struct RightConnector: public FillerConnector{
private:
	HMCont2D::Container<HMCont2D::Contour> left, bot;
	HMCont2D::Contour right, top;
	void BuildInternals() override {
		//1) assemble borders
		HMCont2D::Contour _left = prev->rect->BottomContour();
		left = HMCont2D::Contour::CutByWeight(_left, prev->wend, 1.0); 
		left.ReallyReverse();
		HMCont2D::Contour _bot = next->rect->BottomContour();
		bot = HMCont2D::Contour::CutByWeight(_bot, 0.0, next->wstart); 
		HMCont2D::Contour right = next->LeftContour();
		HMCont2D::Contour top = prev->RightContour();

		//2) assemble connector
		connection_area = FourLineRect::Factory(left, right, bot, top);
		connection_grid.reset(new MappedMesher(connection_area.get(), 0.0, 1.0));

		//3) build a grid there
		auto bot_part=[&](double, double)->vector<double>{
			//This function is called once since its ok to calculate here
			//get a partition from top and stretch it
			vector<double> ret {0};
			auto lens = HMCont2D::Contour::ELengths(top);
			std::copy(lens.begin(), lens.end(), std::back_inserter(ret));
			std::partial_sum(ret.begin(), ret.end(), ret.begin());
			double lentop = std::accumulate(lens.begin(), lens.end(), 0.0);
			double lenbot = bot.length();
			for (auto& v: ret) v = v/lentop*lenbot;
			return ret;
		};
		vector<double> f2out {0};
		auto lens = HMCont2D::Contour::ELengths(right);
		std::copy(lens.begin(), lens.end(), std::back_inserter(f2out));
		std::partial_sum(f2out.begin(), f2out.end(), f2out.begin());
		auto vert_part=[&f2out](double)->vector<double>{
			//calculate before since function is called multiple times
			return f2out;
		};
		connection_grid->Fill(bot_part, vert_part);
	}
public:
	RightConnector(MappedMesher* _prev, MappedMesher* _next): FillerConnector(_prev, _next){
		//find cross point between top lines of previous and
		//next mesh_areas
		HMCont2D::Contour prev_top = prev->rect->TopContour();
		HMCont2D::Contour next_top = next->rect->TopContour();
		std::tuple<bool, Point, double, double> cross =
			HMCont2D::Contour::Cross(prev_top, next_top);
		assert(std::get<0>(cross));
		//cut meshed area by found weights
		prev->wend = std::get<2>(cross);
		next->wstart = std::get<3>(cross);
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
		connection_area = FourLineRect::Factory(left, right, bot, top);
		connection_grid.reset(new MappedMesher(connection_area.get(), 0.0, 1.0));

		//4) build a grid there
		auto bot_part=[&](double, double)->vector<double>{
			vector<double> ret {0};
			auto lens = HMCont2D::Contour::ELengths(bot);
			std::copy(lens.begin(), lens.end(), std::back_inserter(ret));
			std::partial_sum(ret.begin(), ret.end(), ret.begin());
			return ret;
		};
		vector<double> f2out {0};
		auto lens = HMCont2D::Contour::ELengths(left);
		std::copy(lens.begin(), lens.end(), std::back_inserter(f2out));
		std::partial_sum(f2out.begin(), f2out.end(), f2out.begin());
		auto vert_part=[&f2out](double)->vector<double>{
			return f2out;
		};
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
	void Add(shared_ptr<BGrid>& g) override{ _THROW_NOT_IMP_; }
};

shared_ptr<MConnector> MConnector::Build(CornerTp tp, MappedMesher* prev, MappedMesher* next){
	switch (tp){
		case CornerTp::SHARP:  return std::make_shared<SharpConnector>(prev, next);
		case CornerTp::OBTUSE:  return std::make_shared<ObtuseConnector>(prev, next);
		case CornerTp::CORNER: return std::make_shared<RightConnector>(prev, next);
		default: assert(false);
	}
}

}

shared_ptr<BGrid> BGrid::MeshFullPath(const ExtPath& epath){
	//1. divide by angles
	vector<ExtPath> pths = ExtPath::DivideByAngle(epath,
			{CornerTp::CORNER, CornerTp::OBTUSE, CornerTp::SHARP});

	//2. build conform mapping for each subpath
	ShpVector<MappedRect> mps;
	for (auto& p: pths){
		double h = p.largest_depth();
		mps.push_back(MappedCavity::Factory(p.leftbc, p.rightbc, p, h));
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
	if (epath.is_closed()){ _THROW_NOT_IMP_; }


	//5. build rectangular meshes
	for (int i=0; i<mesher4.size(); ++i){
		ExtPath& ipth = pths[i];
		auto bpart = [&ipth](double len1, double len2)->vector<double>{
			return ipth.PathPartition(len1, len2);
		};
		auto wpart = [&ipth](double len1)->vector<double>{
			return ipth.VerticalPartition(len1);
		};
		mesher4[i]->Fill(bpart, wpart);
	}

	//6. Connection procedures
	for (auto& c: connectors) c->Apply();
	
	//7. Gather all resulting meshes
	shared_ptr<BGrid> g;
	for (auto& c: connectors) c->Add(g);
	g->merge_congruent_points();

	return g;
}

shared_ptr<BGrid> BGrid::MeshSequence(vector<Options*>& data){ 
	auto ret = shared_ptr<BGrid>();

	//1) assemble extended path = path + boundaries + angles
	ExtPath fullpath = ExtPath::Assemble(data);

	//2) if corner angle section is very short then
	//   make it sharp or regular depending on adjecent corner types
	ExtPath::ReinterpretCornerTp(fullpath);

	//3) build grid for a path
	ret = BGrid::MeshFullPath(fullpath);

	//4) guarantee no self-intersections
	ret = NoSelfIntersections(ret);


	return ret;
}

shared_ptr<BGrid> BGrid::NoSelfIntersections(shared_ptr<BGrid> g){
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
