#include "bgrid.hpp"

using namespace HMBlay::Impl;

auto HMBlay::Impl::BuildPathPntData(Point* p1, Point* p2, Point* p3, Options* d) -> PathPntData{
	//building of object for 'ext_data' field
	PathPntData ret;
	ret.opt = d;
	if (p1 == 0 || p3 == 0){
		ret.tp = CornerTp::NO;
		ret.angle = M_PI;
		if (p1 == 0) ret.normal = vecRotate((*p3 - *p2), M_PI/2.0);
		if (p3 == 0) ret.normal = vecRotate((*p2 - *p1), M_PI/2.0);
	} else {
		ret.angle = Angle(*p1, *p2, *p3);
		//type
		if (ret.angle<DegToAngle(d->sharp_angle)) ret.tp = CornerTp::SHARP;
		else if (ret.angle<DegToAngle(d->corner_angle)) ret.tp = CornerTp::CORNER;
		else if (ret.angle<DegToAngle(d->regular_angle)) ret.tp = CornerTp::REGULAR;
		else ret.tp = CornerTp::OBTUSE;
		//normal, option
		ret.normal = vecRotate((*p3 - *p2), ret.angle/2.0);
	}
	vecNormalize(ret.normal);
	return ret;
};

ExtPath ExtPath::Assemble(const vector<Options*>& data){
	//data[i]->path are in strict sequence.
	//all data may have different options
	ExtPath ret;
	for (auto d: data){
		auto p = d->get_path();
		//check for ordering
		assert(ret.size() == 0 || ret.last() == p->first());
		auto pnt = p->ordered_points();
		//add temporary values for first and last points
		if (ret.size() == 0) ret.ext_data.push_back({CornerTp::NO, Vect(), d});
		//fill internal
		for (int i=1; i < pnt.size()-1; ++i){
			ret.ext_data.push_back(BuildPathPntData(pnt[i-1], pnt[i], pnt[i+1], d));
		}
		ret.Unite(*p);
	}
	//add temporary values for first and last points
	ret.ext_data.push_back({CornerTp::NO, Vect(), data.back()});
	//check first and last points connection
	if (ret.is_closed()){
		//if closed -> compute last-first connection
		auto last3 = Ed::PointOrder(*ret.data.back(), *ret.edge(0));
		ret.ext_data[0] = BuildPathPntData(last3[0], last3[1], last3[2], ret.ext_data[0].opt);
		ret.ext_data.back() = ret.ext_data[0];
	} else {
		//if open -> place normal at right angle
		auto first3 = Ed::PointOrder(*ret.edge(0), *ret.edge(1));
		ret.ext_data[0] = BuildPathPntData(0, first3[0], first3[1], data[0]);
		auto last3 = Ed::PointOrder(*ret.edge(ret.size()-2), *ret.edge(ret.size()-1));
		ret.ext_data.back() = BuildPathPntData(last3[1], last3[2], 0, data.back());
	}
	return ret;
}

BGrid::TEpMap BGrid::RemoveObtuse(vector<ExtPath*>&& v){
	_DUMMY_FUN_;
	assert(v.size()>0);
	BGrid::TEpMap ret;
	ret[v[0]] = vector<ExtPath>( {ExtPath(*v[0])} );
	return ret;
}

BGrid::TEpMap BGrid::RemoveSharp(vector<ExtPath*>&& v){
	_DUMMY_FUN_;
	assert(v.size()>0);
	BGrid::TEpMap ret;
	ret[v[0]] = vector<ExtPath>( {ExtPath(*v[0])} );
	return ret;
}

BGrid::TEpMap BGrid::RemoveCorner(vector<ExtPath*>&& v){
	_DUMMY_FUN_;
	assert(v.size()>0);
	BGrid::TEpMap ret;
	ret[v[0]] = vector<ExtPath>( {ExtPath(*v[0])} );
	return ret;
}

BGrid BGrid::JoinCornerNodes(vector<std::pair<ExtPath*, BGrid*>>& inp){
	_DUMMY_FUN_;
	assert(inp.size()>0);
	return BGrid(*inp[0].second);
}

BGrid BGrid::JoinSharpNodes(vector<std::pair<ExtPath*, BGrid*>>& inp){
	_DUMMY_FUN_;
	assert(inp.size()>0);
	return BGrid(*inp[0].second);
}

BGrid BGrid::JoinObtuseNodes(vector<std::pair<ExtPath*, BGrid*>>& inp){
	assert(inp.size()>0);
	_DUMMY_FUN_;
	return BGrid(*inp[0].second);
}

SimpleBGrid::SimpleBGrid(ExtPath& pth):BGrid(){
	//Only NO and REGULAR corners in pth
	assert(pth.size() == pth.ext_data.size() - 1);
	assert(std::all_of(pth.ext_data.begin(), pth.ext_data.end(),
		[](PathPntData& d){
			return (d.tp == CornerTp::NO ||
				d.tp == CornerTp::REGULAR);
		}
	));
	is_closed = pth.is_closed();
	//1) recalculate normals to guarantee smooth normal changes from first to last
	_DUMMY_FUN_;
	//2) assemble grid
	//2.1) fill stencil
	auto basepnts = pth.ordered_points();
	for (int i=0; i<pth.size() + (is_closed?0:1); ++i){
		Point* b = basepnts[i];
		Vect n = pth.ext_data[i].normal;
		auto& part = pth.ext_data[i].opt->partition;
		stencil.push_back(vector<GridPoint*>());
		for (int j = 0; j<part.size(); ++j){
			GridPoint* p = aa::add_shared(points, GridPoint(*b + n * part[j]));
			stencil.back().push_back(p);
		}
	}
	//2.2) fill cells
	for (int edind = 0; edind < pth.size(); ++edind){
		int stenind1 = edind;
		int stenind2 = (is_closed && edind == pth.size() - 1) ? 0 : edind+1;
		//number of cells for this edge
		int n_height = std::min(stencil[stenind1].size(), stencil[stenind2].size()) - 1;
		for (int ih = 0; ih < n_height; ++ih){
			Cell* c = aa::add_shared(cells, Cell());
			add_point_to_cell(c, stencil[stenind1][ih]);
			add_point_to_cell(c, stencil[stenind2][ih]);
			add_point_to_cell(c, stencil[stenind2][ih+1]);
			add_point_to_cell(c, stencil[stenind1][ih+1]);
			priority[c] = HIGHEST_PRIORITY - ih;
		}
	}
	set_indicies();
}

ShpVector<BGrid> BGrid::MeshSequence(vector<Options*>& data){ 
	auto ret = ShpVector<BGrid>();
	auto BuildVP = [](TEpMap& inp)->vector<ExtPath*>{
		//Extracts all keys of input map to vector of pointers to keys
		vector<ExtPath*> ret;
		for (auto& kv: inp) for (auto& v: kv.second) ret.push_back(&v);
		return ret;
	};
	auto assemble_path_grid = [](vector<ExtPath>& dt, std::map<ExtPath*, BGrid*>& mp)
			->vector<std::pair<ExtPath*, BGrid*>>{
		//merge pathes from vector with grids from map values
		vector<std::pair<ExtPath*, BGrid*>> inpdata;
		for (auto& v: dt){
			ExtPath* pth = &v;
			BGrid* g = mp[pth];
			inpdata.push_back(std::make_pair(pth, g));
		};
		return inpdata;
	};

	// ============================ Dividing paths
	// Building a tree-like structure of pathes.
	// On each level path is divided by splitting with angles
	// of special type. Sum of pathes on each level equals full path.
	// TEpMap is std::map<ExtPath*, vector<ExtPath>> -- parent to children dictionary.
	//1) assemble extended path
	ExtPath fullpath = ExtPath::Assemble(data);

	//2) Remove all obtuse entries
	TEpMap no_obtuse = RemoveObtuse({&fullpath});
	
	//3) Remove all sharp entries 
	TEpMap no_sharp = RemoveSharp(BuildVP(no_obtuse));

	//4) Remove all corner entries
	TEpMap no_corner = RemoveCorner(BuildVP(no_sharp));

	// ============================ Assembling primitives
	ShpVector<BGrid> primitives;  //grids pool
	std::map<ExtPath*, BGrid*> no_corner_primitives;
	for (auto& kv: no_corner){
		for (auto& v: kv.second){
			primitives.push_back(std::make_shared<BGrid>(
				SimpleBGrid(v)
			));
			no_corner_primitives[&v] = primitives.back().get();
		}
	}

	// ============================ Joining primitives
	// 3) merging no corner
	std::map<ExtPath*, BGrid*> no_sharp_primitives;
	for (auto& kv: no_corner){
		auto inpdata = assemble_path_grid(kv.second, no_corner_primitives);
		primitives.push_back(std::make_shared<BGrid>(
			JoinCornerNodes(inpdata)
		));
		no_sharp_primitives[kv.first] = primitives.back().get();
	}

	// 2) merging no sharp
	std::map<ExtPath*, BGrid*> no_obtuse_primitives;
	for (auto& kv: no_sharp){
		auto inpdata = assemble_path_grid(kv.second, no_sharp_primitives);
		primitives.push_back(std::make_shared<BGrid>(
			JoinSharpNodes(inpdata) 
		));
		no_obtuse_primitives[kv.first] = primitives.back().get();
	}
	
	// 1) merging no obtuse
	std::map<ExtPath*, BGrid*> fingrd;
	for (auto& kv: no_obtuse){
		auto inpdata = assemble_path_grid(kv.second, no_obtuse_primitives);
		primitives.push_back(std::make_shared<BGrid>(
			JoinObtuseNodes(inpdata)
		));
		fingrd[kv.first] = primitives.back().get();
	}
	
	// ============================ Assemble result
	//find a fingrd in a grids pool and return it
	auto full_path_pointer = &fullpath;
	for (auto& g: primitives){
		if (fingrd.find(full_path_pointer) != fingrd.end()){
				ret.push_back(g);
				return ret;
		}
	}
	throw;
}

BGrid BGrid::ImposeBGrids(ShpVector<BGrid>& gg){
	std::cout<<"DUMMY Impose BGrids"<<std::endl;
	return BGrid(*gg[0]);
}

