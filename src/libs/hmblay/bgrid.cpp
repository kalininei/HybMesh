#include "bgrid.hpp"
#include "simple_bgrid.hpp"

using namespace HMBlay::Impl;

void PathPntData::fill(Point* p1, Point* p2, Point* p3){
	//building of object for 'ext_data' field
	if (p1 == 0 || p3 == 0){
		tp = CornerTp::NO;
		angle = M_PI;
		if (p1 == 0) normal = vecRotate((*p3 - *p2), M_PI/2.0);
		if (p3 == 0) normal = vecRotate((*p2 - *p1), M_PI/2.0);
	} else {
		angle = Angle(*p1, *p2, *p3);
		//type
		if (angle<DegToAngle(opt->sharp_angle)) tp = CornerTp::SHARP;
		else if (angle<DegToAngle(opt->corner_angle)) tp = CornerTp::CORNER;
		else if (angle<DegToAngle(opt->regular_angle)) tp = CornerTp::REGULAR;
		else tp = CornerTp::OBTUSE;
		//normal, option
		normal = vecRotate((*p3 - *p2), angle/2.0);
	}
	vecNormalize(normal);
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
		//add empty extended data object.
		//conflicting i=0 edge: set priority to longest partition
		if (ret.size()>0){
			auto L1 = ret.ext_data.back().opt->partition.back();
			auto L2 = d->partition.back();
			if (L2>L1) ret.ext_data.back().opt = d;
		} else ret.ext_data.push_back(PathPntData(d));
		for (int i=1; i<pnt.size(); ++i){
			ret.ext_data.push_back(PathPntData(d));
		}
		//add edges
		ret.Unite(*p);
	}
	//fill extended data
	auto p = ret.ordered_points();
	if (ret.is_closed()){
		//conflicting first-last
		auto L1 = ret.ext_data[0].opt->partition.back(),
		     L2 = ret.ext_data.back().opt->partition.back();
		if (L1>=L2) ret.ext_data.back().opt = ret.ext_data[0].opt;
		else ret.ext_data[0].opt = ret.ext_data.back().opt;
		//take into account last-first connection
		for (int i=0; i<p.size(); ++i){
			Point* pcur = p[i];
			Point* pprev = (i==0) ? p[p.size()-2] : p[i-1];
			Point* pnext = (i==p.size()-1) ? p[1] : p[i+1];
			ret.ext_data[i].fill(pprev, pcur, pnext);
		}
	} else {
		//assemble without additional connections
		for (int i=0; i<p.size(); ++i){
			Point* pcur = p[i];
			Point* pprev = (i==0) ? 0 : p[i-1];
			Point* pnext = (i==p.size()-1) ? 0 :p[i+1];
			ret.ext_data[i].fill(pprev, pcur, pnext);
		}
	}
	return ret;
}

void ExtPath::AdoptEndNormals(const Contour& outer){
	if (is_closed()) return;
	auto need_to_adopt = [](Vect n1, Vect n2){
		if (triarea(Point(0, 0), n1, n2) < 0) return true;
		if (Angle(n2, Point(0, 0), n1) < M_PI/6) return true;
		return false;
	};

	//start point
	Point* p0 = first();
	auto sib0 = outer.point_siblings(p0);
	assert(sib0[1] && sib0[2]);
	if (sib0[0] != 0){
		Vect n = *sib0[0] - *sib0[1];
		if (need_to_adopt(ext_data[0].normal, n)){
			vecNormalize(n);
			ext_data[0].normal = n;
		}
	}
	//end point
	Point* p1 = last();
	auto sib1 = outer.point_siblings(p1);
	assert(sib1[1] && sib1[0]);
	if (sib1[2] != 0){
		Vect n = *sib1[2] - *sib1[1];
		if (need_to_adopt(n, ext_data.back().normal)){
			vecNormalize(n);
			ext_data.back().normal = n;
		}
	}
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
	// basic options
	int Nsmooth = 0;
	for (auto op: data) if (op->smooth_normals_steps>Nsmooth) 
		Nsmooth = op->smooth_normals_steps;

	// ============================ Dividing paths
	// Building a tree-like structure of pathes.
	// On each level path is divided by splitting with angles
	// of special type. Sum of pathes on each level equals full path.
	// TEpMap is std::map<ExtPath*, vector<ExtPath>> -- parent to children dictionary.
	//1) assemble extended path
	ExtPath fullpath = ExtPath::Assemble(data);

	//2) Adopt first/last normal according to full contour if full contour
	//   is closed and meshed contour is open
	fullpath.AdoptEndNormals(*data[0]->get_full_source());

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
				SimpleBGrid(v, Nsmooth)
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
	_DUMMY_FUN_;
	return BGrid(*gg[0]);
}

