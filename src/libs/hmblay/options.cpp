#include "options.hpp"
#include "algos.hpp"
#include "cont_partition.hpp"
#include "cont_assembler.hpp"
#include "treverter2d.hpp"

using namespace HMBlay::Impl;

template <class C>
C string_option(std::string opt,
		std::initializer_list<std::string> opt_val,
		std::initializer_list<C> ret_val){
	std::transform(opt.begin(), opt.end(), opt.begin(), ::toupper);
	auto ival = opt_val.begin();
	auto rval = ret_val.begin();
	while (ival != opt_val.end() || rval != ret_val.end()){
		auto s = *ival;
		std::transform(s.begin(), s.end(), s.begin(), ::toupper);
		if (*ival == opt) return *rval;
		++ival;
		++rval;
	}
	throw std::runtime_error("impossible to parse option " + opt);
}

HMBlay::BndStepMethod HMBlay::MethFromString(const char* str){
	return string_option(str,
			{"NO", "KEEP_ORIGIN", "KEEP_SHAPE", "IGNORE_ALL", "KEEP_ALL", "INCREMENTAL"},
				{BndStepMethod::NO_BND_STEPPING,
				 BndStepMethod::CONST_BND_STEP_KEEP_ALL,
				 BndStepMethod::CONST_BND_STEP_KEEP_SHAPE,
				 BndStepMethod::CONST_BND_STEP,
				 BndStepMethod::CONST_BND_STEP_KEEP_ALL,
				 BndStepMethod::INCREMENTAL});
}

HMBlay::Direction HMBlay::DirectionFromString(const char* str){
	return string_option(str,
		{"INSIDE", "OUTSIDE", "LEFT", "RIGHT", "AROUND", "INNER", "OUTER", "BOTH"},
		{Direction::INNER, Direction::OUTER,
		 Direction::INNER,  Direction::OUTER, Direction::BOTH,
		 Direction::INNER, Direction::OUTER, Direction::BOTH});
}


vector<Options> Options::CreateFromParent(const vector<HMBlay::Input>& par){
	//makes deep copy of parent contours, scales them,
	//writes everything to ret
	vector<Options> ret;
	//1) copy unique edge collections
	std::map<HM2D::EdgeData*, shared_ptr<HM2D::EdgeData>> edges_old_new;
	for (auto& p: par){
		auto er = edges_old_new.emplace(p.edges, std::make_shared<HM2D::EdgeData>());
		if (er.second == true) HM2D::DeepCopy(*p.edges, *er.first->second);
	}
	HM2D::VertexData _av;
	for (auto& m: edges_old_new){
		auto __av = HM2D::AllVertices(*m.second);
		_av.insert(_av.end(), __av.begin(), __av.end());
	}
	aa::no_duplicates(_av);

	//2) non-dimensioning copied edges
	shared_ptr<ScaleBase> sc = std::make_shared<ScaleBase>(HM2D::Scale01(_av));
	
	//3) construct deepcopied input
	for (auto& p: par){
		ret.push_back(Options(p));
		auto& r = ret.back();
		r.__edges_data = edges_old_new[p.edges];
		r.edges = r.__edges_data.get();
		r.scaling = sc;
		r.bnd_step /= sc->L;
		for (auto& b: r.bnd_step_basis){
			sc->scale(b.first);
			b.second /=sc->L;
		}
		sc->scale(r.start);
		sc->scale(r.end);
		for(auto& p: r.partition) p/=sc->L;
	}
	
	//4) Put start, end points to contours
	//   Fill Option.path field
	for (auto& r: ret){
		r.Initialize();
	}
	return ret;
}

void Options::Initialize(){
	//make correct partition of contours, assemble paths
	typedef HM2D::EdgeData ECol;
	typedef HM2D::EdgeData ECont;
	typedef HM2D::Contour::Algos::PartitionTp Ptp;
	typedef HM2D::Contour::Tree Etree;
	//Start and End points
	if (start != end){
		pnt_start = std::get<1>(HM2D::Contour::GuaranteePoint(*edges, start)).get();
		pnt_end = std::get<1>(HM2D::Contour::GuaranteePoint(*edges, end)).get();
	} else {
		auto av = HM2D::AllVertices(*edges);
		auto fnd1 = HM2D::FindClosestNode(av, start);
		auto fnd2 = HM2D::FindClosestNode(av, end);
		pnt_start = av[std::get<0>(fnd1)].get();
		pnt_end = av[std::get<0>(fnd2)].get();
	}
	/*
	auto av = HM2D::AllVertices(*edges);
	auto fnd1 = HM2D::FindClosestNode(av, start);
	auto fnd2 = HM2D::FindClosestNode(av, end);
	pnt_start = av[std::get<0>(fnd1)].get();
	pnt_end = av[std::get<0>(fnd2)].get();

	//if points were not found as edges vertices then place them straight to __edges_data
	//pool so that other option instances could find them and don't create their own copies
	//of these points (only for open contours)
	if (start != end){
		auto split_edge = [&](Point psplit)->Point*{
			auto ce = HM2D::FindClosestEdge(*__edges_data, psplit);
			assert(std::get<0>(ce)!=0);
			HM2D::Edge* eold = (*__edges_data)[std::get<0>(ce)].get();
			if (ISZERO(std::get<2>(ce))) return  eold->first().get();
			if (ISZERO(1 - std::get<2>(ce))) return  eold->last().get();
			auto newp = std::make_shared<HM2D::Vertex>(
					Point::Weigh(*eold->first(), *eold->last(), std::get<2>(ce)));
			auto enew1 = std::make_shared<HM2D::Edge>(eold->first(), newp);
			auto enew2 = std::make_shared<HM2D::Edge>(eold->last(), newp);
			__edges_data->erase(__edges_data->begin() + std::get<0>(ce));
			__edges_data->push_back(enew1);
			__edges_data->push_back(enew2);
			return newp.get();
		};
		if (*pnt_start != start) pnt_start = split_edge(start);
		if (*pnt_end != end) pnt_end = split_edge(end);
	}
	*/

	//assemble a tree -> find contour with pnt_start/end -> cut it
	auto tree = HM2D::Contour::Tree::Assemble(*edges);
	//HM2D::Contour::R::RevertTree::Permanent(tree);
	//HM2D::EdgeData* contour_in_tree = &tree.find_node(pnt_start)->contour;
	//full_source = HM2D::EdgeData(*contour_in_tree);
	auto contour_in_tree = tree.find_node(pnt_start);
	HM2D::DeepCopy(contour_in_tree->contour, full_source, 0);
	if (HM2D::Contour::IsClosed(full_source)){
		double fsarea = HM2D::Contour::Area(full_source);
		if ((contour_in_tree->level % 2 == 0 && fsarea < 0) ||
		    (contour_in_tree->level % 2 == 1 && fsarea > 0)){
			HM2D::Contour::R::ReallyRevert::Permanent(full_source);
		} else {
			HM2D::Contour::R::ReallyDirect::Permanent(full_source);
		}
	}
	assert(HM2D::Contains(full_source, pnt_end));
	//place Start/End points to contour if necessary
	if (!(start==end && HM2D::Contour::IsClosed(full_source))){
		//throw if points coincide
		if (pnt_start == pnt_end)
			throw std::runtime_error("Zero length boundary layer section");
	}

	//assemble path
	path = HM2D::Contour::Assembler::ShrinkContour(full_source, pnt_start, pnt_end);
	//all edges in path are directed according to global direction. May be usefull.
	HM2D::Contour::R::ForceFirst::Permanent(path, *pnt_start);
	//reverse full_source if path is in opposite direction.
	if (HM2D::Contour::PInfo(full_source, pnt_start).pnext !=
	    HM2D::Contour::PInfo(path, pnt_start).pnext){
		HM2D::Contour::Reverse(full_source);
	}
	assert(HM2D::Contour::PInfo(full_source, pnt_start).pnext ==
	       HM2D::Contour::PInfo(path, pnt_start).pnext);
	
	//reverse path if necessary.
	if (direction == Direction::OUTER){
		HM2D::Contour::R::ReallyRevert::Permanent(path);
		HM2D::Contour::Reverse(full_source);
	}
}

vector<vector<Options*>> Options::BuildSequence(vector<Options>& inp){
	//Connect each paths of options into sequences
	auto ret = vector<vector<Options*>>();
	std::set<Options*> setop;
	for (auto& x: inp) setop.insert(&x);
	while (setop.size() > 0){
		std::list<Options*> lst;
		//1) place first
		lst.push_back(*setop.begin()); setop.erase(setop.begin());
		//2) in a loop find if edge->last==lst.first => push_front
		//                  if edge->first==lst.last => push_back
		for (auto it=setop.begin(); it!=setop.end();){
			Point* p0 = HM2D::Contour::First(lst.front()->path).get();
			Point* p1 = HM2D::Contour::Last(lst.back()->path).get();
			//path is already closed
			if (p0 == p1) break;
			//last=first
			if (p1 == HM2D::Contour::First((*it)->path).get()){
				lst.push_back(*it); it = setop.erase(it);
			//first=last
			} else if (p0 == HM2D::Contour::Last((*it)->path).get()){
				lst.push_front(*it); it = setop.erase(it);
			} else ++it;
		}
		//3) add lst to ret and clear it
		ret.push_back(vector<Options*>(lst.begin(), lst.end()));
	}
	return ret;
}

CornerTp Options::CornerType(double a){
	double a1 = DegToAngle(acute_angle);
	double a2 = DegToAngle(right_angle);
	double a3 = DegToAngle(straight_angle);
	double a4 = DegToAngle(reentrant_angle);
	if (ISZERO(a3)) a3 = 2*M_PI+1;
	if (ISZERO(a4)) a4 = 2*M_PI+1;

	CornerTp tp;
	if (fabs(a)<1e-6) tp = CornerTp::ZERO;
	else if (a<a1) tp = CornerTp::ACUTE;
	else if (a<a2) tp = CornerTp::RIGHT;
	else if (a<a3) tp = CornerTp::STRAIGHT;
	else if (a<a4) tp = CornerTp::REENTRANT;
	else tp = CornerTp::ROUND;
	return tp;
};
