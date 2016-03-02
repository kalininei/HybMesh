#include "options.hpp"

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
	typedef HMCont2D::ECollection Edc;
	typedef HMCont2D::Container<Edc> EdCon;
	typedef HMCont2D::Edge Ed;
	vector<Options> ret;
	//1) copy unique edge collections
	shared_ptr<EdCon> all_data(new EdCon());
	std::map<Edc*, shared_ptr<Edc>> edges_old_new;
	for (auto& p: par){
		auto fnd = edges_old_new.find(p.edges);
		if (fnd == edges_old_new.end()){
			Edc& eold = *p.edges;
			Edc& enew = *(edges_old_new[p.edges] = std::make_shared<Edc>(Edc()));
			//Deepcopy data to temp. container and then share points and
			//edge info to add_data container and enew collection
			EdCon tmp; EdCon::DeepCopy(eold, tmp);
			all_data->Unite(tmp);
			enew.Unite(tmp);
		}
	}
	//2) non-dimensioning copied edges
	shared_ptr<ScaleBase> sc = std::make_shared<ScaleBase>(HMCont2D::Scale01(*all_data));
	
	//3) construct deepcopied input
	for (auto& p: par){
		ret.push_back(Options(p));
		auto& r = ret.back();
		r.__all_data = all_data;
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
	for (auto& r: ret) r.Initialize();
	return ret;
}

void Options::Initialize(){
	//make correct partition of contours, assemble paths
	typedef HMCont2D::ECollection ECol;
	typedef HMCont2D::Contour ECont;
	typedef HMCont2D::PartitionTp Ptp;
	typedef HMCont2D::ExtendedTree Etree;
	//Start and End points
	pnt_start = ECol::FindClosestNode(*edges, start);
	pnt_end = ECol::FindClosestNode(*edges, end);

	//if points were not found as edges vertices then place them straight to __edges_data
	//pool so that other option instances could find them and don't create their own copies
	//of these points (only for open contours)
	if (start != end){
		auto split_edge = [&](Point psplit)->Point*{
			auto ce = HMCont2D::ECollection::FindClosestEdge(*__edges_data, psplit);
			assert(std::get<0>(ce)!=0);
			HMCont2D::Edge* eold = std::get<0>(ce);
			if (ISZERO(std::get<2>(ce))) return  eold->pstart;
			if (ISZERO(1 - std::get<2>(ce))) return  eold->pend;
			Point* newp = __all_data->pdata.add_value(Point::Weigh(*eold->pstart,
						*eold->pend, std::get<2>(ce)));
			HMCont2D::Edge enew1(eold->pstart, newp), enew2(eold->pend, newp);
			__edges_data->RemoveAt({std::get<3>(ce)});
			__edges_data->add_value(enew1);
			__edges_data->add_value(enew2);
			return newp;
		};
		if (*pnt_start != start) pnt_start = split_edge(start);
		if (*pnt_end != end) pnt_end = split_edge(end);
	}

	//assemble a tree -> find contour with pnt_start/end -> cut it
	auto tree = HMCont2D::Assembler::ETree(*edges);
	HMCont2D::Contour* contour_in_tree = tree.get_contour(pnt_start);
	assert(contour_in_tree != 0);
	full_source = HMCont2D::Contour(*contour_in_tree);
	assert(full_source.contains_point(pnt_end));
	//place Start/End points to contour if necessary
	if (!(start==end && full_source.is_closed())){
		//pnt_start = std::get<1>(full_source.GuaranteePoint(start, __all_data->pdata));
		//pnt_end   = std::get<1>(full_source.GuaranteePoint(end, __all_data->pdata));
		//throw if points coincide
		if (pnt_start == pnt_end)
			throw std::runtime_error("Zero length boundary layer section");
	}

	//assemble path
	path = HMCont2D::Assembler::Contour1(full_source, pnt_start, pnt_end);
	//all edges in path are directed according to global direction. May be usefull.
	path.DirectEdges();
	//if path has only 1 segment we need to implicitly define its direction
	if (path.size() == 1){ path.edge(0)->pstart = pnt_start; path.edge(0)->pend = pnt_end; }
	//reverse full_source if path is in opposite direction.
	if (full_source.next_point(pnt_start) != path.next_point(pnt_start)) full_source.Reverse();
	assert(full_source.next_point(pnt_start) == path.next_point(pnt_start));
	
	//reverse path if necessary.
	if (direction == Direction::OUTER) { path.ReallyReverse(); full_source.Reverse(); }
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
			Point* p0 = lst.front()->path.first();
			Point* p1 = lst.back()->path.last();
			//path is already closed
			if (p0 == p1) break;
			//last=first
			if (p1 == (*it)->path.first()){
				lst.push_back(*it); it = setop.erase(it);
			//first=last
			} else if (p0 == (*it)->path.last()){
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
