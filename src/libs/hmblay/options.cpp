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
			{"NO", "KEEP_ORIGIN", "KEEP_SHAPE", "IGNORE_ALL"},
				{BndStepMethod::NO_BND_STEPPING,
				 BndStepMethod::CONST_BND_STEP_KEEP_ALL,
				 BndStepMethod::CONST_BND_STEP_KEEP_SHAPE,
				 BndStepMethod::CONST_BND_STEP});
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
		for(auto& p: r.partition) p/=sc->L;
	}
	
	//4) Put start, end points to contours, make a partitions
	//   Fill Option.path field
	for (auto& r: ret) r.Initialize();
	return ret;
}

void Options::Initialize(){
	//make correct partition of contours, assemble pathes
	typedef HMCont2D::ECollection ECol;
	typedef HMCont2D::Contour ECont;
	typedef HMCont2D::PartitionTp Ptp;
	typedef HMCont2D::ExtendedTree Etree;
	//1. Start and End points
	pnt_start = ECol::FindClosestNode(*edges, start);
	pnt_end = ECol::FindClosestNode(*edges, end);
	//assemble a tree -> find contour with pnt_start/end -> cut it
	auto tree = Etree::Assemble(*edges);
	auto cont = tree.get_contour(pnt_start);
	assert(cont->contains_point(pnt_end));
	path = ECont::Assemble(*cont, pnt_start, pnt_end);

	//reverse path if nessesary.
	if (direction == Direction::OUTER) path.Reverse();

	//partition
	switch (bnd_step_method){
		case BndStepMethod::CONST_BND_STEP_KEEP_ALL:
			ECont::Partition(bnd_step, path, __all_data->pdata, Ptp::KEEP_ALL); 
			break;
		case BndStepMethod::CONST_BND_STEP:
			ECont::Partition(bnd_step, path, __all_data->pdata, Ptp::IGNORE_ALL); 
			break;
		case BndStepMethod::CONST_BND_STEP_KEEP_SHAPE:
			ECont::Partition(bnd_step, path, __all_data->pdata, Ptp::KEEP_SHAPE); 
			break;
		case BndStepMethod::NO_BND_STEPPING:
			break;
	};
}

vector<vector<Options*>> Options::BuildSequence(vector<Options>& inp){
	//Connect each pathes of options into sequences
	auto ret = vector<vector<Options*>>();
	std::set<Options*> setop;
	for (auto& x: inp) setop.insert(&x);
	while (setop.size() > 0){
		ret.push_back(vector<Options*>());
		auto& nd = ret.back();
		nd.push_back(*setop.begin()); setop.erase(setop.begin());
		if (nd[0]->path.is_closed()) continue;
		else {
			_THROW_NOT_IMP_;
		}
	}
	return ret;
}
