#include <unordered_map>
#include "hmproject.h"
#include "cport_cont2d.h"
#include "hybmesh_contours2d.hpp"

namespace{
HMCont2D::ECollection* to_ecol(void* g){
	return static_cast<HMCont2D::ECollection*>(g);
}
const HMCont2D::ECollection* to_ecol(const void* g){
	return static_cast<const HMCont2D::ECollection*>(g);
}
}

namespace{

HMCont2D::Container<HMCont2D::ECollection>
const_partition(const HMCont2D::ExtendedTree& etree, double step, const vector<Point*>& keep){
	HMCont2D::PCollection _pdata;
	vector<HMCont2D::Contour> partdata;
	for (int i=0; i<etree.cont_count(); ++i){
		partdata.push_back(
			HMCont2D::Algos::Partition(step, *etree.get_contour(i),
				_pdata, keep));
	}
	//assemble data
	HMCont2D::Container<HMCont2D::ECollection> ret;
	for (int i=0; i<partdata.size(); ++i){
		HMCont2D::Container<HMCont2D::ECollection>::DeepCopy(partdata[i], ret);
	}
	return ret;
}

HMCont2D::Container<HMCont2D::ECollection>
refp_partition(const HMCont2D::Contour& cont,
		const vector<double>& steps, const vector<Point>& bpoints,
		const vector<Point*>& keep){
	if (bpoints.size() == 1){
		HMCont2D::ExtendedTree etree;
		etree.AddContour(cont);
		return const_partition(etree, steps[0], keep);
	}
	//assemble weights
	std::map<double, double> weights;
	for (int i=0; i<steps.size(); ++i){
		auto info = cont.coord_at(bpoints[i]);
		weights[std::get<1>(info)] = steps[i];
	}
	//call partition
	HMCont2D::PCollection _pdata;
	HMCont2D::Contour cret = HMCont2D::Algos::WeightedPartition(weights, cont, _pdata, keep); 
	//assemble return data
	HMCont2D::Container<HMCont2D::ECollection> ret;
	HMCont2D::Container<HMCont2D::ECollection>::DeepCopy(cret, ret);
	return ret;
}

}

// cont: ECOllection pointer
// btypes: size = cont.size(). Boundary feature for each contour edge
// algo: 0 - const step; 1 - refference point step
// n_steps: number of reference points in case of algo=1
// steps: if algo = 0 => [const_step], if algo = 1 => [step0, x0, y0, step1, x1, y1, ...]
// a0: insignificant angle [180-a0, 180 + a0]
// keepbnd: =true if all boundary type changing nodes should be preserved
// n_outbnd - number of output contour edges
// outbnd - boundary feature for each output contour edge
// returns new ECollection pointer or NULL if failed
void* contour_partition(void* cont, int* btypes, int algo,
		int n_steps, double* steps, double a0, int keepbnd,
		int* n_outbnd, int** outbnd){
	typedef HMCont2D::ECollection TCol;
	typedef HMCont2D::Container<TCol> TCont;
	//gather data
	TCont* ret = NULL;
	TCol* ecol = to_ecol(cont);
	vector<Point> basic_points;
	vector<double> basic_steps;
	if (algo == 1){
		for (int i=0; i<n_steps; ++i){
			basic_steps.push_back(steps[3*i]);
			basic_points.emplace_back(steps[3*i+1], steps[3*i+2]);
		}
	} else {
		basic_steps.push_back(steps[0]);
	}
	//scaling
	ScaleBase sc = HMCont2D::ECollection::Scale01(*ecol);
	sc.scale(basic_points.begin(), basic_points.end());
	for (auto& v: basic_steps) v/=sc.L;
	//main procedure
	try{
		//assemble tree from input data
		HMCont2D::ExtendedTree ext = HMCont2D::Assembler::ETree(*ecol);
		//assemble boundary types for edges of input contour
		std::unordered_map<HMCont2D::Edge*, int> btypes_map;
		if (keepbnd){
			int i=0;
			for (auto e: *ecol) btypes_map[e.get()] = btypes[i++]; 
		}
		//assemble significant points
		std::vector<Point*> keep_points;
		for (int i=0; i<ext.cont_count(); ++i){
			HMCont2D::Contour& cont = *ext.get_contour(i);
			for (auto& info: cont.ordered_info()){
				//end points of open contour
				if (info.pprev == NULL || info.pnext == NULL){
					keep_points.push_back(info.p);
					continue;
				}
				//angle between edges
				double angle = (Angle(*info.pprev, *info.p, *info.pnext)-M_PI)/M_PI*180.0;
				if (ISGREATER(fabs(angle), a0)){
					keep_points.push_back(info.p);
					continue;
				}
				//boundary change
				if (keepbnd){
					auto bprev = btypes_map[info.eprev.get()];
					auto bnext = btypes_map[info.enext.get()];
					if (bprev != bnext){
						keep_points.push_back(info.p);
						continue;
					}
				}
			}
		}
		//call partition algorithm
		TCont out;
		if (algo == 0){
			out = const_partition(ext, basic_steps[0], keep_points);
		} else if (algo == 1){
			if (ext.cont_count() != 1){
				throw std::runtime_error("only singly connected contours "
						"are allowed for refference point partition");
			}
			out = refp_partition(*ext.get_contour(0), basic_steps, basic_points, keep_points);
		}
		//boundary assignment
		*n_outbnd = out.size();
		*outbnd = new int[*n_outbnd];
		set_ecollection_bc_force(cont, &out, btypes, *outbnd, 3);
		//unscale and return value
		HMCont2D::ECollection::Unscale(out, sc);
		ret = new TCont(std::move(out));
	} catch (std::runtime_error &e){
		ret = NULL;
		std::cout<<e.what()<<std::endl;
	}
	HMCont2D::ECollection::Unscale(*ecol, sc);
	return ret;
}
