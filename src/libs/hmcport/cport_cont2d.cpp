#include <unordered_map>
#include "hmproject.h"
#include "cport_cont2d.h"
#include "primitives2d.hpp"
#include "contour.hpp"
#include "tree.hpp"
#include "cont_partition.hpp"
#include "partition01.hpp"
#include "cont_assembler.hpp"
#include "algos.hpp"
#include "treverter2d.hpp"
#include "finder2d.hpp"
#include "constructor.hpp"
#include "export2d_hm.hpp"
#include "import2d_hm.hpp"

namespace{
HM2D::EdgeData* to_ecol(void* g){
	return static_cast<HM2D::EdgeData*>(g);
}
const HM2D::EdgeData* to_ecol(const void* g){
	return static_cast<const HM2D::EdgeData*>(g);
}
}

namespace{

HM2D::EdgeData
const_partition(const HM2D::Contour::Tree& etree, double step, const HM2D::VertexData& keep, int ne){
	vector<HM2D::EdgeData> partdata;
	int cc = etree.nodes.size();
	if (ne <= 0) for (int i=0; i<cc; ++i){
		partdata.push_back(
			HM2D::Contour::Algos::Partition(step, etree.nodes[i]->contour, keep));
	} else {
		vector<double> lens;
		for (int i=0; i<cc; ++i) lens.push_back(HM2D::Length(etree.nodes[i]->contour));
		//minimum size includes all keep points
		vector<int> minsize(cc, 0);
		for (int i=0; i<keep.size(); ++i){
			auto c = etree.find_node(keep[i].get());
			if (c == nullptr) continue;
			else for (int j=0; j<cc; ++j){
				if (etree.nodes[j] == c) {
					if (HM2D::Contour::IsClosed(c->contour) ||
						(keep[i] != HM2D::Contour::First(c->contour) &&
						 keep[i] != HM2D::Contour::Last(c->contour)))
							++minsize[j];
					break;
				}
			}
		}
		for (int i=0; i<cc; ++i){
			if (!HM2D::Contour::IsClosed(etree.nodes[i]->contour)) minsize[i]++;
			if (HM2D::Contour::IsClosed(etree.nodes[i]->contour) && minsize[i]<3) minsize[i] = 3;
		}
		//nums double
		vector<double> nums_double;
		double full_len = std::accumulate(lens.begin(), lens.end(), 0.0);
		for (int i=0; i<lens.size(); ++i){
			nums_double.push_back(lens[i]/full_len*ne);
		}
		//nums int
		vector<int> nums = HMMath::RoundVector(nums_double, minsize);
		//building
		std::map<double, double> m; m[0]=step;
		for (int i=0; i<etree.nodes.size(); ++i){
			partdata.push_back(
				HM2D::Contour::Algos::WeightedPartition(m, etree.nodes[i]->contour,
					nums[i], keep));
		}
	}
	//assemble data
	HM2D::EdgeData ret;
	for (int i=0; i<partdata.size(); ++i){
		HM2D::DeepCopy(partdata[i], ret);
	}
	return ret;
}

HM2D::EdgeData
refp_partition(const HM2D::EdgeData& cont,
		const vector<double>& steps, const vector<Point>& bpoints,
		const HM2D::VertexData& keep, int ne){
	if (bpoints.size() == 1){
		HM2D::Contour::Tree etree;
		etree.add_contour(cont);
		return const_partition(etree, steps[0], keep, ne);
	}
	//assemble weights
	std::map<double, double> weights;
	for (int i=0; i<steps.size(); ++i){
		auto info = HM2D::Contour::CoordAt(cont, bpoints[i]);
		weights[std::get<1>(info)] = steps[i];
	}
	//call partition
	HM2D::EdgeData cret = (ne<=0) ? HM2D::Contour::Algos::WeightedPartition(weights, cont, keep)
	                              : HM2D::Contour::Algos::WeightedPartition(weights, cont, ne, keep);
	//assemble return data
	HM2D::EdgeData ret;
	HM2D::DeepCopy(cret, ret);
	return ret;
}

}

// cont: ECOllection pointer
// btypes: size = cont.size(). Boundary feature for each contour edge
// algo: 0 - const step; 1 - refference point step; 2,3 - ref weights/lengths steps
// n_steps: number of reference points in case of algo=1
// steps: if algo = 0 => [const_step], if algo = 1 => [step0, x0, y0, step1, x1, y1, ...]
//        if algo = 2,3=> [step0, s0, step1, s1, ...]
// a0: insignificant angle [180-a0, 180 + a0]
// keepbnd: =true if all boundary type changing nodes should be preserved
// n_outbnd - number of output contour edges
// outbnd - boundary feature for each output contour edge
// returns new ECollection pointer or NULL if failed
void* contour_partition(void* cont, int* btypes, int algo,
		int n_steps, double* steps, double a0, int keepbnd, int nedges,
		int n_crosses, void** crosses,
		int n_keep_pts, double* keep_pts,
		double* start_point, double* end_point,
		int* n_outbnd, int** outbnd){
	typedef HM2D::EdgeData TCol;
	typedef HM2D::EdgeData TCont;
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
	} else if (algo == 0){
		basic_steps.push_back(steps[0]);
	} else if (algo == 2 || algo == 3){
		for (int i=0; i<n_steps; ++i){
			basic_steps.push_back(steps[2*i]);
		}
	}
	vector<TCol*> vcrosses(n_crosses);
	for (int i=0; i<n_crosses; ++i) vcrosses[i] = static_cast<TCol*>(crosses[i]);

	Point start(0, 0), end(0, 0);
	if (start_point){ start.set(start_point[0], start_point[1]); }
	if (end_point){ end.set(end_point[0], end_point[1]); }

	vector<Point> kp(n_keep_pts);
	for (int i=0; i<n_keep_pts; ++i) kp[i] = Point(keep_pts[2*i], keep_pts[2*i+1]);
	//scaling
	ScaleBase sc = HM2D::Scale01(*ecol);
	sc.scale(basic_points.begin(), basic_points.end());
	for (auto& v: basic_steps) v/=sc.L;
	for (auto& c: vcrosses) HM2D::Scale(*c, sc);
	sc.scale(kp.begin(), kp.end());
	sc.scale(start); sc.scale(end);
	//main procedure
	try{
		//assemble tree from input data
		HM2D::Contour::Tree ext = HM2D::Contour::Tree::Assemble(*ecol);
		//cut by start and end point
		HM2D::EdgeData non_processed_edges;
		if (start_point){
			Point *p1, *p2;
			auto _av1 = HM2D::AllVertices(ext.alledges());
			auto _f1 = HM2D::Finder::ClosestPoint(_av1, start);
			p1 = _av1[std::get<0>(_f1)].get();
			HM2D::EdgeData cn = ext.find_node(p1)->contour;
			if (end_point){
				auto _av2 = HM2D::AllVertices(cn);
				auto _f2 = HM2D::Finder::ClosestPoint(_av2, end);
				p2 = _av2[std::get<0>(_f2)].get();
			} else if (HM2D::Contour::IsClosed(cn)) p2 = p1;
			else {
				double d1 = Point::meas(*HM2D::Contour::First(cn), *p1);
				double d2 = Point::meas(*HM2D::Contour::Last(cn), *p1);
				if (d1 < d2) p2 = HM2D::Contour::Last(cn).get();
				else p2 = HM2D::Contour::First(cn).get();
			}
			cn = HM2D::Contour::Assembler::ShrinkContour(cn, p1, p2);
			for (auto e: ext.alledges()){
				if (std::find(cn.begin(), cn.end(), e) == cn.end()){
					non_processed_edges.push_back(e);
				}
			}
			ext = HM2D::Contour::Tree();
			ext.add_contour(cn);
		}

		//assemble significant points
		if (keepbnd){
			int i = 0;
			for (auto e: *ecol) e->id = btypes[i++];
		}
		HM2D::EdgeData simpcol = HM2D::ECol::Algos::Simplified(ext.alledges(), a0, keepbnd);
		HM2D::VertexData keep_points = HM2D::AllVertices(simpcol);
		//keep points
		for (int i=0; i<n_keep_pts; ++i){
			auto _ae = ext.alledges();
			HM2D::Edge* efnd = _ae[std::get<0>(HM2D::Finder::ClosestEdge(_ae, kp[i]))].get();
			HM2D::EdgeData* cont = &ext.find_node(efnd)->contour;
			auto gp = HM2D::Contour::GuaranteePoint(*cont, kp[i]);
			keep_points.push_back(std::get<1>(gp));
		}
		//cross points
		for (auto& c: vcrosses){
			auto cconts = HM2D::Contour::Assembler::AllContours(*c);
			for (auto& cond: cconts)
			for (int i=0; i<ext.nodes.size(); ++i){
				HM2D::EdgeData& cont = ext.nodes[i]->contour;
				auto cres = HM2D::Contour::Finder::CrossAll(cont, cond);
				for (auto cross: cres){
					auto gp = HM2D::Contour::GuaranteePoint(cont, std::get<1>(cross));
					keep_points.push_back(std::get<1>(gp));
				}
			}
		}
		//checks
		if (algo != 0 && ext.nodes.size() != 1){
			throw std::runtime_error("only singly connected contours "
					"are allowed for refference point partition");
		}
		if (algo > 1 && !start_point){
			throw std::runtime_error("define start point "
				"for refference partition");
		}
		if (ext.nodes.size() == 0 || ext.nodes[0]->contour.size() == 0){
			throw std::runtime_error("zero length contour can not be parted");
		}
		//call partition algorithm
		TCont out;
		auto& c0 = ext.nodes[0]->contour;
		if (algo == 0){
			out = const_partition(ext, basic_steps[0], keep_points, nedges);
		} else if (algo == 1){
			out = refp_partition(c0, basic_steps, basic_points, keep_points, nedges);
		} else if (algo == 2){
			HM2D::Contour::R::ForceFirst ff(c0, start);
			for (int i=0; i<n_steps; ++i){
				basic_points.push_back(HM2D::Contour::WeightPoint(
					c0, steps[2*i+1]));
			}
			out = refp_partition(c0, basic_steps, basic_points, keep_points, nedges);
		} else if (algo == 3){
			HM2D::Contour::R::ForceFirst ff(c0, start);
			double len = HM2D::Length(c0);
			for (int i=0; i<n_steps; ++i){
				double s = steps[2*i+1]/sc.L/len;
				if (s < 0) s = 1+s;
				basic_points.push_back(HM2D::Contour::WeightPoint(c0, s));
			}
			out = refp_partition(c0, basic_steps, basic_points, keep_points, nedges);
		}
		//add non-processed
		HM2D::EdgeData outcol = out;
		if (non_processed_edges.size() > 0){
			outcol.insert(outcol.end(), non_processed_edges.begin(),
					non_processed_edges.end());
			HM2D::ECol::Algos::MergePoints(outcol);
		}
		//boundary assignment
		*n_outbnd = outcol.size();
		*outbnd = new int[*n_outbnd];
		set_ecollection_bc_force(cont, &outcol, btypes, *outbnd, 3);
		//unscale and return value
		HM2D::Unscale(outcol, sc);
		TCont contret;
		DeepCopy(outcol, contret);
		ret = new TCont(std::move(contret));
	} catch (std::runtime_error &e){
		ret = NULL;
		std::cout<<e.what()<<std::endl;
	}
	HM2D::Unscale(*ecol, sc);
	for (auto& c: vcrosses) HM2D::Unscale(*c, sc);
	return ret;
}

void* spline(int npnt, double* pnts, int nbtypes, int* btypes, int nedges,
		int* n_outbnd, int** outbnd){
	HM2D::EdgeData* ret = NULL;
	try{
		//build points
		vector<Point> p;
		for (int i=0; i<2*npnt; i+=2) p.push_back(Point(pnts[i], pnts[i+1]));
		auto sc = ScaleBase::doscale(p.begin(), p.end());
		//create contour
		HM2D::EdgeData spline = HM2D::Contour::Constructor::Spline(p, nedges);
		//assign boundary types
		vector<int> vbt(btypes, btypes + nbtypes);
		vbt.resize(p.size()-1, vbt.back());
		auto op = HM2D::Contour::OrderedPoints(spline);
		vector<int> basis_index;
		for (int i=0; i<p.size(); ++i){
			for (int j=0; j<op.size(); ++j){
				if (p[i] == *op[j]){
					basis_index.push_back(j);
					break;
				}
			}
		}
		basis_index.back() = op.size()-1;
		*n_outbnd = spline.size();
		*outbnd = new int[spline.size()];
		for (int i=0; i<basis_index.size() - 1; ++i){
			int i1 = basis_index[i];
			int i2 = basis_index[i+1];
			for (int j=i1; j<i2; ++j) (*outbnd)[j] = vbt[i];
		}
		// Unscale and return
		HM2D::Unscale(spline, sc);
		ret = new HM2D::EdgeData();
		HM2D::DeepCopy(spline, *ret);
	} catch (std::runtime_error &e){
		if (ret != 0) delete ret; 
		ret = NULL;
		std::cout<<e.what()<<std::endl;
	}
	return ret;
}

void* matched_partition(void* cont, int ncond, void** conds, int npts, double* pcond, double step,
		double influence_dist, double pw, double a0){
	typedef HM2D::EdgeData ContEC;
	ContEC* ret = new ContEC();
	vector<HM2D::EdgeData> input;
	{
		auto ac = HM2D::Contour::Assembler::AllContours(
				*static_cast<ContEC*>(cont));
		for (auto& c: ac) input.push_back(c);
	}
	vector<HM2D::EdgeData> conditions;
	for (int i=0; i<ncond; ++i){
		auto ac = HM2D::Contour::Assembler::AllContours(
				*static_cast<ContEC*>(conds[i]));
		for (auto& c: ac) conditions.push_back(c);
	}
	vector<std::pair<Point, double>> pconditions;
	for (int i=0; i<npts; ++i){
		Point p1(pcond[3*i+1], pcond[3*i+2]);
		pconditions.emplace(pconditions.begin(), p1, pcond[3*i]);
	}

	ScaleBase sc = HM2D::Scale01(*static_cast<ContEC*>(cont));
	for (auto& c: conditions) HM2D::Scale(c, sc);
	step/=sc.L;
	influence_dist/=sc.L;
	for (auto& x: pconditions){ sc.scale(x.first); x.second/=sc.L; }

	try{
		for (auto& icont: input){
			//set angle points
			HM2D::VertexData fixpoints;
			if (a0>0) for (auto& info: HM2D::Contour::OrderedInfo(icont)){
				//end points of open contour
				if (info.pprev == NULL || info.pnext == NULL){
					fixpoints.push_back(info.p);
					continue;
				}
				//angle between edges
				double angle = (Angle(*info.pprev, *info.p, *info.pnext)-M_PI)/M_PI*180.0;
				if (ISGREATER(fabs(angle), a0)){
					fixpoints.push_back(info.p);
					continue;
				}
			} else {
				fixpoints = HM2D::AllVertices(icont);
			}
			//set cross points
			for (auto& cond: conditions){
				auto cr = HM2D::Contour::Finder::CrossAll(icont, cond);
				for (auto& icr: cr){
					Point p = std::get<1>(icr);
					auto p1 = std::get<1>(HM2D::Contour::GuaranteePoint(icont, p));
					fixpoints.push_back(p1);
				}
			}
			//build partition
			auto r = HM2D::Contour::Algos::ConditionalPartition(icont, step, influence_dist,
					conditions, pconditions, pw, fixpoints);
			//add to answer
			HM2D::DeepCopy(r, *ret);
		}
	} catch (std::runtime_error &e){
		if (ret!=0) delete ret;
		ret = NULL;
		std::cout<<e.what()<<std::endl;
	}
	HM2D::Unscale(*static_cast<ContEC*>(cont), sc);
	for (auto& c: conditions) HM2D::Unscale(c, sc);
	if (ret!=NULL) HM2D::Unscale(*static_cast<ContEC*>(ret), sc);
	return ret;

}

int segment_part(double start, double end, double h0, double h1,
		int n_internals, double* h_internals, int* nout, double** hout){
	try{
		if (start>=end) throw std::runtime_error("start > end");
		//scaling
		double sc = end - start;
	
		//assembling conditions data
		std::map<double, double> conds;
		conds[0] = h0/sc;
		conds[1] = h1/sc;
		for (int i=0; i<n_internals; ++i){
			if (h_internals[2*i]>start && h_internals[2*i]<end){
				conds[(h_internals[2*i]-start)/sc] = h_internals[2*i+1]/sc;
			}
		}
		//make partition
		auto icont = HM2D::Contour::Constructor::FromPoints({Point(0, 0), Point(1, 0)});
		auto ocont = HM2D::Contour::Algos::WeightedPartition(conds, icont);
		auto opoints = HM2D::Contour::OrderedPoints(ocont);
		if (opoints.size()<2) throw std::runtime_error("failed to build partition");

		*nout = opoints.size();
		*hout = new double[*nout];
		for (int i=0; i<*nout; ++i) (*hout)[i] = sc*opoints[i]->x+start;
		return 1;
	} catch (std::runtime_error &e){
		std::cout<<e.what()<<std::endl;
		return 0;
	}
}

int extract_contour(void* source, int npnts, double* pnts, const char* method, void*** ret){
	int ok;
	vector<Point> p0;
	for (int i=0; i<npnts; ++i) p0.emplace_back(pnts[2*i], pnts[2*i+1]);
	std::string m(method);
	HM2D::EdgeData* ss = static_cast<HM2D::EdgeData*>(source);
	ScaleBase sc = HM2D::Scale01(*ss);
	sc.scale(p0.begin(), p0.end());
	try{
		vector<Point*> p1(p0.size(), nullptr);
		if (p0.size()<2) throw std::runtime_error("insufficient number of base points");
		vector<HM2D::EdgeData> et = HM2D::Contour::Assembler::SimpleContours(*ss);
		//1) find source contour
		HM2D::EdgeData* src=nullptr;
		double mind = 1e32;
		for (auto& c: et){
			auto ce1 = HM2D::Finder::ClosestEdge(*ss, p0[0]);
			if (std::get<0>(ce1) >= 0)
			if (std::get<1>(ce1) < mind){
				mind = std::get<1>(ce1);
				src = &c;
			}
		}
		if (src == nullptr) throw std::runtime_error("source contour was not found");
		//2) project base points
		if (m == "corner"){
			auto cp = HM2D::Contour::CornerPoints1(*src);
			for (int i=0; i<p0.size(); ++i){
				double d0 = 1e32;
				for (int j=0; j<cp.size(); ++j){
					double m = Point::meas(p0[i], *cp[j]);
					if (m < d0){
						d0 = m;
						p0[i].set(*cp[j]);
						p1[i] = cp[j].get();
					}
				}
			}
		} else if (m == "line"){
			for (int i=0; i<p0.size(); ++i){
				auto gp0 = HM2D::Contour::GuaranteePoint(*src, p0[i]);
				p0[i].set(*std::get<1>(gp0));
				p1[i] = std::get<1>(gp0).get();
			}
		} else {
			auto sv = HM2D::AllVertices(*src);
			for (int i=0; i<p0.size(); ++i){
				auto fndc = HM2D::Finder::ClosestPoint(sv, p0[i]);
				p1[i] = sv[std::get<0>(fndc)].get();
			}
		}
		//3) Reverse p0 if needed
		bool reversed_order = false;
		if (HM2D::Contour::IsClosed(*src)){
			if (HM2D::Contour::Area(*src) < 0) HM2D::Contour::Reverse(*src);
			if (p0.size() > 2){
				double w1 = std::get<1>(HM2D::Contour::CoordAt(*src, p0[0]));
				double w2 = std::get<1>(HM2D::Contour::CoordAt(*src, p0[1]));
				double w3 = std::get<1>(HM2D::Contour::CoordAt(*src, p0[2]));
				if (w2 <= w1) w2 += 1.0;
				if (w3 <= w1) w3 += 1.0;
				if (w3 < w2) reversed_order = true;
			}
		}
		//4) assemble
		*ret = new void*[npnts-1];
		for (int i=0; i<npnts-1; ++i){
			int i1 = i, i2 = i+1;
			if (reversed_order) std::swap(i1, i2);
			//auto c = HM2D::Contour::Constructor::CutContour(*src, p0[i1], p0[i2]);
			auto c = HM2D::Contour::Assembler::ShrinkContour(*src, p1[i1], p1[i2]);
			auto econt = new HM2D::EdgeData(); 
			HM2D::DeepCopy(c, *econt);
			HM2D::Unscale(*econt, sc);
			(*ret)[i] = econt;
		}
		ok = 1;
	} catch (std::runtime_error &e){
		std::cout<<e.what()<<std::endl;
		ok = 0;
	}
	HM2D::Unscale(*ss, sc);
	return ok;
}

void* cwriter_create(const char* cname, void* cont, void* awriter, void* subnode, const char* fmt){
	try{
		HMXML::ReaderA* wr = static_cast<HMXML::ReaderA*>(awriter);
		HMXML::Reader* sn = static_cast<HMXML::ReaderA*>(subnode);
		HM2D::EdgeData* c = static_cast<HM2D::EdgeData*>(cont);
		return new HM2D::Export::EColWriter(*c, wr, sn, cname, fmt);
	} catch (std::runtime_error &e){
		std::cout<<e.what()<<std::endl;
		return 0;
	}
}
void cwriter_free(void* cwriter){
	try{
		delete static_cast<HM2D::Export::EColWriter*>(cwriter);
	} catch (std::runtime_error &e){
		std::cout<<e.what()<<std::endl;
	}
}
int cwriter_add_edge_field(void* cwriter, const char* fieldname, void* field, int fsize, const char* type){
	try{
		auto cw = static_cast<HM2D::Export::EColWriter*>(cwriter);
		std::string fn(fieldname);
		std::string tpstr(type);
		if (tpstr == "int"){
			int* intfield = static_cast<int*>(field);
			std::vector<int> data(intfield, intfield+fsize);
			cw->AddEdgeData(fn, data, cw->is_binary<int>());
		} else if (tpstr == "char"){
			char* chrfield = static_cast<char*>(field);
			std::vector<char> data(chrfield, chrfield+fsize);
			cw->AddEdgeData(fn, data, cw->is_binary<char>());
		} else if (tpstr == "double"){
			double* chrfield = static_cast<double*>(field);
			std::vector<double> data(chrfield, chrfield+fsize);
			cw->AddEdgeData(fn, data, cw->is_binary<double>());
		} else if (tpstr == "float"){
			float* chrfield = static_cast<float*>(field);
			std::vector<float> data(chrfield, chrfield+fsize);
			cw->AddEdgeData(fn, data, cw->is_binary<float>());
		} else {
			throw std::runtime_error("unknown data type "+tpstr);
		}
		return 1;
	} catch (std::runtime_error &e){
		std::cout<<e.what()<<std::endl;
		return 0;
	}
}
void* creader_create(void* awriter, void* subnode, char* outname){
	try{
		auto wr = static_cast<HMXML::ReaderA*>(awriter);
		auto sn = static_cast<HMXML::Reader*>(subnode);
		//name
		std::string nm = sn->attribute(".", "name");
		if (nm.size()>1000) throw std::runtime_error("grid name is too long: " + nm);
		strcpy(outname, nm.c_str());
		//reader
		HM2D::Import::EColReader* ret = new HM2D::Import::EColReader(wr, sn);
		return ret;
	} catch (std::runtime_error &e){
		std::cout<<e.what()<<std::endl;
		return 0;
	}
}
void* creader_getresult(void* rd){
	try{
		auto reader = static_cast<HM2D::Import::EColReader*>(rd);
		return reader->result.release();
	} catch (std::runtime_error &e){
		std::cout<<e.what()<<std::endl;
		return 0;
	}
}
void* creader_read_edge_field(void* rd, const char* fieldname, const char* type){
	try{
		auto reader = static_cast<HM2D::Import::EColReader*>(rd);
		std::string fname(fieldname);
		std::string tpname(type);
		void* ret = NULL;
		if (tpname == "int"){
			auto outv = reader->read_edges_field<int>(fieldname);
			if (outv.size() != reader->Ne) throw std::runtime_error("field size doesn't match");
			ret = new int[outv.size()];
			std::copy(outv.begin(), outv.end(), (int*)ret);
		} else if (tpname == "char"){
			auto outv = reader->read_edges_field<char>(fieldname);
			if (outv.size() != reader->Ne) throw std::runtime_error("field size doesn't match");
			ret = new char[outv.size()];
			std::copy(outv.begin(), outv.end(), (char*)ret);
		} else if (tpname == "float"){
			auto outv = reader->read_edges_field<float>(fieldname);
			if (outv.size() != reader->Ne) throw std::runtime_error("field size doesn't match");
			ret = new float[outv.size()];
			std::copy(outv.begin(), outv.end(), (float*)ret);
		} else if (tpname == "double"){
			auto outv = reader->read_edges_field<double>(fieldname);
			if (outv.size() != reader->Ne) throw std::runtime_error("field size doesn't match");
			ret = new double[outv.size()];
			std::copy(outv.begin(), outv.end(), (double*)ret);
		}
		return ret;
	} catch (std::runtime_error &e){
		return 0;
	}
}
void creader_free(void* creader){
	try{
		delete (HM2D::Import::EColReader*)creader;
	} catch (std::runtime_error &e){
		std::cout<<e.what()<<std::endl;
	}
}

int unite_contours(int ncont, void** conts, void** retcont, int* Nlinks, int** links){
	typedef HM2D::EdgeData TECol;
	typedef HM2D::EdgeData TECont;
	TECol inpall;
	TECont ret;
	for (int i=0; i<ncont; ++i){
		auto c = static_cast<TECol*>(conts[i]);
		inpall.insert(inpall.end(), c->begin(), c->end());
	}
	ScaleBase sc = HM2D::Scale01(inpall);
	int r = 0;
	try{
		HM2D::DeepCopy(inpall, ret);
		for (int i=0; i<ret.size(); ++i) ret[i]->id = i;
		HM2D::ECol::Algos::MergePoints(ret);
		*Nlinks = ret.size();
		*links = new int[*Nlinks];
		for (int i=0; i<ret.size(); ++i){
			(*links)[i] = ret[i]->id;
		}
		HM2D::Unscale(ret, sc);
		(*retcont) = new TECont(std::move(ret));
		r = 1;
	} catch (std::exception &e){
		std::cout<<e.what()<<std::endl;
	}
	HM2D::Unscale(inpall, sc);
	return r;
}
int simplify_contour(void* cont, double angle, int* btypes, void** ret_cont, int* Nretb, int** retb){
	typedef HM2D::EdgeData TECol;
	typedef HM2D::EdgeData TECont;
	TECol* inp = static_cast<TECol*>(cont);
	ScaleBase sc = HM2D::Scale01(*inp);
	int r=0;
	try{
		for (int i=0; i<inp->size(); ++i) (*inp)[i]->id = btypes[i];
		TECol ret1 = HM2D::ECol::Algos::Simplified(*inp, angle, true);
		TECont ret_cont1;
		HM2D::DeepCopy(ret1, ret_cont1);
		*Nretb = ret_cont1.size();
		*retb = new int[*Nretb];
		for (int i=0; i<*Nretb; ++i) (*retb)[i] = ret_cont1[i]->id;
		HM2D::Unscale(ret_cont1, sc);
		*ret_cont = new TECont(std::move(ret_cont1));
		r = 1;
	} catch (std::exception &e){
		std::cout<<e.what()<<std::endl;
	}
	HM2D::Unscale(*inp, sc);
	return r;
}
int separate_contour(void* cont, int* btypes, int* Nretc, void*** retc, int* Nretb, int** retb){
	HM2D::EdgeData* inp = static_cast<HM2D::EdgeData*>(cont);
	ScaleBase sc = HM2D::Scale01(*inp);
	int r = 0;
	try{
		for (int i=0; i<inp->size(); ++i) (*inp)[i]->id = btypes[i];
		vector<HM2D::EdgeData> sep = HM2D::Contour::Constructor::ExtendedSeparate(*inp);
		vector<HM2D::EdgeData> sepr;
		for (auto& s: sep){
			sepr.emplace_back();
			HM2D::DeepCopy(s, sepr.back());
		}
		*Nretc = sepr.size();
		*Nretb = 0;
		for (auto& s: sepr) *Nretb += s.size();
		*retb = new int[*Nretb];
		int k=0;
		for (auto& s: sepr)
		for (int i=0; i<s.size(); ++i){
			(*retb)[k++] = s[i]->id;
		}
		(*retc) = new void*[sepr.size()];
		for (int i=0; i<sepr.size(); ++i){
			HM2D::Unscale(sepr[i], sc);
			(*retc)[i] = new HM2D::EdgeData(std::move(sepr[i]));
		}
		r = 1;
	} catch (std::exception &e){
		std::cout<<e.what()<<std::endl;
	}
	HM2D::Unscale(*inp, sc);
	return r;
}

int quick_separate_contour(void* cont, int* btypes, int* Nretc, void*** retc, int* Nretb, int** retb){
	HM2D::EdgeData* inp = static_cast<HM2D::EdgeData*>(cont);
	int r = 0;
	try{
		for (int i=0; i<inp->size(); ++i) (*inp)[i]->id = btypes[i];
		vector<HM2D::EdgeData> sep = HM2D::SplitData(*inp);
		vector<HM2D::EdgeData> sepr(sep.size());
		for (int i=0; i<sep.size(); ++i){
			HM2D::DeepCopy(sep[i], sepr[i]);
		}
		*Nretc = sep.size();
		*Nretb = 0;
		for (auto& s: sep) *Nretb += s.size();
		*retb = new int[*Nretb];
		int k=0;
		for (auto& s: sep)
		for (int i=0; i<s.size(); ++i){
			(*retb)[k++] = s[i]->id;
		}
		(*retc) = new void*[sepr.size()];
		for (int i=0; i<sepr.size(); ++i){
			(*retc)[i] = new HM2D::EdgeData(std::move(sepr[i]));
		}
		r = 1;
	} catch (std::exception &e){
		std::cout<<e.what()<<std::endl;
	}
	return r;
}

int connect_subcontours(int nconts, void** conts, int nfix, int* fix,  const char* close, int shiftnext, void** cret){
	int ok;
	vector<HM2D::EdgeData*> ecols(nconts);
	int k=0;
	for (int i=0; i<nconts; ++i){
		ecols[i] = static_cast<HM2D::EdgeData*>(conts[i]);
		//assign indices for sorting edges in return collection
		for (int j=0; j<ecols[i]->size(); ++j) (*ecols[i])[j]->id = k++;
	}
	HM2D::VertexData allpoints;
	for (auto e: ecols)
	for (auto p: HM2D::AllVertices(*e))
		allpoints.push_back(p);
	
	ScaleBase sc = ScaleBase::p_doscale(allpoints.begin(), allpoints.end());
	std::string close_method(close);
	try{
		//assemble contours
		vector<HM2D::EdgeData> vconts;
		for (int i=0; i<nconts; ++i) if (ecols[i]->size()>0){
			vconts.push_back(HM2D::Contour::Assembler::Contour1(*ecols[i], (*ecols[i])[0]->first().get()));
		}
		//construct new contour
		std::set<int> fixset(fix, fix+nfix);
		HM2D::EdgeData ret = HM2D::Contour::Constructor::FromContours(
			vconts, close_method=="yes", shiftnext==1, fixset);
		if (close_method == "force" && HM2D::Contour::IsOpen(ret)){
			HM2D::Contour::AddLastPoint(ret, HM2D::Contour::First(ret));
		}
		HM2D::Unscale(ret, sc);
		//return value
		auto mm = new HM2D::EdgeData();
		HM2D::DeepCopy(ret, *mm);
		std::sort(mm->begin(), mm->end(),
				[](shared_ptr<HM2D::Edge> e1, shared_ptr<HM2D::Edge> e2)->bool{
					return e1->id < e2->id;
				});
		*cret = mm;
		ok = 1;
	} catch (std::exception &e){
		std::cout<<e.what()<<std::endl;
		ok = 0;
	}
	sc.p_unscale(allpoints.begin(), allpoints.end());
	return ok;
}
