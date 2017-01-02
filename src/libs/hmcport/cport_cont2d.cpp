#include <unordered_map>
#include "hmproject.h"
#include "cport_cont2d.h"
#include "hybmesh_contours2d.hpp"
#include "hmc_imex.hpp"

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
const_partition(const HMCont2D::ExtendedTree& etree, double step, const vector<Point*>& keep, int ne){
	HMCont2D::PCollection _pdata;
	vector<HMCont2D::Contour> partdata;
	int cc = etree.cont_count();
	if (ne <= 0) for (int i=0; i<cc; ++i){
		partdata.push_back(
			HMCont2D::Algos::Partition(step, *etree.get_contour(i),
				_pdata, keep));
	} else {
		vector<double> lens;
		for (int i=0; i<cc; ++i) lens.push_back(etree.get_contour(i)->length());
		//minimum size includes all keep points
		vector<int> minsize(cc, 0);
		for (int i=0; i<keep.size(); ++i){
			auto c = etree.get_contour(keep[i]);
			if (c == NULL) continue;
			else for (int j=0; j<cc; ++j){
				if (etree.get_contour(j) == c) {
					if (c->is_closed() ||
						(keep[i] != c->first() && keep[i] != c->last()))
							++minsize[j];
					break;
				}
			}
		}
		for (int i=0; i<cc; ++i){
			if (!etree.get_contour(i)->is_closed()) minsize[i]++;
			if (etree.get_contour(i)->is_closed() && minsize[i]<3) minsize[i] = 3;
		}
		//nums double
		vector<double> nums_double;
		double full_len = std::accumulate(lens.begin(), lens.end(), 0.0);
		for (int i=0; i<lens.size(); ++i){
			nums_double.push_back(lens[i]/full_len*ne);
		}
		//nums int
		vector<int> nums = HMCont2D::Algos::RoundVector(nums_double, minsize);
		//building
		std::map<double, double> m; m[0]=step;
		for (int i=0; i<etree.cont_count(); ++i){
			partdata.push_back(
				HMCont2D::Algos::WeightedPartition(m, *etree.get_contour(i),
					_pdata, nums[i], keep));
		}
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
		const vector<Point*>& keep, int ne){
	if (bpoints.size() == 1){
		HMCont2D::ExtendedTree etree;
		etree.AddContour(cont);
		return const_partition(etree, steps[0], keep, ne);
	}
	//assemble weights
	std::map<double, double> weights;
	for (int i=0; i<steps.size(); ++i){
		auto info = cont.coord_at(bpoints[i]);
		weights[std::get<1>(info)] = steps[i];
	}
	//call partition
	HMCont2D::PCollection _pdata;
	HMCont2D::Contour cret = (ne<=0) ? HMCont2D::Algos::WeightedPartition(weights, cont, _pdata, keep)
	                                 : HMCont2D::Algos::WeightedPartition(weights, cont, _pdata, ne, keep);
	//assemble return data
	HMCont2D::Container<HMCont2D::ECollection> ret;
	HMCont2D::Container<HMCont2D::ECollection>::DeepCopy(cret, ret);
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
	ScaleBase sc = HMCont2D::ECollection::Scale01(*ecol);
	sc.scale(basic_points.begin(), basic_points.end());
	for (auto& v: basic_steps) v/=sc.L;
	for (auto& c: vcrosses) HMCont2D::ECollection::Scale(*c, sc);
	sc.scale(kp.begin(), kp.end());
	sc.scale(start); sc.scale(end);
	//main procedure
	try{
		//assemble tree from input data
		HMCont2D::ExtendedTree ext = HMCont2D::Assembler::ETree(*ecol);
		//cut by start and end point
		HMCont2D::ECollection non_processed_edges;
		if (start_point){
			Point *p1, *p2;
			p1 = HMCont2D::ECollection::FindClosestNode(ext, start);
			HMCont2D::Contour cn = *ext.get_contour(p1);
			if (end_point){
				p2 = HMCont2D::ECollection::FindClosestNode(cn, end);
			} else if (cn.is_closed()) p2 = p1;
			else {
				double d1 = Point::meas(*cn.first(), *p1);
				double d2 = Point::meas(*cn.last(), *p1);
				if (d1 < d2) p2 = cn.last();
				else p2 = cn.first();
			}
			cn = HMCont2D::Assembler::Contour1(cn, p1, p2);
			for (auto e: ext.data){
				if (std::find(cn.begin(), cn.end(), e) == cn.end()){
					non_processed_edges.add_value(*e);
				}
			}
			ext = HMCont2D::ExtendedTree();
			ext.AddContour(cn);
		}

		//assemble significant points
		if (keepbnd){
			int i = 0;
			for (auto e: *ecol) e->id = btypes[i++];
		}
		HMCont2D::ECollection simpcol = HMCont2D::Algos::Simplified(ext, a0, keepbnd);
		std::vector<Point*> keep_points = simpcol.all_points();
		HMCont2D::PCollection pcol;
		//keep points
		for (int i=0; i<n_keep_pts; ++i){
			HMCont2D::Edge* efnd = std::get<0>(HMCont2D::ECollection::FindClosestEdge(ext, kp[i]));
			HMCont2D::Contour* cont = ext.get_contour(efnd);
			auto gp = cont->GuaranteePoint(kp[i], pcol);
			keep_points.push_back(std::get<1>(gp));
		}
		//cross points
		for (auto& c: vcrosses){
			auto cconts = HMCont2D::Assembler::AllContours(*c);
			for (auto& cond: cconts)
			for (int i=0; i<ext.cont_count(); ++i){
				HMCont2D::Contour& cont = *ext.get_contour(i);
				auto cres = HMCont2D::Algos::CrossAll(cont, cond);
				for (auto cross: cres){
					auto gp = cont.GuaranteePoint(std::get<1>(cross), pcol);
					keep_points.push_back(std::get<1>(gp));
				}
			}
		}
		//checks
		if (algo != 0 && ext.cont_count() != 1){
			throw std::runtime_error("only singly connected contours "
					"are allowed for refference point partition");
		}
		if (algo > 1 && !start_point){
			throw std::runtime_error("define start point "
				"for refference partition");
		}
		if (ext.cont_count() == 0 || ext.get_contour(0)->size() == 0){
			throw std::runtime_error("zero length contour can not be parted");
		}
		//call partition algorithm
		TCont out;
		auto& c0 = *ext.get_contour(0);
		if (algo == 0){
			out = const_partition(ext, basic_steps[0], keep_points, nedges);
		} else if (algo == 1){
			out = refp_partition(c0, basic_steps, basic_points, keep_points, nedges);
		} else if (algo == 2){
			c0.StartFrom(start);
			for (int i=0; i<n_steps; ++i){
				basic_points.push_back(HMCont2D::Contour::WeightPoint(
					c0, steps[2*i+1]));
			}
			out = refp_partition(c0, basic_steps, basic_points, keep_points, nedges);
		} else if (algo == 3){
			c0.StartFrom(start);
			double len = c0.length();
			for (int i=0; i<n_steps; ++i){
				double s = steps[2*i+1]/sc.L/len;
				if (s < 0) s = 1+s;
				basic_points.push_back(HMCont2D::Contour::WeightPoint(c0, s));
			}
			out = refp_partition(c0, basic_steps, basic_points, keep_points, nedges);
		}
		//add non-processed
		HMCont2D::ECollection outcol = out;
		if (non_processed_edges.size() > 0){
			outcol.Unite(non_processed_edges);
			HMCont2D::Algos::MergePoints(outcol);
		}
		//boundary assignment
		*n_outbnd = outcol.size();
		*outbnd = new int[*n_outbnd];
		set_ecollection_bc_force(cont, &outcol, btypes, *outbnd, 3);
		//unscale and return value
		HMCont2D::ECollection::Unscale(outcol, sc);
		TCont contret;
		TCont::DeepCopy(outcol, contret);
		ret = new TCont(std::move(contret));
	} catch (std::runtime_error &e){
		ret = NULL;
		std::cout<<e.what()<<std::endl;
	}
	HMCont2D::ECollection::Unscale(*ecol, sc);
	for (auto& c: vcrosses) HMCont2D::ECollection::Unscale(*c, sc);
	return ret;
}

void* spline(int npnt, double* pnts, int nbtypes, int* btypes, int nedges,
		int* n_outbnd, int** outbnd){
	HMCont2D::Container<HMCont2D::Contour>* ret = NULL;
	try{
		//build points
		vector<Point> p;
		for (int i=0; i<2*npnt; i+=2) p.push_back(Point(pnts[i], pnts[i+1]));
		vector<Point*> pvec;
		for (int i=0; i<p.size(); ++i) pvec.push_back(&p[i]);
		ScaleBase sc = ScaleBase::p_doscale(pvec.begin(), pvec.end());
		if (*pvec[0] == *pvec.back()) pvec.back() = pvec[0];
		//create contour
		HMCont2D::PCollection pcol;
		HMCont2D::Contour spline = HMCont2D::Constructor::Spline(pvec, pcol, nedges);
		//assign boundary types
		vector<int> vbt(btypes, btypes + nbtypes);
		vbt.resize(pvec.size()-1, vbt.back());
		vector<Point*> op = spline.ordered_points();
		vector<int> basis_index;
		for (int i=0; i<pvec.size(); ++i){
			for (int j=0; j<op.size(); ++j){
				if (pvec[i] == op[j]){
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
		HMCont2D::ECollection::Unscale(spline, sc);
		ret = new HMCont2D::Container<HMCont2D::Contour>();
		HMCont2D::Container<HMCont2D::Contour>::DeepCopy(spline, *ret);
	} catch (std::runtime_error &e){
		if (ret != 0) delete ret; 
		ret = NULL;
		std::cout<<e.what()<<std::endl;
	}
	return ret;
}

void* matched_partition(void* cont, int ncond, void** conds, int npts, double* pcond, double step,
		double influence_dist, double pw, double a0){
	typedef HMCont2D::Container<HMCont2D::ECollection> ContEC;
	ContEC* ret = new ContEC();
	vector<HMCont2D::Contour> input;
	{
		auto ac = HMCont2D::Assembler::AllContours(
				*static_cast<ContEC*>(cont));
		for (auto& c: ac) input.push_back(c);
	}
	vector<HMCont2D::Contour> conditions;
	for (int i=0; i<ncond; ++i){
		auto ac = HMCont2D::Assembler::AllContours(
				*static_cast<ContEC*>(conds[i]));
		for (auto& c: ac) conditions.push_back(c);
	}
	vector<std::pair<Point, double>> pconditions;
	for (int i=0; i<npts; ++i){
		Point p1(pcond[3*i+1], pcond[3*i+2]);
		pconditions.emplace(pconditions.begin(), p1, pcond[3*i]);
	}

	ScaleBase sc = HMCont2D::ECollection::Scale01(*static_cast<ContEC*>(cont));
	for (auto& c: conditions) HMCont2D::ECollection::Scale(c, sc);
	step/=sc.L;
	influence_dist/=sc.L;
	for (auto& x: pconditions){ sc.scale(x.first); x.second/=sc.L; }

	try{
		HMCont2D::PCollection pcol;
		for (auto& icont: input){
			//set angle points
			vector<Point*> fixpoints;
			if (a0>0) for (auto& info: icont.ordered_info()){
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
				fixpoints = icont.all_points();
			}
			//set cross points
			for (auto& cond: conditions){
				auto cr = HMCont2D::Algos::CrossAll(icont, cond);
				for (auto& icr: cr){
					Point p = std::get<1>(icr);
					Point* p1 = std::get<1>(icont.GuaranteePoint(p, pcol));
					fixpoints.push_back(p1);
				}
			}
			//build partition
			auto r = HMCont2D::Algos::ConditionalPartition(icont, step, influence_dist,
					conditions, pconditions, pw, pcol, fixpoints);
			//add to answer
			ContEC::DeepCopy(r, *ret);
		}
	} catch (std::runtime_error &e){
		if (ret!=0) delete ret;
		ret = NULL;
		std::cout<<e.what()<<std::endl;
	}
	HMCont2D::ECollection::Unscale(*static_cast<ContEC*>(cont), sc);
	for (auto& c: conditions) HMCont2D::ECollection::Unscale(c, sc);
	if (ret!=NULL) HMCont2D::ECollection::Unscale(*static_cast<ContEC*>(ret), sc);
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
		HMCont2D::PCollection pstore;
		auto icont = HMCont2D::Constructor::ContourFromPoints({Point(0, 0), Point(1, 0)});
		auto ocont = HMCont2D::Algos::WeightedPartition(conds, icont, pstore);
		auto opoints = ocont.ordered_points();
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
	HMCont2D::ECollection* ss = static_cast<HMCont2D::ECollection*>(source);
	ScaleBase sc = HMCont2D::ECollection::Scale01(*ss);
	sc.scale(p0.begin(), p0.end());
	try{
		if (p0.size()<2) throw std::runtime_error("insufficient number of base points");
		vector<HMCont2D::Contour> et = HMCont2D::Assembler::SimpleContours(*ss);
		//1) find source contour
		HMCont2D::Contour* src=nullptr;
		double mind = 1e32;
		for (auto& c: et){
			auto ce1 = HMCont2D::ECollection::FindClosestEdge(*ss, p0[0]);
			if (std::get<0>(ce1) != nullptr)
			if (std::get<1>(ce1) < mind){
				mind = std::get<1>(ce1);
				src = &c;
			}
		}
		if (src == nullptr) throw std::runtime_error("source contour was not found");
		//2) project base points
		HMCont2D::PCollection pcol;
		if (m == "corner"){
			auto cp = src->corner_points1();
			for (int i=0; i<p0.size(); ++i){
				double d0 = 1e32;
				for (int j=0; j<cp.size(); ++j){
					double m = Point::meas(p0[i], *cp[j]);
					if (m < d0){
						d0 = m;
						p0[i].set(*cp[j]);
					}
				}
			}
		} else if (m == "line"){
			for (int i=0; i<p0.size(); ++i){
				auto gp0 = src->GuaranteePoint(p0[i], pcol);
				p0[i].set(*std::get<1>(gp0));
			}
		}
		//3) Reverse p0 if needed
		bool reversed_order = false;
		if (src->is_closed()){
			if (HMCont2D::Area(*src) < 0) src->Reverse();
			if (p0.size() > 2){
				double w1 = std::get<1>(src->coord_at(p0[0]));
				double w2 = std::get<1>(src->coord_at(p0[1]));
				double w3 = std::get<1>(src->coord_at(p0[2]));
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
			auto c = HMCont2D::Assembler::Contour1(*src, p0[i1], p0[i2]);
			auto econt = new HMCont2D::Container<HMCont2D::ECollection>(); 
			HMCont2D::Container<HMCont2D::ECollection>::DeepCopy(c, *econt);
			HMCont2D::ECollection::Unscale(*econt, sc);
			(*ret)[i] = econt;
		}
		ok = 1;
	} catch (std::runtime_error &e){
		std::cout<<e.what()<<std::endl;
		ok = 0;
	}
	HMCont2D::ECollection::Unscale(*ss, sc);
	return ok;
}

void* cwriter_create(const char* cname, void* cont, void* awriter, void* subnode, const char* fmt){
	try{
		HMXML::ReaderA* wr = static_cast<HMXML::ReaderA*>(awriter);
		HMXML::Reader* sn = static_cast<HMXML::ReaderA*>(subnode);
		HMCont2D::ECollection* c = static_cast<HMCont2D::ECollection*>(cont);
		return new HMCont2D::Export::EColWriter(*c, wr, sn, cname, fmt);
	} catch (std::runtime_error &e){
		std::cout<<e.what()<<std::endl;
		return 0;
	}
}
void cwriter_free(void* cwriter){
	try{
		delete static_cast<HMCont2D::Export::EColWriter*>(cwriter);
	} catch (std::runtime_error &e){
		std::cout<<e.what()<<std::endl;
	}
}
int cwriter_add_edge_field(void* cwriter, const char* fieldname, void* field, int fsize, const char* type){
	try{
		auto cw = static_cast<HMCont2D::Export::EColWriter*>(cwriter);
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
		HMCont2D::Import::EColReader* ret = new HMCont2D::Import::EColReader(wr, sn);
		return ret;
	} catch (std::runtime_error &e){
		std::cout<<e.what()<<std::endl;
		return 0;
	}
}
void* creader_getresult(void* rd){
	try{
		auto reader = static_cast<HMCont2D::Import::EColReader*>(rd);
		return reader->result.release();
	} catch (std::runtime_error &e){
		std::cout<<e.what()<<std::endl;
		return 0;
	}
}
void* creader_read_edge_field(void* rd, const char* fieldname, const char* type){
	try{
		auto reader = static_cast<HMCont2D::Import::EColReader*>(rd);
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
		delete (HMCont2D::Import::EColReader*)creader;
	} catch (std::runtime_error &e){
		std::cout<<e.what()<<std::endl;
	}
}

int unite_contours(int ncont, void** conts, void** retcont, int* Nlinks, int** links){
	typedef HMCont2D::ECollection TECol;
	typedef HMCont2D::Container<HMCont2D::ECollection> TECont;
	TECol inpall;
	TECont ret;
	for (int i=0; i<ncont; ++i){
		auto c = static_cast<TECol*>(conts[i]);
		inpall.Unite(*c);
	}
	ScaleBase sc = HMCont2D::ECollection::Scale01(inpall);
	int r = 0;
	try{
		TECont::DeepCopy(inpall, ret);
		for (int i=0; i<ret.size(); ++i) ret.data[i]->id = i;
		HMCont2D::Algos::MergePoints(ret);
		HMCont2D::Algos::DeleteUnusedPoints(ret);
		*Nlinks = ret.size();
		*links = new int[*Nlinks];
		for (int i=0; i<ret.size(); ++i){
			(*links)[i] = ret.data[i]->id;
		}
		HMCont2D::ECollection::Unscale(ret, sc);
		(*retcont) = new TECont(std::move(ret));
		r = 1;
	} catch (std::exception &e){
		std::cout<<e.what()<<std::endl;
	}
	HMCont2D::ECollection::Unscale(inpall, sc);
	return r;
}
int simplify_contour(void* cont, double angle, int* btypes, void** ret_cont, int* Nretb, int** retb){
	typedef HMCont2D::ECollection TECol;
	typedef HMCont2D::Container<HMCont2D::ECollection> TECont;
	TECol* inp = static_cast<TECol*>(cont);
	ScaleBase sc = TECol::Scale01(*inp);
	int r=0;
	try{
		for (int i=0; i<inp->size(); ++i) inp->data[i]->id = btypes[i];
		TECol ret1 = HMCont2D::Algos::Simplified(*inp, angle, true);
		TECont ret_cont1;
		TECont::DeepCopy(ret1, ret_cont1);
		*Nretb = ret_cont1.size();
		*retb = new int[*Nretb];
		for (int i=0; i<*Nretb; ++i) (*retb)[i] = ret_cont1.data[i]->id;
		TECol::Unscale(ret_cont1, sc);
		*ret_cont = new TECont(std::move(ret_cont1));
		r = 1;
	} catch (std::exception &e){
		std::cout<<e.what()<<std::endl;
	}
	TECol::Unscale(*inp, sc);
	return r;
}
int separate_contour(void* cont, int* btypes, int* Nretc, void*** retc, int* Nretb, int** retb){
	typedef HMCont2D::ECollection TECol;
	typedef HMCont2D::Container<HMCont2D::ECollection> TECont;
	TECol* inp = static_cast<TECol*>(cont);
	ScaleBase sc = TECol::Scale01(*inp);
	int r = 0;
	try{
		for (int i=0; i<inp->size(); ++i) inp->data[i]->id = btypes[i];
		HMCont2D::PCollection pcol;
		vector<TECol> sep = HMCont2D::Assembler::ExtendedSeparate(*inp, pcol);
		vector<TECont> sepr(sep.size());
		for (int i=0; i<sep.size(); ++i){
			TECont::DeepCopy(sep[i], sepr[i]);
		}
		*Nretc = sep.size();
		*Nretb = 0;
		for (auto& s: sep) *Nretb += s.size();
		*retb = new int[*Nretb];
		int k=0;
		for (auto& s: sep)
		for (int i=0; i<s.size(); ++i){
			(*retb)[k++] = s.data[i]->id;
		}
		(*retc) = new void*[sepr.size()];
		for (int i=0; i<sepr.size(); ++i){
			TECol::Unscale(sepr[i], sc);
			(*retc)[i] = new TECont(std::move(sepr[i]));
		}
		r = 1;
	} catch (std::exception &e){
		std::cout<<e.what()<<std::endl;
	}
	TECol::Unscale(*inp, sc);
	return r;
}

int quick_separate_contour(void* cont, int* btypes, int* Nretc, void*** retc, int* Nretb, int** retb){
	typedef HMCont2D::ECollection TECol;
	typedef HMCont2D::Container<HMCont2D::ECollection> TECont;
	TECol* inp = static_cast<TECol*>(cont);
	int r = 0;
	try{
		for (int i=0; i<inp->size(); ++i) inp->data[i]->id = btypes[i];
		vector<TECol> sep = HMCont2D::Assembler::QuickSeparate(*inp);
		vector<TECont> sepr(sep.size());
		for (int i=0; i<sep.size(); ++i){
			TECont::DeepCopy(sep[i], sepr[i]);
		}
		*Nretc = sep.size();
		*Nretb = 0;
		for (auto& s: sep) *Nretb += s.size();
		*retb = new int[*Nretb];
		int k=0;
		for (auto& s: sep)
		for (int i=0; i<s.size(); ++i){
			(*retb)[k++] = s.data[i]->id;
		}
		(*retc) = new void*[sepr.size()];
		for (int i=0; i<sepr.size(); ++i){
			(*retc)[i] = new TECont(std::move(sepr[i]));
		}
		r = 1;
	} catch (std::exception &e){
		std::cout<<e.what()<<std::endl;
	}
	return r;
}

int connect_subcontours(int nconts, void** conts, int nfix, int* fix,  const char* close, int shiftnext, void** cret){
	int ok;
	vector<HMCont2D::ECollection*> ecols(nconts);
	int k=0;
	for (int i=0; i<nconts; ++i){
		ecols[i] = static_cast<HMCont2D::ECollection*>(conts[i]);
		//assign indices for sorting edges in return collection
		for (int j=0; j<ecols[i]->size(); ++j) ecols[i]->data[j]->id = k++;
	}
	vector<Point*> allpoints;
	for (auto e: ecols)
	for (auto p: e->all_points())
		allpoints.push_back(p);
	
	ScaleBase sc = ScaleBase::p_doscale(allpoints.begin(), allpoints.end());
	std::string close_method(close);
	try{
		//assemble contours
		vector<HMCont2D::Contour> vconts;
		for (int i=0; i<nconts; ++i) if (ecols[i]->size()>0){
			vconts.push_back(HMCont2D::Assembler::Contour1(*ecols[i], ecols[i]->data[0]->pstart));
		}
		//construct new contour
		std::set<int> fixset(fix, fix+nfix);
		HMCont2D::Container<HMCont2D::Contour> ret = HMCont2D::Constructor::ContourFromContours(
			vconts, close_method=="yes", shiftnext==1, fixset);
		if (close_method == "force" && ret.is_open()){
			ret.add_value(HMCont2D::Edge(ret.last(), ret.first()));
		}
		HMCont2D::Unscale(ret, sc);
		//return value
		auto mm = new HMCont2D::Container<HMCont2D::ECollection>();
		HMCont2D::Container<HMCont2D::ECollection>::DeepCopy(ret, *mm);
		std::sort(mm->data.begin(), mm->data.end(),
				[](shared_ptr<HMCont2D::Edge> e1, shared_ptr<HMCont2D::Edge> e2)->bool{
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
