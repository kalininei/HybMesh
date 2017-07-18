#include "sizefun.hpp"
#include "femassembly.hpp"
#include "assemble2d.hpp"
#include "piecewise.hpp"
#include "finder2d.hpp"
#include "treverter2d.hpp"
#include "hmfem.hpp"
#include "partcont.hpp"
#include "modcont.hpp"
//#######
#include "export2d_vtk.hpp"

using namespace HM2D;
using namespace HM2D::Grid;

vector<double> SizeFun::sz(const vector<Point>& p) const{
	vector<double> ret(p.size());
	for (int i=0; i<ret.size(); ++i) {
		try{ 
			ret[i] = sz(p[i]);
		} catch (std::runtime_error){
			ret[i] = -1;
		}
	}
	return ret;
}

namespace{

//This exception is used to detect whether is site segment belongs could be
//restored from aux grid edges or not.
//Method sz_boundary_segment throws this exception on fail.
struct _PointNotOnGridEdge: public std::runtime_error{
	_PointNotOnGridEdge(): std::runtime_error("point is not on grid edge"){}
};

//Size function for multiple non-crossed trees
struct SumSizeFun: public SizeFun{
	vector<shared_ptr<SizeFun>> szfuns;
	vector<Contour::Tree> bounds;

	SumSizeFun(const vector<shared_ptr<SizeFun>>& szf,
			const vector<Contour::Tree>& bnds){
		assert(szf.size() == bnds.size());
		szfuns = szf;
		bounds.resize(szf.size());
		for (int i=0; i<szf.size(); ++i){
			bounds[i] = Contour::Tree::DeepCopy(bnds[i], 3);
		}
	}

	double minstep() const override{
		double ret = szfuns[0]->minstep();
		for (int i=1; i<szfuns.size(); ++i){
			double m = szfuns[i]->minstep();
			if (m<ret) ret = m;
		}
		return ret;
	}
	double maxstep() const override{
		double ret = szfuns[0]->maxstep();
		for (int i=1; i<szfuns.size(); ++i){
			double m = szfuns[i]->maxstep();
			if (m>ret) ret = m;
		}
		return ret;
	}

	double sz(const Point& p) const override{
		for (int i=0; i<szfuns.size(); ++i){
			if (bounds[i].whereis(p) != OUTSIDE) return szfuns[i]->sz(p);
		}
		throw std::runtime_error("error detecting size outside of the area");
	}
	double sz_proj(const Point& p) const override{
		try{
			return sz(p);
		} catch (std::runtime_error){
			vector<Point> proj(szfuns.size());
			for (int i=0; i<szfuns.size(); ++i){
				proj[i] = HM2D::Finder::ClosestEPoint(bounds[i].alledges(), p);
			}
			int ibest = std::min_element(proj.begin(), proj.end(),
				[&p](const Point& a, const Point& b){
					return Point::meas(a, p)<Point::meas(b,p);
				})-proj.begin();
			return szfuns[ibest]->sz(proj[ibest]);
		}
	}

	vector<double> sz(const vector<Point>& p) const override{
		vector<int> pos(p.size(), -1);
		for (int i=0; i<szfuns.size(); ++i){
			vector<Point> p2;
			for (int j=0; j<p.size(); ++j) if (pos[j] == -1){
				p2.push_back(p[j]);
			}
			vector<int> wh = Contour::Finder::SortOutPoints(bounds[i], p2);
			auto whit = wh.begin();
			for (int j=0; j<p.size(); ++j) if (pos[j] == -1){
				if (*whit++ != OUTSIDE) pos[j] = i;
			}
		}

		vector<double> ret(p.size(), -1);
		for (int i=0; i<szfuns.size(); ++i){
			vector<Point> p2;
			for (int j=0; j<p.size(); ++j) if (pos[j] == i){
				p2.push_back(p[j]);
			}
			vector<double> ret2 = szfuns[i]->sz(p2);
			auto rit = ret2.begin();
			for (int j=0; j<p.size(); ++j) if (pos[j] == i){
				ret[j] = *rit++;
			}
			
		}

		return ret;
	}

	HMMath::LinearPiecewise sz_boundary_segment(const HM2D::EdgeData& pos) const override{
		Point p = *HM2D::Contour::First(pos);
		double ksi;
		for (int i=0; i<szfuns.size(); ++i)
		for (int j=0; j<bounds[i].nodes.size(); ++j)
		for (int k=0; k<bounds[i].nodes[j]->contour.size(); ++k){
			auto& e = bounds[i].nodes[j]->contour[k];
			if (isOnSection(p, *e->pfirst(), *e->plast(), ksi)){
				return szfuns[i]->sz_boundary_segment(pos);
			}
		}
		throw _PointNotOnGridEdge();
	}
};

//size function with constant value
struct ConstSizeFun: public SizeFun{
	const double val;

	ConstSizeFun(double v): val(v){};

	double minstep() const override{ return val; }
	double maxstep() const override{ return val; }

	double sz(const Point& p) const override{ return val; }
	double sz_proj(const Point& p) const override{ return val; }

	HMMath::LinearPiecewise sz_boundary_segment(const HM2D::EdgeData& pos) const override{
		HMMath::LinearPiecewise ret;
		ret.add_point(0, val);
		ret.add_point(1, val);
		return ret;
	}
};

//size function defined on a triangle grid
struct TriBasedSizeFun: public SizeFun{
	shared_ptr<HM2D::GridData> grid;
	mutable shared_ptr<HMFem::Grid43::Approximator> approx;
	VertexData known_pts;
	vector<double> known_sz;
	vector<double> common_size;

	double minstep() const override{
		return *std::min_element(known_sz.begin(), known_sz.end());
	}
	double maxstep() const override{
		return *std::max_element(known_sz.begin(), known_sz.end());
	}

	//detached contours of source will be ignored
	void set_zone(const Contour::Tree& source,
			const vector<HM2D::EdgeData>& constraints, double gridh){
		grid = std::make_shared<HM2D::GridData>(
			HMFem::AuxGrid3(source, constraints, 2000, 20000, gridh));
		approx.reset(new HMFem::Grid43::Approximator(grid.get()));
	}

	void compute(){
		HMFem::LaplaceProblem lprob(*grid);
		std::map<const Vertex*, double> mp;
		std::vector<int> bp(known_pts.size());
		aa::enumerate_ids_pvec(grid->vvert);
		for (int i=0; i<known_pts.size(); ++i) {
			mp[known_pts[i].get()] = known_sz[i];
			bp[i] = known_pts[i]->id;
		}
		lprob.SetDirichlet(bp, [mp](const Vertex* v)->double{ return mp.find(v)->second; });
		common_size.resize(grid->vvert.size(), 0);
		lprob.Solve(common_size);

		/*
		//resort points so that all known points was at the end of point list
		aa::constant_ids_pvec(grid->vvert, -1);
		aa::enumerate_ids_pvec(known_pts);
		std::sort(grid->vvert.begin(), grid->vvert.end(),
			[](const shared_ptr<Vertex>& p1, const shared_ptr<Vertex>& p2){
				return p1->id < p2->id;
			});
		// create gradient matrix
		auto dxmat = HMFem::Assemble::DDx(*grid);
		auto dymat = HMFem::Assemble::DDy(*grid);
		dxmat->data.resize(2*grid->vvert.size());
		std::copy(dymat->data.begin(), dymat->data.end(),
				dxmat->data.begin() + grid->vvert.size());
		auto& gradmat = *dxmat;
		// create rhs
		vector<double> rhs(gradmat.data.size(), 0);
		// apply known values
		int start_known = grid->vvert.size() - known_sz.size();
		for (int i=0; i<gradmat.data.size(); ++i){
			auto jt = gradmat.data[i].lower_bound(start_known);
			while (jt != gradmat.data[i].end()){
				rhs[i] -= jt->second * known_sz[jt->first-start_known];
				gradmat.data[i].erase(jt++);
			}
		}

		// QR solution
		common_size.reserve(grid->vvert.size());
		common_size.resize(start_known, 1e8);
		HMMath::SuiteSparseQRSolver slv(gradmat, common_size.size());
		slv.Solve(rhs, common_size);
		//add known values
		common_size.resize(grid->vvert.size());
		std::copy(known_sz.begin(), known_sz.end(), common_size.begin()+start_known);
		*/
		//#############################
		HM2D::Export::GridVDataVTK(*grid, common_size, "g2.vtk");
	}

	double sz(const Point& p) const override{
		if (!approx){
			approx.reset(new HMFem::Grid43::Approximator(grid.get()));
		}
		return approx->Val(p, common_size);
	}
	double sz_proj(const Point& p) const override{
		try {
			return sz(p);
		} catch (HMFem::Grid43::Approximator::EOutOfArea){
			return sz(HM2D::Finder::ClosestEPoint(
				ECol::Assembler::GridBoundary(*grid), p));
		}
	}

	HMMath::LinearPiecewise sz_boundary_segment(const HM2D::EdgeData& site) const override{
		HMMath::LinearPiecewise ret;
		HM2D::EdgeData gedges = cont_along(site);
		aa::constant_ids_pvec(HM2D::AllVertices(gedges), -1);
		aa::enumerate_ids_pvec(grid->vvert);

		double lencoord = 0;
		Point* plast = 0;
		for (auto& p: HM2D::Contour::OrderedPoints(gedges)){
			if (plast != nullptr) lencoord += Point::dist(*p, *plast);
			plast = p.get();
			double val = (p->id == -1) ? sz(*p) : common_size[p->id];
			ret.add_point(lencoord, val);
		}
		return ret;
	}
	
	HM2D::EdgeData cont_along(const HM2D::EdgeData& site) const{
		shared_ptr<Vertex> before1, after2;
		shared_ptr<Vertex> v1, v2;

		auto treat_end =[&](Point pfirst, shared_ptr<Vertex>& v1, shared_ptr<Vertex>& before1){
			auto pos1 = approx->FindNodePos(pfirst);
			double ksi;
			switch (pos1.pos){
				case HMFem::Grid43::Approximator::NodePos::InternalVertex:
				case HMFem::Grid43::Approximator::NodePos::BndVertex:
					v1 = grid->vvert[pos1.nvert];
					return;
				case HMFem::Grid43::Approximator::NodePos::BndEdge:
				case HMFem::Grid43::Approximator::NodePos::InternalEdge:
				{
					before1.reset(new HM2D::Vertex(pfirst));
					auto cand1 = grid->vvert[pos1.nve1];
					auto cand2 = grid->vvert[pos1.nve2];
					for (auto it: site){
						if (isOnSection(*cand1, *it->pfirst(), *it->plast(), ksi)){
							v1 = cand1; return;
						}
						if (isOnSection(*cand2, *it->pfirst(), *it->plast(), ksi)){
							v1 = cand2; return;
						}
					}
				}
			}
			throw _PointNotOnGridEdge();
		};


		treat_end(*HM2D::Contour::First(site), v1, before1);
		if (HM2D::Contour::IsClosed(site)){ after2 = before1; v2 = v1; }
		else{ treat_end(*HM2D::Contour::Last(site), v2, after2); }

		HM2D::EdgeData ret = AssembleContourAlong(*grid, site, v1, v2);
		if (ret.size() == 0 && before1!=nullptr && after2 != nullptr){
			ret.emplace_back(new Edge(before1, after2));
		} if (before1 != nullptr){
			HM2D::Contour::Algos::AddFirstPoint(ret, before1);
		} else if (after2 != nullptr){
			HM2D::Contour::Algos::AddFirstPoint(ret, after2);
		}
		return ret;
	}

	//v1, v2 should belong to grid.
	//src should lie strictly along grid edges
	static HM2D::EdgeData AssembleContourAlong(const GridData& grid, const EdgeData& src,
			shared_ptr<Vertex> v1, shared_ptr<Vertex> v2){
		HM2D::EdgeData ret;
		//1. get source contour using grid primitives
		auto ve = HM2D::Connectivity::VertexEdge(grid.vedges, grid.vvert);
		aa::enumerate_ids_pvec(grid.vvert);
		int veit = v1->id;
		auto bb = BBox(src);
		BoundingBoxFinder bbfinder(bb, bb.maxlen()/20.);
		for (auto& e: src) bbfinder.addentry(BoundingBox(*e->pfirst(), *e->plast()));
		vector<bool> usededges(grid.vedges.size(), false);

		double ksi;
		int ntry = 0;
		while(1){
			for (auto& ei: ve[veit].eind) if (!usededges[ei]){
				usededges[ei] = true;
				Point cnt = grid.vedges[ei]->center();
				for (auto& ei2: bbfinder.suspects(cnt)){
					if (isOnSection(cnt, *src[ei2]->pfirst(), *src[ei2]->plast(), ksi)){
						ret.push_back(grid.vedges[ei]);
						veit = ret.back()->sibling(ve[veit].v.get())->id;
						goto END_WHILE;
					}
				}
			}
			// we must reach v2 anyway. If we are here then we've failed. 
			// Try one more time in the opposite direction and then raise an exception
			if (ntry++ != 0) throw _PointNotOnGridEdge();
			ret.clear();
			veit = v1->id;
			continue;
	END_WHILE:
			if (veit == v1->id || veit == v2->id) break;
		}
		//2. sort gridedges according to src
		bool dorevert = false;
		if (HM2D::Contour::IsClosed(ret) && HM2D::Contour::Area(ret)*HM2D::Contour::Area(src)<0){
			dorevert = true;
		}
		if (dorevert){
			HM2D::Contour::R::ReallyRevert::Permanent(ret);
		} else {
			HM2D::Contour::R::ReallyDirect::Permanent(ret);
		}
		return ret;
	}
};

void extract_vital_id1(const Contour::Tree& source, EdgeData& vital_edges, VertexData& vital_vertices){
	for (auto n: source.nodes)
	for (auto e: n->contour) if (e->id == 1) {
		vital_edges.push_back(e);
		//we don't need to explicitly define
		//edge end points as needed
		e->pfirst()->id = 0;
		e->plast()->id = 0;
	}
	for (auto n: source.nodes)
	for (auto e: n->contour) if (e->id != 1) {
		if (e->pfirst()->id == 1) {
			vital_vertices.push_back(e->vertices[0]);
			e->pfirst()->id = 0;
		}
		if (e->plast()->id == 1) {
			vital_vertices.push_back(e->vertices[1]);
			e->plast()->id = 0;
		}
	}
}

VertexData get_source_vertices(const EdgeData& cnt){
	VertexData ret = HM2D::Contour::OrderedPoints(cnt);
	// move end points half edge backward because
	// their size is defined on the edge middle point.
	if (ret[0]!=ret.back()){
		ret[0] = std::make_shared<Vertex>((*ret[0] + *ret[1])/2.0);
		if (ret.size()>2){
			ret.back() = std::make_shared<Vertex>((*ret.back() + *ret[ret.size()-2])/2.0);
		} else {
			ret.resize(1);
		}
	}
	return ret;
}

int set_first_last(HM2D::GridData& grid, HMFem::Grid43::Approximator& approx,
		const VertexData& vrt, shared_ptr<Vertex>& first, shared_ptr<Vertex>& last){
	HMFem::Grid43::Approximator::NodePos firstpos;
	first = HMFem::Grid43::GuaranteePoint(grid, *vrt[0], approx, firstpos);
	if (firstpos.pos == firstpos.Out) return -1;
	last  = (vrt[0]==vrt.back()) ? first
	                             : HMFem::Grid43::GuaranteePoint(grid, *vrt.back(), approx);

	if (vrt.size() == 1){
		if (firstpos.pos==firstpos.Internal ||
			firstpos.pos == firstpos.InternalEdge ||
			firstpos.pos == firstpos.Internal) return 0;
		else return 1;
	} else if (firstpos.pos == firstpos.BndVertex || firstpos.pos == firstpos.BndEdge){
		return (first == last) ? 2 : 3;
	} else {
		return 4;
	}
}

shared_ptr<Vertex> set_internal_point(HM2D::GridData& grid, HMFem::Grid43::Approximator& approx, Point p){
	return HMFem::Grid43::GuaranteePoint(grid, p, approx);
}

void calculate_sizes_open(VertexData& opts, vector<double>& osz,
		const EdgeData& src, const EdgeData& gcont){
	//src and gcont are already really directed in the same direction
	HMMath::LinearPiecewise lenfun;
	//fill lenfun
	std::vector<double> cntlens = HM2D::Contour::ELengths(src);
	//!!! src is half section longer in both sides then gcont
	std::vector<double> cntfuns = cntlens;
	cntlens[0]/=2.; cntlens.back()/=2.;
	double cntlen = std::accumulate(cntlens.begin(), cntlens.end(), 0.0);
	lenfun.add_point(0, cntfuns[0]);
	lenfun.add_point(1, cntfuns.back());
	double partlen = cntlens[0];
	for (int i=1; i<src.size(); ++i){
		lenfun.add_point(partlen/cntlen, (cntfuns[i-1] + cntfuns[i])/2.0);
		partlen += cntlens[i];
	}
	//fill return values using lenfun
	std::vector<double> griddist(gcont.size()+1, 0);
	for (int i=1; i<gcont.size()+1; ++i)
		griddist[i] = griddist[i-1] + gcont[i-1]->length();
	for (auto& v: griddist) v/=griddist.back();
	for (int i=0; i<gcont.size(); ++i){
		opts.push_back(gcont[i]->first());
		osz.push_back(lenfun(griddist[i]));
	}
	opts.push_back(gcont.back()->last());
	osz.push_back(lenfun(griddist.back()));
}

void calculate_sizes_closed(VertexData& opts, vector<double>& osz,
		const EdgeData& src, const EdgeData& gcont){
	HMMath::LinearPiecewise lenfun;
	//fill lenfun
	std::vector<double> cntlens = HM2D::Contour::ELengths(src);
	double cntlen = std::accumulate(cntlens.begin(), cntlens.end(), 0.0);
	lenfun.add_point(0, (cntlens[0]+cntlens.back())/2.0);
	lenfun.add_point(1, (cntlens[0]+cntlens.back())/2.0);
	double partlen = cntlens[0];
	for (int i=1; i<src.size(); ++i){
		lenfun.add_point(partlen/cntlen, (cntlens[i-1] + cntlens[i])/2.0);
		partlen += cntlens[i];
	}
	//fill return values using lenfun
	std::vector<double> griddist(gcont.size()+1, 0);
	for (int i=1; i<gcont.size()+1; ++i)
		griddist[i] = griddist[i-1] + gcont[i-1]->length();
	for (auto& v: griddist) v/=griddist.back();
	for (int i=0; i<gcont.size(); ++i){
		opts.push_back(gcont[i]->first());
		osz.push_back(lenfun(griddist[i]));
	}
}



void treat_internal_point(VertexData& opts, vector<double>& osz, shared_ptr<Vertex>& v, double sz){
	opts.push_back(v);
	osz.push_back(sz);
}
void treat_boundary_point(VertexData& opts, vector<double>& osz, shared_ptr<Vertex>& v, double sz){
	opts.push_back(v);
	osz.push_back(sz);
}

void treat_closed_bnd_cont(VertexData& opts, vector<double>& osz,
		shared_ptr<Vertex> v1, const HM2D::EdgeData& src, const HM2D::EdgeData& gridbnd){
	//src is in counterclockwise direction
	HM2D::EdgeData gcont = HM2D::Contour::Assembler::Contour1(gridbnd, v1.get());
	shared_ptr<HM2D::Contour::R::Reverter> rev;
	if (HM2D::Contour::Area(gcont) < 0){
		rev.reset(new HM2D::Contour::R::ReallyRevert(gcont));
	} else {
		rev.reset(new HM2D::Contour::R::ReallyDirect(gcont));
	}
	calculate_sizes_closed(opts, osz, src, gcont);
}

void treat_open_bnd_cont(VertexData& opts, vector<double>& osz,
		shared_ptr<Vertex> v1, shared_ptr<Vertex> v2,
		const EdgeData& src, const EdgeData& gridbnd){
	//src is a part of a closed contour in counterclockwise direction
	HM2D::EdgeData g1 = HM2D::Contour::Assembler::Contour1(gridbnd, v1.get());
	shared_ptr<HM2D::Contour::R::Reverter> rev;
	if (HM2D::Contour::Area(g1) > 0){
		rev.reset(new HM2D::Contour::R::ReallyDirect(g1));
	} else {
		rev.reset(new HM2D::Contour::R::ReallyRevert(g1));
	}
	HM2D::EdgeData gcont = HM2D::Contour::Assembler::ShrinkContour(g1, v1.get(), v2.get());
	calculate_sizes_open(opts, osz, src, gcont);
}


void treat_internal_cont(VertexData& opts, vector<double>& osz,
		shared_ptr<Vertex> v1, shared_ptr<Vertex> v2,
		const EdgeData& src, const GridData& grid){
	//1. get source contour using grid primitives
	HM2D::EdgeData gridedges = TriBasedSizeFun::AssembleContourAlong(grid, src, v1, v2);
	//2. fill opts, sz
	if (HM2D::Contour::IsClosed(gridedges)){
		calculate_sizes_closed(opts, osz, src, gridedges);
	} else {
		calculate_sizes_open(opts, osz, src, gridedges);
	}
}

void remove_bound_cross_edges(EdgeData& vital_detached, const Contour::Tree& source){
	if (vital_detached.size() == 0) return;
	EdgeData ae = vital_detached;
	for (auto& e: source.alledges()) ae.push_back(e);
	BoundingBox bb = BBox(AllVertices(ae));
	BoundingBoxFinder bbfinder(bb, bb.maxlen()/20);
	for (auto& e: ae) bbfinder.addentry(BoundingBox(*e->pfirst(), *e->plast()));
	std::vector<bool> to_delete(vital_detached.size(), false);
	double ksieta[2];
	for (int i=0; i<to_delete.size(); ++i){
		auto e = ae[i];
		for (int j: bbfinder.suspects(BoundingBox(*e->pfirst(), *e->plast()))){
			auto e2 = ae[j];
			if (j<=i || e->connected_to(*e2)) continue;
			SectCross(*e->pfirst(), *e->plast(), *e2->pfirst(), *e2->plast(), ksieta);
			if (ISIN_NN(ksieta[0], 0, 1)) {
				to_delete[i] = true;
				break;
			}
		}
	}
	aa::constant_ids_pvec(vital_detached, 0);
	for (int i=0; i<to_delete.size(); ++i)
		if (to_delete[i]) vital_detached[i]->id=1;
	aa::keep_by_id(vital_detached, 0);
}

void remove_outer_edges(EdgeData& vital_detached, const Contour::Tree& source){
	auto cc = HM2D::Contour::Assembler::SimpleContours(vital_detached);
	vector<Point> cv;
	for (auto& c: cc) cv.push_back(*c[0]->pfirst());
	auto srt = HM2D::Contour::Finder::SortOutPoints(source, cv);
	aa::constant_ids_pvec(vital_detached, 0);
	for (int i=0; i<srt.size(); ++i) if (srt[i] != INSIDE){
		aa::constant_ids_pvec(cc[i], 1);
	}
	aa::remove_by_id(vital_detached, 1);
}

shared_ptr<SizeFun> tri_based_size_function(const Contour::Tree& source, const EdgeData& cbound,
		const EdgeData& cdetached, vector<std::pair<Point, double>> psrc, double h){
	//!!! Here all detached and all psrc lies strictly within the source
	// initialize return value
	shared_ptr<SizeFun> ret = shared_ptr<SizeFun>(new TriBasedSizeFun());
	TriBasedSizeFun* sf = static_cast<TriBasedSizeFun*>(ret.get());
	
	// assemble connected groups of source data
	auto simpcont = HM2D::Contour::Assembler::SimpleContours(cbound);
	auto simpcont_det = HM2D::Contour::Assembler::SimpleContours(cdetached);
	for (auto& c: simpcont_det){
		if (c.size() == 1) psrc.emplace_back(c[0]->center(), c[0]->length());
		else simpcont.push_back(c);
	}

	// building aux grid
	sf->set_zone(source, simpcont_det, h);

	// force all simpcont be in counterclockwise direction
	vector<shared_ptr<HM2D::Contour::R::Clockwise>> rev;
	for (auto& n: source.bound_contours()){
		rev.emplace_back(new HM2D::Contour::R::Clockwise(n->contour, false));
	}
	for (auto& c: simpcont) if (c.size()>1){
		if (!HM2D::Contour::CorrectlyDirectedEdge(c, 0)){
			std::reverse(c.begin(), c.end());
		}
	}
	
	// extract source ordered vertices
	vector<VertexData> vvrt(simpcont.size());
	for (int icont=0; icont<simpcont.size(); ++icont)
		vvrt[icont] = get_source_vertices(simpcont[icont]);

	// first and last points of vital contours to aux grid
	// and types: 0 - internal point
	//            1 - bound point
	//            2 - closed bound contour
	//            3 - open bound contour
	//            4 - internal contour
	//            -1 - outside
	VertexData first_gp(simpcont.size()), last_gp(simpcont.size());
	vector<int> tps(simpcont.size());
	for (int icont=0; icont<simpcont.size(); ++icont){
		tps[icont] = set_first_last(*sf->grid, *sf->approx, vvrt[icont],
				first_gp[icont], last_gp[icont]);
	}

	// set internal points conditions at the end of sources list
	for (auto& p: psrc){
		auto gp = set_internal_point(*sf->grid, *sf->approx, p.first);
		first_gp.push_back(gp);
		last_gp.push_back(gp);
		tps.push_back(0);
	}

	//set size values for conditional sources
	int psrc_it = 0;
	auto gbnd = HM2D::ECol::Assembler::GridBoundary(*sf->grid);
	for (int i=0; i<tps.size(); ++i){
		switch (tps[i]){
		case 0: treat_internal_point(sf->known_pts, sf->known_sz,
				first_gp[i], psrc[psrc_it++].second);
			break;
		case 1: treat_boundary_point(sf->known_pts, sf->known_sz,
				first_gp[i], simpcont[i][0]->length());
			break;
		case 2: treat_closed_bnd_cont(sf->known_pts, sf->known_sz,
				first_gp[i], simpcont[i], gbnd);
			break;
		case 3: treat_open_bnd_cont(sf->known_pts, sf->known_sz,
				first_gp[i], last_gp[i], simpcont[i], gbnd);
			break;
		case 4: treat_internal_cont(sf->known_pts, sf->known_sz,
				first_gp[i], last_gp[i], simpcont[i], *sf->grid);
			break;
		}
	}
	//proceed with computation
	sf->compute();

	return ret;
}

shared_ptr<SizeFun> build_size_function1(const Contour::Tree& source, const EdgeData& vital_edges,
		const vector<std::pair<Point, double>>& psrc){
	//sort out vital edges and psrc
	//so that they lie inside source
	HM2D::EdgeData ae = source.alledges();
	aa::constant_ids_pvec(ae, 0);
	aa::constant_ids_pvec(vital_edges, 1);
	HM2D::EdgeData vital_detached;
	HM2D::EdgeData vital_bound;
	for (auto& n: source.bound_contours())
	for (auto& e: n->contour) if (e->id == 1) vital_bound.push_back(e);
	for (auto& n: source.detached_contours())
	for (auto& e: n->contour) if (e->id == 1) vital_detached.push_back(e);
	VertexData vital_detached_pts = HM2D::AllVertices(vital_detached);
	vector<Point> to_sort;
	for (auto& p: psrc) to_sort.push_back(p.first);
	for (auto& p: vital_detached_pts) to_sort.push_back(*p);
	auto sorted = HM2D::Contour::Finder::SortOutPoints(source, to_sort);
	vector<std::pair<Point, double>> psrc2;
	for (int i=0; i<psrc.size(); ++i)
		if (sorted[i] == INSIDE) psrc2.push_back(psrc[i]);
	for (int i=psrc.size(); i<sorted.size(); ++i)
		vital_detached_pts[i-psrc.size()]->id = (sorted[i] != OUTSIDE) ? 1 : 0;
	for (auto& e: vital_detached)
		e->id =  (e->pfirst()->id == 1 && e->plast()->id == 1) ? 1 : 0;
	aa::keep_by_id(vital_detached, 1);
	remove_bound_cross_edges(vital_detached, source);
	remove_outer_edges(vital_detached, source);

	//if there are no size constraints in input data use source bound
	//segmentation as a size source.
	int n = (vital_detached.size() + vital_bound.size() + psrc2.size());
	if (n == 0){
		for (auto& n: source.bound_contours())
		for (auto& e: n->contour) vital_bound.push_back(e);
		n = vital_bound.size();
	}

	//get minimum and maximum steps
	double minval=1e100, maxval=-1, aver=0;
	auto found_a_range = [&minval, &maxval, &aver](double len){
		if (len < minval) minval = len;
		if (len > maxval) maxval = len;
		aver += 1/len;
	};
	for (auto& e: vital_detached) found_a_range(e->length());
	for (auto& e: vital_bound) found_a_range(e->length());
	for (auto& p: psrc2) found_a_range(p.second);
	aver = n/aver;

	//check if we can use constant size function
	if (maxval/minval < 1.1)
		return shared_ptr<SizeFun>(new ConstSizeFun(aver));
	else{
		return tri_based_size_function(source, vital_bound, vital_detached, psrc2, aver);
	}
}

shared_ptr<SizeFun> build_size_function(const Contour::Tree& source, const EdgeData& vital_edges,
		const vector<std::pair<Point, double>>& psrc){
	auto ae = source.alledges();
	//1) divide source
	vector<Contour::Tree> subtrees = Contour::Tree::CropLevel01(source);

	//2) assemble answer
	if (subtrees.size() < 2){
		return build_size_function1(source, vital_edges, psrc);
	} else {
		vector<shared_ptr<SizeFun>> ret(subtrees.size());
		for (int i=0; i<subtrees.size(); ++i){
			for (auto& n: source.detached_contours()) subtrees[i].nodes.push_back(n);
			ret[i] = build_size_function1(subtrees[i], vital_edges, psrc);
		}
		return std::shared_ptr<SizeFun>(new SumSizeFun(ret, subtrees));
	}
}

void check_edges_vitality(shared_ptr<SizeFun> sfun, 
		const HM2D::EdgeData& src, const HM2D::EdgeData& reparted,
		const HM2D::VertexData& vv,
		vector<HM2D::EdgeData>& ret_src, vector<HM2D::EdgeData>& ret_to){
	shared_ptr<HM2D::Contour::R::Reverter> rev;
	if (HM2D::Contour::IsClosed(src)){
		rev.reset(new  HM2D::Contour::R::ForceFirst(src, *HM2D::Contour::First(reparted)));
	}
	aa::constant_ids_pvec(HM2D::AllVertices(reparted), 0);
	aa::constant_ids_pvec(vv, 1);
	vector<int> new_vital;
	for (int i=0; i<reparted.size(); ++i){
		if (reparted[i]->pfirst()->id == 1 && reparted[i]->plast()->id == 1){
			double true_len = reparted[i]->length();
			double sfun_len = sfun->sz_proj(reparted[i]->center());
			if (sfun_len / true_len >= 2.) new_vital.push_back(i);
		}
	}
	if (new_vital.size() == 0) return; 
	auto op1 = HM2D::Contour::OrderedPoints(src);
	auto op2 = HM2D::Contour::OrderedPoints(reparted);
	auto itop1 = op1.begin();
	for (int k: new_vital){
		auto p1 = op2[k], p2 = op2[k+1];
		while (itop1 != op1.end() && *itop1 != p1) ++itop1;
		assert(itop1 != op1.end());
		int ifirst = itop1 - op1.begin();
		while (itop1 != op1.end() && *itop1 != p2) ++itop1;
		assert(itop1 != op1.end());
		int isecond = itop1 - op1.begin();
		ret_to.push_back(HM2D::EdgeData {reparted[k]});
		ret_src.emplace_back();
		for (int i=ifirst; i<isecond; ++i) ret_src.back().push_back(src[i]);
	}
}

//first makes partition of the source using given sfun.
//then makes a check if edge length between vital_vertices (if any) satisfies sfun condition.
//If everything is fine builds weighted partition of source and returns true,
//else inserts edges between vital_vertices, adds them to vital_edge and returns false
//so that sfun can be rebuilt using new source and new vital_edge list.
bool boundary_repart(const shared_ptr<SizeFun>& sfun, Contour::Tree& source,
		EdgeData& vital_edges, const VertexData& vital_vertices){
	EdgeData partlist = source.alledges();
	//get list of contours for repartition
	aa::constant_ids_pvec(partlist, 0);
	aa::constant_ids_pvec(vital_edges, 1);
	aa::keep_by_id(partlist, 0);
	vector<HM2D::EdgeData> conts;
	for (auto& n: source.nodes){
		auto c2 = HM2D::Contour::Assembler::SimpleContours(
				aa::copy_by_id(n->contour, 0));
		conts.insert(conts.end(), c2.begin(), c2.end());
	}
	//auto conts = HM2D::Contour::Assembler::SimpleContours(partlist);
	//assign conts to source contours
	for (int i=0; i<source.nodes.size(); ++i) aa::constant_ids_pvec(source.nodes[i]->contour, i);
	vector<int> iconts;
	for (auto& cc: conts) iconts.push_back(cc[0]->id);
	//map vital_vertices to those contours
	aa::constant_ids_pvec(HM2D::AllVertices(partlist), 0);
	aa::constant_ids_pvec(vital_vertices, 1);
	vector<VertexData> vv(conts.size());
	for (int i=0; i<conts.size(); ++i)
	for (auto& e: conts[i]){
		if (e->pfirst()->id == 1){
			vv[i].push_back(e->first());
			e->pfirst()->id = 0;
		}
		if (e->plast()->id == 1){
			vv[i].push_back(e->last());
			e->plast()->id = 0;
		}
	}
	vector<EdgeData> conts2(conts.size());
	//repart each contour
	for (int i=0; i<conts.size(); ++i){
		std::map<double, double> bmap;
		try {
			//first try boundary procedures
			HMMath::LinearPiecewise bfun = sfun->sz_boundary_segment(conts[i]);
			std::map<double, double> bmap1 = bfun.to_map();
			double v0 = bmap1.begin()->first;
			double nrm = bmap1.rbegin()->first - v0;
			for (auto& n: bmap1){
				bmap.emplace((n.first-v0)/nrm, n.second);
			}
		} catch (_PointNotOnGridEdge){
			//if conts[i] doesn't lie along the set of grid edges the
			//do simple linear interpolation
			double len = HM2D::Contour::Length(conts[i]);
			double step = sfun->minstep();
			int N = std::min(50, int(len/step));
			vector<double> w(N+1);
			for (int i=0; i<N+1; ++i) w[i] = (double)i/(N+1);
			vector<Point> pp = HM2D::Contour::WeightPoints(conts[i], w);
			vector<double> svals = sfun->sz(pp);
			for (int i=0; i<N+1; ++i){
				if (svals[i]>0) bmap.emplace(w[i], svals[i]);
			}
		}

		if (bmap.size() > 0) conts2[i] = HM2D::Contour::Algos::WeightedPartition(bmap, conts[i], vv[i]);
		else HM2D::DeepCopy(conts[i], conts2[i], 0);
	}

	//make a check for resultin edges bounded by vital vertices.
	//from, to is a list of conts/conts2 respective edges so that
	//each of conts2[i] contains single edge bounded by vital vertices.
	bool ret = true;
	vector<HM2D::EdgeData> from, to;
	vector<int> new_iconts;
	for (int i=0; i<conts.size(); ++i) if (vv[i].size() > 1){
		int szold = from.size();
		check_edges_vitality(sfun, conts[i], conts2[i], vv[i], from, to);
		for (int j=szold; j<from.size(); ++j) new_iconts.push_back(iconts[i]);
	}
	if (from.size() > 0){
		for (auto& it: to) vital_edges.push_back(it[0]);
		std::swap(from, conts);
		std::swap(to, conts2);
		std::swap(iconts, new_iconts);
		ret = false;
	}
	
	//substitute conts2 to source
	aa::constant_ids_pvec(source.alledges(), 0);
	for (auto& ic: conts) for (auto& ie: ic) ie->id = 1;
	for (int i=0; i<source.nodes.size(); ++i){
		aa::remove_by_id(source.nodes[i]->contour, 1);
	}
	for (int i=0; i<conts.size(); ++i){
		auto& cnt = source.nodes[iconts[i]]->contour;
		cnt.insert(cnt.end(), conts2[i].begin(), conts2[i].end());
	}
	for (auto& n: source.nodes)
		n->contour = HM2D::Contour::Assembler::Contour1(n->contour);
	return ret;
}

}

shared_ptr<SizeFun> Grid::BuildSizeFunction(const Contour::Tree& source,
		const vector<std::pair<Point, double>>& psrc){
	EdgeData vital_edges;
	for (auto& n: source.nodes)
	for (auto& e: n->contour) if (e->id == 1) vital_edges.push_back(e);
	return build_size_function(source, vital_edges, psrc);
}

shared_ptr<SizeFun> Grid::ApplySizeFunction(Contour::Tree& source,
		const vector<std::pair<Point, double>>& psrc, bool force_sizefun){
	// extract primitives which should be saved
	EdgeData vital_edges;
	VertexData vital_vertices;
	extract_vital_id1(source, vital_edges, vital_vertices);
	bool need_repart = vital_edges.size() < source.alledges().size();

	// if no need for repartition
	if (!force_sizefun && !need_repart) return shared_ptr<SizeFun>();

	// construct size function
	auto sfun = build_size_function(source, vital_edges, psrc);

	// do repartition
	if (need_repart) {
		int tries = 0;
		while (tries < 10 && !boundary_repart(sfun, source, vital_edges, vital_vertices)){
			sfun = build_size_function(source, vital_edges, psrc);
			++tries;
		}
		if (tries >= 10) throw std::runtime_error(
			"Failed to satisfy 'keep vertex' condition "
			"while building size function");
	}

	return sfun;
}
