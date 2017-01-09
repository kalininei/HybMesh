#include "wireframegrid.h"
#include "addalgo.hpp"
#include <algorithm>
#include <numeric>
#include "assert.h"
#include "intrusion.h"
#include "procgrid.h"

//constructors
PtsGraph::PtsGraph(const GridGeom& g2){
	//nodes
	nodes.reserve(g2.n_points());
	for (int i=0; i<g2.n_points(); ++i){
		nodes.push_back(Point(*g2.get_point(i)));
	}
	//edges
	auto eds = g2.get_edges();
	lines.reserve(eds.size());
	for (auto e: eds){
		lines.push_back(GraphLine(e.p1, e.p2));
	}
}

PtsGraph::PtsGraph(const ContoursCollection& cc){
	for (auto c: cc.contours_list()){
		for (int i=0; i<c.n_points(); ++i){
			nodes.push_back(Point(*c.get_point(i)));
			if (i<c.n_points()-1){
				lines.push_back(GraphLine(nodes.size()-1, nodes.size()));
			} else {
				lines.push_back(GraphLine(nodes.size()-1, nodes.size()-c.n_points()));
			}
		}
	}
}

PtsGraph::PtsGraph(const PContour& c){
	for (int i=0; i<c.n_points(); ++i){
		nodes.push_back(Point(*c.get_point(i)));
		if (i<c.n_points()-1){
			lines.push_back(GraphLine(nodes.size()-1, nodes.size()));
		} else {
			lines.push_back(GraphLine(nodes.size()-1, nodes.size()-c.n_points()));
		}
	}
}

PtsGraph::PtsGraph(const HM2D::EdgeData& cc){
	auto av = HM2D::AllVertices(cc);
	aa::enumerate_ids_pvec(av);
	for (auto p: av) nodes.push_back(*p);
	for (auto e: cc) lines.push_back(GraphLine(e->first()->id, e->last()->id));
}

auto PtsGraph::_impose_impl(const PtsGraph& main_graph, const PtsGraph& imp_graph, double eps)
		-> impResT {
	auto ret = impResT(main_graph, vector<int>(), vector<double>());
	auto& G=std::get<0>(ret);
	auto& ig_nodes=std::get<1>(ret);
	auto& cross_nodes=std::get<2>(ret);

	// === Acceleration initialization
	////this guarantees that points in G will not change their addresses
	////on push_back invocation. This is important due to ndfinder implementation details.
	//G.nodes.reserve(main_graph.Nnodes()+imp_graph.Nnodes());
	//find rectangle which contains all nodes
	Point top1 = Point::GetTop(main_graph.nodes.begin(), main_graph.nodes.end());
	Point bot1 = Point::GetBot(main_graph.nodes.begin(), main_graph.nodes.end());
	Point top2 = Point::GetTop(imp_graph.nodes.begin(), imp_graph.nodes.end());
	Point bot2 = Point::GetBot(imp_graph.nodes.begin(), imp_graph.nodes.end());
	Point maxp = Point(std::max(top1.x, top2.x), std::max(top1.y, top2.y));
	Point minp = Point(std::min(bot1.x, bot2.x), std::min(bot1.y, bot2.y));
	PtsGraphAccel Accel(G, minp, maxp, eps);

	vector<int> crnodes_ind;
	//add points from g; fill ig_nodes array
	ig_nodes.reserve(imp_graph.Nnodes());
	for (auto& p: imp_graph.nodes){
		auto ares = Accel.add_point(p);
		ig_nodes.push_back(std::get<1>(ares));
		if (std::get<0>(ares)==2) crnodes_ind.push_back(G.nodes.size()-1);
	}

	//all nodes which were added after last_node are cross nodes
	int last_node=G.Nnodes();
	//add connections from g
	for (auto& line: imp_graph.lines){
		int i0=ig_nodes[line.i0];
		int i1=ig_nodes[line.i1];
		Accel.add_connection(i0, i1);
	}
	for (int i=last_node; i<G.Nnodes(); ++i) crnodes_ind.push_back(i);

	//calculate cross_nodes coordinates
	cross_nodes=vector<double>(G.Nnodes(), -1);
	const PtsGraphAccel AccelMainG(const_cast<PtsGraph&>(main_graph),eps);
	for (int i: crnodes_ind){
		cross_nodes[i] = AccelMainG.find_gline(G.nodes[i]);
	}
	return ret;
}

auto PtsGraph::impose(const PtsGraph& main_graph, const PtsGraph& imp_graph, double eps) -> impResT{
	//TODO: parallel implementation
	return _impose_impl(main_graph, imp_graph, eps);
}

//PtsGrapth::togrid specific routines and classes
namespace{
struct tgPoint;
struct tgHalfEdge;
struct tgCell;

struct geps_lower{
	bool operator()(double a, double b){ return ISLOWER(a,b); }
};

struct tgPoint: public Point{
	tgPoint(const Point& p, int oind): Point(p), origindex(oind){}
	std::map<double, tgHalfEdge*, geps_lower> he;
	int NumEd() const {return he.size(); }
	tgHalfEdge* popfirst(){
		tgHalfEdge* ret=he.begin()->second;
		he.erase(he.begin());
		return ret;
	}
	tgHalfEdge* popnext(double angle){
		assert(he.size() > 0);
		angle=AngleAdd(angle, M_PI, geps); //revert angle of incoming edge
		auto it=he.lower_bound(angle);
		if (it==he.begin()) it=he.end();
		auto ret=(--it)->second;
		he.erase(it);
		return ret;
	}
	int origindex;
};
struct tgHalfEdge{
	tgPoint *first, *last;
	double angle;
	bool direction;
	int origindex;
};
struct tgCell{
	std::list<tgHalfEdge*> he;
	bool is_close() const { return (*he.begin())->first==(*he.rbegin())->last; }
	void add(tgHalfEdge* h) { he.push_back(h); }
	tgHalfEdge* lasthe() const{ return *he.rbegin(); };
	double Area() const{
		PContour cnt;
		for (auto e: he) cnt.add_point(e->first);
		return cnt.area();
	}
};

void build_tg(const vector<Point>& pts, const vector<GraphLine>& lines,
		vector<tgPoint>& P, vector<tgHalfEdge>& HE
){
	//Points
	P.reserve(pts.size());
	for (size_t i=0; i<pts.size(); ++i) P.push_back(tgPoint(pts[i], i));
	//HalfEdges
	HE.reserve(2*lines.size());
	int ind=0;
	for (auto& line: lines){
		tgPoint* p0=&P[line.i0];
		tgPoint* p1=&P[line.i1];
		double angle=ToAngle(atan2(p1->y-p0->y, p1->x-p0->x), geps);
		HE.push_back(tgHalfEdge()); auto& v1=HE.back();
		HE.push_back(tgHalfEdge()); auto& v2=HE.back();
		//direct
		v1.first=p0; v1.last=p1;
		v1.origindex=ind; v1.direction=true;
		v1.angle=angle;
		//backward
		v2.first=p1; v2.last=p0;
		v2.origindex=ind; v2.direction=false;
		v2.angle=AngleAdd(angle, M_PI, geps);
		//points connectivity
		p0->he[v1.angle]=&v1;
		p1->he[v2.angle]=&v2;
		++ind;
	}
}

GridGeom formgrid(const vector<tgPoint>& P, const vector<tgHalfEdge>& HE,
		const vector<tgCell*>& cell_inner, const vector<tgCell*>& cell_outer)
{
	vector<Edge> ed2(HE.size()/2);
	//edges (points connectivity)
	for (auto& e: HE) if (e.direction){
		int ind1=e.first->origindex, ind2=e.last->origindex;
		ed2[e.origindex]=Edge(ind1, ind2);
	}
	//edges cells connectivity
	//inner cells
	int cell_index=0;
	for (auto& c: cell_inner){
		for (auto& he: c->he){
			if (he->direction) ed2[he->origindex].cell_left=cell_index;
			else ed2[he->origindex].cell_right=cell_index;
		}
		++cell_index;
	}
	//edges to cells->nodes
	vector<std::list<int>> cells_edges(cell_inner.size());
	for (auto& e: ed2){
		if (e.cell_left>=0){
			cells_edges[e.cell_left].push_back(e.p1);
			cells_edges[e.cell_left].push_back(e.p2);
		}
		if (e.cell_right>=0){
			cells_edges[e.cell_right].push_back(e.p2);
			cells_edges[e.cell_right].push_back(e.p1);
		}
	}
	vector<int> cells_nodes;
	for (auto& c2: cells_edges){
		cells_nodes.push_back(c2.size()/2);
		cells_nodes.push_back(c2.front()); c2.pop_front();
		int nextnode = c2.front(); c2.pop_front();
		while (c2.size()>0){
			bool noerror = false;
			for (auto it=c2.begin(); it!=c2.end(); std::advance(it,2)){
				if (*it==nextnode){
					cells_nodes.push_back(*it);
					auto it2 = it++;
					nextnode = *(it);
					c2.erase(it2, ++it);  //erases two entries
					noerror = true;
					break;
				}
			}
			if (!noerror) throw std::runtime_error("Failed to assemble a grid from the wireframe");
		}
	}
	vector<double> raw_pnt;
	for (auto& p: P){ raw_pnt.push_back(p.x); raw_pnt.push_back(p.y); }
	return GridGeom(raw_pnt.size()/2, cell_inner.size(), &raw_pnt[0], &cells_nodes[0]);
}
}//namespace


GridGeom PtsGraph::togrid() const{
	//assemble subgraphs each of which can be
	//processed to a single connected grid
	auto subs = SubGraph::build(*this);

	//assemble single connected grids for each sub graph
	vector<GridGeom> grids; grids.reserve(subs.size());
	for (auto& s: subs) grids.push_back(s.togrid());

	//if only a single grid exists return it
	if (grids.size() == 1) return GridGeom::sum(grids);

	//assembling system of outer contours
	HM2D::Contour::Tree cont;
	for (int i=0; i<grids.size(); ++i){
		auto c1 = GGeom::Info::Contour(grids[i]);
		for (int j=0; j<c1.nodes.size(); ++j){
			cont.add_contour(c1.nodes[j]->contour);
		}
	}

	//if no inner contours make a summation
	if (cont.roots().size() == cont.nodes.size()) return GridGeom::sum(grids);

	//some grids lay within cells of parent grids. We need the intrusion algo.
	return PtsGraph::intrusion_algo(grids);
}

GridGeom PtsGraph::intrusion_algo(const vector<GridGeom>& g){
	//1 divide grids into subgrids
	vector<GridGeom> grids; grids.reserve(10*g.size());
	for (auto& it: g){
		auto sg = GGeom::Modify::SubGrids(it);
		if (sg.size() == 1){ 
			grids.push_back(std::move(sg[0]));
		} else {
			for (int i=0; i<sg.size(); ++i){
				GGeom::Modify::ReallocatePrimitives(sg[i]);
				grids.push_back(std::move(sg[i]));
			}
		}
	}
	//2 assemble tree for each singly connected grid
	HM2D::Contour::Tree concol;
	for (int i=0; i<grids.size(); ++i){
		auto c1 = GGeom::Info::Contour(grids[i]);
		for (int j=0; j<c1.nodes.size(); ++j){
			concol.add_contour(c1.nodes[j]->contour);
		}
	}
	
	//set of level, parent_index, child_index
	std::set<std::tuple<int, int, int>> lpc;
	for (int i=0; i<concol.nodes.size(); ++i){
		int level = concol.nodes[i]->level;
		if (level == 0) continue;
		auto parent_cont = concol.nodes[i]->parent;
		int parent = -1;
		for (int j=0; j<concol.nodes.size(); ++j){
			if (concol.nodes[j] == parent_cont.lock()){
				parent = j;
				break;
			}
		}
		assert(parent != -1);
		lpc.insert(std::make_tuple(level, parent, i));
	}
	//moving from the highest level of nesting
	int parent=-1;
	vector<int> childs;
	auto intrude = [&](){
		std::vector<GridGeom*> cg;
		for (auto i: childs) cg.push_back(&grids[i]);
		GridForIntrusion::intrude_grids(grids[parent], cg);
	};
	for (auto it=lpc.rbegin(); it!=lpc.rend(); ++it){
		auto itnext = it; ++itnext;
		parent = std::get<1>(*it);
		childs.push_back(std::get<2>(*it));
		if (itnext == lpc.rend() || std::get<1>(*itnext) != parent){
			intrude();
			parent=-1; childs.clear();
		}
	}

	//assemble first level grids
	vector<GridGeom*> lev0;
	for (int i=0; i<concol.nodes.size(); ++i) if (concol.nodes[i]->level == 0){
		lev0.push_back(&grids[i]);
	}
	return GridGeom::sum(lev0);
}

HM2D::EdgeData PtsGraph::toecollection() const{
	HM2D::EdgeData ret;
	HM2D::VertexData pcol;
	for (auto& p: nodes) pcol.push_back(std::make_shared<HM2D::Vertex>(p));

	for (int i=0; i<Nlines(); ++i){
		auto pr = get_line(i);
		ret.emplace_back(new HM2D::Edge(pcol[pr.first], pcol[pr.second]));
	}
	return ret;
}

vector<vector<int>> PtsGraph::lines_lines_tab() const{
	vector<vector<int>> pts_lines(nodes.size());
	for (int i=0; i<lines.size(); ++i){
		pts_lines[lines[i].i0].push_back(i);
		pts_lines[lines[i].i1].push_back(i);
	}
	vector<vector<int>> ret(lines.size());
	for (auto& pl: pts_lines){
		for (int i=0; i<pl.size(); ++i){
			for (int j=i+1; j<pl.size(); ++j){
				ret[pl[i]].push_back(pl[j]);
				ret[pl[j]].push_back(pl[i]);
			}
		}
	}
	return ret;
}

// =============================== SubGraph
SubGraph::SubGraph(const std::list<int>& used_lines, const PtsGraph* p):parent(p), lines(used_lines){}

vector<SubGraph> SubGraph::build(const PtsGraph& p){
	vector<SubGraph> ret;
	std::vector<int> lines_usage(p.lines.size(), 0);
	auto contab = p.lines_lines_tab();
	while (1){
		auto fnd0 = std::find(lines_usage.begin(), lines_usage.end(), 0);
		if (fnd0 == lines_usage.end()) break;
		std::list<int> used_lines;
		lines_sort_out(fnd0 - lines_usage.begin(), lines_usage, used_lines, contab);
		ret.push_back(SubGraph(used_lines, &p));
	};
	return ret;
}

//Recursive algorithm fails with stack overflow on debug mode in Windows (MinGW 64 at least).
//Don't think it's a big problem but it's better to avoid recursion here.
void SubGraph::lines_sort_out(int iline, std::vector<int>& lines_usage, std::list<int>& result,
			const vector<std::vector<int>>& contab){
	//add sibling lambda: pushes indexes to result if they are not already there
	auto add_siblings = [&](int par){
		for (auto j: contab[par]){
			if (lines_usage[j] == 0){
				lines_usage[j] = 1;
				result.push_back(j);
			}
		}
	};
	//1) add first line siblings
	result.push_back(iline); add_siblings(iline);
	//2) add all other siblings until its possible
	auto it = result.begin(); ++it;
	while (it != result.end()) { add_siblings(*it); ++it; }
}

GridGeom SubGraph::togrid() const{
	vector<tgPoint> P;
	vector<tgHalfEdge> HEdges;
	vector<tgCell> VC;
	
	auto nl = rebuild_nodes_lines();
	auto nodes = std::get<0>(nl);
	auto lines = std::get<1>(nl);

	//1) Build tgHalfEdge, tgPoint vectors.
	build_tg(nodes, lines, P, HEdges);

	//2) Fill VC vector
	for (auto& p: P){
		while (p.NumEd()!=0){
			VC.push_back(tgCell());
			tgCell& c=VC.back();
			c.add(p.popfirst());
			while (!c.is_close()){
				auto h=c.lasthe();
				auto hnext=h->last->popnext(h->angle);
				c.add(hnext);
			}
		}
	}

	//3) Sort inner and outer cells using sign of their areas.
	vector<tgCell*> inner_cells, outer_cells;
	for (auto& c: VC){
		if (c.Area()>0) inner_cells.push_back(&c);
		else outer_cells.push_back(&c);
	}
	//4) build Grid2D
	return formgrid(P, HEdges, inner_cells, outer_cells);
}

std::tuple<
	vector<Point>,
	vector<GraphLine>
> SubGraph::rebuild_nodes_lines() const{
	std::tuple<vector<Point>, vector<GraphLine>> ret;
	auto& pts = std::get<0>(ret);
	auto& lns = std::get<1>(ret);

	std::set<int> used_points;
	for (auto& oln: lines){
		used_points.insert(parent->lines[oln].i0);
		used_points.insert(parent->lines[oln].i1);
	}

	std::map<int, int> oldnew_p;
	for (auto up: used_points){
		oldnew_p[up] = pts.size();
		pts.push_back(parent->nodes[up]);
	}
	
	for (auto& nl: lines){
		lns.push_back(parent->lines[nl]);
		lns.back().i0 = oldnew_p[lns.back().i0];
		lns.back().i1 = oldnew_p[lns.back().i1];
	}

	return ret;
};

// =============================== PtsGraphAccel
PtsGraphAccel::PtsGraphAccel(PtsGraph& g, const Point& pmin, const Point& pmax, double e) noexcept:
		eps(e), epsPnt(e,e), ip(Tind2Proc(Naux, Naux)), G(&g), nodecomp(e), ndfinder(nodecomp)
{
	init(pmin, pmax);
}
PtsGraphAccel::PtsGraphAccel(PtsGraph& g,  double e) noexcept:
		eps(e), epsPnt(e,e), ip(Tind2Proc(Naux,Naux)), G(&g), nodecomp(e), ndfinder(nodecomp)
{
	Point pmax = Point::GetTop(G->nodes.begin(), G->nodes.end());
	Point pmin = Point::GetBot(G->nodes.begin(), G->nodes.end());
	init(pmin, pmax);
}
void PtsGraphAccel::init(Point pmin, Point pmax) noexcept{
	pmin-=epsPnt; pmax+=epsPnt;
	p0 = pmin;
	hx = (pmax.x - pmin.x)/ip.sizex();
	hy = (pmax.y - pmin.y)/ip.sizey();
	edmap.resize(ip.size()); ndmap.resize(ip.size());
	//1) fill points
	for (int i = 0; i<G->Nnodes(); ++i) add_node_to_map(i);
	//2) fill edges
	for (int i = 0; i<G->Nlines(); ++i) add_edge_to_map(i);
}

std::tuple<int, int>  PtsGraphAccel::add_point(const Point& p) noexcept{
	//find equal point
	auto fnd = ndfinder.find(p);
	if (fnd != ndfinder.end()){
		return std::make_tuple(1, fnd->second);
	}
	//find equal edge
	double ksi;
	for (auto i: candidates_edges(p)){
		if (isOnSection(p, G->ednode0(i), G->ednode1(i), ksi, eps)){
			//return std::make_tuple(2, break_graph_line(i, p));
			return std::make_tuple(2, break_graph_line(i, ksi));
		}
	}
	//add new node to Graph
	return std::make_tuple(0, add_graph_node(p));
}

int PtsGraphAccel::add_graph_node(const Point& p) noexcept{
	G->nodes.push_back(p);
	add_node_to_map(G->Nnodes()-1);
	return G->Nnodes()-1;
}

int PtsGraphAccel::break_graph_line(int iline, double ksi) noexcept{
	int reti = add_graph_node(Point::Weigh(G->ednode0(iline), G->ednode1(iline), ksi));
	add_graph_edge(reti, G->lines[iline].i1);
	reset_graph_edge(G->lines[iline].i0, reti, iline);
	return reti;
}

int PtsGraphAccel::break_graph_line(int iline, const Point& p) noexcept{
	int reti = add_graph_node(p);
	add_graph_edge(reti, G->lines[iline].i1);
	reset_graph_edge(G->lines[iline].i0, reti, iline);
	return reti;
}

void PtsGraphAccel::add_graph_edge(int i0, int i1) noexcept{
	G->lines.push_back(GraphLine(i0,i1));
	add_edge_to_map(G->Nlines()-1);
}

void PtsGraphAccel::reset_graph_edge(int i0, int i1, int iline) noexcept{
	delete_edge_from_map(iline);
	G->lines[iline].set(i0, i1);
	add_edge_to_map(iline);
}

void PtsGraphAccel::add_connection(int i0, int i1) noexcept{
	if (i1<i0) std::swap(i0,i1);
	else if (i0==i1) return;
	//find all points which lie on i0, i1 line
	std::map<double, int, geps_lower> linepts;
	linepts[0]=i0; linepts[1]=i1;
	double ksi;
	for (auto i: candidates_points(G->nodes[i0], G->nodes[i1])) {
		if (i0==i || i1==i) continue;
		if (isOnSection(G->nodes[i], G->nodes[i0], G->nodes[i1], ksi, eps)){
			linepts[ksi]=i;
		}
	}
	//add all segments
	auto it=linepts.begin();
	while (1){
		int p0ind=(it++)->second;
		if (it==linepts.end()) break;
		int p1ind=it->second;
		add_clear_connection(p0ind, p1ind);
	}
}

void PtsGraphAccel::add_clear_connection(int i0, int i1) noexcept{
	//find intersections
	Point &p0 = G->nodes[i0], &p1 = G->nodes[i1];
	auto ce = candidates_edges(G->nodes[i0], G->nodes[i1]);
	//check if edge already presents -> do nothing
	for (auto iline: ce){
		GraphLine& line = G->lines[iline];
		if (line.has_node(i0) && line.has_node(i1)) return;
	}
	//check for crosses
	double ksieta[2];
	for (auto iline: ce){
		GraphLine& line = G->lines[iline];
		//if edge starts with i0 or i1 point then it is not a candidate for intersection
		if (line.has_node(i0) || line.has_node(i1)) continue;
		//check for crosses
		if (SectCross(p0, p1, G->ednode0(iline), G->ednode1(iline), ksieta) ){
			int an = break_graph_line(iline, ksieta[1]);
			add_clear_connection(i0, an);
			add_clear_connection(an, i1);
			return;
		}
	}
	//clear conection
	add_graph_edge(i0, i1);
}

double PtsGraphAccel::find_gline(const Point& p) const noexcept{
	double ksi;
	for (auto ind: candidates_edges(p)){
		if (isOnSection(p, G->ednode0(ind), G->ednode1(ind), ksi, eps)) return ind+ksi;
	}
	return -1;
}

std::vector<Point> PtsGraph::center_line_points() const{
	std::vector<Point> ret; ret.reserve(lines.size());
	for (auto& ln: lines){
		ret.push_back( (nodes[ln.i0] + nodes[ln.i1])/2.0 );
	}
	return ret;
}

void PtsGraph::delete_unused_points(){
	std::vector<int> usage(nodes.size(), 0);
	for (auto& ln: lines){
		++usage[ln.i0];
		++usage[ln.i1];
	}
	std::set<int> bad_pts;
	for (size_t i = 0; i<usage.size(); ++i)
		if (usage[i] == 0) bad_pts.insert(i);

	//return if all points are used
	if (bad_pts.size() == 0) return;
	//old index -> new index dictionary
	std::vector<int> old_new(nodes.size());
	std::iota(old_new.begin(), old_new.end(), 0);
	for (auto it = bad_pts.rbegin(); it != bad_pts.rend(); ++it){
		std::for_each(old_new.begin() + *it, old_new.end(), [](int& i){ --i; });
	}
	//remove points
	aa::remove_entries(nodes, bad_pts);
	//renumber line entries
	for (auto& ln: lines){
		ln.i0 = old_new[ln.i0];
		ln.i1 = old_new[ln.i1];
	}
}

PtsGraph PtsGraph::cut(const PtsGraph& wmain, const ContoursCollection& conts, int dir){
	auto ccut = PtsGraph(conts);
	//get the imposition: wmain + conts
	auto ires = impose(wmain, ccut, geps);
	PtsGraph& pg = std::get<0>(ires);
	//filter lines which lie within bad region
	//using the position of line center point
	/* bad algorithm because of the short edges problem
	std::vector<Point> line_cnt = pg.center_line_points();
	auto flt = conts.filter_points_i(line_cnt);
	std::vector<int>& badi = (dir==INSIDE) ? std::get<0>(flt) : std::get<2>(flt);
	std::set<int> bad_lines(badi.begin(), badi.end());
	*/

	//TODO: need better algorithm
	std::vector<bool> is_on_conts(pg.nodes.size(), false);
	for (auto a: std::get<1>(ires)) is_on_conts[a] = true;
	for (int i=0; i<pg.nodes.size(); ++i) if (std::get<2>(ires)[i] >= 0.0) is_on_conts[i] = true;

	std::vector<Point> line_cnt = pg.center_line_points();
	auto flt = conts.filter_points_i(line_cnt);
	std::vector<int> ptypes(line_cnt.size(), 1); //1 - good, 2 - bad, 3 - bnd
	for(auto a: std::get<1>(flt)) ptypes[a] = 3;
	if (dir == INSIDE) for (auto a: std::get<0>(flt)) ptypes[a] = 2;
	else for (auto a: std::get<2>(flt)) ptypes[a] = 2;
	std::set<int> bad_lines;
	for (int i=0; i<line_cnt.size(); ++i){
		if (ptypes[i]==1) continue;
		if (ptypes[i]==2) {bad_lines.insert(i); continue;}
		int p1 = pg.lines[i].i0;
		int p2 = pg.lines[i].i1;
		int p1_type, p2_type;
		if (is_on_conts[p1]) p1_type = 3;
		else{
			int r1 = conts.is_inside(pg.nodes[p1]);
			if (r1 == INSIDE) p1_type =  (dir == INSIDE) ? 2 : 1;
			else if (r1 == OUTSIDE) p1_type = (dir == INSIDE) ? 1 : 2;
			else {p1_type = 3;}
		}
		if (is_on_conts[p2]) p2_type = 3;
		else{
			int r1 = conts.is_inside(pg.nodes[p2]);
			if (r1 == INSIDE) p2_type =  (dir == INSIDE) ? 2 : 1;
			else if (r1 == OUTSIDE) p2_type = (dir == INSIDE) ? 1 : 2;
			else {p2_type = 3;}
		}

		if (p1_type == 2 || p2_type == 2) bad_lines.insert(i);
	}
	
	aa::remove_entries(pg.lines, bad_lines);

	pg.delete_unused_points();

	return pg;
}

void PtsGraph::add_edges(const PContour& c){
	PtsGraph pgr(c);
	auto ires = std::get<0>(impose(*this, pgr, geps));
	nodes = ires.nodes;
	lines = ires.lines;
}

void PtsGraph::exclude_area(const ContoursCollection& cont, int dir){
	std::vector<Point> line_cnt = center_line_points();
	auto flt = cont.filter_points_i(line_cnt);
	std::vector<int>& badi = (dir==INSIDE) ? std::get<0>(flt) : std::get<2>(flt);
	std::set<int> bad_lines(badi.begin(), badi.end());
	exclude_lines(bad_lines);
};

void PtsGraph::exclude_lines(const std::set<int>& exlines){
	aa::remove_entries(lines, exlines);
	delete_unused_points();
}

PtsGraph PtsGraph::overlay(const PtsGraph& wmain, const PtsGraph& wsec){
	auto ret = impose(wmain, wsec, geps);
	return std::get<0>(ret);
}
