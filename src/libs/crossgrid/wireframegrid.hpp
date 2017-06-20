#ifndef WIREFRAMEGRID_CROSSGRID_H
#define WIREFRAMEGRID_CROSSGRID_H
//procedures with a grid built by connected edges without cells definition

#include <tuple>
#include <array>
#include <list>
#include <map>
#include "primitives2d.hpp"
#include "contour_tree.hpp"

namespace HM2D{ namespace Grid{ namespace Impl{

struct GraphLine{
	GraphLine(int a, int b):i0(a), i1(b) { if (i0>i1) std::swap(i0,i1); }
	int i0, i1;
	void set(int a, int b) { i0=a; i1=b; if (i0>i1) std::swap(i0, i1); }
	bool has_node(int i) const  { return (i0==i || i1 ==i); }
};
inline bool operator<(const GraphLine& x, const GraphLine& y){
	return (x.i0!=y.i0) ? (x.i0<y.i0) : (x.i1<y.i1);
}

//Represents wireframe of the 2D grid: points and edges
struct PtsGraph{
	PtsGraph(){};
	// ==== create wireframe from grid
	explicit PtsGraph(const GridData& tg);
	explicit PtsGraph(const HM2D::EdgeData& cc);

	//get data
	int Nnodes() const {return (int)nodes.size();}
	int Nlines() const {return (int)lines.size();}
	const Point* get_point(int i) const { return &nodes[i]; }
	std::pair<int, int> get_line(int i) const { return std::make_pair(lines[i].i0, lines[i].i1); }

	//generate data
	std::vector<Point> center_line_points() const;

	//create 2D grid on the basis of current wireframe
	//all boundary edges has an adjacent cell with index=-1.
	GridData togrid() const;
	EdgeData toedges() const;

	//imposes contour edges.
	void add_edges(const EdgeData&);

	//deletes edges which centers lie outside (dir = OUTSIDE) or inside (dir = INSIDE) contour
	void exclude_area(const Contour::Tree& cont, int dir);
	void exclude_lines(const std::set<int>& exlines);

	//cuts wmain with the contours collection internals
	//all collection contours lines will be present in the resulting graph
	//dir = INSIDE leave only inner zone, dir = OUTSIDE leave outer zone
	static PtsGraph cut(const PtsGraph& wmain, const Contour::Tree& conts, int dir);
	
	//overalay two graphs
	static PtsGraph overlay(const PtsGraph& wmain, const PtsGraph& wsec);
private:
	// ==== main data
	vector<Point> nodes;
	vector<GraphLine> lines;
	
	// ==== info
	Point& ednode0(int i)  { return nodes[lines[i].i0]; }
	Point& ednode1(int i)  { return nodes[lines[i].i1]; }

	// ===== algorithms
	//Imposes wireframe @imp_graph upon @main_graph. 
	//  @eps -- is radius at which points of @main_graph and @imp_graph are considered equal
	//It is guaranteed that:
	//-- resulting wireframe sections will not intersect any sections of old ones,
	//-- points of @main_graph will not be deleted or moved, their order will remain the same,
	//-- nodes in @imp_grid can be translated within @eps radius to node or line of @main_graph.
	//Returns tuple of:
	//  <0>  -- Resulting graph,
	//  <1>  -- vector[@imp_graph.size()] which represents the position of secondary graph points
	//          within resulting graph points list
	//  <2>  -- vector[@resulting_graph.nodes.size()] of @main_graph edge coordinates 
	//          of points which were placed on its edges  or -1
	typedef std::tuple<
		PtsGraph,      // resulting graph
		vector<int>,   // ig_nodes
		vector<double> // crossnodes
	> impResT;
	static impResT impose(const PtsGraph& main_graph, const PtsGraph& imp_graph);
	
	// subprocedure for togrid() function
	static GridData intrusion_algo(const vector<GridData>& grids);
private:
	static impResT _impose_impl(const PtsGraph& main_graph, const PtsGraph& imp_graph);

	friend struct PtsGraphAccel;
	friend struct SubGraph;

	//data management
	void delete_unused_points();
	//lines->lines connectivity
	vector<vector<int>> lines_lines_tab() const;
};

//single connected part of the ptsgraph
struct SubGraph{
	const PtsGraph* parent;
	std::list<int> lines;

	static vector<SubGraph> build(const PtsGraph& p);
	GridData togrid() const;
	HM2D::EdgeData tocontour() const;
private:
	SubGraph(const std::list<int>& used_lines, const PtsGraph* p);

	//recursive algo for subgraph assembling
	static void lines_sort_out(int iline, std::vector<int>& lines_usage, std::list<int>& result,
			const vector<std::vector<int>>& contab);

	//build Graph base data
	std::tuple<
		vector<Point>,
		vector<GraphLine>
	> rebuild_nodes_lines() const;
};

typedef std::array<int, 2> Tind2;
struct Tind2Proc{
	Tind2Proc(int nx, int ny) : Nx(nx), Ny(ny), N(nx*ny), indicies(Ny){
		for (int iy=0; iy<Ny; ++iy){
			indicies[iy].resize(Nx);
			std::iota(indicies[iy].begin(), indicies[iy].end(), iy*Nx);
		}
	}
	int size()  const  { return N; }
	int sizex() const  { return Nx; }
	int sizey() const  { return Ny; }
	//get plain index procedures
	int gindex(Tind2 ind) const  { return indicies[ind[1]][ind[0]]; }
	std::list<int> gindex(Tind2 ind1, Tind2 ind2) const  {
		if (ind1[0]>ind2[0]) std::swap(ind1[0], ind2[0]);
		if (ind1[1]>ind2[1]) std::swap(ind1[1], ind2[1]);
		std::list<int> ret;
		for (int i=ind1[1]; i<=ind2[1]; ++i){
			std::copy(indicies[i].begin()+ind1[0], 
			          indicies[i].begin()+ind2[0]+1, 
			          std::back_inserter(ret));
		}
		return ret;
	}
	//get index
	Tind2 get_tind(double ix, double iy) const  {
		Tind2 ret = {{int(ix), int(iy)}};
		if (ret[0]<0) ret[0]=0; else if (ret[0]>=Nx) ret[0]=Nx-1;
		if (ret[1]<0) ret[1]=0; else if (ret[1]>=Ny) ret[1]=Ny-1;
		return ret;
	}
private:
	const int Nx;
	const int Ny;
	const int N;
	vector<vector<int>> indicies;

};

struct PtsGraphAccel{
	//Auxilliary grid partition Naux x Naux
	constexpr static const int Naux = 100;
	//constructors
	PtsGraphAccel(PtsGraph& g, const Point& pmin, const Point& pmax) ;
	PtsGraphAccel(PtsGraph& g) ;
	//Adds the point to the node vector if it lies further then @eps from existing node.
	//Returns a tuple which describes added point position:
	//<0> -- type of addition:
	//       = 0 -- node has been added to the end of point list
	//       = 1 -- node was equal to existing one and was not added (see <1>)
	//       = 2 -- node was moved to the edge of a graph and added to the @nodes end
	//<1> -- array index of added node. Equals last array index if <0>==0, 2
	//       and equals equivalent node index if @p was close to existing node (<0>==1).  
	std::tuple<int, int>  add_point(const Point& p) ;
	//Adds new graph connection.
	void add_connection(int i1, int i2) ;
	
	//finds graph line with contains point &p.
	//Returns lineIndex+lineWeight of the node or -1 if there is no such graph line.
	double find_gline(const Point& p) const ;
private:
	//Data
	const double eps;
	const Point epsPnt;
	const Tind2Proc ip;
	PtsGraph* G;
	Point p0;
	double hx, hy;
	//construction
	void init(Point pmin, Point pmax) ;

	//additional functionality for PtsGraph data access
	Point ednodep(int i)  { 
	       Point& p0 = G->ednode0(i); Point& p1 = G->ednode1(i);
	       return Point(std::max(p0.x, p1.x)+eps,
	       	     std::max(p0.y, p1.y)+eps);
	}
	Point ednodem(int i)  { 
		Point& p0 = G->ednode0(i); Point& p1 = G->ednode1(i);
		return Point(std::min(p0.x, p1.x)-eps,
			     std::min(p0.y, p1.y)-eps);
	}

	Tind2 point_sqind(const Point& p) const {
		return  ip.get_tind((p.x-p0.x)/hx, (p.y-p0.y)/hy);
	}
	std::set<int> candidates_edges(const Point& p0, const Point& p1) const {
		std::set<int> ret;
		auto gind = ip.gindex(point_sqind(p0), point_sqind(p1));
		for (auto i: gind) ret.insert(edmap[i].begin(), edmap[i].end());
		return ret;
	}
	std::set<int> candidates_edges(const Point& p0) const {
		return edmap[ ip.gindex(point_sqind(p0)) ];
	}
	std::set<int> candidates_points(Point p0, Point p1) const {
		//add epsilons to widen square if p0, p1 line is parallel to x or y axis
		if (fabs(p0.x-p1.x)<eps) { p0.x+=eps; p1.x-=eps; }
		if (fabs(p0.y-p1.y)<eps) { p0.y+=eps; p1.y-=eps; }
		//build return set
		std::set<int> ret;
		auto gind = ip.gindex(point_sqind(p0), point_sqind(p1));
		for (auto i: gind) ret.insert(ndmap[i].begin(), ndmap[i].end());
		return ret;
	}

	//add new data to original graph and search data
	//returns the index of newly added node
	int add_graph_node(const Point& p) ;
	int break_graph_line(int iline, double ksi) ;
	int break_graph_line(int iline, const Point& p) ;
	//adds edge to data and graph which can intersect existing one
	void add_clear_connection(int i1, int i2) ;
	//adds non-intersecting edge to data and graph
	void add_graph_edge(int i0, int i1) ;
	//resets existing edge
	void reset_graph_edge(int i0, int i1, int iline) ;
	
	// ====================== data for accelerated search procedures
	//add data which already exists in G to search data
	void add_node_to_map(int inode)  { 
		ndfinder.add(G->nodes[inode], inode);
		ndmap[ ip.gindex(point_sqind(G->nodes[inode])) ].insert(inode);
	}
	void add_edge_to_map(int iline) {
		Tind2 i0 = point_sqind(ednodem(iline)), i1 = point_sqind(ednodep(iline));
		for (auto k: ip.gindex(i0,i1)) edmap[k].insert(iline);
	}
	void delete_edge_from_map(int iline) {
		Tind2 i0 = point_sqind(ednodem(iline)), i1 = point_sqind(ednodep(iline));
		for (auto k: ip.gindex(i0,i1)) edmap[k].erase(iline);
	}
	//edges collection: square -> edge set
	vector<std::set<int>> edmap;
	vector<std::set<int>> ndmap;
	CoordinateMap2D<int> ndfinder;
};

}}}
#endif
