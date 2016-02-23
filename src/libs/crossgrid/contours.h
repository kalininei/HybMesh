#ifndef CROSSGRID_CONTOURS_H
#define CROSSGRID_CONTOURS_H
#include "crossgrid.h"
#include "bgeom2d.h"
#include <list>

//Contour basic class
class PContour{
	std::vector<Point*> pts;
	//points which lies within contour or on its edge
	vector<const Point*> find_inner(const vector<const Point*>& pts) const;
	//meas from point to contour
	double meas_to_point(const Point& p) const;
public:
	//construction and modification
	PContour(const vector<Point*>& _pts=vector<Point*>()): pts(_pts){}
	template<class C>
	static PContour build(C& _pts){
		PContour ret = PContour();
		ret.pts = aa::Cfill_vector(_pts, [](Point* x){ return (Point*)x; });
		return ret;
	}

	virtual void add_point(Point* p){ pts.push_back(p); }
	void delete_by_index(const std::set<int>& badind);

	//build reversed contour
	PContour reverse() const;
	//reverses itself
	void reverse_self();

	//build simplified contour
	PContour simplify() const;
	//simplify itself
	void simplify_self();

	//contour data access
	int n_points() const { return pts.size(); }
	const Point* get_point(int i) const {
		if (i<0) return get_point(i+n_points());
		else if (i>=n_points()) return get_point(i-n_points());
		else return pts[i];
	}

	//returns squared distances to points. Sign depends on whether points lies
	//outside(-) or inside(+) the contour with respect ot its sence of rotation
	vector<double> meas_points(const vector<const Point*>& pts) const;

	//additional procedures
	//return point position with respect to the contour:
	//    OUTSIDE if outside, INSIDE if inside, BOUND if lies on the contour
	//contour direction matters
	int is_inside(const Point& p, const bool* inner_hint = 0) const;
	//filter inner, outer and contour points from point list
	//Contour direction matters
	std::tuple<
		vector<const Point*>,  //internal points
		vector<const Point*>,  //points on contour
		vector<const Point*>   //outer points
	> filter_points(const vector<const Point*>& points) const;

	//contour area
	double area() const;
	//length of each section
	vector<double> section_lenghts() const;
	//return distances which are covered by each node  = (hleft+hright)/2.0
	vector<double> chdist() const;
	//return distances which are covered by each node  = max(hleft,hright)
	vector<double> max_chdist() const;
	//is the i-th point lies on the section between i-1 and i+1 point
	bool is_corner_point(int i) const;

	//find internal point: cross algorithm
	//direction is not taken into account
	Point inside_point_ca() const;
	

	//finds intersections with another point.
	//for each intersection return contour coordinate (edge + [0,1]) and the point
	vector<std::pair<double, Point>> intersections(const PContour& p) const;
	vector<std::pair<double, Point>> intersections(const vector<PContour>& p) const;
};


//Contour which owns all its points
class Contour: public PContour{
	ShpVector<Point> pdata;
public:
	Contour(const std::vector<Point>& _pts=std::vector<Point>());
	Contour(const PContour& c);

	void add_point(Point* p);
	void add_point(const Point& p);
	void add_point(double x, double y){ add_point(Point(x,y)); }
};


//collection of contours.
//Builds a contour tree.
//Automatically reverses contours in the way
//that all first level contours be inner
struct ContoursCollection: public Cont{
	ContoursCollection(){};
	explicit ContoursCollection(const vector<PContour>& cnts);
	virtual void add_contour(const PContour& cnt);
	void remove_contour(const PContour* cnt);
	vector<PContour> contours_list() const;
	PContour contour(int i) const { return PContour(*contours[i]); }
	//-> BOUND if point lies on contour,
	//-> INSIDE if point lies inside, OUTSIDE if outside
	int is_inside(const Point& p) const;
	int n_cont() const { return contours.size(); }
	int n_inner_cont() const {
		int ret = 0;
		for (auto c: entries) if (c->is_inner) ++ret;
		return ret;
	}

	//contours management
	const PContour* get_contour(int i) const { return entries[i]->data; }
	bool is_inner(int i) const { return entries[i]->is_inner; }
	const PContour* get_parent(int i) const {
		return (entries[i]->upper == 0) ? 0: entries[i]->upper->data;
	}
	int get_parent_index(int i) const;
	int get_level(int i) const { return entries[i]->get_level();}

	//returns contours collection which contains i-th contour and its
	//first level children
	ContoursCollection level_01(int i) const;
	//return collection which contains only [ist, iend] levels of current
	ContoursCollection cut_by_level(int ist, int iend) const;

	std::list<const PContour*> get_childs(int i) const {
		std::list<const PContour*> ret;
		for (auto c: entries[i]->lower) ret.push_back(c->data);
		return ret;
	}
	int num_childs(int i) const { return entries[i]->lower.size(); }
	
	//contours geometry procedures
	std::tuple<
		vector<int>,  //internal points
		vector<int>,  //points on contour
		vector<int>   //outer points
	> filter_points_i(const vector<Point>& points) const;
	std::tuple<
		vector<int>,  //internal points
		vector<int>,  //points on contour
		vector<int>   //outer points
	> filter_points_i(const vector<const Point*>& points) const;

	//total area
	double area() const;

protected:
	struct _entry{
		_entry* upper;
		_entry(PContour* d): upper(0), is_inner(d->area()>0), data(d){}
		std::list<_entry*> lower;
		bool is_inner;
		mutable PContour* data;
		int geom_inside(const Point& p) const;
		void set_nesting(bool inner);
		//finds the entry wich contains p among this and this->lower
		const _entry* find(const Point& p) const;
		int get_level() const{
			return (upper==0) ? 0 : 1 + upper->get_level();
		}
	};
	ShpVector<PContour>  contours;
	ShpVector<_entry>  entries;
	std::list<_entry*> top_level;
	void set_nesting();
	const _entry* efind(const Point& p) const;

	//simplify/unsimplify entries: 
	//modifies entires.data to simplified contours
	mutable ShpVector<PContour> _simpcont;
	mutable vector<PContour*> _origcont;
	void simplify_entries() const;
	void unsimplify_entries() const;
};

//Contour collection which owns all its points
class PointsContoursCollection: public ContoursCollection{
	ShpVector<Point> pdata;
	void build(const vector<Point>& pts, const vector<int>& eds);
	//edges management
	struct Edge{
		Edge(const PointsContoursCollection* par, int p0, int p1, int gi): 
			parent(par), i0(p0), i1(p1), index(gi), _angle(-1){}
		const PointsContoursCollection* parent;
		int i0, i1;  //point indicies
		int index;   //edge index
		mutable double _angle; //angle [0, pi) between edge and 0x axis
		//proc
		double get_angle() const;
		const Point* pnt0() const { return parent->pdata[i0].get(); }
		const Point* pnt1() const { return parent->pdata[i1].get(); }
		const Point center() const { return Point::Weigh(*pnt0(), *pnt1(), 0.5); }
	};
	vector<Edge> edges;
public:
	void add_contour(const PContour& cnt){
		throw std::runtime_error(
			"No way of the new contours addition "\
		        "to PointsContoursCollection at runtime"
		);
	}

	PointsContoursCollection(const vector<double>& pts, const vector<int>& eds);
	PointsContoursCollection(const vector<Point>& pts, const vector<int>& eds);
	PointsContoursCollection(const ContoursCollection& col);

	int n_edges() const { return edges.size();}
	int n_total_points() const { return pdata.size(); }
	const Point* get_point(int i) const { return pdata[i].get(); }
	std::pair<int, int> get_edge(int i) const { return std::make_pair(edges[i].i0, edges[i].i1); }

	virtual ~PointsContoursCollection(){}
	ContoursCollection shallow_copy();

	void do_scale(const ScaleBase& sc);
	void undo_scale(const ScaleBase& sc);

	//for each source edge returns corresponding target edge
	//or -1 if no edge was found
	static vector<int> edge_correlation(
			const PointsContoursCollection& src,
			const PointsContoursCollection& tar);
};

class CGBoundingBox: public BoundingBox{
//crossgrid copy of bounding box
//TODO: place to conours2d library
public:
	CGBoundingBox(double x0, double y0, double x1, double y1):BoundingBox(x0, y0, x1, y1){}
	CGBoundingBox(const vector<CGBoundingBox>&, double e=0.0);
	CGBoundingBox(const PContour& cont, double e=0.0);
	CGBoundingBox(const ContoursCollection& col, double e=0.0);
	CGBoundingBox(const Point& p1, const Point& p2, double e=0.0);

	Contour get_contour() const;
};


#endif
