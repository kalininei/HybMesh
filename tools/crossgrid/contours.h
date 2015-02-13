#ifndef CROSSGRID_CONTOURS_H
#define CROSSGRID_CONTOURS_H
#include "bgeom.h"

//Contour basic class
class PContour{
	std::vector<Point*> pts;
	std::pair<Point, Point> rectangle_bnd() const;
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

	//build reversed contour
	PContour reverse() const;
	//reverses itself
	void reverse_self();

	//contour data access
	int n_points() const { return pts.size(); }
	const Point* get_point(int i) const { return pts[i]; }

	//contour geometry procedures
	void select_points(const vector<Point*>& pts, 
			vector<Point*>& inner, vector<Point*>& outer) const;

	//returns squared distances to points. Sign depends on whether points lies 
	//outside(-) or inside(+) the contour
	vector<double> meas_points(const vector<const Point*>& pts) const;

	//additional procedures
	//return point position with respect to the contour:
	//    -1 if outside, 1 if inside, 0 if lies on the contour
	int is_inside(const Point& p, const bool* inner_hint = 0) const;
	//filter inner, outer and contour points from point list
	std::tuple<
		vector<Point*>,  //internal points
		vector<Point*>,  //points on contour
		vector<Point*>   //outer points
	> filter_points(const vector<Point*>& points) const;
	//contour area
	double area() const;
	//length of each section
	vector<double> section_lenghts() const;
	//return distances which are covered by each node  = (hleft+hright)/2.0
	vector<double> chdist() const;
};

//collection of contours.
//Builds a contour tree.
//Automatically reverses contours in the way
//that all first level contours be inner
struct ContoursCollection{
	ContoursCollection(){};
	ContoursCollection(const vector<PContour>& cnts);
	void add_contour(const PContour& cnt);
	vector<PContour> contours_list() const;
	int is_inside(const Point& p) const;
	int n_cont() const { return contours.size(); }
private:
	struct _entry{
		_entry* upper;
		_entry(PContour* d): upper(0), is_inner(d->area()>0), data(d){}
		vector<_entry*> lower;
		bool is_inner;
		PContour* data;
		int geom_inside(const Point& p) const;
		void set_nesting(bool inner);
		const _entry* find(const Point& p) const;
	};
	shp_vector<PContour>  contours;
	shp_vector<_entry>  entries;
	vector<_entry*> top_level;
	void set_nesting();
	const _entry* efind(const Point& p) const;
};

//Contour which owns all its points
class Contour: public PContour{
	shp_vector<Point> pdata;
public:
	Contour(const std::vector<Point>& _pts=std::vector<Point>());
	Contour(const PContour& c);

	void add_point(Point* p);
	void add_point(const Point& p);
	void add_point(double x, double y){ add_point(Point(x,y)); }
};



#endif
