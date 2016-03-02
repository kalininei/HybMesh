#ifndef HMCONT2D_CONTOUR_HPP
#define HMCONT2D_CONTOUR_HPP
#include "collections.hpp"
#include "containers.hpp"

namespace HMCont2D{

// ============== Contour: collection of connected edges
struct Contour: public ECollection{
	Contour(const Contour& col): ECollection(col){}
	Contour(): ECollection(){};
	//get info
	//endpoints
	Point* first() const;
	Point* last() const;
	//first point == last point for closed paths
	vector<Point*> ordered_points() const;
	//like ordered points but removes zero length sections
	vector<Point*> unique_points() const;
	//list of all corner points in correct order without doubling
	//the last point for closed contours
	//for open contours includes first and last points by default
	vector<Point*> corner_points() const;
	//same as corner_points but doubles end points for open contours
	vector<Point*> corner_points1() const;
	//returns [point previous, point currant, point next]
	//or NULLS if no such points
	std::array<Point*, 3> point_siblings(Point* p) const;
	std::array<Point*, 3> point_siblings(int i) const;
	//return next point or null
	Point* next_point(Point* p) const{ return std::get<2>(point_siblings(p));}
	Point* prev_point(Point* p) const{ return std::get<0>(point_siblings(p));}
	//properties
	bool is_closed() const { return first() == last(); }
	bool is_straight() const { return (is_closed()) ? false : corner_points().size() == 2; }
	bool correctly_directed_edge(int i) const;

	struct PInfo{
		Point *p, *pprev, *pnext;
		shared_ptr<Edge>  eprev, enext;
		int index;
	};
	//returns vector of length (size()+1)
	//detailed information about each node connection.
	vector<PInfo> ordered_info() const;
	PInfo pinfo(Point* p) const;

	//gives
	//<0>  contour length coordinate, 
	//<1>  contour weight coordinate,
	//<2>  index of edge
	//<3>  edge weight coordinate
	//<4>  distance to contour
	//of the contour point closest to given one.
	std::tuple<double, double, int, double, double>
	coord_at(const Point& p) const;
	
	// ====== overridden from Collection
	//only if TTarget is also a contour.
	//only if endpoints are connected.
	//if impossible -> throws HMCont2D::GeomError
	template<class TTarget, class = Tpp::IsBase<Contour, TTarget>>
	void Unite(const TTarget& c);
	
	//Methods
	//reverse edge order
	void Reverse(){std::reverse(data.begin(), data.end());}
	//directs internal edges order according to points direction in contour
	void DirectEdges();
	//reverse edge order and internal points order within edges
	void ReallyReverse();
	//adds edge (last(), p) to the end
	void AddLastPoint(Point* p);
	//deletes edge[i] end point
	void RemoveEdge(int i);

	//removes edge next to deleted point,
	//moves previous edge start or end point
	//closed/open contours stay close/open or become zero sized;
	//removed edge (if it is shared by someone else) has NULL data entries
	void RemovePoint(const Point* p);
	void RemovePoints(const vector<const Point*>& p);

	//force positive(true)/negative(false) direction.
	//->true if contour was reversed
	bool ForceDirection(bool dir);
	//Find a point on a path, closest to p and place it there by splitting edge.
	//If p already exists do nothing.
	//Args: 
	//    p - point coordinate, pcol - collection for placing created point
	//Returns:
	//    <0> if point was placed
	//    <1> pointer to a new added point or existed one which equals p
	std::tuple<bool, Point*>
	GuaranteePoint(const Point& p, PCollection& pcol);

	//find coordinates of closest contour point
	Point ClosestPoint(const Point& p) const;

	//Returns true if point lies strictly within/without closed contour
	//Direction is not taken into account
	bool IsWithin(const Point& p) const;
	bool IsWithout(const Point& p) const;
	bool AllWithin(const vector<Point>& p) const;
	bool AllWithout(const vector<Point>& p) const;
	//-1 - point is on polygon
	// 0 - point is outside polygon
	// 1 - point is in polygon  -> this is not reliable
	int WhereIs(const Point& p) const;
	//Return arbitrary inner point for a closed contour
	//Direction is not considered
	Point InnerPoint() const;

	// ======= Algorithms
	//positive/negative value for inner/outer contours
	static double Area(const Contour& c);

	//weigth points by [0,1] weights of full contour length.
	//Returns point in sorted order from start to end point of c
	static PCollection WeightPoints(const Contour& c, vector<double> w);
	static Point WeightPoint(const Contour& c, double w){
		return *WeightPoints(c, {w}).point(0);
	}
	//weigth points by lenght
	static PCollection WeightPointsByLen(const Contour& c, vector<double> lens);
	static Point WeightPointByLen(const Contour& c, double len){
		return *WeightPointsByLen(c, {len}).point(0);
	}

	//return weights of points as given by c.ordered_points()
	//first is 0, last is always 1.
	static vector<double> EWeights(const Contour& c);

};


template<class TTarget, class>
void Contour::Unite(const TTarget& c){
	auto self0 = first(), self1 = last();
	auto target0 = c.first(), target1 = c.last();
	//choosing option for unition
	if (size() == 0 || c.size() == 0 ) goto COPY12;
	else if (self0 == self1 || target0 == target1) goto THROW;
	else if (c.size() == 1 && c.edge(0)->contains(self1)) goto COPY12;
	else if (c.size() == 1 && c.edge(0)->contains(self0)) goto COPY03;
	//try to add new contour to the end of current
	else if (self1 == target0) goto COPY12;
	else if (self1 == target1) goto NEED_REVERSE;
	//if failed try to add before the start
	else if (self0 == target1) goto COPY03;
	else if (self0 == target0) goto NEED_REVERSE;
	else goto THROW;

	COPY03:{
		ShpVector<Tvalue> dt;
		std::copy(c.begin(), c.end(), std::back_inserter(dt));
		std::copy(begin(), end(), std::back_inserter(dt));
		data = dt;
		return;
	}
	COPY12:{
		std::copy(c.begin(), c.end(), std::back_inserter(data));
		return;
	}
	NEED_REVERSE:{
		TTarget tmp(c);
		tmp.Reverse();
		return Unite(tmp);
	}
	THROW:{
		throw GeomError("Impossible to unite non-connected contours");
	}
}



}

#endif
