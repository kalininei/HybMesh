#ifndef HMCONT2D_CONTOUR_HPP
#define HMCONT2D_CONTOUR_HPP
#include "collections.hpp"
#include "containers.hpp"
#include "algos.hpp"

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
	//list of all corner points in correct order without doubling
	//the last point for closed contours
	//for closed contours may not include first and last points
	vector<Point*> corner_points() const;
	//returns [point previous, point currant, point next]
	//or NULLS if no such points
	std::array<Point*, 3> point_siblings(Point* p) const;
	//return next point or null
	Point* next_point(Point* p) const{ return std::get<2>(point_siblings(p));}
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

	//Returns true if point lies strictly within/without closed contour
	bool IsWithin(const Point& p) const;
	bool IsWithout(const Point& p) const;

	// ======= Algorithms
	//calculates vector representing direction of contour at point p smoother by lengh len
	//if p doesn't lie on c -> project it to c and calculate
	static Vect SmoothedDirection(const Contour& c, Point* p, int direction, double len);

	//finds first cross (with respect to length of c1) of contours c1, c2.
	//returns <0>: if cross was found
	//        <1>: cross point
	//        <2,3>: normalized length coordinate of intersection
	static std::tuple<bool, Point, double, double> Cross(const Contour& c1, const Contour& c2);
		
	static double Area(const Contour& c);
	//weigth points by [0,1] weights of full contour length.
	//Returns point in sorted order from start to end point of c
	static PCollection WeightPoints(const Contour& c, vector<double> w);
	static Point WeightPoint(const Contour& c, double len){
		return *WeightPoints(c, {len}).point(0);
	}
	//weigth points by lenght
	static PCollection WeightPointsByLen(const Contour& c, vector<double> lens);
	static Point WeightPointByLen(const Contour& c, double len){
		return *WeightPointsByLen(c, {len}).point(0);
	}

	//return weights of points as given by c.ordered_points()
	//first is 0, last is always 1.
	static vector<double> EWeights(const Contour& c);

	//Offset
	static Container<ContourTree> Offset(const Contour& source, double delta, OffsetTp tp);
	static Container<Contour> OffsetOuter(const Contour& source, double delta);
	//Cut contour
	static Container<Contour> CutByWeight(const Contour& source, double w1, double w2);
	static Container<Contour> CutByLen(const Contour& source, double len1, double len2);

	//Assemble from shattered edges
	static Contour Assemble(const ECollection& col, const Point* pnt_start, const Point* pnt_end);
	static Contour Assemble(const ECollection& col, const Point* pnt_start);
	//Assemble from another contour
	static Contour Assemble(const Contour& col, const Point* pnt_start, const Point* pnt_end);
	//assemles for pnt_start in the direction (+-1) til the length of
	//resulting contour will be more then len
	static Contour Assemble(const Contour& col, const Point* pnt_start, int direction, double len);

	//Partition contour.
	//step - step of partitioning
	//pstore - point collection where to put new generated points with defined ShpGenerator
	//tp: keep all old points/keep only shaped points/ignore all old points except end points
	static Contour Partition(double step, const Contour& contour, PCollection& pstore, PartitionTp tp);

	//make a partition keeping points defined in keepit
	static Contour Partition(double step, const Contour& contour, PCollection& pstore,
			const std::vector<Point*>& keepit = {});

	//this can be called with Container<Contour> object.
	template<class TContainer,
		class = Tpp::IsBase<Contour, typename TContainer::TParent>>
	static TContainer Partition(double step, const TContainer& container, PartitionTp tp){
		PCollection pc;
		auto x = Partition(step, container, pc, tp);
		TContainer ret;
		TContainer::DeepCopy(x, ret);
		return ret;
	}

	template<class TContainer,
		class = Tpp::IsBase<Contour, typename TContainer::TParent>>
	static TContainer Partition(double step, const TContainer& container,
			const vector<Point*>& keepit = {}){
		return Partition(step, container, container.pdata, keepit);
	}
};


template<class TTarget, class>
void Contour::Unite(const TTarget& c){
	auto self0 = first(), self1 = last();
	auto target0 = c.first(), target1 = c.last();
	//choosing option for unition
	if (size() == 0 || c.size() == 0 ) goto COPY12;
	else if (self0 == self1 || target0 == target1) goto THROW;
	else if (c.size() == 1 && c.edge(0)->contains(self0)) goto COPY03;
	else if (c.size() == 1 && c.edge(0)->contains(self1)) goto COPY12;
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
