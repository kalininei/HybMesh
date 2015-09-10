#ifndef HMCONT2D_CONTOUR_HPP
#define HMCONT2D_CONTOUR_HPP
#include "collections.hpp"
#include "containers.hpp"
#include "algos.hpp"

namespace HMCont2D{

// ============== Contour: collection of connected edges
struct Contour: public ECollection{
	Contour(ECollection& col) = delete;
	Contour(): ECollection(){};
	//get info
	//endpoints
	Point* first() const;
	Point* last() const;
	bool is_closed() const { return first() == last(); }
	//first point == last point for closed paths
	vector<Point*> ordered_points() const;
	//list of all non-corner points
	vector<Point*> corner_points() const;
	//returns [point previous, point currant, point next]
	//or NULLS if no such points
	std::array<Point*, 3> point_siblings(Point* p) const;
	//return next point or null
	Point* next_point(Point* p) const{ return std::get<2>(point_siblings(p));}
	
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

	//Returns true if point lies strictly within closed contour
	bool IsWithin(const Point& p) const;

	// ======= Algorithms
	static double Area(const Contour& c);
	static PCollection WeightPoints(const Contour& c, vector<double> w);

	//Offset
	static Container<ContourTree> Offset(const Contour& source, double delta, OffsetTp tp);

	//Assemble from shattered edges
	static Contour Assemble(const ECollection& col, Point* pnt_start, Point* pnt_end);
	static Contour Assemble(const ECollection& col, Point* pnt_start);
	//Assemble from another contour
	static Contour Assemble(const Contour& col, Point* pnt_start, Point* pnt_end);

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
	//try to unite to the end first
	else if (self1 == target0) goto COPY12;
	else if (self1 == target1) goto NEED_REVERSE;
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
