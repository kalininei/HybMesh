#ifndef HMCONT2D_CONTOUR_HPP
#define HMCONT2D_CONTOUR_HPP
#include "collections.hpp"
#include "containers.hpp"
#include "algos.hpp"

namespace HMCont2D{

// ============== Contour: collection of connected edges
struct Contour: public ECollection{
	//get info
	//endpoints
	Point* first() const;
	Point* last() const;
	bool is_closed() const { return first() == last(); }
	//first point == last point for closed paths
	vector<Point*> ordered_points() const;
	//list of all non-corner points
	vector<Point*> corner_points() const;
	

	//Methods
	void Reverse(){ std::reverse(data.begin(), data.end()); }
	//force positive(true)/negative(false) direction.
	//->true if conotour was reversed
	bool ForceDirection(bool dir);

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





}

#endif
