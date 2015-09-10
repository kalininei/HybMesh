#ifndef HYBMESH_SIMPLE_BGRID_HPP
#define HYBMESH_SIMPLE_BGRID_HPP
#include "bgrid.hpp"

namespace HMBlay{ namespace Impl{

//grid with rectangular structure
class SimpleBGrid: public BGrid{
	//noraml smoothing options: number and weight
	const int Nsmooth;
	static constexpr double Wsmooth = 0.0;

	//each stencil entry is a line, perpendicular to source
	class TStencil: public vector<ShpVector<GridPoint>>{
		//properties
		bool is_closed;
		//inclination angles of segments
		vector<double> _angles;
		void fill_angles();
		vector<double> _basedist;
		void fill_basedist();
		double basedist(int i, int inext) const;
	public:
		//construction procedures
		//Only builds along normals. No checks.
		void Fill(ExtPath& pth);

		//tries to remove intersections in stencil.
		//No guarantee. All bad elements will be removed
		//during grids impositions.
		void NoIntersections(); 

		//change vectors angles with respect to adjecent vectors
		//newpos(i) = w*oldpos(i) + (1-w)/2*(k1*oldpos(i-1)+k2*oldpos(i+1));
		//doesn't change first/last if unclosed.
		//Internal weights k1, k2 are computing according to distances
		//between stencil base points
		void SmoothStep(double w); 

		//methods
		//Returns all points
		ShpVector<GridPoint> GetAll(); 

		//inclination angles of segments
		double get_angle(int i) const { return _angles[i]; }
		void set_angle(int, double);
	} stencil;

	//properties
	const bool _is_closed;  //is source path closed
	const int _hor_size;    //number of segments in a source path

	// ===== constructing procedures
	// writes data to stencil field
	// points are stored in GridGeom::points ShpVector
	void FillStencil(ExtPath& pth);
	// assemble cells
	void AssembleFromStencil();

public:
	//construct from extended path
	SimpleBGrid(ExtPath& pth, int _Nsmooth);

	//get properties
	bool is_closed() const { return _is_closed; }
	int hor_size() const {return _hor_size;}
};


}}
#endif
