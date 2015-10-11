#ifndef HYBMESH_CONNECTORS_HPP
#define HYBMESH_CONNECTORS_HPP
#include "bgrid.hpp"
#include "canonic_bgrid.hpp"
namespace HMBlay{namespace Impl{

struct MConnector{
protected:
	virtual BGrid* ConnectionGrid(){ return 0;}
	virtual void ModifyAdjacents(){};
	virtual void BuildInternals(){};
	MappedMesher *prev, *next;
	MConnector(MappedMesher* _prev, MappedMesher* _next): prev(_prev), next(_next){}
public:
	static shared_ptr<MConnector>
	Build(CornerTp tp, MappedMesher* prev, MappedMesher* next);
	
	//modifies prev, next meshes,
	//build connection mesh if nessessary.
	void Apply(){ ModifyAdjacents(); BuildInternals();}

	//adds next grid along with connection section to g
	void Add(shared_ptr<BGrid>& g, bool with_prev, bool with_next);
};

struct PlainConnector: public MConnector{
//connects grids which have congruent left|right nodes
private:
	void ModifyAdjacents() override;
public:
	PlainConnector(MappedMesher* prev, MappedMesher* next): MConnector(prev, next){};
};

struct AcuteConnector: public MConnector{
private:
	shared_ptr<BGrid> filler;
	BGrid* ConnectionGrid() override { return filler.get();}
	void BuildInternals() override;
	void ModifyAdjacents() override;
public:
	AcuteConnector(MappedMesher* _prev, MappedMesher* _next): MConnector(_prev, _next){}
};

struct RoundConnector: public MConnector{
private:
	shared_ptr<BGrid> filler;
	BGrid* ConnectionGrid() override { return filler.get();}
	void BuildInternals() override;
	void AssembleGrid(vector<ShpVector<GridPoint>>& pallpts, shared_ptr<GridPoint> cpoint);
public:
	RoundConnector(MappedMesher* _prev, MappedMesher* _next): MConnector(_prev, _next){};
};

struct RightConnector: public MConnector{
private:
	shared_ptr<MappedRect> connection_area;
	shared_ptr<MappedMesher> connection_grid;
	HMCont2D::Container<HMCont2D::Contour> left, bot;
	HMCont2D::Contour right, top;
	BGrid* ConnectionGrid() override { return &connection_grid->result;}
	void BuildInternals() override;
public:
	RightConnector(MappedMesher* _prev, MappedMesher* _next);
};

struct ReentrantConnector: public MConnector{
private:
	shared_ptr<MappedRect> connection_area;
	shared_ptr<MappedMesher> connection_grid;
	Point top_right_pnt;
	HMCont2D::Contour left, bot, right, top;
	BGrid* ConnectionGrid() override { return &connection_grid->result;}
	void BuildInternals() override;
public:
	ReentrantConnector(MappedMesher* _prev, MappedMesher* _next): MConnector(_prev, _next){}
};

}}
#endif

