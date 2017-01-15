#ifndef HYBMESH_CONNECTORS_HPP
#define HYBMESH_CONNECTORS_HPP
#include "bgrid.hpp"
#include "canonic_bgrid.hpp"
namespace HMBlay{namespace Impl{

struct MConnector{
protected:
	virtual BGrid* ConnectionGrid(){ return 0;}
	// modifies prev and next meshes
	virtual void ModifyAdjacents(){};
	// builds grid in a connection area itself
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
	void Add(BGrid& g, bool with_prev, bool with_next);
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
	void AssembleGrid(vector<ShpVector<HM2D::Vertex>>& pallpts, shared_ptr<HM2D::Vertex> cpoint);
public:
	RoundConnector(MappedMesher* _prev, MappedMesher* _next): MConnector(_prev, _next){};
};

struct RightConnector: public MConnector{
private:
	shared_ptr<MappedRect> connection_area;
	shared_ptr<MappedMesher> connection_grid;
	HM2D::EdgeData left, bot;
	HM2D::EdgeData right, top;
	BGrid* ConnectionGrid() override { return &connection_grid->result;}
	void BuildInternals() override;
public:
	RightConnector(MappedMesher* _prev, MappedMesher* _next);
};

struct ReentrantConnector: public MConnector{
private:
	shared_ptr<MappedRect> connection_area;
	shared_ptr<MappedMesher> connection_grid;
	shared_ptr<HM2D::Vertex> top_right_pnt;
	HM2D::EdgeData left, bot, right, top;
	BGrid* ConnectionGrid() override { return &connection_grid->result;}
	void BuildInternals() override;
public:
	ReentrantConnector(MappedMesher* _prev, MappedMesher* _next): MConnector(_prev, _next){}
};

}}
#endif

