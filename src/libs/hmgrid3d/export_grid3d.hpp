#ifndef EXPORT_GRID_3D_HPP
#define EXPORT_GRID_3D_HPP
#include "hmgrid3d.hpp"

namespace HMGrid3D{namespace Export{


std::tuple<
	std::vector<double>,   //points
	std::vector<int>       //
>
Serialize(const HMGrid3D::Grid& g);


// ===== VTK format
void GridVTK(const HMGrid3D::Grid& g, const char* fn);
void BoundaryVTK(const HMGrid3D::Grid& g, const char* fn);

// ====== Fluent msh format
struct PeriodicDataEntry{
	int bt, bt_shadow;
	Vertex v, v_shadow;
	bool reversed;
};
struct PeriodicData{
	PeriodicData(){};
	std::vector<PeriodicData> data;
	void add_condition(int b1, int b2, Vertex v1, Vertex v2, bool is_reversed);
	int size() const { return data.size(); }

	//assemble procedure
	HMGrid3D::Grid assemble(const HMGrid3D::Grid& g,
			std::map<Face*, Face*>& outmap);
};
//Main function
void GridMSH(const HMGrid3D::Grid& g, const char* fn,
		std::function<std::string(int)> btype_name,
		PeriodicData periodic);
//Simplified forms: no periodic, all boundaries are called "boundary1", 2, 3 etc.
void GridMSH(const HMGrid3D::Grid& g, const char* fn);

//Simplified forms: no periodic
void GridMSH(const HMGrid3D::Grid& g, const char* fn,
		std::function<std::string(int)> btype_name);





}}
#endif
