#ifndef VTK_EXPORT_GRID_3D_HPP
#define VTK_EXPORT_GRID_3D_HPP

#include "hmproject.h"
#include "hmcallback.hpp"
#include "primitives_grid3d.hpp"
#include "serialize_grid3d.hpp"

namespace HMGrid3D{namespace Export{

//save grid (only tetra, hex, wedge cells
//signature void GridVTK(const HMGrid3D::Grid& g, std::string fn);
struct TGridVTK: public HMCallback::ExecutorBase{
	HMCB_SET_PROCNAME("Exporting 3d grid to *.vtk");

	HMCB_SET_DURATION(80, const SGrid&, std::string);
	void _run(const SGrid& g, std::string fn);

};
extern HMCallback::FunctionWithCallback<TGridVTK> GridVTK;


//signature void BoundaryVTK(const HMGrid3D::Grid& g, std::string* fn);
struct TBoundaryVTK: public HMCallback::ExecutorBase{
	HMCB_SET_PROCNAME("Exporting 3d grid surface to *.vtk");

	HMCB_SET_DURATION(50, SGrid, std::string);
	void _run(const SGrid&, std::string);

};
extern HMCallback::FunctionWithCallback<TBoundaryVTK> BoundaryVTK;

//boundary + grid to 2 separate files
struct TAllVTK: public HMCallback::ExecutorBase{
	HMCB_SET_PROCNAME("Exporting 3d data to *.vtk");
	HMCB_SET_DEFAULT_DURATION(
		HMCB_DURATION(TBoundaryVTK, SGrid, std::string)+ 
		HMCB_DURATION(TGridVTK, SGrid, std::string)
	);

	void _run(const SGrid& g, std::string fngrid, std::string fnbnd);
};
extern HMCallback::FunctionWithCallback<TAllVTK> AllVTK;


//structure which tries to treat 3D cell
//as a tetrahedron/hexahedron/prism/wedge.
//Used in vtk and gmsh exporting routines
class vtkcell_expression{
	bool _false_return();

	static int find_opposite(std::vector<std::vector<int>>& data, int cind, int pind);
	static int is_opposite(std::vector<std::vector<int>>& data, int i1, int i2);
	static std::pair<int, int> get_opposite(std::vector<std::vector<int>>& data, int ic);
	bool try_tetrahedron(std::vector<std::vector<int>>& data);
	bool try_hexahedron(std::vector<std::vector<int>>& data);
	bool try_pyramid(std::vector<std::vector<int>>& data);
	bool try_wedge(std::vector<std::vector<int>>& data);
	static vtkcell_expression build(std::vector<std::vector<int>>& cint);
public:
	//holds points indicies in vtk type orders
	vector<int> pts;
	int celltype; //10 - tetrahedron; 12 - hexahedron; 14 - pyramid; 13 - wedge

	//tries to build expressions for all cells in ser.
	//aface is face_vertex connectivity table
	static vector<vtkcell_expression> cell_assembler(const SGrid& ser,
			const vector<vector<int>>& aface);
	virtual std::string to_string() const;

	int wsize() const;  //number of points + 1
	std::string stype() const; //string(celltype)
};

}}

#endif
