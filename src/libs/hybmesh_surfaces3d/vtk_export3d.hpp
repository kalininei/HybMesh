#ifndef VTK_EXPORT_GRID_3D_HPP
#define VTK_EXPORT_GRID_3D_HPP

#include "hmproject.h"
#include "hmcallback.hpp"
#include "primitives3d.hpp"
#include "serialize3d.hpp"

namespace HM3D{namespace Export{

//save grid (only tetra, hex, wedge cells
//signature void GridVTK(const HMGrid3D::Grid& g, std::string fn);
struct TGridVTK: public HMCallback::ExecutorBase{
	HMCB_SET_PROCNAME("Exporting 3d grid to *.vtk");
	HMCB_SET_DEFAULT_DURATION(80);

	void _run(const Ser::Grid& g, std::string fn);
	void _run(const GridData& g, std::string fn);

};
extern HMCallback::FunctionWithCallback<TGridVTK> GridVTK;


//signature void BoundaryVTK(const HMGrid3D::Grid& g, std::string* fn);
struct TBoundaryVTK: public HMCallback::ExecutorBase{
	HMCB_SET_PROCNAME("Exporting 3d grid surface to vtk");
	HMCB_SET_DEFAULT_DURATION(50);

	void _run(const Ser::Grid&, std::string);
	void _run(const GridData& g, std::string fn);

};
extern HMCallback::FunctionWithCallback<TBoundaryVTK> BoundaryVTK;

//boundary + grid to 2 separate files
struct TAllVTK: public HMCallback::ExecutorBase{
	HMCB_SET_PROCNAME("Exporting 3d data to vtk");
	HMCB_SET_DEFAULT_DURATION(
		HMCB_DURATION(TBoundaryVTK, Ser::Grid, std::string)+ 
		HMCB_DURATION(TGridVTK, Ser::Grid, std::string)
	);

	void _run(const Ser::Grid& g, std::string fngrid, std::string fnbnd);
	void _run(const GridData& g, std::string, std::string fnbnd);
};
extern HMCallback::FunctionWithCallback<TAllVTK> AllVTK;

struct TSurfaceVTK: public HMCallback::ExecutorBase{
	HMCB_SET_PROCNAME("Exporting surface to vtk");
	HMCB_SET_DURATION(80, const Ser::Surface&, std::string);
	HMCB_SET_DURATION(100, const FaceData&, std::string);

	void _run(const Ser::Surface& s, std::string fn);
	void _run(const FaceData& s, std::string fn);
};
extern HMCallback::FunctionWithCallback<TSurfaceVTK> SurfaceVTK;


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
	bool try_polygon(std::vector<std::vector<int>>& data);
	static vtkcell_expression build(std::vector<std::vector<int>>& cint);
public:
	//holds points indicies in vtk type orders
	vector<int> pts;
	int celltype; //10 - tetrahedron; 12 - hexahedron; 14 - pyramid; 13 - wedge

	//tries to build expressions for all cells in ser.
	//aface is face_vertex connectivity table
	static vector<vtkcell_expression> cell_assembler(const Ser::Grid& ser,
			const vector<vector<int>>& aface, bool ignore_errors=false);
	virtual std::string to_string() const;

	int wsize() const;  //number of points + 1
	std::string stype() const; //string(celltype)
};

}}

#endif
