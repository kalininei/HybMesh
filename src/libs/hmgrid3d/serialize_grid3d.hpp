#ifndef SERIALIZE_GRID3D_HPP
#define SERIALIZE_GRID3D_HPP
#include "primitives_grid3d.hpp"
#include <type_traits>
#include "hmcallback.hpp"

namespace HMGrid3D{

//only linear data is saved here
struct SimpleSerialize{
	//primitives sizes
	int n_vert, n_faces, n_edges, n_cells;

	//data
	std::vector<double> vert;  //x0, y0, z0, x1, y1, z1, ...
	std::vector<int> edges;    //edge0_start, edge0_end, edge1_start, edge1_end, ...
	std::vector<int> faces;    //face0_num_edges, face0_edge0, ..., face0_leftcell, face0_rightcell, ...
	std::vector<int> cells;    //cell0_num_faces, cell0_face0, cell0_face1, ....
	std::vector<int> bnd;      //btype0 num_faces_in_btype0 btype0_face0, ....  #only non-zero boundaries: WHY???

	//convert functions

	//possible arguments: (const Grid&),  (const Grid&, Grid::Talldata&)
	//returns: SimpleSerialize
	struct ConvertExe: public HMCallback::ExecutorBase{
		typedef SimpleSerialize TRet1;
		HMCB_SET_PROCNAME("Serialization");
		HMCB_SET_DEFAULT_DURATION(100);

		TRet1 _run(const Grid& g);
		TRet1 _run(const Grid& g, Grid::Talldata& dt);
	};
	static HMCallback::FunctionWithCallback<ConvertExe> Convert;
protected:
	//uses 100 units of HMCallback::Caller2
	void fill_from_grid(const Grid& g, Grid::Talldata& dt, HMCallback::Caller2&);
};


struct ExtendedSimpleSerialize: public SimpleSerialize{
	//==== additional data
	//get<0> - ShpVector<Vertex>,
	//get<1> - ShpVector<Edge>,
	//get<2> - ShpVector<Face>,
	//get<3> - ShpVector<Cell>
	Grid::Talldata _alldata;
	vector<int> icell, iface;

	//==== data access
	ShpVector<Vertex>& data_vertex(){ return std::get<0>(_alldata); }

	//==== additional methods
	//face as a sequence on vertices
	vector<int> face_assembler(int nface) const;
	vector<vector<int>> face_assembler() const;
	Grid to_grid() const;
	
	//convert functions
	//possible arguments -- (const Grid& g)
	//returns ExtendedSimpleSerialize
	struct ConvertExe: public HMCallback::ExecutorBase{
		typedef ExtendedSimpleSerialize TRet1;
		HMCB_SET_PROCNAME("Extended serialization");
		HMCB_SET_DEFAULT_DURATION(
			10+
			HMCB_DURATION(SimpleSerialize::ConvertExe, Grid, Grid::Talldata)
		);

		TRet1 _run(const Grid& g);
	};
	static HMCallback::FunctionWithCallback<ConvertExe> Convert;
};

typedef SimpleSerialize SS;
typedef ExtendedSimpleSerialize ESS;

}

#endif
