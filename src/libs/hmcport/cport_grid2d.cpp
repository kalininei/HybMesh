#include "hmproject.h"
#include "cport_grid2d.h"
#include "grid.h"
#include "fluent_export_grid2d.hpp"


namespace{
GridGeom* togrid(void* g){
	return static_cast<GridGeom*>(g);
}
const GridGeom* togrid(const void* g){
	return static_cast<const GridGeom*>(g);
}
}

Grid2DBoundaryStruct*
set_grid2_boundary_types(int n, int* e1, int* e2, int* bt){
	Grid2DBoundaryStruct* ret = new Grid2DBoundaryStruct();
	ret->n = n;
	ret->edge_start_nodes = new int[n];
	ret->edge_end_nodes = new int[n];
	ret->btypes = new int[n];
	for (int i=0; i<n; ++i){
		ret->edge_start_nodes[i] = e1[i];
		ret->edge_end_nodes[i] = e2[i];
		ret->btypes[i] = bt[i];
	}
	return ret;
}
void free_grid2_boundary_types(Grid2DBoundaryStruct* p){
	delete[] p->edge_start_nodes;
	delete[] p->edge_end_nodes;
	delete[] p->btypes;
	delete p;
}

int export_msh_grid(const Grid* grid, const char* fname,
		const Grid2DBoundaryStruct* bstr,
		const BoundaryNamesStruct* bnames,
		int n_periodic,
		int* data_periodic){
	const GridGeom* g = togrid(grid);
	try {
		//1) construct boundary index vector
		std::set<Edge> eds_set = g->get_edges();
		std::vector<Edge> eds(eds_set.begin(), eds_set.end());
		std::vector<int> bindex(eds.size(), -1);
		for (int i=0; i<bstr->n; ++i){
			int p1 = bstr->edge_start_nodes[i];
			int p2 = bstr->edge_end_nodes[i];
			Edge e(p1, p2);
			auto fnd = std::find(eds.begin(), eds.end(), e);
			assert(fnd != eds.end());
			int index = fnd - eds.begin();
			bindex[index] = bstr->btypes[i];
		}
		//2) construct boundary names function
		std::map<int, std::string> bnames_map;
		for (int i=0; i<bnames->n; ++i){
			bnames_map[bnames->values[i]] = std::string(bnames->names[i]);
		}
		auto fnames = [&bnames_map](int i)->std::string{
				auto fnd = bnames_map.find(i);
				if (fnd == bnames_map.end()) return "boundary" + std::to_string(i);
				else return fnd->second;
			};
		//3 construct periodic data
		GGeom::Export::PeriodicData pd;
		for (int i=0; i<n_periodic; ++i){
			int b1 = *data_periodic++;
			int b2 = *data_periodic++;
			bool is_rev = (bool)(*data_periodic++);
			pd.add_data(b1, b2, is_rev);
		}
		//4 call function
		GGeom::Export::GridMSH(*g, fname, bindex, fnames, pd);
		return 0;
	} catch (std::runtime_error &e){
		std::cout<<e.what()<<std::endl;
		return 1;
	}
}
