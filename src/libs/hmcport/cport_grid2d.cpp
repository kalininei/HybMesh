#include "hmproject.h"
#include "cport_grid2d.h"
#include "grid.h"
#include "hybmesh_contours2d.hpp"
#include "fluent_export_grid2d.hpp"
#include "tecplot_export_grid2d.hpp"
#include "procgrid.h"
#include "hmmapping.hpp"


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

namespace{
vector<int> construct_bindex(const Grid2DBoundaryStruct* bstr, const GridGeom& g){
	std::set<Edge> eds_set = g.get_edges();
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
	return bindex;
}

GGeom::Export::BFun construct_bnames(const BoundaryNamesStruct* bnames){
	std::map<int, std::string> bnames_map;
	if (bnames != NULL) for (int i=0; i<bnames->n; ++i){
		bnames_map[bnames->values[i]] = std::string(bnames->names[i]);
	}
	auto fnames = [bnames_map](int i)->std::string{
			auto fnd = bnames_map.find(i);
			if (fnd == bnames_map.end()) return "boundary" + std::to_string(i);
			else return fnd->second;
		};
	return fnames;
}
}

int export_msh_grid(const Grid* grid, const char* fname,
		const Grid2DBoundaryStruct* bstr,
		const BoundaryNamesStruct* bnames,
		int n_periodic,
		int* data_periodic){
	const GridGeom* g = togrid(grid);
	try {
		vector<int> bindex = construct_bindex(bstr, *g);
		auto fnames = construct_bnames(bnames);
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

int export_tecplot_grid(const Grid* grid, const char* fname,
		const Grid2DBoundaryStruct* bstr,
		const BoundaryNamesStruct* bnames){
	const GridGeom* g = togrid(grid);
	try{
		vector<int> bindex = construct_bindex(bstr, *g);
		auto fnames = construct_bnames(bnames);
		//call function
		GGeom::Export::GridTecplot(*g, fname, bindex, fnames);
		return 0;
	} catch (std::runtime_error &e){
		std::cout<<e.what()<<std::endl;
		return 1;
	}
}

void* custom_rectangular_grid(int algo, void* left, void* bot,
		void* right, void* top, double* her_w, hmcport_callback cb){
	//gather data
	HMCont2D::Contour* left1 = static_cast<HMCont2D::Contour*>(left);
	HMCont2D::Contour* bot1 = static_cast<HMCont2D::Contour*>(bot);
	HMCont2D::Contour* right1 = static_cast<HMCont2D::Contour*>(right);
	HMCont2D::Contour* top1 = static_cast<HMCont2D::Contour*>(top);
	GridGeom* ret = 0;
	//scale
	ScaleBase sc = HMCont2D::ECollection::Scale01(*left1);
	HMCont2D::ECollection::Scale(*bot1, sc);
	HMCont2D::ECollection::Scale(*right1, sc);
	HMCont2D::ECollection::Scale(*top1, sc);
	try{
		//assemble grid
		if (algo == 0){
			ret = new GridGeom(HMGMap::LinearRectGrid(*left1, *bot1, *right1, *top1));
		} else if (algo == 1){
			ret = new GridGeom(HMGMap::LaplaceRectGrid.WithCallback(
				cb, *left1, *bot1, *right1, *top1, "inverse-laplace"));
		} else if (algo == 2){
			ret = new GridGeom(HMGMap::LaplaceRectGrid.WithCallback(
				cb, *left1, *bot1, *right1, *top1, "direct-laplace"));
		} else if (algo == 3){
			ret = new GridGeom(HMGMap::OrthogonalRectGrid.WithCallback(
				cb, *left1, *bot1, *right1, *top1));
		} else if (algo == 4){
			ret = new GridGeom(HMGMap::LinearTFIRectGrid(
				*left1, *bot1, *right1, *top1));
		} else if (algo == 5){
			ret = new GridGeom(HMGMap::CubicTFIRectGrid(
				*left1, *bot1, *right1, *top1, {her_w[0], her_w[1], her_w[2], her_w[3]}));
		} else throw std::runtime_error("unknown algorithm");
		ret->undo_scale(sc);
	} catch (std::runtime_error &e){
		std::cout<<e.what()<<std::endl;
		ret = 0;
	}
	//unscale
	HMCont2D::ECollection::Unscale(*left1, sc);
	HMCont2D::ECollection::Unscale(*bot1, sc);
	HMCont2D::ECollection::Unscale(*right1, sc);
	HMCont2D::ECollection::Unscale(*top1, sc);
	return ret;
}

Grid* circ4grid(int algo, double* center, double rad, double step, double sqrside, double outer_refinement){
	GridGeom* ret = NULL;
	try{
		double n = 2*M_PI*rad/step;
		int n1 = round(n/8.0);
		std::string stralgo;
		switch (algo){
			case 0: stralgo="linear"; break;
			case 1: stralgo="laplace"; break;
			case 2: stralgo="orthogonal-circ"; break;
			case 3: stralgo="orthogonal-rect"; break;
			default: throw std::runtime_error("unknown algorithm");
		};
		ret = new GridGeom(HMGMap::Circ4Prototype(Point(0, 0), 1.0, 8*n1,
			stralgo, sqrside, outer_refinement));
		//unscale
		ScaleBase sc(center[0], center[1], rad);
		ret->undo_scale(sc);
		return ret;
	} catch (std::runtime_error &e){
		std::cout<<e.what()<<std::endl;
		if (ret!=0) delete ret;
		return NULL;
	}
}
