#include "hmproject.h"
#include "cport_grid3d.h"
#include "grid.h"
#include "hmgrid3d.hpp"
#include "vtk_export_grid3d.hpp"
#include "tecplot_export_grid3d.hpp"

namespace{
HMGrid3D::Grid& hmgrid(CPortGrid3D* g){
	return *static_cast<HMGrid3D::Grid*>(g->grid);
}
const HMGrid3D::Grid& hmgrid(const CPortGrid3D* g){
	return *static_cast<const HMGrid3D::Grid*>(g->grid);
}
HMGrid3D::ExtendedSimpleSerialize& hmser(CPortGrid3D* g){
	return *static_cast<HMGrid3D::ExtendedSimpleSerialize*>(g->serialized);
}
const HMGrid3D::ExtendedSimpleSerialize& hmser(const CPortGrid3D* g){
	return *static_cast<const HMGrid3D::ExtendedSimpleSerialize*>(g->serialized);
}
}

void free_grid3d(CPortGrid3D* g){
	if (g->grid != NULL){
		delete static_cast<HMGrid3D::Grid*>(g->grid);
	}
	if (g->serialized != NULL){
		delete static_cast<HMGrid3D::ExtendedSimpleSerialize*>(g->serialized);
	}
	delete g;
}

namespace{
struct SideBndFunctor{
	std::vector<int> edge_bc;

	SideBndFunctor(const GridGeom& g, const Grid2DBoundaryStruct& bnd,
			int algo, int sidebnd){
		std::set<Edge> set_eds = g.get_edges();
		edge_bc.resize(set_eds.size(), sidebnd);
		if (algo == 0){
			return;
		} else if (algo == 1){
			std::vector<Edge> eds(set_eds.begin(), set_eds.end());
			for (int i=0; i<bnd.n; ++i){
				int p1 = bnd.edge_start_nodes[i];
				int p2 = bnd.edge_end_nodes[i];
				Edge e(p1, p2);
				auto gindex_fnd = std::find(eds.begin(), eds.end(), e);
				assert(gindex_fnd != eds.end());
				int gindex = gindex_fnd - eds.begin();
				edge_bc[gindex] = bnd.btypes[i];
			}
		} else {
			throw std::runtime_error("invalid side bc algorithm");
		}
	}
	int operator()(int i){
		assert(i>=0);
		return edge_bc[i];
	}
};
}
CPortGrid3D* grid2_sweep_z(const Grid* g, const Grid2DBoundaryStruct* bc,
		int nz, double* zvals,
		int algo_top, int* btop,
		int algo_bot, int* bbot,
		int algo_side, int bside){
	const GridGeom* gg = static_cast<const GridGeom*>(g);
	auto ret = new CPortGrid3D({0, 0});
	try{
		vector<double> z(zvals, zvals+nz);
		std::function<int(int)> botfun, topfun;
		switch (algo_top){
			case 0: topfun=[&btop](int){ return btop[0]; }; break;
			case 1: topfun=[&btop](int i){ return btop[i]; }; break;
			default: throw std::runtime_error("Invalid algo_top");
		}
		switch (algo_bot){
			case 0: botfun=[&bbot](int){ return bbot[0]; }; break;
			case 1: botfun=[&bbot](int i){ return bbot[i]; }; break;
			default: throw std::runtime_error("Invalid algo_bot");
		}
		auto sidefun = SideBndFunctor(*gg, *bc, algo_side, bside);
		ret->grid = new HMGrid3D::Grid(
			HMGrid3D::Constructor::SweepGrid2D(
				*gg, z, botfun, topfun, sidefun)
			);
		return ret;
	} catch (const std::runtime_error &e) {
		delete ret;
		std::cout<<e.what()<<std::endl;
		return NULL;
	}
}

CPortGrid3D* grid2_revolve(Grid* g, double* vec, int n_phi, double* phi,
		Grid2DBoundaryStruct* bc,
		int b1, int b2, int is_trian){
	CPortGrid3D* ret = NULL;
	//Scale
	GridGeom* grid = static_cast<GridGeom*>(g);
	ScaleBase sc = grid->do_scale();
	Point pstart(vec[0], vec[1]), pend(vec[2], vec[3]);
	sc.scale(pstart); sc.scale(pend);
	try{
		std::vector<double> vphi(phi, phi + n_phi);
		auto sidefun = SideBndFunctor(*grid, *bc, 1, 0);
		HMGrid3D::ESS ser = HMGrid3D::Constructor::RevolveGrid2D_S(
			*grid, vphi, pstart, pend, is_trian!=0, sidefun,
			[b1](int){ return b1; }, [b2](int){ return b2; });
		HMGrid3D::Vertex::Unscale2D(ser.data_vertex(), sc);
		HMGrid3D::Grid ans = ser.to_grid();
		ret = new CPortGrid3D;
		ret->serialized = new HMGrid3D::ESS(std::move(ser));
		ret->grid = new HMGrid3D::Grid(std::move(ans));
	} catch (const std::runtime_error &e){
		std::cout<<e.what()<<std::endl;
		if (ret != NULL) {free_grid3d(ret); ret = NULL; }
	}
	//Unscale
	grid->undo_scale(sc);
	//return
	return ret;
}

int export_vtk_grid3(const CPortGrid3D* grid, const char* fname, hmcport_callback f2){
	try{
		if (!grid->serialized){
			HMGrid3D::Export::GridVTK.WithCallback(f2, hmgrid(grid), fname);
			grid->serialized = new HMGrid3D::ExtendedSimpleSerialize(std::move(
						*HMGrid3D::Export::GridVTK.functor().last_ser_result
						));
		} else {
			HMGrid3D::Export::GridVTK.WithCallback(f2, hmser(grid), fname);
		}
		return 0;
	} catch (const std::runtime_error &e) {
		std::cout<<e.what()<<std::endl;
		return 1;
	}
}

int export_surface_vtk_grid3(const CPortGrid3D* grid, const char* fname, hmcport_callback f2){
	try{
		if (!grid->serialized){
			HMGrid3D::Export::BoundaryVTK.WithCallback(f2, hmgrid(grid), fname);
			grid->serialized = new HMGrid3D::ExtendedSimpleSerialize(std::move(
						*HMGrid3D::Export::BoundaryVTK.functor().last_ser_result
						));
		} else {
			HMGrid3D::Export::BoundaryVTK.WithCallback(f2, hmser(grid), fname);
		}
		return 0;
	} catch (const std::runtime_error &e) {
		std::cout<<e.what()<<std::endl;
		return 1;
	}

}

namespace{

HMGrid3D::Export::BFun construct_bnames(const BoundaryNamesStruct* bnames){
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

};

int export_msh_grid3(const CPortGrid3D* grid, const char* fname, const BoundaryNamesStruct* bnames,
		int n_periodic, double* data_periodic, hmcport_callback f2){
	try{
		// name function
		//std::map<int, std::string> bnames_map;
		//if (bnames != NULL) for (int i=0; i<bnames->n; ++i){
		//        bnames_map[bnames->values[i]] = std::string(bnames->names[i]);
		//}
		//auto nmfunc = [&bnames_map](int i)->std::string{
		//        auto fnd = bnames_map.find(i);
		//        if (fnd != bnames_map.end()) return fnd->second;
		//        else return "boundary" + std::to_string(i);
		//};
		auto nmfunc = construct_bnames(bnames);
		// building periodic
		HMGrid3D::Export::PeriodicData pd;
		for (int i=0; i<n_periodic; ++i){
			int b1 = (int)(*data_periodic++);
			int b2 = (int)(*data_periodic++);
			double x1 = *data_periodic++;
			double y1 = *data_periodic++;
			double z1 = *data_periodic++;
			double x2 = *data_periodic++;
			double y2 = *data_periodic++;
			double z2 = *data_periodic++;
			pd.add_condition(b1, b2, HMGrid3D::Vertex(x1, y1, z1), HMGrid3D::Vertex(x2, y2, z2), true);
		}
		//call function
		if (pd.size() == 0){
			HMGrid3D::Export::GridMSH.WithCallback(f2, hmgrid(grid), fname, nmfunc);
		} else {
			HMGrid3D::Export::GridMSH.WithCallback(f2, hmgrid(grid), fname, nmfunc, pd);
		}
		return 0;
	} catch (const std::runtime_error &e){
		std::cout<<e.what()<<std::endl;
		return 1;
	}
}

int export_tecplot_grid3(const CPortGrid3D* grid, const char* fname, const BoundaryNamesStruct* bnames,
		hmcport_callback f2){
	try{
		auto nmfunc = construct_bnames(bnames);
		HMGrid3D::Export::GridTecplot.WithCallback(f2, hmgrid(grid), fname, nmfunc);
		return 0;
	} catch (const std::runtime_error &e){
		std::cout<<e.what()<<std::endl;
		return 1;
	}
}
