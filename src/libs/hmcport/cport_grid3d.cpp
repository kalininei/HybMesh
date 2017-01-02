#include "hmproject.h"
#include "cport_grid3d.h"
#include "grid.h"
#include "hmgrid3d.hpp"
#include "bgeom3d.h"

using HMGrid3D::SGrid;

namespace{
HMGrid3D::SGrid& hmgrid(CPortGrid3D* g){
	return *static_cast<HMGrid3D::SGrid*>(g);
}
const HMGrid3D::SGrid& hmgrid(const CPortGrid3D* g){
	return *static_cast<const HMGrid3D::SGrid*>(g);
}
}

void free_grid3d(CPortGrid3D* g){
	if (g != NULL){
		delete static_cast<HMGrid3D::SGrid*>(g);
	}
}

void grid3_dims(const CPortGrid3D* g, int* nums){
	nums[0] = hmgrid(g).n_vert;
	nums[1] = hmgrid(g).n_edges;
	nums[2] = hmgrid(g).n_faces;
	nums[3] = hmgrid(g).n_cells;
}

void grid3_surface_dims(const CPortGrid3D* g, int* dims){
	dims[0] = hmgrid(g).bvert().size();
	dims[1] = hmgrid(g).bedges().size();
	dims[2] = hmgrid(g).bfaces().size();
}

int grid3_surface_btypes(const CPortGrid3D* g3, int* ret){
	try{
		auto g = hmgrid(g3);
		const auto& bf = g.bfaces();
		for (size_t i=0; i<bf.size(); ++i){
			ret[i] = g.btypes[bf[i]];
		}
	} catch (std::exception& e){
		std::cout<<e.what()<<std::endl;
		return 0;
	}
	return 1;
}

int grid3_volume(const CPortGrid3D* g, double* ret){
	try{
		auto& g3 = hmgrid(g);
		auto srf = HMGrid3D::Surface::GridSurface(g3);
		auto rv = HMGrid3D::GridSurfaceReverter(srf, true);
		*ret = HMGrid3D::Surface::Volume(srf);
	} catch (std::exception& e){
		std::cout<<e.what()<<std::endl;
		return 0;
	}
	return 1;
}

int grid3_merge(const CPortGrid3D* _g1, const CPortGrid3D* _g2, CPortGrid3D** gret){
	auto& g1 = hmgrid(_g1);
	auto& g2 = hmgrid(_g2);
	HMGrid3D::VertexData av(g1.vvert);
	std::copy(g2.vvert.begin(), g2.vvert.end(),
			std::back_inserter(av));
	ScaleBase3 sc = ScaleBase3::doscale(av);

	int res = 0;
	try{
		HMGrid3D::GridData g = HMGrid3D::MergeGrids(g1, g2);
		sc.unscale(g.vvert.begin(), g.vvert.end());
		*gret = new SGrid(std::move(g));
		res = 1;
	} catch (std::exception& e){
		std::cout<<e.what()<<std::endl;
	}

	sc.unscale(av.begin(), av.end());
	return res;
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
	HMGrid3D::SGrid* ret = NULL;
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
		HMGrid3D::GridData sg = HMGrid3D::Constructor::SweepGrid2D(
				*gg, z, botfun, topfun, sidefun);
		ret = new HMGrid3D::SGrid(std::move(sg));
		return ret;
	} catch (const std::exception &e) {
		if (ret != NULL) delete ret;
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
		HMGrid3D::GridData ser = HMGrid3D::Constructor::RevolveGrid2D(
			*grid, vphi, pstart, pend, is_trian!=0, sidefun,
			[b1](int){ return b1; }, [b2](int){ return b2; });
		HMGrid3D::Vertex::Unscale2D(ser.vvert, sc);
		ret = new HMGrid3D::SGrid(std::move(ser));
	} catch (const std::exception &e){
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
		HMGrid3D::Export::GridVTK.WithCallback(f2, hmgrid(grid), fname);
		return 0;
	} catch (const std::exception &e) {
		std::cout<<e.what()<<std::endl;
		return 1;
	}
}

int export_surface_vtk_grid3(const CPortGrid3D* grid, const char* fname, hmcport_callback f2){
	try{
		HMGrid3D::Export::BoundaryVTK.WithCallback(f2, hmgrid(grid), fname);
		return 0;
	} catch (const std::exception &e) {
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
	} catch (const std::exception &e){
		std::cout<<e.what()<<std::endl;
		return 1;
	}
}
int export_gmsh_grid3(const CPortGrid3D* grid, const char* fname, const BoundaryNamesStruct* bnames,
		hmcport_callback f2){
	try{
		// name function
		auto nmfunc = construct_bnames(bnames);
		//call function
		HMGrid3D::Export::GridGMSH.WithCallback(f2, hmgrid(grid), fname, nmfunc);
		return 0;
	} catch (const std::exception &e){
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
	} catch (const std::exception &e){
		std::cout<<e.what()<<std::endl;
		return 1;
	}
}

void* g3writer_create(const char* gname, void* grid, void* awriter, void* subnode,
		const char* fmt){
	try{
		HMXML::ReaderA* wr = static_cast<HMXML::ReaderA*>(awriter);
		HMXML::Reader* sn = static_cast<HMXML::Reader*>(subnode);
		SGrid* g = static_cast<SGrid*>(grid);
		auto ret = new HMGrid3D::Export::GridWriter(*g, wr, sn, gname, fmt);
		return ret;
	} catch (const std::exception &e){
		std::cout<<e.what()<<std::endl;
		return 0;
	}
}
void g3writer_free(void* gwriter){
	try{
		delete static_cast<HMGrid3D::Export::GridWriter*>(gwriter);
	} catch (const std::exception &e){
		std::cout<<e.what()<<std::endl;
	}
}
int g3writer_add_defined_field(void* gwriter, const char* field){
	try{
		auto gw = static_cast<HMGrid3D::Export::GridWriter*>(gwriter);
		std::string f(field);
		if (f == "face_vertices" || f == "face-vertices") gw->AddFaceVertexConnectivity();
		else if (f == "cell_faces" || f == "cell-faces") gw->AddCellFaceConnectivity();
		else if (f == "cell_vertices" || f == "cell-vertices") gw->AddCellVertexConnectivity();
		else if (f == "linfem") gw->AddLinFemConnectivity();
		else throw std::runtime_error("unknown field "+f);
		return 1;
	} catch (const std::exception &e){
		std::cout<<e.what()<<std::endl;
		return 0;
	}
}
void* g3reader_create(void* awriter, void* subnode, char* outname, hmcport_callback f2){
	try{
		auto wr = static_cast<HMXML::ReaderA*>(awriter);
		auto sn = static_cast<HMXML::Reader*>(subnode);
		//name
		std::string nm = sn->attribute(".", "name");
		if (nm.size()>1000) throw std::runtime_error("grid name is too long: " + nm);
		strcpy(outname, nm.c_str());
		//reader
		auto ret = HMGrid3D::Import::ReadHMG.WithCallback(f2, wr, sn);
		return ret.release();
	} catch (const std::exception &e){
		std::cout<<e.what()<<std::endl;
		return 0;
	}
}
CPortGrid3D* g3reader_getresult(void* rd){
	try{
		auto reader = static_cast<HMGrid3D::Import::GridReader*>(rd);
		return reader->result.release();
	} catch (const std::exception &e){
		std::cout<<e.what()<<std::endl;
		return 0;
	}
}
void g3reader_free(void* greader){
	try{
		delete (HMGrid3D::Import::GridReader*)greader;
	} catch (const std::exception &e){
		std::cout<<e.what()<<std::endl;
	}
}

int tetrahedral_fill(int nsurf, void** surf,
		int nconstr, void** constr,
		int npts, double* pcoords, double* psizes,
		void** ret,
		hmcport_callback cb){
	int r = 0;
	//copy constraint points
	HMGrid3D::VertexData cp;
	for (int i=0; i<npts; ++i){
		cp.emplace_back(new HMGrid3D::Vertex(pcoords[3*i], pcoords[3*i+1],
					pcoords[3*i+2]));
	}
	//shallow copy source surfaces
	HMGrid3D::Surface source;
	for (int i=0; i<nsurf; ++i){
		auto s = static_cast<HMGrid3D::SSurface*>(surf[i]);
		std::copy(s->faces.begin(), s->faces.end(),
				std::back_inserter(source.faces));
	}
	//shallow copy constraint surfaces
	HMGrid3D::Surface cs;
	for (int i=0; i<nconstr; ++i){
		auto s = static_cast<HMGrid3D::Surface*>(constr[i]);
		std::copy(s->faces.begin(), s->faces.end(),
				std::back_inserter(cs.faces));
	}
	//scaling
	HMGrid3D::VertexData srcpoints = source.allvertices();
	HMGrid3D::VertexData cspoints = cs.allvertices();
	ScaleBase3 sc = ScaleBase3::doscale(srcpoints);
	sc.scale(cspoints.begin(), cspoints.end());
	sc.scale(cp.begin(), cp.end());
	vector<double> ps(npts);
	for (int i=0; i<npts; ++i) ps[i] = psizes[i] / sc.L;
	try{
		HMGrid3D::GridData ans = HMGrid3D::Mesher::UnstructuredTetrahedral.WithCallback(
				cb, source, cs, cp, ps);
		sc.unscale(ans.vvert.begin(), ans.vvert.end());
		*ret = new SGrid(std::move(ans));
		r = 1;
	} catch (const std::exception &e){
		std::cout<<e.what()<<std::endl;
		r = 0;
	}
	//unscale input and return
	sc.unscale(srcpoints.begin(), srcpoints.end());
	sc.unscale(cspoints.begin(), cspoints.end());
	return r;
}
