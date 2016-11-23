#include "hmproject.h"
#include "cport_grid2d.h"
#include "grid.h"
#include "hybmesh_contours2d.hpp"
#include "fluent_export_grid2d.hpp"
#include "tecplot_export_grid2d.hpp"
#include "procgrid.h"
#include "hmmapping.hpp"
#include "trigrid.h"
#include "hmblay.hpp"
#include "pebi.h"
#include "hmg_export_grid2d.hpp"
#include "hmg_import_grid2d.hpp"


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
			ret = new GridGeom(HMMap::LinearRectGrid(*left1, *bot1, *right1, *top1));
		} else if (algo == 1){
			ret = new GridGeom(HMMap::LaplaceRectGrid.WithCallback(
				cb, *left1, *bot1, *right1, *top1, "inverse-laplace"));
		} else if (algo == 2){
			ret = new GridGeom(HMMap::LaplaceRectGrid.WithCallback(
				cb, *left1, *bot1, *right1, *top1, "direct-laplace"));
		} else if (algo == 3){
			ret = new GridGeom(HMMap::OrthogonalRectGrid.WithCallback(
				cb, *left1, *bot1, *right1, *top1));
		} else if (algo == 4){
			ret = new GridGeom(HMMap::LinearTFIRectGrid(
				*left1, *bot1, *right1, *top1));
		} else if (algo == 5){
			ret = new GridGeom(HMMap::CubicTFIRectGrid(
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
		ret = new GridGeom(HMMap::Circ4Prototype(Point(0, 0), 1.0, 8*n1,
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

void* triangulate_domain(void* domain, void* constr, int nemb, double* emb, int algo){
	GridGeom* ret = NULL;
	HMCont2D::ECollection* dom = static_cast<HMCont2D::ECollection*>(domain);
	HMCont2D::ECollection* con = (constr==0)?0:static_cast<HMCont2D::ECollection*>(constr);

	ScaleBase sc = HMCont2D::ECollection::Scale01(*dom);
	if (con!=0) HMCont2D::ECollection::Scale(*con, sc);
	std::vector<double> ep;
	for (int i=0; i<nemb; ++i){
		ep.push_back((emb[3*i+0]-sc.p0.x)/sc.L);
		ep.push_back((emb[3*i+1]-sc.p0.y)/sc.L);
		ep.push_back(emb[3*i+2]/sc.L);
	}
	try{
		auto etree = HMCont2D::Assembler::ETree(*dom);
		auto tree = etree.ExtractTree(etree);
		if (tree.nodes.size() == 0)
			throw std::runtime_error("Failed to find bounding contour");

		ShpVector<HMCont2D::Contour> cc;
		vector<HMCont2D::Contour> cas;
		if (con!=0) cas = HMCont2D::Assembler::AllContours(*con);
		for (int i=0; i<cas.size(); ++i) aa::add_shared(cc, cas[i]);

		if (algo == 0) ret = new TriGrid(TriGrid(tree, cc, ep));
		else if (algo == 1) ret = new GridGeom(QuadGrid(tree, cc, ep));

		ret->undo_scale(sc);
	} catch (std::runtime_error &e){
		std::cout<<e.what()<<std::endl;
		if (ret!=0) delete ret;
		ret = NULL;
	}
	HMCont2D::ECollection::Unscale(*dom, sc);
	if (con!=0) HMCont2D::ECollection::Unscale(*con, sc);
	return ret;
}

void* pebi_fill(void* domain, void* constr, int nemb, double* emb){
	GridGeom* ret = NULL;
	HMCont2D::ECollection* dom = static_cast<HMCont2D::ECollection*>(domain);
	HMCont2D::ECollection* con = (constr==0)?0:static_cast<HMCont2D::ECollection*>(constr);

	ScaleBase sc = HMCont2D::ECollection::Scale01(*dom);
	if (con!=0) HMCont2D::ECollection::Scale(*con, sc);
	std::vector<double> ep;
	for (int i=0; i<nemb; ++i){
		ep.push_back((emb[3*i+0]-sc.p0.x)/sc.L);
		ep.push_back((emb[3*i+1]-sc.p0.y)/sc.L);
		ep.push_back(emb[3*i+2]/sc.L);
	}
	try{
		auto etree = HMCont2D::Assembler::ETree(*dom);
		auto tree = etree.ExtractTree(etree);

		if (tree.nodes.size() == 0)
			throw std::runtime_error("Failed to find bounding contour");

		ShpVector<HMCont2D::Contour> cc;
		vector<HMCont2D::Contour> cas;
		if (con!=0) cas = HMCont2D::Assembler::AllContours(*con);
		for (int i=0; i<cas.size(); ++i) aa::add_shared(cc, cas[i]);

		TriGrid g3 = TriGrid(tree, cc, ep);
		ret = new GridGeom(TriToPebi(g3));

		ret->undo_scale(sc);
	} catch (std::runtime_error &e){
		std::cout<<e.what()<<std::endl;
		if (ret!=0) delete ret;
		ret = NULL;
	}
	HMCont2D::ECollection::Unscale(*dom, sc);
	if (con!=0) HMCont2D::ECollection::Unscale(*con, sc);
	return ret;
}

void* convex_cells(void* input_grid, double an){
	GridGeom* ret = NULL;
	try{
		GridGeom* inp = static_cast<GridGeom*>(input_grid);
		ret = new GridGeom(GGeom::Constructor::DeepCopy(*inp));
		ScaleBase sc = ret->do_scale();
		GGeom::Repair::NoConcaveCells(*ret, an/180.0*M_PI);
		ret->undo_scale(sc);
		return ret;
	} catch (std::runtime_error &e){
		std::cout<<e.what()<<std::endl;
		if (ret!=0) delete ret;
		return NULL;
	}
}

void* stripe_grid(void* input_contour, int npart, double* part, int tip_algo, 
		void** bot, void** left, void** top, void** right, hmcport_callback cb){
	GridGeom* ret = NULL;
	auto ecol = static_cast<HMCont2D::ECollection*>(input_contour);
	ScaleBase sc = HMCont2D::ECollection::Scale01(*ecol);
	HMCont2D::Contour cbot, cright, cleft, ctop;
	try{
		vector<double> spart(part, part+npart);
		for (auto& x: spart) x/=sc.L;
		Point bl, br, tr, tl;
		HMCont2D::Contour ic = HMCont2D::Assembler::Contour1(*ecol, ecol->data[0]->pstart);
		ret = new GridGeom(HMBlay::BuildStripeGrid.WithCallback(cb, ic, spart, tip_algo, bl, br, tr, tl));

		if (bl != br){
			//if open
			auto gc = GGeom::Info::Contour1(*ret);
			cbot = HMCont2D::Assembler::Contour1(gc, bl, br);
			ctop = HMCont2D::Assembler::Contour1(gc, tr, tl);
			cleft = HMCont2D::Assembler::Contour1(gc, tl, bl);
			cright = HMCont2D::Assembler::Contour1(gc, br, tr);
		} else {
			auto gc = GGeom::Info::Contour(*ret);
			auto gout = gc.roots()[0];
			auto gin = gout->children[0];
			cbot = HMCont2D::Assembler::Contour1(*gin, bl, br);
			ctop = HMCont2D::Assembler::Contour1(*gout, tr, tl);
		}
	} catch (std::runtime_error &e){
		std::cout<<e.what()<<std::endl;
		if (ret!=0) delete ret;
		return NULL;
	}
	ret->undo_scale(sc);

	auto s1 = new HMCont2D::Container<HMCont2D::ECollection>();
	HMCont2D::Container<HMCont2D::ECollection>::DeepCopy(cbot, *s1);
	auto s2 = new HMCont2D::Container<HMCont2D::ECollection>();
	HMCont2D::Container<HMCont2D::ECollection>::DeepCopy(cright, *s2);
	auto s3 = new HMCont2D::Container<HMCont2D::ECollection>();
	HMCont2D::Container<HMCont2D::ECollection>::DeepCopy(ctop, *s3);
	auto s4 = new HMCont2D::Container<HMCont2D::ECollection>();
	HMCont2D::Container<HMCont2D::ECollection>::DeepCopy(cleft, *s4);

	*bot = s1;
	*right =s2;
	*top = s3;
	*left = s4;
		
	return ret;
}

void* regular_hex_grid(double* area, int area_type, double cell_rad, int strict_area){
	try{
		ScaleBase sc;
		GridGeom* ret = NULL;
		if (area_type == 0){
			sc = BoundingBox(area[0], area[1], area[2], area[3]).to_scale();
			ret = new GridGeom(RegularHexagonal(
				Point(0, 0),
				Point((area[2]-area[0])/sc.L, (area[3]-area[1])/sc.L),
				cell_rad/sc.L, (bool)strict_area));
		} else if (area_type == 1){
			sc = ScaleBase(area[0], area[1], area[2]);
			ret = new GridGeom(RegularHexagonal(
				Point(0, 0),
				area[2]/sc.L,
				cell_rad/sc.L,
				(bool)strict_area));
		} else {throw std::runtime_error("unknown algo");}
		ret->undo_scale(sc);
		return ret;
	} catch (std::runtime_error& e){
		std::cout<<e.what()<<std::endl;
		return NULL;
	}
}

void* gwriter_create(const char* gname, void* grid, void* awriter, void* subnode, const char* fmt){
	try{
		HMXML::ReaderA* wr = static_cast<HMXML::ReaderA*>(awriter);
		HMXML::Reader* sn = static_cast<HMXML::ReaderA*>(subnode);
		GridGeom* g = static_cast<GridGeom*>(grid);
		auto ret = new GGeom::Export::GridWriter(*g, wr, sn, gname, fmt);
		return ret;
	} catch (std::runtime_error& e){
		std::cout<<e.what()<<std::endl;
		return NULL;
	}
}

void gwriter_free(void* gwriter){
	delete static_cast<GGeom::Export::GridWriter*>(gwriter);
}

int gwriter_add_defined_field(void* gwriter, const char* field){
	try{
		auto gw = static_cast<GGeom::Export::GridWriter*>(gwriter);
		std::string f(field);
		if (f == "cell_edges" || f == "cell-edges") gw->AddCellEdgeConnectivity();
		else if (f == "cell_vertices" || f == "cell-vertices") gw->AddCellVertexConnectivity();
		else throw std::runtime_error("unknown field "+f);
		return 1;
	} catch (std::runtime_error& e){
		std::cout<<e.what()<<std::endl;
		return 0;
	}
}
int gwriter_add_edge_field(void* gwriter, const char* fieldname, void* field, int fsize, const char* type){
	try{
		auto gw = static_cast<GGeom::Export::GridWriter*>(gwriter);
		std::string fn(fieldname);
		std::string tpstr(type);
		if (tpstr == "int"){
			int* intfield = static_cast<int*>(field);
			std::vector<int> data(intfield, intfield+fsize);
			gw->AddEdgeData(fn, data, gw->is_binary<int>());
		} else if (tpstr == "char"){
			char* chrfield = static_cast<char*>(field);
			std::vector<char> data(chrfield, chrfield+fsize);
			gw->AddEdgeData(fn, data, gw->is_binary<char>());
		} else if (tpstr == "double"){
			double* chrfield = static_cast<double*>(field);
			std::vector<double> data(chrfield, chrfield+fsize);
			gw->AddEdgeData(fn, data, gw->is_binary<double>());
		} else if (tpstr == "float"){
			float* chrfield = static_cast<float*>(field);
			std::vector<float> data(chrfield, chrfield+fsize);
			gw->AddEdgeData(fn, data, gw->is_binary<float>());
		} else {
			throw std::runtime_error("unknown data type "+tpstr);
		}
		return 1;
	} catch (std::runtime_error& e){
		std::cout<<e.what()<<std::endl;
		return 0;
	}
}

void* greader_create(void* awriter, void* subnode, char* outname){
	try{
		auto wr = static_cast<HMXML::ReaderA*>(awriter);
		auto sn = static_cast<HMXML::Reader*>(subnode);
		//name
		std::string nm = sn->attribute(".", "name");
		if (nm.size()>1000) throw std::runtime_error("grid name is too long: " + nm);
		strcpy(outname, nm.c_str());
		//reader
		GGeom::Import::GridReader* ret = new GGeom::Import::GridReader(wr, sn);
		return ret;
	} catch (std::runtime_error& e){
		std::cout<<e.what()<<std::endl;
		return 0;
	}
}

void* greader_getresult(void* rd){
	auto reader = static_cast<GGeom::Import::GridReader*>(rd);
	return reader->result.release();
}

void* greader_read_edge_field(void* rd, const char* fieldname, const char* type){
	try{
		auto reader = static_cast<GGeom::Import::GridReader*>(rd);
		std::string fname(fieldname);
		std::string tpname(type);
		void* ret = NULL;
		if (tpname == "int"){
			auto outv = reader->read_edges_field<int>(fieldname);
			if (outv.size() != reader->Ne) throw std::runtime_error("field size doesn't match");
			ret = new int[outv.size()];
			std::copy(outv.begin(), outv.end(), (int*)ret);
		} else if (tpname == "char"){
			auto outv = reader->read_edges_field<char>(fieldname);
			if (outv.size() != reader->Ne) throw std::runtime_error("field size doesn't match");
			ret = new char[outv.size()];
			std::copy(outv.begin(), outv.end(), (char*)ret);
		} else if (tpname == "float"){
			auto outv = reader->read_edges_field<float>(fieldname);
			if (outv.size() != reader->Ne) throw std::runtime_error("field size doesn't match");
			ret = new float[outv.size()];
			std::copy(outv.begin(), outv.end(), (float*)ret);
		} else if (tpname == "double"){
			auto outv = reader->read_edges_field<double>(fieldname);
			if (outv.size() != reader->Ne) throw std::runtime_error("field size doesn't match");
			ret = new double[outv.size()];
			std::copy(outv.begin(), outv.end(), (double*)ret);
		}
		return ret;
	} catch (std::runtime_error& e){
		return NULL;
	}
}

void greader_free(void* greader){
	delete (GGeom::Import::GridReader*)greader;
}
