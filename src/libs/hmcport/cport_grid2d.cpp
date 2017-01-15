#include "hmproject.h"
#include "cport_grid2d.h"
#include "hmmapping.hpp"
#include "trigrid.hpp"
#include "pebi.hpp"
#include "hmblay.hpp"
#include "cont_assembler.hpp"
#include "contabs2d.hpp"
#include "healgrid.hpp"
#include "finder2d.hpp"

#include "export2d_fluent.hpp"
#include "export2d_tecplot.hpp"
#include "export2d_hm.hpp"
#include "import2d_hm.hpp"


namespace{
HM2D::GridData* togrid(void* g){
	return static_cast<HM2D::GridData*>(g);
}
const HM2D::GridData* togrid(const void* g){
	return static_cast<const HM2D::GridData*>(g);
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
vector<int> construct_bindex(const Grid2DBoundaryStruct* bstr, const HM2D::GridData& g){
	auto vertedges = HM2D::Connectivity::VertexEdge(g.vedges);
	int k=0;
	for (auto& ve: vertedges) ve.v->id = k++;
	std::vector<int> bindex(g.vedges.size(), -1);
	aa::enumerate_ids_pvec(g.vedges);
	for (int i=0; i<bstr->n; ++i){
		int p1 = g.vvert[bstr->edge_start_nodes[i]]->id;
		int p2 = g.vvert[bstr->edge_end_nodes[i]]->id;
		HM2D::Edge* efnd = 0;
		for (auto ei: vertedges[p1].eind){
			auto e = g.vedges[ei];
			if ((e->first()->id == p1 && e->last()->id == p2) ||
			    (e->last()->id == p1 && e->first()->id == p2)){
				efnd = e.get(); break;
			}
		}
		if (efnd == 0) throw std::runtime_error("edge was not found");
		int index = efnd->id;
		bindex[index] = bstr->btypes[i];
	}
	return bindex;
}

HM2D::Export::BNamesFun construct_bnames(const BoundaryNamesStruct* bnames){
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
	const HM2D::GridData* g = togrid(grid);
	try {
		vector<int> bindex = construct_bindex(bstr, *g);
		auto fnames = construct_bnames(bnames);
		//3 construct periodic data
		HM2D::Export::PeriodicData pd;
		for (int i=0; i<n_periodic; ++i){
			int b1 = *data_periodic++;
			int b2 = *data_periodic++;
			bool is_rev = (bool)(*data_periodic++);
			pd.add_data(b1, b2, is_rev);
		}
		//4 call function
		HM2D::Export::GridMSH(*g, fname, bindex, fnames, pd);
		return 0;
	} catch (std::runtime_error &e){
		std::cout<<e.what()<<std::endl;
		return 1;
	}
}

int export_tecplot_grid(const Grid* grid, const char* fname,
		const Grid2DBoundaryStruct* bstr,
		const BoundaryNamesStruct* bnames){
	const HM2D::GridData* g = togrid(grid);
	try{
		vector<int> bindex = construct_bindex(bstr, *g);
		auto fnames = construct_bnames(bnames);
		//call function
		HM2D::Export::GridTecplot(*g, fname, bindex, fnames);
		return 0;
	} catch (std::runtime_error &e){
		std::cout<<e.what()<<std::endl;
		return 1;
	}
}

void* custom_rectangular_grid(int algo, void* left, void* bot,
		void* right, void* top, double* her_w, int return_invalid, hmcport_callback cb){
	//gather data
	HM2D::EdgeData* _left1 = static_cast<HM2D::EdgeData*>(left);
	HM2D::EdgeData* _bot1 = static_cast<HM2D::EdgeData*>(bot);
	HM2D::EdgeData* _right1 = static_cast<HM2D::EdgeData*>(right);
	HM2D::EdgeData* _top1 = static_cast<HM2D::EdgeData*>(top);
	HM2D::GridData* ret = 0;
	//scale
	ScaleBase sc = HM2D::Scale01(*_left1);
	HM2D::Scale(*_bot1, sc);
	HM2D::Scale(*_right1, sc);
	HM2D::Scale(*_top1, sc);
	try{
		HM2D::EdgeData left1 = HM2D::Contour::Assembler::SimpleContours(*_left1)[0];
		HM2D::EdgeData bot1 = HM2D::Contour::Assembler::SimpleContours(*_bot1)[0];
		HM2D::EdgeData right1 = HM2D::Contour::Assembler::SimpleContours(*_right1)[0];
		HM2D::EdgeData top1 = HM2D::Contour::Assembler::SimpleContours(*_top1)[0];
		//assemble grid
		if (algo == 0){
			ret = new HM2D::GridData(HMMap::LinearRectGrid(left1, bot1, right1, top1));
		} else if (algo == 1){
			ret = new HM2D::GridData(HMMap::LaplaceRectGrid.WithCallback(
				cb, left1, bot1, right1, top1, "inverse-laplace"));
		} else if (algo == 2){
			ret = new HM2D::GridData(HMMap::LaplaceRectGrid.WithCallback(
				cb, left1, bot1, right1, top1, "direct-laplace"));
		} else if (algo == 3){
			ret = new HM2D::GridData(HMMap::OrthogonalRectGrid.WithCallback(
				cb, left1, bot1, right1, top1));
		} else if (algo == 4){
			ret = new HM2D::GridData(HMMap::LinearTFIRectGrid(
				left1, bot1, right1, top1));
		} else if (algo == 5){
			ret = new HM2D::GridData(HMMap::CubicTFIRectGrid(
				left1, bot1, right1, top1, {her_w[0], her_w[1], her_w[2], her_w[3]}));
		} else throw std::runtime_error("unknown algorithm");
		HM2D::Unscale(ret->vvert, sc);
	} catch (HMMap::EInvalidGrid &e){
		if (!return_invalid){
			std::cout<<e.what()<<std::endl;
			ret = 0;
		} else{
			std::cout<<"Ignored error: "<<e.what()<<std::endl;
			HM2D::Unscale(e.invalid_grid.vvert, sc);
			ret = new HM2D::GridData(std::move(e.invalid_grid));
		}
	} catch (std::runtime_error &e){
		std::cout<<e.what()<<std::endl;
		ret = 0;
	}
	//unscale
	HM2D::Unscale(*_left1, sc);
	HM2D::Unscale(*_bot1, sc);
	HM2D::Unscale(*_right1, sc);
	HM2D::Unscale(*_top1, sc);
	return ret;
}

Grid* circ4grid(int algo, double* center, double rad, double step, double sqrside, double outer_refinement){
	HM2D::GridData* ret = NULL;
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
		ret = new HM2D::GridData(HMMap::Circ4Prototype(Point(0, 0), 1.0, 8*n1,
			stralgo, sqrside, outer_refinement));
		//unscale
		ScaleBase sc(center[0], center[1], rad);
		HM2D::Unscale(ret->vvert, sc);
		return ret;
	} catch (std::runtime_error &e){
		std::cout<<e.what()<<std::endl;
		if (ret!=0) delete ret;
		return NULL;
	}
}

void* triangulate_domain(void* domain, void* constr, int nemb, double* emb, int algo){
	using namespace HM2D::Mesher;
	HM2D::GridData* ret = NULL;
	HM2D::EdgeData* dom = static_cast<HM2D::EdgeData*>(domain);
	HM2D::EdgeData* con = static_cast<HM2D::EdgeData*>(constr);

	ScaleBase sc = HM2D::Scale01(*dom);
	if (con!=0) HM2D::Scale(*con, sc);
	std::map<Point, double> ep;
	for (int i=0; i<nemb; ++i){
		Point p(emb[3*i], emb[3*i+1]);
		sc.scale(p);
		ep.emplace(p, emb[3*i+2]/sc.L);
	}
	try{
		auto etree = HM2D::Contour::Tree::Assemble(*dom);
		if (etree.roots().size() == 0)
			throw std::runtime_error("Failed to find bounding contour");

		vector<HM2D::EdgeData> cas;
		if (con!=0) cas = HM2D::Contour::Assembler::AllContours(*con);
		for (auto& c: cas){ etree.add_detached_contour(c); }

		HM2D::GridData gret;
		switch (algo){
			case 0: gret = UnstructuredTriangle(etree, ep); break;
			case 1: gret = UnstructuredTriangleRecomb(etree, ep); break;
			default: throw std::runtime_error("unknown unstructured fill algorithm");
				
		};
		HM2D::Unscale(gret.vvert, sc);
		ret = new HM2D::GridData(std::move(gret));
	} catch (std::runtime_error &e){
		std::cout<<e.what()<<std::endl;
		if (ret!=0) delete ret;
		ret = NULL;
	}
	HM2D::Unscale(*dom, sc);
	if (con!=0) HM2D::Unscale(*con, sc);
	return ret;
}

void* pebi_fill(void* domain, void* constr, int nemb, double* emb){
	HM2D::GridData* ret = NULL;
	HM2D::EdgeData* dom = static_cast<HM2D::EdgeData*>(domain);
	HM2D::EdgeData* con = (constr==0)?0:static_cast<HM2D::EdgeData*>(constr);

	ScaleBase sc = HM2D::Scale01(*dom);
	if (con!=0) HM2D::Scale(*con, sc);
	std::map<Point, double> ep;
	for (int i=0; i<nemb; ++i){
		Point p(emb[3*i], emb[3*i+1]);
		sc.scale(p);
		ep.emplace(p, emb[3*i+2]/sc.L);
	}
	try{
		auto tree = HM2D::Contour::Tree::Assemble(*dom);

		if (tree.roots().size() == 0)
			throw std::runtime_error("Failed to find bounding contour");

		HM2D::GridData gret = HM2D::Mesher::UnstructuredPebi(tree, ep);
		HM2D::Unscale(gret.vvert, sc);
		ret = new HM2D::GridData(std::move(gret));
	} catch (std::runtime_error &e){
		std::cout<<e.what()<<std::endl;
		if (ret!=0) delete ret;
		ret = NULL;
	}
	HM2D::Unscale(*dom, sc);
	if (con!=0) HM2D::Unscale(*con, sc);
	return ret;
}

void* convex_cells(void* input_grid, double an){
	HM2D::GridData* ret=0;
	try{
		HM2D::GridData* inp = static_cast<HM2D::GridData*>(input_grid);
		ret = new HM2D::GridData();
		HM2D::DeepCopy(*inp, *ret);
		ScaleBase sc = HM2D::Scale01(ret->vvert);
		HM2D::Grid::Algos::NoConcaveCells(*ret, an);
		HM2D::Unscale(ret->vvert, sc);
		return ret;
	} catch (std::runtime_error &e){
		std::cout<<e.what()<<std::endl;
		if (ret!=0) delete ret;
		return NULL;
	}
}

void* stripe_grid(void* input_contour, int npart, double* part, int tip_algo, 
		void** bot, void** left, void** top, void** right, hmcport_callback cb){
	HM2D::GridData* ret = NULL;
	auto ecol = static_cast<HM2D::EdgeData*>(input_contour);
	ScaleBase sc = HM2D::Scale01(*ecol);
	HM2D::EdgeData cbot, cright, cleft, ctop;
	try{
		vector<double> spart(part, part+npart);
		for (auto& x: spart) x/=sc.L;
		Point bl, br, tr, tl;
		HM2D::EdgeData ic = HM2D::Contour::Assembler::Contour1(*ecol, ecol->at(0)->first().get());
		ret = new HM2D::GridData(HMBlay::BuildStripeGrid.WithCallback(cb, ic, spart, tip_algo, bl, br, tr, tl));

		auto gc = HM2D::Contour::Tree::GridBoundary(*ret);
		auto av = HM2D::AllVertices(gc.alledges());
		auto bl2 = av[std::get<0>(HM2D::Finder::ClosestPoint(av, bl))].get();
		auto br2 = av[std::get<0>(HM2D::Finder::ClosestPoint(av, br))].get();
		auto tl2 = av[std::get<0>(HM2D::Finder::ClosestPoint(av, tl))].get();
		auto tr2 = av[std::get<0>(HM2D::Finder::ClosestPoint(av, tr))].get();

		if (bl != br){
			auto& cont = gc.nodes[0]->contour;
			//if open
			cbot = HM2D::Contour::Assembler::Contour1(cont, bl2, br2);
			ctop = HM2D::Contour::Assembler::Contour1(cont, tr2, tl2);
			cleft = HM2D::Contour::Assembler::Contour1(cont, tl2, bl2);
			cright = HM2D::Contour::Assembler::Contour1(cont, br2, tr2);
		} else {
			auto& gout = gc.roots()[0]->contour;
			auto& gin = gc.roots()[0]->children[0].lock()->contour;
			cbot = HM2D::Contour::Assembler::Contour1(gin, bl2, br2);
			ctop = HM2D::Contour::Assembler::Contour1(gout, tr2, tl2);
		}
	} catch (std::runtime_error &e){
		std::cout<<e.what()<<std::endl;
		if (ret!=0) delete ret;
		return NULL;
	}
	HM2D::Unscale(ret->vvert, sc);

	auto s1 = new HM2D::EdgeData();
	HM2D::DeepCopy(cbot, *s1);
	auto s2 = new HM2D::EdgeData();
	HM2D::DeepCopy(cright, *s2);
	auto s3 = new HM2D::EdgeData();
	HM2D::DeepCopy(ctop, *s3);
	auto s4 = new HM2D::EdgeData();
	HM2D::DeepCopy(cleft, *s4);

	*bot = s1;
	*right =s2;
	*top = s3;
	*left = s4;
		
	return ret;
}

void* regular_hex_grid(double* area, int area_type, double cell_rad, int strict_area){
	using namespace HM2D::Grid::Constructor;
	try{
		ScaleBase sc;
		HM2D::GridData* ret = NULL;
		if (area_type == 0){
			sc = BoundingBox(area[0], area[1], area[2], area[3]).to_scale();
			ret = new HM2D::GridData(RegularHexagonal(
				Point(0, 0),
				Point((area[2]-area[0])/sc.L, (area[3]-area[1])/sc.L),
				cell_rad/sc.L, (bool)strict_area));
		} else if (area_type == 1){
			sc = ScaleBase(area[0], area[1], area[2]);
			ret = new HM2D::GridData(RegularHexagonal(
				Point(0, 0),
				area[2]/sc.L,
				cell_rad/sc.L,
				(bool)strict_area));
		} else {throw std::runtime_error("unknown algo");}
		HM2D::Unscale(ret->vvert, sc);
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
		HM2D::GridData* g = static_cast<HM2D::GridData*>(grid);
		auto ret = new HM2D::Export::GridWriter(*g, wr, sn, gname, fmt);
		return ret;
	} catch (std::runtime_error& e){
		std::cout<<e.what()<<std::endl;
		return NULL;
	}
}

void gwriter_free(void* gwriter){
	delete static_cast<HM2D::Export::GridWriter*>(gwriter);
}

int gwriter_add_defined_field(void* gwriter, const char* field){
	try{
		auto gw = static_cast<HM2D::Export::GridWriter*>(gwriter);
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
		auto gw = static_cast<HM2D::Export::GridWriter*>(gwriter);
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
		HM2D::Import::GridReader* ret = new HM2D::Import::GridReader(wr, sn);
		return ret;
	} catch (std::runtime_error& e){
		std::cout<<e.what()<<std::endl;
		return 0;
	}
}

void* greader_getresult(void* rd){
	auto reader = static_cast<HM2D::Import::GridReader*>(rd);
	return reader->result.release();
}

void* greader_read_edge_field(void* rd, const char* fieldname, const char* type){
	try{
		auto reader = static_cast<HM2D::Import::GridReader*>(rd);
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
	delete (HM2D::Import::GridReader*)greader;
}
