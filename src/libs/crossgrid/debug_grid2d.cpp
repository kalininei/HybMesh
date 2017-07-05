#ifndef NDEBUG

#include "debug_grid2d.hpp"
#include "debug2d.hpp"
#include "buildgrid.hpp"
#include <fstream>
#include "nodes_compare.h"
#include "infogrid.hpp"

using namespace HM2D;

void HM2D::Grid::Debug::report_grid_problems(const GridData& g){
	return report_grid_problems(g, 2);
}

void HM2D::Grid::Debug::report_grid_problems(const GridData& g, double maxskew){
	const int DEF = -90943212;
	const int DEF2 = -90983212;
	aa::RestoreIds<CellData> r1(g.vcells);
	aa::RestoreIds<EdgeData> r2(g.vedges);
	aa::RestoreIds<VertexData> r3(g.vvert);

	std::cout<<"Analyzing grid: ";
	int Nv = g.vvert.size();
	int Ne = g.vedges.size();
	int Nc = g.vcells.size();
	std::cout<<Nc<<" cells, ";
	std::cout<<Ne<<" edges, ";
	std::cout<<Nv<<" vertices."<<std::endl;
	if (Nc == 0) return;

	//================== invalid pointer problems
	std::cout<<"--- Validity check."<<std::endl;
	for (int i=0; i<Nc; ++i){
		if (!g.vcells[i]){
			std::cout<<i<<" cell is nullptr"<<std::endl;
		} else for (int j=0; j<g.vcells[i]->edges.size(); ++j){
			if (!g.vcells[i]->edges[j]){
				std::cout<<i<<" cell "<<j<<" edge is nullptr"<<std::endl;
			}
		}
	}
	for (int i=0; i<Ne; ++i){
		if (!g.vedges[i]){
			std::cout<<i<<" edge is nullptr"<<std::endl;
		} else {
			if (!g.vedges[i]->first())
				std::cout<<i<<" edge first vertex is nullptr"<<std::endl;
			if (!g.vedges[i]->last())
				std::cout<<i<<" edge last vertex is nullptr"<<std::endl;
		}
	}

	//================== Cells dimension
	std::cout<<"--- Cells dimension check."<<std::endl;
	vector<int> cdims(Nc);
	for (int i=0; i<Nc; ++i) {
		cdims[i] = g.vcells[i]->edges.size();
		if (cdims[i]<3){
			std::cout<<i<<" cell contains "<<cdims[i]<<" edges"<<std::endl;
		}
	}
	
	//================== complexity
	std::cout<<"--- Complexity check."<<std::endl;
	aa::constant_ids_pvec(g.vcells, DEF);
	aa::constant_ids_pvec(g.vedges, DEF);
	aa::constant_ids_pvec(g.vvert, DEF);
	//cell->edges
	for (int i=0; i<Nc; ++i){
		auto c = g.vcells[i];
		for (int j=0; j<cdims[i]; ++j){
			auto e = c->edges[j];
			if (e->id != DEF){
				std::cout<<i<<" cell "<<j<<" edge is not in vedges."<<std::endl;
			}
		}
	}
	//edges->vertices, cells
	for (int i=0; i<Ne; ++i){
		auto e = g.vedges[i];
		if (e->pfirst()->id != DEF){
			std::cout<<"first vertex of "<<i<<" edge is not in vvert"<<std::endl;
		}
		if (e->plast()->id != DEF){
			std::cout<<"second vertex of "<<i<<" edge is not in vvert"<<std::endl;
		}
		if (e->has_left_cell() && e->left.lock()->id !=DEF){
			std::cout<<"left cell of "<<i<<" edge is not in vcells"<<std::endl;
		}
		if (e->has_right_cell() && e->right.lock()->id !=DEF){
			std::cout<<"right cell of "<<i<<" edge is not in vcells"<<std::endl;
		}
	}
	for (auto c: g.vcells){
		for (auto e: c->edges){
			e->id = DEF2;
			e->pfirst()->id = DEF2;
			e->plast()->id = DEF2;
		}
	}
	for (int i=0; i<Ne; ++i) if (g.vedges[i]->id == DEF){
		std::cout<<i<<" edge does not present in cell connections"<<std::endl;
	}
	for (int i=0; i<Nv; ++i) if (g.vvert[i]->id == DEF){
		std::cout<<i<<" vertex does not present in cell connections"<<std::endl;
	}

	//doubled vertices
	std::cout<<"--- Doubled vertices."<<std::endl;
	CoordinateMap2D<int> mp;
	for (int i=0; i<Nv; ++i){
		bool r = mp.add(*g.vvert[i], i);
		if (!r){
			int i2=mp.find(*g.vvert[i]).data();
			std::cout<<i<<" vertex equals "<<i2<<" vertex"<<std::endl;
		}
	}
	
	//edges length
	std::cout<<"--- Edges length."<<std::endl;
	std::vector<double> EL(Ne);
	for (int i=0; i<Ne; ++i){
		EL[i] = g.vedges[i]->length();
	}
	double maxlen = *max_element(EL.begin(), EL.end());
	for (int i=0; i<Ne; ++i){
		if (EL[i]/maxlen < 10*geps){
			std::cout<<"edge "<<i<<" is too short: "<<EL[i]<<" vs. "<<maxlen<<std::endl;
		}
	}

	//crossed edges
	std::cout<<"--- Crossed edges."<<std::endl;
	BoundingBox bb = HM2D::BBox(g.vedges);
	BoundingBoxFinder bbfinder(bb, bb.maxlen()/100);
	for (int i=0; i<Ne; ++i){
		auto e0 = g.vedges[i];
		auto bb2 = BoundingBox(*e0->pfirst(), *e0->plast());
		for (int j: bbfinder.suspects(bb2)){
			auto e1 = g.vedges[j];
			auto cr = SectCrossGeps(*e0->pfirst(), *e0->plast(),
					*e1->pfirst(), *e1->plast());
			if (cr.inner_cross()){
				std::cout<<"edges "<<i<<" and "<<j<<" have an intersection"<<std::endl;
			}
			if (cr.segment_cross()){
				std::cout<<"edges "<<i<<" and "<<j<<" are overlaid"<<std::endl;
			}
		}
		bbfinder.addentry(bb2);
	}
	
	//cell areas and rotation
	std::cout<<"--- Cell area and rotation."<<std::endl;
	vector<double> car(Nc, 0);
	for (int i=0; i<Nc; ++i){
		bool gc = HM2D::Contour::IsContour(g.vcells[i]->edges) &&
			HM2D::Contour::IsClosed(g.vcells[i]->edges);
		if (!gc) std::cout<<"cell "<<i<<" is not a closed contour"<<std::endl;
		else {
			car[i] = HM2D::Contour::Area(g.vcells[i]->edges);
			if (car[i]<geps*geps)
				std::cout<<"cell "<<i<<" area equals "<<car[i]<<std::endl;
		}
	}
	
	//edge->cell connectivity
	std::cout<<"--- Edge cell connectivity."<<std::endl;
	vector<vector<bool>> pres(Nc);
	for (int i=0; i<Nc; ++i) pres[i].resize(cdims[i], false);
	aa::enumerate_ids_pvec(g.vcells);
	for (int i=0; i<Ne; ++i){
		auto e0=g.vedges[i];
		if (e0->has_left_cell()){
			auto c = e0->left.lock();
			auto fnd = std::find(c->edges.begin(), c->edges.end(), e0);
			if (fnd == c->edges.end()){
				std::cout<<"left cell("<<c->id<<") of "<<i<<" edge doesn't contain it"<<std::endl;
			} else {
				int ifnd = fnd-c->edges.begin();
				pres[c->id][ifnd] = true;
				bool cd = HM2D::Contour::CorrectlyDirectedEdge(c->edges, ifnd);
				if ((car[c->id]>0 && !cd) || (car[c->id]<0 && cd)){
					std::cout<<i<<" edge right cell is marked as left"<<std::endl;
				}
			}
		}
		if (e0->has_right_cell()){
			auto c = e0->right.lock();
			auto fnd = std::find(c->edges.begin(), c->edges.end(), e0);
			if (fnd == c->edges.end()){
				std::cout<<"right cell("<<c->id<<") of "<<i<<" edge doesn't contain it"<<std::endl;
			} else {
				int ifnd = fnd-c->edges.begin();
				pres[c->id][ifnd] = true;
				bool cd = HM2D::Contour::CorrectlyDirectedEdge(c->edges, ifnd);
				if ((car[c->id]>0 && cd) || (car[c->id]<0 && !cd)){
					std::cout<<i<<" edge left cell is marked as right"<<std::endl;
				}
			}
		}
	}
	for (int i=0; i<Nc; ++i)
	for (int j=0; j<cdims[i]; ++j) if (!pres[i][j]){
		std::cout<<i<<" cell "<<j<<" edge doesn't present amoung its left or right connections"<<std::endl;
	}
	
	//skewness
	if (maxskew > 1) return;
	std::cout<<"--- Skewness"<<std::endl;
	vector<double> sk = HM2D::Grid::Skewness(g);
	for (int i=0; i<sk.size(); ++i) if (sk[i]>=maxskew){
		std::cout<<i<<" cell skewness equals "<<sk[i]<<std::endl;
	}
}



void HM2D::Grid::Debug::save_wf_vtk(const Grid::Impl::PtsGraph& c){
	auto ecol = c.toedges();
	HM2D::Debug::save_edges_vtk(ecol);

}

void HM2D::Grid::Debug::save_bbfinder_vtk(const BoundingBoxFinder& bf){
	Point p1 = bf.pmin();
	Point p2 = bf.pmax();
	auto rg = Grid::Constructor::RectGrid(p1, p2, bf.nx(), bf.ny());
	HM2D::Debug::save_grid_vtk(rg);
	std::ofstream fs("_dbgout.vtk", std::ios::app);
	fs<<"CELL_DATA "<<bf.nsqr()<<std::endl;
	fs<<"SCALARS entries_num int 1"<<std::endl;
	fs<<"LOOKUP_TABLE default"<<std::endl;
	for (int i=0; i<bf.nsqr(); ++i) fs<<bf.sqr_entries(i).size()<<std::endl;
}
void HM2D::Grid::Debug::save_bbfinder_vtk(const BoundingBoxFinder& bf, const vector<int>& dt){
	save_bbfinder_vtk(bf);
	std::ofstream fs("_dbgout.vtk", std::ios::app);
	fs<<"SCALARS data int 1"<<std::endl;
	fs<<"LOOKUP_TABLE default"<<std::endl;
	for (auto a: dt) fs<<a<<std::endl;
}

void HM2D::Grid::Debug::save_raster_vtk(const Finder::RasterizeEdges& rs){
	vector<double> x(rs.bbfinder().nx()+1);
	vector<double> y(rs.bbfinder().ny()+1);
	for (int i=0; i<x.size(); ++i){
		x[i] = rs.bbfinder().pmin().x + i*rs.bbfinder().stepx();
	}
	for (int i=0; i<y.size(); ++i){
		y[i] = rs.bbfinder().pmin().y + i*rs.bbfinder().stepy();
	}
	GridData g = HM2D::Grid::Constructor::RectGrid(x, y);

	vector<int> dt = rs.colour_squares(true);
	HM2D::Debug::save_grid_vtk(g);
	std::ofstream fs("_dbgout.vtk", std::ios::app);
	fs<<"CELL_DATA "<<g.vcells.size()<<std::endl;
	fs<<"SCALARS id int 1"<<std::endl;
	fs<<"LOOKUP_TABLE default"<<std::endl;
	for (int i=0; i<g.vcells.size(); ++i) fs<<dt[i]<<std::endl;
}

void HM2D::Grid::Debug::save_raster_vtk(shared_ptr<Finder::RasterizeEdges> rs){
	return save_raster_vtk(*rs);
}



#endif
