#include "vtk_export2d.hpp"
#include <fstream>

using namespace HM2D;

void Export::ContourVTK(const EdgeData& ec, std::string fn){
	std::ofstream fs(fn);
	fs<<"# vtk DataFile Version 3.0"<<std::endl;
	fs<<"EdgesData"<<std::endl;
	fs<<"ASCII"<<std::endl;
	VertexData pts = AllVertices(ec);
	aa::enumerate_ids_pvec(pts);
	//Points
	fs<<"DATASET UNSTRUCTURED_GRID"<<std::endl;
	fs<<"POINTS "<<pts.size()<< " float"<<std::endl;
	for (int i=0; i<pts.size(); ++i){
		fs<<float(pts[i]->x)<<" "<<float(pts[i]->y)<<" 0"<<std::endl;
	}
	//Cells
	fs<<"CELLS  "<<ec.size()<<"   "<<3*ec.size()<<std::endl;
	for (int i=0; i<ec.size(); ++i){
		fs<<2<<" "<<ec[i]->first()->id<<" "<<ec[i]->last()->id<<std::endl;
	}
	fs<<"CELL_TYPES  "<<ec.size()<<std::endl;
	for (int i=0;i<ec.size();++i) fs<<3<<std::endl;
	fs.close();
}

void Export::VerticesVTK(const VertexData& data, std::string fn){
	std::ofstream fs(fn);
	fs<<"# vtk DataFile Version 3.0"<<std::endl;
	fs<<"VertexData"<<std::endl;
	fs<<"ASCII"<<std::endl;
	//Points
	fs<<"DATASET UNSTRUCTURED_GRID"<<std::endl;
	fs<<"POINTS "<<data.size()<< " float"<<std::endl;
	for (int i=0; i<data.size(); ++i){
		fs<<float(data[i]->x)<<" "<<float(data[i]->y)<<" 0"<<std::endl;
	}
	//Cells
	fs<<"CELLS  "<<data.size()<<"   "<<2*data.size()<<std::endl;
	for (int i=0; i<data.size(); ++i){
		fs<<"1 "<<i<<std::endl;
	}
	fs<<"CELL_TYPES  "<<data.size()<<std::endl;
	for (int i=0;i<data.size();++i) fs<<1<<std::endl;
	fs.close();
}
