#include "hybmesh_contours2d.hpp"
#include <fstream>

using namespace HMCont2D;

void HMCont2D::SaveVtk(const ContourTree& tree, const char* fname){
	std::ofstream fs(fname);
	fs<<"# vtk DataFile Version 3.0"<<std::endl;
	fs<<"Contour 2D"<<std::endl;
	fs<<"ASCII"<<std::endl;
	//Points
	int numpts = tree.NumPoints();
	fs<<"DATASET UNSTRUCTURED_GRID"<<std::endl;
	fs<<"POINTS "<<numpts<< " float"<<std::endl;
	std::map<const Point*, int> point_ind;
	int n=0;
	for (int i=0; i<tree.NumContours(); ++i){
		auto c = tree.Cont(i);
		for (int j=0; j<c->NumPoints();++j){
			auto p = c->Pnt(j);
			fs<<p->x<<" "<<p->y<<" 0"<<std::endl;
			point_ind[p] = n++;
		}
	}
	//Edges
	int numedges = tree.NumEdges();
	//Cells
	fs<<"CELLS  "<<numedges<<"   "<<3*numedges<<std::endl;
	vector<int> edges_shift; edges_shift.push_back(0);
	for (int i=0; i<tree.NumContours(); ++i){
		auto c = tree.Cont(i);
		for (int j=0; j<c->NumEdges(); ++j){
			auto ed = c->Edge(j);
			fs<<2<<" "<<point_ind[std::get<0>(ed)]<<" "
				<<point_ind[std::get<1>(ed)]<<std::endl;
		}
	}
	fs<<"CELL_TYPES  "<<numedges<<std::endl;
	for (int i=0;i<numedges;++i) fs<<3<<std::endl;
	fs.close();
}
