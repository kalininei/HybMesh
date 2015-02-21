#ifndef CROSSGRID_FILEPROC_H
#define CROSSGRID_FILEPROC_H
#include "bgeom.h"
#include "grid.h"
#include "wireframegrid.h"

void save_vtk(const GridGeom* g, const char* fn);
void save_vtk(const PContour* c, const char* fn);
void save_vtk(const vector<PContour>& c, const char* fn);
void save_vtk(const vector<PContour>& c, const vector<double>& data, const char* fn);
void save_vtk(const PtsGraph* g, const char* fn);

inline void save_vtk(const ContoursCollection& c, const char* fn){ save_vtk(c.contours_list(), fn); }
inline void save_vtk(const ContoursCollection* c, const char* fn){ save_vtk(c->contours_list(), fn); }




#endif
