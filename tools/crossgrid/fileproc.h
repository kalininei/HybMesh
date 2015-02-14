#ifndef CROSSGRID_FILEPROC_H
#define CROSSGRID_FILEPROC_H
#include "bgeom.h"
#include "grid.h"
#include "wireframegrid.h"

void save_vtk(const GridGeom* g, const char* fn);
void save_vtk(const PContour* c, const char* fn);
void save_vtk(const PtsGraph* g, const char* fn);




#endif
