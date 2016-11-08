#ifndef CROSSGRID_PEBI_H
#define CROSSGRID_PEBI_H

#include "trigrid.h"

GridGeom TriToPebi(const TriGrid& g);
GridGeom RegularHexagonal(Point cnt, double area_rad, double cell_rad, bool strict_area=false);
GridGeom RegularHexagonal(Point cnt1, Point cnt2, double cell_rad, bool strict_area=false);




#endif
