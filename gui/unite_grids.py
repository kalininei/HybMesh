' grid unification algorithm '
import os
import ctypes as ct
import grid2


def _clib():
    dn = os.path.dirname(__file__)
    libfn = os.path.join(dn, '../tools/bin/libcrossgrid.so')
    return ct.cdll.LoadLibrary(libfn)


def _grid_to_c(g):
    " return c pointer to a python grid "
    lib_fa = _clib()
    npt = g.n_points()
    ncls = g.n_cells()
    #fill points
    c_pnt = (ct.c_double * (2 * npt))()
    for i, p in enumerate(g.points):
        c_pnt[2 * i] = p.x
        c_pnt[2 * i + 1] = p.y
    #fill cells
    cls_a = g.cells_nodes_connect()
    cls_a_num = len(cls_a)
    for ce in cls_a:
        cls_a_num += len(ce)
    c_cls = (ct.c_int * cls_a_num)()
    ind = 0
    for ce in cls_a:
        c_cls[ind] = len(ce)
        for i in range(1, len(ce) + 1):
            c_cls[ind + i] = ce[i - 1]
        ind += (len(ce) + 1)
    return lib_fa.grid_construct(ct.c_int(npt), ct.c_int(ncls), c_pnt, c_cls)


def _grid_from_c(c_gr):
    " builds a grid from a c object"
    lib_fa = _clib()
    npt = lib_fa.grid_npoints(c_gr)
    ncl = lib_fa.grid_ncells(c_gr)
    ncldim = lib_fa.grid_cellsdim(c_gr)
    c_pnt = (ct.c_double * (2 * npt))()
    c_cls = (ct.c_int * (ncldim + ncl))()
    lib_fa.grid_get_points_cells(c_gr, c_pnt, c_cls)
    pnt, cls = [], []
    for i in range(len(c_pnt)):
        pnt.append(c_pnt[i])
    for i in range(len(c_cls)):
        cls.append(c_cls[i])
    return grid2.Grid2.from_points_cells(npt, ncl, pnt, cls)


def unite_grids(g1, g2, buf, density):
    'adds g2 to g1. Returns new grid'
    lib_fa = _clib()
    c_g1 = _grid_to_c(g1)
    c_g2 = _grid_to_c(g2)
    c_cross = lib_fa.cross_grids(c_g1, c_g2, ct.c_double(buf))

    ret = _grid_from_c(c_cross)

    #lib_fa.grid_save_vtk(c_cross, "union_grid.vtk")

    #free c memory
    lib_fa.grid_free(c_g1)
    lib_fa.grid_free(c_g2)
    lib_fa.grid_free(c_cross)

    return ret
