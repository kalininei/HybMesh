' grid unification algorithm '
import os
import ctypes as ct
import globvars
import grid2
import dlgs
import bgeom


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
    #c types
    ct_pint = ct.POINTER(ct.c_int)
    ct_ppint = ct.POINTER(ct.POINTER(ct.c_int))
    ct_pd = ct.POINTER(ct.c_double)
    ct_ppd = ct.POINTER(ct.POINTER(ct.c_double))
    #c data allocation
    npt, ned, ncl = ct.c_int(), ct.c_int(), ct.c_int()
    pt, ed, cdims, ced = ct_pd(), ct_pint(), ct_pint(), ct_pint()
    #call c function
    lib_fa.grid_get_edges_info.argtypes = [ct.c_void_p, ct_pint, ct_pint,
            ct_pint, ct_ppd, ct_ppint, ct_ppint, ct_ppint]
    lib_fa.grid_get_edges_info(c_gr, ct.byref(npt), ct.byref(ned),
            ct.byref(ncl), ct.byref(pt), ct.byref(ed),
            ct.byref(cdims), ct.byref(ced))
    # ---- construct grid
    ret = grid2.Grid2()
    #points
    it = iter(pt)
    for i in range(npt.value):
        ret.points.append(bgeom.Point2(next(it), next(it)))
    #edges
    it = iter(ed)
    for i in range(ned.value):
        ret.edges.append([next(it), next(it)])
    #cells
    it1, it2 = iter(cdims), iter(ced)
    for i in range(ncl.value):
        ied = [next(it2) for j in range(next(it1))]
        ret.cells.append(ied)
    #free c memory
    lib_fa.grid_free_edges_info.argtypes = [ct_ppd, ct_ppint, ct_ppint,
            ct_ppint]
    lib_fa.grid_free_edges_info(ct.byref(pt), ct.byref(ed), ct.byref(cdims),
            ct.byref(ced))

    return ret

def unite_grids(g1, g2, buf, density, fix_bnd):
    'adds g2 to g1. Returns new grid'
    lib_fa = _clib()
    c_g1 = _grid_to_c(g1)
    c_g2 = _grid_to_c(g2)
    c_buf = ct.c_double(buf)
    c_den = ct.c_double(density / 10.0)
    c_fix = ct.c_int(1) if fix_bnd else ct.c_int(0)
    args = (c_g1, c_g2, c_buf, c_den, c_fix)

    lib_fa.cross_grids_wcb.restype = ct.c_void_p
    e = dlgs.ProgressProcedureDlg(lib_fa.cross_grids_wcb, args,
            globvars.mainWindow)
    e.exec_()
    c_cross = ct.c_void_p(e.get_result())

    #if result was obtained (no errors, no cancel)
    if (c_cross.value is not None):
        ret = _grid_from_c(c_cross)
        #lib_fa.grid_save_vtk(c_cross, "union_grid.vtk")
        lib_fa.grid_free(c_cross)
    else:
        ret = None

    lib_fa.grid_free(c_g1)
    lib_fa.grid_free(c_g2)

    return ret
