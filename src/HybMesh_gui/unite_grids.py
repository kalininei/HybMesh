' grid unification algorithm '
import ctypes as ct
import globvars
import grid2
import dlgs
import bgeom


def _clib():
    libfn = globvars.prog_options.lib_crossgrid
    return ct.cdll.LoadLibrary(libfn)


def _grid_to_c(g):
    " return c pointer to a python grid"
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


def _cont_to_c(c):
    "return c pointer to a contours2.Contour2 object"
    lib_fa = _clib()
    #fill points
    c_pnt = (ct.c_double * (2 * c.n_points()))()
    for i, p in enumerate(c.points):
        c_pnt[2 * i] = p.x
        c_pnt[2 * i + 1] = p.y
    #fill edges
    edg = c.edges_points()
    c_edg = (ct.c_int * (2 * c.n_edges()))()
    for i, ce in enumerate(edg):
        c_edg[2 * i] = ce[0]
        c_edg[2 * i + 1] = ce[1]

    return lib_fa.contour_construct(ct.c_int(c.n_points()),
            ct.c_int(c.n_edges()), c_pnt, c_edg)


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


def add_bc_from_cont(tar_cont, src_cont, c_src=None, c_tar=None):
    """ (contour2.Contour2, contour2.Contour2) -> None

        Adds boundary condition to tar_cont from src_cont
        if grid boundary edge lies on the contour edge.
        Changes grd.cont data.

        Helper c_src and c_tar are the c representation of
        target and source contours of ct.c_void_p type.
    """
    lib_fa = _clib()

    #create c representations of contours
    need_src_clear, need_tar_clear = False, False
    if c_src is None:
        need_src_clear = True
        c_src = _cont_to_c(src_cont)

    if c_tar is None:
        need_tar_clear = True
        c_tar = _cont_to_c(tar_cont)

    #create contours boundary array
    c_src_bnd = (ct.c_int * src_cont.n_edges())()  # input
    c_tar_bnd = (ct.c_int * tar_cont.n_edges())()  # output
    for i in range(src_cont.n_edges()):
        c_src_bnd[i] = src_cont.edge_bnd(i)

    #call lib_fa function
    lib_fa.add_contour_bc(c_src, c_tar, c_src_bnd, c_tar_bnd, ct.c_int(-1))

    #write data to target contour
    new_bnd = {}
    for i in range(tar_cont.n_edges()):
        b = c_tar_bnd[i]
        if b >= 0:
            new_bnd[i] = b
    tar_cont.add_edge_bnd(new_bnd)

    #clear c data
    if need_src_clear:
        lib_fa.cont_free(c_src)
    if need_tar_clear:
        lib_fa.cont_free(c_tar)


def setbc_from_conts(tar_cont, src):
    """ (contour2.Contour2, [contour2.Contour2]) -> None
        Sets boundary condition to tar_cont from list of source contours.
    """
    lib_fa = _clib()
    #set default zero to all target boundaries
    tar_cont.set_edge_bnd({i: 0 for i in range(tar_cont.n_edges())})

    #loop over each source contours
    c_tar = _cont_to_c(tar_cont)
    for cont in src:
        c_src = _cont_to_c(cont)
        add_bc_from_cont(tar_cont, cont, c_src, c_tar)
        lib_fa.cont_free(c_src)
    lib_fa.cont_free(c_tar)


def unite_grids(g1, g2, buf, fix_bnd, empty_holes):
    'adds g2 to g1. Returns new grid'
    lib_fa = _clib()
    c_g1 = _grid_to_c(g1)
    c_g2 = _grid_to_c(g2)
    c_buf = ct.c_double(buf)
    c_fix = ct.c_int(1) if fix_bnd else ct.c_int(0)
    c_eh = ct.c_int(1) if empty_holes else ct.c_int(0)
    args = (c_g1, c_g2, c_buf, c_fix, c_eh)

    lib_fa.cross_grids_wcb.restype = ct.c_void_p
    e = dlgs.ProgressProcedureDlg(lib_fa.cross_grids_wcb, args,
            True, globvars.mainWindow)
    e.exec_()
    c_cross = ct.c_void_p(e.get_result())

    #if result was obtained (no errors, no cancel)
    if (c_cross.value is not None):
        ret = _grid_from_c(c_cross)
        lib_fa.grid_free(c_cross)
    else:
        ret = None

    lib_fa.grid_free(c_g1)
    lib_fa.grid_free(c_g2)

    return ret


def grid_excl_cont(grd, cnt, is_inner):
    """ ->grid2.Grid2
        Returns a grid with excluded contour area (inner or outer)
    """
    lib_fa = _clib()
    c_g = _grid_to_c(grd)
    c_c = _cont_to_c(cnt)
    c_isinner = ct.c_int(1 if is_inner else 0)
    args = (c_g, c_c, c_isinner)
    lib_fa.grid_exclude_cont_wcb.restype = ct.c_void_p
    e = dlgs.ProgressProcedureDlg(lib_fa.grid_exclude_cont_wcb, args,
            True, globvars.mainWindow)
    e.exec_()
    res = ct.c_void_p(e.get_result())
    if res.value is not None:
        newg = _grid_from_c(res)
        newg.build_contour()

        #assign boundary conditions
        bs = cnt.bnd_types().union(grd.cont.bnd_types()).difference(set([0]))
        if len(bs) > 0:
            res_cont = _cont_to_c(newg.cont)
            #1. from source grid
            add_bc_from_cont(newg.cont, grd.cont, c_tar=res_cont)
            #2. from contour
            add_bc_from_cont(newg.cont, cnt, c_tar=res_cont, c_src=c_c)
            #3. free contour memory
            lib_fa.cont_free(res_cont)

        #free grid memory
        lib_fa.grid_free(res)
    else:
        newg = None

    lib_fa.cont_free(c_c)
    lib_fa.grid_free(c_g)

    return newg
