' grid unification algorithm '
import ctypes as ct
import cobj
import basic.interf


def add_bc_from_cont(tar_cont, src_cont, c_src=None, c_tar=None):
    """ (contour2.Contour2, contour2.Contour2) -> None

        Adds boundary condition to tar_cont from src_cont
        if grid boundary edge lies on the contour edge.
        Changes grd.cont data.

        Helper c_src and c_tar are the c representation of
        target and source contours of ct.c_void_p type.
    """
    lib_fa = cobj.crossgrid_lib()

    #create c representations of contours
    need_src_clear, need_tar_clear = False, False
    if c_src is None:
        need_src_clear = True
        c_src = cobj.cont_to_c(src_cont)

    if c_tar is None:
        need_tar_clear = True
        c_tar = cobj.cont_to_c(tar_cont)

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
    lib_fa = cobj.crossgrid_lib()
    #set default zero to all target boundaries
    tar_cont.set_edge_bnd({i: 0 for i in range(tar_cont.n_edges())})

    #loop over each source contours
    c_tar = cobj.cont_to_c(tar_cont)
    for cont in src:
        c_src = cobj.cont_to_c(cont)
        add_bc_from_cont(tar_cont, cont, c_src, c_tar)
        lib_fa.cont_free(c_src)
    lib_fa.cont_free(c_tar)


def unite_grids(g1, g2, buf, fix_bnd, empty_holes, cb=None):
    """ adds g2 to g1. Returns new grid. cb -- callback object
        inherited from inter.SilentCallback2PG'
    """
    lib_fa = cobj.crossgrid_lib()
    c_g1 = cobj.grid_to_c(g1)
    c_g2 = cobj.grid_to_c(g2)
    c_buf = ct.c_double(buf)
    c_fix = ct.c_int(1) if fix_bnd else ct.c_int(0)
    c_eh = ct.c_int(1) if empty_holes else ct.c_int(0)
    args = (c_g1, c_g2, c_buf, c_fix, c_eh)

    if cb is None:
        cb = basic.interf.SilentCallbackCancel2PG()
    lib_fa.cross_grids_wcb.restype = ct.c_void_p
    cb.initialize(lib_fa.cross_grids_wcb, args)
    cb.execute_command()
    c_cross = ct.c_void_p(cb.get_result())

    #if result was obtained (no errors, no cancel)
    if (c_cross.value is not None):
        ret = cobj.grid_from_c(c_cross)
        lib_fa.grid_free(c_cross)
    else:
        ret = None

    lib_fa.grid_free(c_g1)
    lib_fa.grid_free(c_g2)

    return ret


def grid_excl_cont(grd, cnt, is_inner, cb=None):
    """ ->grid2.Grid2
        Returns a grid with excluded contour area (inner or outer)
    """
    lib_fa = cobj.crossgrid_lib()
    c_g = cobj.grid_to_c(grd)
    c_c = cobj.cont_to_c(cnt)
    c_isinner = ct.c_int(1 if is_inner else 0)
    args = (c_g, c_c, c_isinner)
    lib_fa.grid_exclude_cont_wcb.restype = ct.c_void_p

    # e = dlgs.ProgressProcedureDlg(lib_fa.grid_exclude_cont_wcb, args,
    #         True, globvars.mainWindow)
    # e.exec_()
    # res = ct.c_void_p(e.get_result())
    if cb is None:
        cb = basic.interf.SilentCallbackCancel2PG()
    cb.initialize(lib_fa.grid_exclude_cont_wcb, args)
    cb.execute_command()
    res = ct.c_void_p(cb.get_result())

    if res.value is not None:
        newg = cobj.grid_from_c(res)
        newg.build_contour()

        #assign boundary conditions
        bs = cnt.bnd_types().union(grd.cont.bnd_types()).difference(set([0]))
        if len(bs) > 0:
            res_cont = cobj.cont_to_c(newg.cont)
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


def boundary_layer_grid(cont, opt):
    """-> Grid2 or None. Builds a grid around contour.
    cont - AbstractContour2
    opt - dlgs._default_odata() object
    """
    lib_cg = cobj.crossgrid_lib()
    lib_c2 = cobj.cont2_lib()
    # parse arguments
    c_cont = cobj.cont2_to_c(cont)
    c_isinner = ct.c_int(1 if opt['direct'] == 'inner' else 0)
    c_npart = ct.c_int(len(opt['partition']))
    c_part = (ct.c_double * len(opt['partition']))()
    for i, p in enumerate(opt['partition']):
        c_part[i] = p
    _t = {'No': 0,
            'Mesh. Ignore all nodes': 1,
            'Mesh. Keep shape': 2,
            'Mesh. Keep origin': 3,
          }[opt['mesh_cont']]
    c_mesh_cont_type = ct.c_int(_t)
    c_mesh_cont_step = ct.c_double(opt['mesh_cont_step'])
    c_round_off = ct.c_int(1 if opt['rnd'] else 0)
    c_minsharp = ct.c_double(opt['minsharp'])
    lib_cg.boundary_layer_grid.restype = ct.c_void_p
    res = lib_cg.boundary_layer_grid(c_cont, c_isinner,
            c_npart, c_part, c_mesh_cont_type, c_mesh_cont_step,
            c_round_off, c_minsharp)
    if res is not None:
        newg = cobj.grid_from_c(res)
        newg.build_contour()
        #assign boundary conditions
        bs = cont.bnd_types().difference(set([0]))
        if len(bs) > 0:
            add_bc_from_cont(newg.cont, cont, c_src=c_cont)

        #free grid memory
        lib_cg.grid_free(res)
    else:
        newg = None

    lib_c2.free_contour_tree(c_cont)
    return newg
