' grid unification algorithm '
import ctypes as ct
import cobj


def add_bc_from_cont(tar_cont, src_cont, c_tar=None, c_src=None, force=1):
    """ (contour2.AbstractContour, contour2.AbstractContour) -> None

    Adds boundary condition to tar_cont from src_cont
    if grid boundary edge lies on the contour edge.
    Changes tar_cont data.

    Helper c_src and c_tar are the c representation of
    target and source contours of ct.c_void_p type.

    force = 1: only edges which lie on src will be assigned
    force = 2: edges which end points lie on src will be assigned
    force = 3: all edges will be assigned
    """
    lib_fa = cobj.cport_lib()

    #create c representations of contours
    need_src_clear, need_tar_clear = False, False
    if c_src is None:
        need_src_clear = True
        c_src = cobj.cont2_to_c(src_cont)

    if c_tar is None:
        need_tar_clear = True
        c_tar = cobj.cont2_to_c(tar_cont)

    #create contours boundary array
    c_src_bnd = (ct.c_int * src_cont.n_edges())()  # input
    c_tar_bnd = (ct.c_int * tar_cont.n_edges())()  # output
    for i in range(src_cont.n_edges()):
        c_src_bnd[i] = src_cont.edge_bnd(i)

    #call lib_fa function
    if force == 1:
        out = lib_fa.set_ecollection_bc(
            c_src, c_tar, ct.c_int(-1), c_src_bnd, c_tar_bnd)
    else:
        out = lib_fa.set_ecollection_bc_force(
            c_src, c_tar, c_src_bnd, c_tar_bnd, ct.c_int(force))

    #write data to target contour
    new_bnd = {}
    for i in range(tar_cont.n_edges()):
        b = c_tar_bnd[i]
        if b >= 0:
            new_bnd[i] = b
    tar_cont.add_edge_bnd(new_bnd)

    #clear c data
    if need_src_clear:
        lib_fa.free_ecollection_container(c_src)
    if need_tar_clear:
        lib_fa.free_ecollection_container(c_tar)

    if int(out) != 0:
        raise Exception("Error at boundary assignment")


def setbc_from_conts(tar_cont, src, force=1):
    """ (contour2.AbstractContour, [contour2.AbstractContour]) -> None
        Sets boundary condition to tar_cont from list of source contours.

        force is the same as in add_bc_from_cont
    """
    lib_fa = cobj.cport_lib()
    #set default zero to all target boundaries
    tar_cont.set_edge_bnd({i: 0 for i in range(tar_cont.n_edges())})

    #loop over each source contours
    c_tar = cobj.cont2_to_c(tar_cont)
    for cont in src:
        c_src = cobj.cont2_to_c(cont)
        add_bc_from_cont(tar_cont, cont, c_tar, c_src, force)
        lib_fa.free_ecollection_container(c_src)
    lib_fa.free_ecollection_container(c_tar)


def unite_grids(g1, g2, buf, fix_bnd, empty_holes, cb):
    """ adds g2 to g1. Returns new grid.
        cb -- Callback.CB_CANCEL2 callback object
    """
    lib_fa = cobj.cport_lib()
    c_g1 = cobj.grid_to_c(g1)
    c_g2 = cobj.grid_to_c(g2)
    c_buf = ct.c_double(buf)
    c_fix = ct.c_int(1) if fix_bnd else ct.c_int(0)
    c_eh = ct.c_int(1) if empty_holes else ct.c_int(0)
    args = (c_g1, c_g2, c_buf, c_fix, c_eh)

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


def grid_excl_cont(grd, cnt, is_inner, cb):
    """ ->grid2.Grid2
        Returns a grid with excluded contour area (inner or outer)
    """
    lib_fa = cobj.cport_lib()
    c_g = cobj.grid_to_c(grd)
    c_c = cobj.cont_to_c(cnt)
    c_isinner = ct.c_int(1 if is_inner else 0)
    args = (c_g, c_c, c_isinner)
    lib_fa.grid_exclude_cont_wcb.restype = ct.c_void_p

    cb.initialize(lib_fa.grid_exclude_cont_wcb, args)
    cb.execute_command()
    res = ct.c_void_p(cb.get_result())

    if res.value is not None:
        newg = cobj.grid_from_c(res)
        newg.build_contour()

        #assign boundary conditions
        bs = cnt.bnd_types().union(grd.cont.bnd_types()).difference(set([0]))
        if len(bs) > 0:
            res_cont = cobj.cont2_to_c(newg.cont)
            #1. from source grid
            add_bc_from_cont(newg.cont, grd.cont, c_tar=res_cont)
            #2. from contour
            add_bc_from_cont(newg.cont, cnt, c_tar=res_cont)
            #3. free contour memory
            lib_fa.free_ecollection_container(res_cont)

        #free grid memory
        lib_fa.grid_free(res)
    else:
        newg = None

    lib_fa.cont_free(c_c)
    lib_fa.grid_free(c_g)

    return newg


def boundary_layer_grid(opt, cb):
    """ ->grid2.
        opt - options dictionary object.
              Same as gridcom.BuildBoundaryGrid options['opt'] dictionary
              except for 'source' field is a proper Contour2 object.
              All fields must be filled
        cb - callback of CB_CANCEL2 type
    """
    # prepare c input data
    class COptStruct(ct.Structure):
        usedconts = {}

        def __init__(self, opt_entry):
            self.cont = self.usedconts[opt_entry['source']]
            n = len(opt_entry['partition'])
            self.Npartition = ct.c_int(n)
            self.partition = (ct.c_double * n)(*opt_entry['partition'])
            self.direction = ct.c_char_p(
                {1: 'LEFT', -1: 'RIGHT'}[opt_entry['direction']])
            self.mesh_cont = ct.c_char_p({
                0: 'NO', 1: 'KEEP_ORIGIN', 2: 'KEEP_SHAPE', 3: 'IGNORE_ALL',
                4: "INCREMENTAL",
            }[opt_entry['mesh_cont']])
            self.mesh_cont_step = ct.c_double(opt_entry['mesh_cont_step'])
            self.start = (ct.c_double * 2)(opt_entry['start'].x,
                                           opt_entry['start'].y)
            self.end = (ct.c_double * 2)(opt_entry['end'].x,
                                         opt_entry['end'].y)
            self.force_conformal = ct.c_int(
                1 if opt_entry['force_conf'] else 0)
            self.angle_range = (ct.c_double * 4)(
                opt_entry['algo_acute'],
                opt_entry['algo_right'],
                opt_entry['algo_straight'],
                opt_entry['algo_reentr'])
            self.step_start = ct.c_double(opt_entry['step_start'])
            self.step_end = ct.c_double(opt_entry['step_end'])

        _fields_ = [
            ('cont', ct.c_void_p),
            ('Npartition', ct.c_int),
            ('partition', ct.POINTER(ct.c_double)),
            ('direction', ct.c_char_p),
            ('mesh_cont', ct.c_char_p),
            ('mesh_cont_step', ct.c_double),
            ('start', ct.c_double * 2),
            ('end', ct.c_double * 2),
            ('force_conformal', ct.c_int),
            ('angle_range', ct.c_double * 4),
            ('step_start', ct.c_double),
            ('step_end', ct.c_double),
        ]

    lib_fa = cobj.cport_lib()

    # 1) get contour pointers
    for co in opt:
        pycont = co['source']
        if pycont not in COptStruct.usedconts:
            COptStruct.usedconts[pycont] = cobj.cont2_to_c(pycont)

    # 2) build array of c structures
    c_opt_type = COptStruct * len(opt)
    c_opt = c_opt_type()
    for i, co in enumerate(opt):
        c_opt[i] = COptStruct(co)

    # 3) call through callback object
    args = (len(opt), c_opt)
    lib_fa.boundary_layer_grid_wcb.restype = ct.c_void_p
    cb.initialize(lib_fa.boundary_layer_grid_wcb, args)
    cb.execute_command()
    cres = ct.c_void_p(cb.get_result())

    # 4) take result
    if cres.value is not None:
        ret = cobj.grid_from_c(cres)
    elif cb._proceed:
        # if no result and no cancel -> there was an error
        raise Exception("Boundary grid builder error")
    else:
        # if cancel -> return None
        ret = None

    # 5) free data
    lib_fa.grid_free(cres)
    for c in COptStruct.usedconts.values():
        lib_fa.free_ecollection_container(c)
    COptStruct.usedconts = {}

    return ret
