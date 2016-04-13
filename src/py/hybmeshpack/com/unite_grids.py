' grid unification algorithm '
import ctypes as ct
from hybmeshpack.hmcore import libhmcport
from hybmeshpack.hmcore import c2 as c2core
from hybmeshpack.hmcore import g2 as g2core


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
    #create c representations of contours
    need_src_clear, need_tar_clear = False, False
    if c_src is None:
        need_src_clear = True
        c_src = c2core.cont2_to_c(src_cont)

    if c_tar is None:
        need_tar_clear = True
        c_tar = c2core.cont2_to_c(tar_cont)

    #create contours boundary array
    c_src_bnd = (ct.c_int * src_cont.n_edges())()  # input
    c_tar_bnd = (ct.c_int * tar_cont.n_edges())()  # output
    for i in range(src_cont.n_edges()):
        c_src_bnd[i] = src_cont.edge_bnd(i)
    for i in range(tar_cont.n_edges()):
        c_tar_bnd[i] = tar_cont.edge_bnd(i)

    #call libhmcport function
    if force == 1:
        out = libhmcport.set_ecollection_bc(
            c_src, c_tar, ct.c_int(-1), c_src_bnd, c_tar_bnd)
    else:
        out = libhmcport.set_ecollection_bc_force(
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
        libhmcport.free_ecollection_container(c_src)
    if need_tar_clear:
        libhmcport.free_ecollection_container(c_tar)

    if int(out) != 0:
        raise Exception("Error at boundary assignment")


def setbc_from_conts(tar_cont, src, force=1):
    """ (contour2.AbstractContour, [contour2.AbstractContour]) -> None
        Sets boundary condition to tar_cont from list of source contours.

        force is the same as in add_bc_from_cont
    """
    #set default zero to all target boundaries
    tar_cont.set_edge_bnd({i: 0 for i in range(tar_cont.n_edges())})

    #loop over each source contours
    c_tar = c2core.cont2_to_c(tar_cont)
    for cont in src:
        c_src = c2core.cont2_to_c(cont)
        add_bc_from_cont(tar_cont, cont, c_tar, c_src, force)
        libhmcport.free_ecollection_container(c_src)
    libhmcport.free_ecollection_container(c_tar)


def unite_grids(g1, g2, buf, fix_bnd, empty_holes, an0, cb):
    """ adds g2 to g1. Returns new grid.
        cb -- Callback.CB_CANCEL2 callback object
    """
    c_g1 = g2core.grid_to_c(g1)
    c_g2 = g2core.grid_to_c(g2)
    c_buf = ct.c_double(buf)
    c_an0 = ct.c_double(an0)
    c_fix = ct.c_int(1) if fix_bnd else ct.c_int(0)
    c_eh = ct.c_int(1) if empty_holes else ct.c_int(0)
    args = (c_g1, c_g2, c_buf, c_fix, c_eh, c_an0)

    libhmcport.cross_grids_wcb.restype = ct.c_void_p
    cb.initialize(libhmcport.cross_grids_wcb, args)
    cb.execute_command()
    c_cross = ct.c_void_p(cb.get_result())

    #if result was obtained (no errors, no cancel)
    if (c_cross.value is not None):
        ret = g2core.grid_from_c(c_cross)
        libhmcport.grid_free(c_cross)
    else:
        ret = None

    libhmcport.grid_free(c_g1)
    libhmcport.grid_free(c_g2)

    return ret


def grid_excl_cont(grd, cnt, is_inner, cb):
    """ ->grid2.Grid2
        Returns a grid with excluded contour area (inner or outer)
    """
    c_g = g2core.grid_to_c(grd)
    c_c = c2core.cont_to_c(cnt)
    c_isinner = ct.c_int(1 if is_inner else 0)
    args = (c_g, c_c, c_isinner)
    libhmcport.grid_exclude_cont_wcb.restype = ct.c_void_p

    cb.initialize(libhmcport.grid_exclude_cont_wcb, args)
    cb.execute_command()
    res = ct.c_void_p(cb.get_result())

    if res.value is not None:
        newg = g2core.grid_from_c(res)
        newg.build_contour()

        #assign boundary conditions
        bs = cnt.bnd_types().union(grd.cont.bnd_types()).difference(set([0]))
        if len(bs) > 0:
            res_cont = c2core.cont2_to_c(newg.cont)
            #1. from source grid
            add_bc_from_cont(newg.cont, grd.cont, c_tar=res_cont)
            #2. from contour
            add_bc_from_cont(newg.cont, cnt, c_tar=res_cont)
            #3. free contour memory
            libhmcport.free_ecollection_container(res_cont)

        #free grid memory
        libhmcport.grid_free(res)
    else:
        newg = None

    libhmcport.cont_free(c_c)
    libhmcport.grid_free(c_g)

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

    # 1) get contour pointers
    for co in opt:
        pycont = co['source']
        if pycont not in COptStruct.usedconts:
            COptStruct.usedconts[pycont] = c2core.cont2_to_c(pycont)

    # 2) build array of c structures
    c_opt_type = COptStruct * len(opt)
    c_opt = c_opt_type()
    for i, co in enumerate(opt):
        c_opt[i] = COptStruct(co)

    # 3) call through callback object
    args = (len(opt), c_opt)
    libhmcport.boundary_layer_grid_wcb.restype = ct.c_void_p
    cb.initialize(libhmcport.boundary_layer_grid_wcb, args)
    cb.execute_command()
    cres = ct.c_void_p(cb.get_result())

    # 4) take result
    if cres.value is not None:
        ret = g2core.grid_from_c(cres)
    elif cb._proceed:
        # if no result and no cancel -> there was an error
        raise Exception("Boundary grid builder error")
    else:
        # if cancel -> return None
        ret = None

    # 5) free data
    libhmcport.grid_free(cres)
    for c in COptStruct.usedconts.values():
        libhmcport.free_ecollection_container(c)
    COptStruct.usedconts = {}

    return ret


def map_grid(grid, cont, gpoints, cpoints, snap, bt):
    """maps grid on cont using gpoints, cpoints as basis.
       snap = "no", "add_vertices", "shift_vertices"
       bt = "from_grid", "from_contour"
    """
    # build c grid, contour
    cgrid = g2core.grid_to_c(grid)
    ccont = c2core.cont2_to_c(cont)

    # build points array
    n = ct.c_int(min(len(gpoints), len(cpoints)))
    p1 = (ct.c_double * (2 * n.value))()
    p2 = (ct.c_double * (2 * n.value))()
    for i in range(n.value):
        p1[2 * i] = gpoints[i].x
        p1[2 * i + 1] = gpoints[i].y
        p2[2 * i] = cpoints[i].x
        p2[2 * i + 1] = cpoints[i].y

    # snap
    if snap == 'no':
        s = ct.c_int(1)
    elif snap == "add_vertices":
        s = ct.c_int(2)
    elif snap == "shift_vertices":
        s = ct.c_int(3)
    else:
        raise ValueError("Unknown snap = " + str(snap))

    # call c function
    cret = libhmcport.build_grid_mapping(cgrid, ccont, n, p1, p2, s)
    if cret == 0:
        g2core.free_c_grid(cgrid)
        c2core.free_cont2(ccont)
        raise Exception("Failed to map a grid")

    # copy to hm grid
    ret = g2core.grid_from_c(cret)
    ret.build_contour()

    # treat boundaries
    if bt == "from_contour":
        add_bc_from_cont(ret.cont, cont, c_src=ccont, force=3)
    elif bt == "from_grid":
        # if snap != "add_vertices":
            # # topology is same hence we can simply copy bc from source
            # ret.bt = copy.deepcopy(grid.bt)
        # else:
            # all ret points which don't have their siblings in grid
            # were added at the end of ret.points array. Hence we can
            # restore edge btypes from points indices

            # 1. assemble boundary edges: edge index -> bnd feature
            ret_full_bnd = ret.boundary_contours()
            ret_bnd_edges = []
            for ecol in ret_full_bnd:
                ret_bnd_edges.extend(ecol)
            grid_bnd_edges = []
            for ecol in grid.boundary_contours():
                grid_bnd_edges.extend(ecol)
            # 2. whether boundary edge has its sibling in grid
            is_old_edge = []
            oldptsnum = grid.n_points()
            for ei in ret_bnd_edges:
                e = ret.edges[ei]
                is_old_edge.append(e[0] < oldptsnum and e[1] < oldptsnum)

            # 3. build dictionary for boundary edges in grid
            def etostr(e):
                return str(min(e[0], e[1])) + '-' + str(max(e[0], e[1]))
            old_ebt = {}
            for ei in grid_bnd_edges:
                bt = 0
                if ei in grid.bt:
                    bt = grid.bt[ei]
                old_ebt[etostr(grid.edges[ei])] = bt

            # 4. restore bfeatures for old edges
            for i, ei in enumerate(ret_bnd_edges):
                if is_old_edge[i]:
                    bt = old_ebt[etostr(ret.edges[ei])]
                    if bt > 0:
                        ret.bt[ei] = bt

            # 5. boundary types for new edges
            for c in ret_full_bnd:
                # find any old edge
                for istart, ei in enumerate(c):
                    if is_old_edge[ret_bnd_edges.index(ei)]:
                        break
                else:
                    # no old edges in current contour
                    # this should not happen
                    continue
                # starting from old edge fill all new ones
                bcur = ret.get_edge_bnd(c[istart])
                for i in range(len(c)):
                    icur = i + istart
                    if icur >= len(c):
                        icur -= len(c)
                    ecur = c[icur]
                    if is_old_edge[ret_bnd_edges.index(ecur)]:
                        bcur = ret.get_edge_bnd(c[icur])
                    else:
                        if bcur > 0:
                            ret.bt[c[icur]] = bcur
    else:
        raise ValueError("Unknown btypes = " + str(bt))

    # free c memory
    g2core.free_c_grid(cgrid)
    c2core.free_cont2(ccont)
    g2core.free_c_grid(cret)

    # return
    return ret
