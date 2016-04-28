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


