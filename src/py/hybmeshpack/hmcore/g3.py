import ctypes as ct
from . import libhmcport
from hybmeshpack.gdata.grid3 import Grid3


def extrude_call(c_grid, c_btypes, c_zvals, c_bbot, c_btop, bsides):
    """ makes extrusion:
    c_grid: pointer to grid structure get by grid_to_c
    c_btypes: boundary types for c_grid get by boundary_types_to_c
    c_zvals: ctypes.c_double array of z values for extrusion
    c_bbot, c_btop: ctypes.c_int array of resulting boundary types
                    at z=min, z=max surface.
                    Their length could equal number of grid cells or 1.
                    The latter gives same boundary type for whole surface
    c_bsides: integer btype  at side surface or None
              int -> takes bc from c_btypes
              None -> constant value for whole surface
    returns Grid3 c pointer or None if failed
    """
    bot_algo = 0 if len(c_bbot) == 1 else 1
    top_algo = 0 if len(c_btop) == 1 else 1
    if bsides is None:
        side_algo = 1
        bsides = 1
    else:
        side_algo = 0

    ret = libhmcport.grid2_sweep_z(c_grid, c_btypes,
                                   ct.c_int(len(c_zvals)), c_zvals,
                                   ct.c_int(top_algo), c_btop,
                                   ct.c_int(bot_algo), c_bbot,
                                   ct.c_int(side_algo), ct.c_int(bsides))
    if ret == 0:
        return None
    else:
        return ret


def grid3_from_c(cdata):
    return Grid3.from_cdata(cdata)


def free_g3(c_g3):
    libhmcport.free_grid3d(c_g3)


def to_vtk(c_g, fname, cb):
    c_fname = fname.encode('utf-8')
    args = (c_g, c_fname)

    cb.initialize(libhmcport.export_vtk_grid3, args)
    cb.execute_command()

    res = ct.c_int(cb.get_result())

    if res.value != 0:
        raise Exception("VTK export failed")


def surface_to_vtk(c_g, fname, cb):
    c_fname = fname.encode('utf-8')
    args = (c_g, c_fname)

    cb.initialize(libhmcport.export_surface_vtk_grid3, args)
    cb.execute_command()

    res = ct.c_int(cb.get_result())

    if res.value != 0:
        raise Exception("VTK surface export failed")


def to_msh(c_g, fname, c_bnames, c_periodic, cb):
    c_fname = fname.encode('utf-8')
    if c_periodic is not None:
        n_periodic = ct.c_int(len(c_periodic) / 8)
    else:
        n_periodic = ct.c_int(0)

    args = (c_g, c_fname, c_bnames, n_periodic, c_periodic)
    cb.initialize(libhmcport.export_msh_grid3, args)
    cb.execute_command()

    res = cb.get_result()

    if res != 0:
        raise Exception("msh 3d grid export failed")


def to_tecplot(c_g, fname, c_bnames, cb):
    c_fname = fname.encode('utf-8')
    args = (c_g, c_fname, c_bnames)
    cb.initialize(libhmcport.export_tecplot_grid3, args)
    cb.execute_command()
    res = cb.get_result()
    if res != 0:
        raise Exception("tecplot 3d grid export failed")
