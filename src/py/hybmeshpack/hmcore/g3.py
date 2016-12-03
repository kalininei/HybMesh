import ctypes as ct
from . import libhmcport, list_to_c


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
    from hybmeshpack.gdata.grid3 import Grid3
    return Grid3(cdata)


def free_g3(c_g3):
    libhmcport.free_grid3d(c_g3)


def dims(c_g3):
    " -> [nvertices, nedges, nfaces, ncells] "
    ret = (ct.c_int * 4)()
    libhmcport.grid3_dims(c_g3, ret)
    return [ret[0], ret[1], ret[2], ret[3]]


def grid_from_hmxml(reader, subnode, cb):
    greader, c_g = 0, 0
    try:
        name = ct.create_string_buffer(1000)

        args = (reader, subnode, name)
        cb.initialize(libhmcport.g3reader_create, args)
        cb.execute_command()
        greader = ct.c_void_p(cb.get_result())

        if not greader:
            raise Exception("Failed to assemble a grid")

        # get grid
        c_g = libhmcport.g3reader_getresult(greader)
        g = grid3_from_c(c_g)

        return g, str(name.value)
    except:
        raise
    finally:
        libhmcport.g3reader_free(greader) if greader != 0 else None


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


def revolve_call(c_g, c_vec, c_phi, c_btypes, b1, b2, tri_center):
    """ Makes 2d grid revolution:
    c_g - 2d grid on c side
    c_vec - [x0, y0, x1, y1] c-side array which defines vector of rotation
    c_phi - c-side array of angular partition [0]==[-1] for full revolution
    c_btypes - boundary types of 2d grid on c-side
    b1, b2 - integers defining boundary types of c_phi[0], c_phi[-1] surfaces
    tri_ceneter  - boolean, defining wheher to triangulate center
    returns 3d grid as a c-side pointer
    """
    args = (c_g, c_vec, ct.c_int(len(c_phi)), c_phi, c_btypes,
            ct.c_int(b1), ct.c_int(b2),
            ct.c_int(1 if tri_center else 0))
    res = libhmcport.grid2_revolve(*args)
    if res == 0:
        raise Exception("Grid revolution failed")
    return res


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


def to_gmsh(c_g, fname, c_bnames, cb):
    c_fname = fname.encode('utf-8')

    args = (c_g, c_fname, c_bnames)
    cb.initialize(libhmcport.export_gmsh_grid3, args)
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


def gwriter_create(gridname, c_g, c_writer, c_sub, fmt):
    ret = libhmcport.g3writer_create(gridname, c_g, c_writer, c_sub, fmt)
    if ret == 0:
        raise Exception("Error writing grid3d to hmg")
    else:
        return ret


def gwriter_add_field(c_gwriter, f):
    ret = libhmcport.g3writer_add_defined_field(c_gwriter, f)
    if ret == 0:
        raise Exception("Error writing " + f + "field to a grid 3d")
    else:
        return ret


def free_gwriter(c_gwriter):
    libhmcport.g3writer_free(c_gwriter)


def tetrahedral_fill(slist, constrlist, p, psize):
    plen = max(len(p) / 3, len(psize))
    c_p = list_to_c(p, float)
    c_ps = list_to_c(psize, float)
    c_slist = list_to_c(slist, "void*")
    c_constrlist = list_to_c(constrlist, "void*")

    libhmcport.tetrahedral_fill.restype = ct.c_void_p
    r = libhmcport.tetrahedral_fill(
        ct.c_int(len(slist)), c_slist,
        ct.c_int(len(constrlist)), c_constrlist,
        ct.c_int(plen), c_p, c_ps)

    if not r:
        raise Exception("Tetrahedral fill failed")
    else:
        return r
