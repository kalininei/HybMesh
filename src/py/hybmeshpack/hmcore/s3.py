import ctypes as ct
from . import libhmcport, free_cside_array


def surface_from_c(cdata):
    from hybmeshpack.gdata.srf3 import Surface3
    return Surface3(cdata)


def free_srf3(c_s3):
    libhmcport.free_srf3(c_s3)


def dims(surf=None, grid=None):
    " -> [nvertices, nedges, nfaces] "
    ret = (ct.c_int * 3)()
    if surf is not None:
        libhmcport.surf3_dims(surf, ret)
    elif grid is not None:
        libhmcport.grid3_surface_dims(grid, ret)
    else:
        return [0, 0, 0]
    return [ret[0], ret[1], ret[2]]


def btypes(surf=None, grid=None):
    " -> [btype for each surface] "
    cret = (ct.c_int * dims(surf, grid)[2])()
    r = 0
    if surf is not None:
        r = libhmcport.surf3_btypes(surf, cret)
    elif grid is not None:
        r = libhmcport.grid3_surface_btypes(grid, cret)
    if r == 0:
        raise Exception("failed to read boundary types of surface")
    ret = []
    for i in range(len(cret)):
        ret.append(cret[i])
    return ret


def extract_grid3_surface(c_g3):
    "-> c pointer to deep copied surface"
    ret = libhmcport.extract_grid3_surface(c_g3)
    if ret == 0:
        raise Exception("Failed to extract grid surface")
    return ret


def extract_subsurfaces(c_s3):
    "-> c pointers to shallow copied surfaces"
    c_ret = ct.POINTER(ct.c_void_p)()
    c_num = ct.c_int()
    r = libhmcport.extract_subsurfaces(c_s3, ct.byref(c_num), ct.byref(c_ret))
    if r == 0:
        raise Exception("Failed to extract subsurfaces")
    ret = []
    for i in range(c_num.value):
        ret.append(c_ret[i])
    free_cside_array(c_ret, "void*")
    return ret


def surface_from_hmxml(reader, subnode):
    sreader, c_s = 0, 0
    try:
        name = ct.create_string_buffer(1000)

        libhmcport.s3reader_create.restype = ct.c_void_p
        sreader = libhmcport.s3reader_create(reader, subnode, name)
        if not sreader:
            raise Exception("Failed to read a surface")
        c_s = libhmcport.s3reader_getresult(sreader)
        s = surface_from_c(c_s)
        return s, str(name.value)
    except:
        raise
    finally:
        libhmcport.s3reader_free(sreader) if sreader != 0 else None


def swriter_create(sname, c_s, c_writer, c_sub, fmt):
    ret = libhmcport.s3writer_create(sname, c_s, c_writer, c_sub, fmt)
    if ret == 0:
        raise Exception("Error writing surface to hmc")
    else:
        return ret


def swriter_add_field(c_swriter, f):
    ret = libhmcport.g3writer_add_defined_field(c_swriter, f)
    if ret == 0:
        raise Exception("Error writing " + f + "field to a surface")
    else:
        return ret


def free_swriter(c_swriter):
    libhmcport.s3writer_free(c_swriter)
