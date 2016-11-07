' procedures for communication with c libraries'
import ctypes as ct
import os
from hybmeshpack import progdata


_libfn = progdata.get_lib_fn('cport')
libhmcport = ct.cdll.LoadLibrary(os.path.abspath(_libfn))


def list_to_c(lst, tp):
    """ returns (c_tp * len(lst)()
    where tp = ct.c_double or ct.c_int depending on tp = float/int
    if lst is None or zero length returns NULL ptr
    """
    if lst is None:
        return None
    if isinstance(lst, list):
        d = len(lst)
        if d == 0:
            return None
        if tp == int:
            ret = (ct.c_int * len(lst))()
            for i in range(d):
                ret[i] = lst[i]
            return ret
        if tp == float:
            ret = (ct.c_double * len(lst))()
            for i in range(d):
                ret[i] = lst[i]
            return ret
        if tp == str:
            ret = (ct.c_char_p * len(lst))()
            for i in range(d):
                ret[i] = lst[i]
            return ret
        if tp == 'void*':
            ret = (ct.c_void_p * len(lst))()
            for i in range(d):
                ret[i] = lst[i]
            return ret
    raise ValueError(str(tp))


def boundary_names_to_c(bnames):
    nameslist = []
    valslist = []
    for b in bnames._data:
        valslist.append(b.index)
        if b.index != 0:
            nameslist.append(b.name)
        else:
            nameslist.append("default_boundary")
    c_names = list_to_c(nameslist, str)
    c_vals = list_to_c(valslist, int)
    if c_names is not None:
        c_n = ct.c_int(len(c_names))
    else:
        c_n = ct.c_int(0)

    libhmcport.set_boundary_names.restype = ct.c_void_p
    return libhmcport.set_boundary_names(c_n, c_names, c_vals)


def free_boundary_names(c_bnames):
    libhmcport.free_boundary_names(c_bnames)
