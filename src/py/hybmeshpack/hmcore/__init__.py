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
        if tp == 'char':
            s = ''.join([chr(x) for x in lst])
            ret = s
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


def hmxml_new():
    ret = libhmcport.new_writer()
    if ret == 0:
        raise Exception("Failed creating a xml document")
    return ret


def hmxml_read(fname):
    ret = libhmcport.new_reader(fname)
    if ret == 0:
        raise Exception("Failed parsing file " + fname)
    return ret


def hmxml_finalize(writer, fname):
    """ write to file and close"""
    ret = libhmcport.finalize_writer(writer, fname)
    if ret == 0:
        raise Exception("Failed writing a xml document")


def hmxml_free_node(node):
    """ frees and deletes xml substructure if root """
    libhmcport.free_hmxml_node(node)


def hmxml_query(reader, query, required="no"):
    """ required = "no", '>0', '=1', '=0'
        returns hmxml nodes list
    """
    num = ct.c_int(0)
    ans = (ct.c_void_p * 1000)()
    ret = libhmcport.hmxml_query(reader, query, ct.byref(num), ans)
    num = num.value
    if ret == 0:
        raise Exception("Query failed " + query)
    try:
        if required == '>0' and num == 0:
            raise
        if required == '=1' and num != 1:
            raise
        if required == '=0' and num != 0:
            raise
    except:
        for i in range(num):
            hmxml_free_node(ans[i])
        raise Exception("Improper number of entries (=" + str(num) + ") for " +
                        query)
    ret = []
    for i in range(num):
        ret.append(ans[i])
    return ret
