import ctypes as ct
from hybmeshpack.basic.interf import SilentCallbackCancel2
from . import CppLibError, cport


def ccall_cb(cb, func, *args):
    if cb is None:
        cb = SilentCallbackCancel2()
    cb.initialize(func, args)
    cb.execute_command()
    if not cb.get_result():
        raise CppLibError()


def ccall(func, *args):
    ok = func(*args)
    if not ok:
        raise CppLibError()


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
        if tp == int or tp == 'int':
            ret = (ct.c_int * len(lst))()
            for i in range(d):
                ret[i] = lst[i]
            return ret
        if tp == float or tp == 'float' or tp == 'double':
            ret = (ct.c_double * len(lst))()
            for i in range(d):
                ret[i] = lst[i]
            return ret
        if tp == str or tp == 'str':
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


def free_cside_array(a, tp):
    if not tp:
        return
    if tp == int or tp == "int":
        cport.free_int_array(a)
    elif tp == float or tp == "double" or tp == 'float':
        cport.free_double_array(a)
    elif tp == "char":
        cport.free_char_array(a)
    elif tp == "void*":
        cport.free_voidp_array(a)
    else:
        raise ValueError(str(tp))


def move_to_static(n, dyn, tp):
    ctp = None
    if tp == int or tp == 'int':
        ctp = ct.c_int
    elif tp == float or tp == 'double' or tp == 'float':
        ctp = ct.c_double
    elif tp == "char":
        ctp = ct.c_char
    elif tp == "void*":
        ctp = ct.c_void_p
    else:
        raise ValueError

    ret = (ctp * n)()
    ct.memmove(ret, dyn, ct.sizeof(ret))
    free_cside_array(dyn, tp)
    return ret


class CBoundaryNames(ct.Structure):
    def __init__(self, bdict):
        " bdict is {index: name} "
        self.n = ct.c_int(len(bdict))
        self.index = (ct.c_int * len(bdict))()
        self.name = (ct.c_char_p * len(bdict))()
        for i, k in enumerate(sorted(bdict.keys())):
            self.index[i] = k
            self.name[i] = bdict[k].encode('utf-8')

    _fields_ = [('n', ct.c_int),
                ('index', ct.POINTER(ct.c_int)),
                ('name', ct.POINTER(ct.c_char_p))]


def supplement(bnd, sz):
    """ supplements list bnd to fit size *sz* with bnd[-1] value.
        if None -> [0, 0, ....]
        if not list -> [bnd, bnd, ....]
    """

    if not bnd:
        bnd = [0]
    elif not isinstance(bnd, list):
        bnd = [bnd]
    bnd = bnd[:sz]
    while len(bnd) < sz:
        bnd.append(bnd[-1])
    return bnd


def concat(a):
    """ [[a, b, c], [d, e], ...] -> [a, b, c, d, ...] """
    b = []
    for it in a:
        b.extend(it)
    return b
