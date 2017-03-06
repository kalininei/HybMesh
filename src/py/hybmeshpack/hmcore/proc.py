import ctypes as ct
from hybmeshpack.basic.cb import SilentCallbackCancel2, UserInterrupt
from . import CppLibError, cport


def get_last_cside_error():
    ret = ct.c_char_p()
    cport.get_last_error_message(ct.byref(ret))
    out = str(ret.value)
    free_cside_array(ret, "char")
    return out


def ccall_cb(func, cb, *args):
    if cb is None:
        cb = SilentCallbackCancel2()
    cb.initialize(func, args)
    cb.execute_command()
    if not cb.get_result():
        msg = get_last_cside_error()
        if msg.startswith("User interrupt"):
            raise UserInterrupt()
        else:
            raise CppLibError(msg)


def ccall(func, *args):
    ok = func(*args)
    if not ok:
        msg = get_last_cside_error()
        raise CppLibError(msg)


def list_to_c(lst, tp):
    """ returns (c_tp * len(lst)()
    where tp = ct.c_double or ct.c_int depending on tp = float/int
    if lst is None or zero length returns NULL ptr
    """
    if lst is None:
        return None
    if isinstance(lst, list) or isinstance(lst, tuple):
        d = len(lst)
        if d == 0:
            return None
        if tp == int or tp == 'int':
            ret = (ct.c_int * len(lst))()
            for i in range(d):
                ret[i] = int(lst[i])
            return ret
        if tp == float or tp == 'float' or tp == 'double':
            ret = (ct.c_double * len(lst))()
            for i in range(d):
                ret[i] = float(lst[i])
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


class BndTypesDifference(object):
    def __init__(self):
        # data is:
        #    total number of btypes,
        #    whole_btype, 1, -1      # if whole is needed
        #    btype0, number of ind0_i, ind0_0, ind0_1, .....
        #    btype1, number of ind1_i, ind1_0, ind1_1, .....
        super(BndTypesDifference, self).__init__()
        self.data = None
        self.need_free = False

    @classmethod
    def from_pydata(cls, whole, typesdict):
        """ whole is None or btype for all boundary primitives
            typesdict is {btype: [primitives indicies]}
        """
        ret = cls()
        dtsize = 1
        if whole is not None:
            dtsize += 3
        for v in typesdict.values():
            dtsize += 2 + len(v)
        ret.data = (ct.c_int * dtsize)()
        ret.data[0] = len(typesdict)
        if whole is not None:
            ret.data[0] += 1
            ret.data[1], ret.data[2], ret.data[3] = whole, 1, -1
            icur = 4
        else:
            icur = 1
        for k, v in typesdict.iteritems():
            ret.data[icur], ret.data[icur+1] = k, len(v)
            icur += 2
            for t in v:
                ret.data[icur] = t
                icur += 1
        return ret

    @classmethod
    def from_cdata(cls, cd):
        ret = cls()
        if isinstance(cd, ct.POINTER(ct.c_int)):
            ret.need_free = True
        ret.data = cd

    def __del__(self):
        # if it was allocated on c-side
        if self.need_free:
            free_cside_array(self.data, "int")
