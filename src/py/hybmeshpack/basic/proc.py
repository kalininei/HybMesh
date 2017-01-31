#!/usr/bin/env python
import time
import copy
import sys


class AbstractSender(object):
    'Base class which send messages'

    def __init__(self):
        #subscribers list
        self._subscribers = []

    def add_subscriber(self, fun):
        'fun is function(int tp, AbstractSender sender, ....)'
        self._subscribers.append(fun)

    def remove_subscriber(self, fun):
        self._subscribers.remove(fun)

    def _send_message(self, tp, **kwargs):
        "sent a message for subscribers"
        for s in self._subscribers:
            s(tp, self, **kwargs)


class EmbException(Exception):
    """ Exception which remembers the traceback
        of previous Exception. Can be embedded multiple times.
        Example
        try:
            try:
                raise Exception("First")
            except Exception as e:
                raise EmbException("Second", e)
        except EmbException as e:
            EmbException.print_estack(e)
            print e

        will print trace from exception First
        and message from exception Second
    """
    def __init__(self, message, prevexc=None):
        import traceback
        super(EmbException, self).__init__(message)
        self.message = message
        self.__stack = None
        if prevexc is None:
            self.__stack = traceback.format_stack()[:-1]
        elif hasattr(prevexc, '_SmartException__stack'):
            self.__stack = prevexc._SmartException__stack
        else:
            self.__stack = traceback.format_exception(*sys.exc_info())

    def get_stack(self):
        return copy.deepcopy(self.__stack)

    @staticmethod
    def estack(exc):
        if hasattr(exc, '_EmbException__stack'):
            return exc.get_stack()
        else:
            import traceback
            return copy.deepcopy(traceback.format_exc())

    @staticmethod
    def print_estack(exc):
        print ''.join(EmbException.estack(exc))


def exectime(method):
    "decorator for execution time measuring"
    def timed(*args, **kw):
        ts = time.time()
        result = method(*args, **kw)
        te = time.time()
        print '%r:  %2.6f sec' % \
              (method.__name__, te - ts)
        return result
    return timed


def dict_to_plain(d):
    '{k1: v1, k2: v2, k3: v3, ...} -> [k1, v1, k2, v2, ...]'
    ret = []
    for k, v in d.items():
        ret.append(k)
        ret.append(v)
    return ret


def plain_to_dict(p):
    '[k1, v1, k2, v2, ...] -> {k1: v1, k2: v2, k3: v3, ...}'
    ret = {}
    it = iter(p)
    for k, v in zip(it, it):
        ret[k] = v
    return ret


def multifield_list(nd, dim):
    """ ([list], int) -> [[sublist1], [sublist2], ...]

    rebuild plain array into 2d array
    nd - plain array
    dim - dimension of sublist or "auto"
    if "auto" then dimension is set before the start of sublist:
        [1,1, 3,2,2,1, 5,3,3,2,1,1] -> [[1], [2,2,5], [3,3,2,1,1]]
    """
    if dim == 2:
        it = iter(nd)
        return [list(k) for k in zip(it, it)]
    elif dim == 3:
        it = iter(nd)
        return [list(k) for k in zip(it, it, it)]
    elif dim == "auto":
        lst, ret = nd, []
        while lst:
            l1 = lst.pop(0)
            ret.append(lst[:l1])
            lst = lst[l1:]
        return ret
    else:
        lst, ret = nd, []
        while lst:
            ret.append(lst[:dim])
            lst = lst[dim:]
        return ret


def _compress_core(lst, cmpfun, outfun):
    # finds increasing values
    ret = []
    stack = []

    def stack_to_ret(ret, stack):
        if len(stack) < 2:
            ret.extend(stack)
        elif len(stack) > 0:
            ret.append(stack)
        return []

    for v in lst:
        if isinstance(v, str):
            stack = stack_to_ret(ret, stack)
            ret.append(v)
        else:
            if len(stack) > 0 and not cmpfun(v, stack[-1]):
                stack = stack_to_ret(ret, stack)
            stack.append(v)
    stack_to_ret(ret, stack)

    for i in range(len(ret)):
        if isinstance(ret[i], list):
            ret[i] = outfun(ret[i])

    return ret


def _decompress_core(s, symbol, func):
    [a, b] = s.split(symbol)
    return func(int(a), int(b))


def compress_int_list(lst):
    """ changes repeating values with i*k
        changes increasing/decreasing values with i:k
        [18, 17, 16, 1, 2, 3, 4, 4, 4, 4, 4, 2, 1] ->
        "18:16 1:4 4*4 2 1"
    """
    if len(lst) == 0:
        return ""
    ret = lst

    # increasing values
    ret = _compress_core(ret,
                         lambda x, y: x == y + 1,
                         lambda r: str(r[0]) + ":" + str(r[-1]))
    # decreasing values
    ret = _compress_core(ret,
                         lambda x, y: x == y - 1,
                         lambda r: str(r[0]) + ":" + str(r[-1]))
    # repeating values
    ret = _compress_core(ret,
                         lambda x, y: x == y,
                         lambda r: str(r[0]) + "*" + str(len(r)))
    # other values
    ret = map(str, ret)

    # join and return
    return " ".join(ret)


def int_list_from_compress(s):
    """ reciprocal to compress_int_list """
    if len(s) == 0:
        return []
    ret = []
    for a in s.split():
        # separate ints
        try:
            ret.append(int(a))
            continue
        except:
            pass
        # decompress repeating values
        try:
            ret.extend(_decompress_core(a, '*', lambda x, y: [x] * y))
            continue
        except:
            pass
        # decompress increasing/decreasing
        try:
            ret.extend(_decompress_core(a, ':', lambda x, y:
                       range(x, y + 1) if x < y else range(x, y - 1, -1)))
            continue
        except:
            pass
        raise ValueError(a)
    return ret


def dict_readbool(d, key, defval):
    """ -> bool or defval
        read a value from dictionary which is supposed to be str(bool).
        If key exists returns (value == 'True') else defval
    """
    if key in d:
        return d[key] == 'True'
    else:
        return defval


def set_if_no(dic, key, val):
    'set dic[key] if it is empty'
    if key not in dic:
        dic[key] = val
