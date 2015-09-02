#!/usr/bin/env python
import time
import copy
import sys
from collections import OrderedDict


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




class DefaultCallBack(object):
    """ basic implementation of callback
        Flushes progress to stdout
    """
    pass


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
    dim - dimension of sublist.
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


def unique_name(stem, nms):
    """ ->str. Returns unique name on the basis of stem which
        doesn't present in nms list
    """
    if stem not in nms:
        return stem
    else:
        #find last non digit
        index = len(stem) - 1
        while stem[index].isdigit() and index >= 0:
            index -= 1
        if index < 0:
            raise Exception("Invalid name")
        if index < len(stem) - 1:
            i = int(stem[index + 1:])
            stem = stem[:index + 1]
        else:
            i = 0
        i += 1
        #find the unique name
        while stem + str(i) in nms:
            i += 1
        return stem + str(i)


class NamedList(OrderedDict):
    """ Presents modified collection.OrderedDict class.
        It changes the value of the string key (adds number to the end)
        if it was already presented in key list. F.e.
            a = NamedList([("A",1)])
            a["A"]=2; a["A"]=3
        Leads to {"A":1, "A1":2, "A2":3} dictionary
    """
    def __setitem__(self, key, value):
        ' Sets NamedList[key] = val behavior '
        key = unique_name(key, self.keys())
        OrderedDict.__setitem__(self, key, value)

    def change_value(self, key, value):
        """ Changes value for specified key.
        Default behaviour of OrderedDict.__setitem__
        """
        OrderedDict.__setitem__(self, key, value)

    def change_key(self, key_old, key_new):
        ' Changes the key string without breaking the order '
        if (key_new == key_old or key_old not in self):
            return
        key_new = unique_name(key_new, self.keys())
        restored = []
        while len(self) > 0:
            (k, v) = self.popitem()
            if (k == key_old):
                OrderedDict.__setitem__(self, key_new, v)
                break
            else:
                restored.append((k, v))
        for (k, v) in reversed(restored):
            OrderedDict.__setitem__(self, k, v)

    def get_by_index(self, ind):
        '-> (key, value). Get kv pair for certain index '
        return self.items()[ind]

    def get_by_value(self, val):
        ' -> (index, key) for certain value . '
        ind = self.values().index(val)
        return ind, self.keys()[ind]

    def get_by_key(self, key):
        ' -> (index, value) for certain key '
        ind = self.keys().index(key)
        return ind, self.values()[ind]

    def get(self, ind=None, key=None, val=None):
        '-> (index, key, value). Get full entry information'
        if ind is not None:
            key, val = self.get_by_index(ind)
        elif key is not None:
            ind, val = self.get_by_key(key)
        elif val is not None:
            ind, key = self.get_by_value(val)
        return ind, key, val

    #places (key, value) in a certain position
    def insert(self, ind, key, value):
        'inserts key, value at the specified position'
        key = unique_name(key, self.keys())
        restored = []
        k = len(self) - ind
        #remove all entries after ind
        for i in range(k):
            restored.append(self.popitem())
        #add current item
        OrderedDict.__setitem__(self, key, value)
        #add all removed items
        for (k, v) in reversed(restored):
            OrderedDict.__setitem__(self, k, v)


def compress_int_list(a):
    """ ([i0, i1, i2, ....] -> "i0-i5 i6 i7-i8"

        Converts integer array into string with
        eliminated values which go in a row.
        The procedure provides no internal sorting.
    """
    res = []
    i, n = 0, len(a)
    while i < n:
        istart = i
        while i < n - 1 and a[i + 1] - a[i] == 1:
            i += 1
        res.append([a[istart], a[i]])
        i += 1
    res2 = []
    for r in res:
        if r[0] == r[1]:
            res2.append(str(r[0]))
        else:
            res2.append("%i-%i" % (r[0], r[1]))
    return " ".join(res2)


def int_list_from_compress(s):
    'reciprocal to compress_int_list'
    st = s.split()
    ret = []
    for v in st:
        v = v.split('-', 2)
        if len(v) == 1:
            ret.append(int(v[0]))
        else:
            ist, ien = int(v[0]), int(v[1])
            ret.extend(range(ist, ien + 1))
    return ret


def find(f, seq):
    """Return first item in sequence where f(item) == True or None"""
    for item in seq:
        if f(item):
            return item
    else:
        return None


def dict_readbool(d, key, defval):
    """ -> bool or defval
        read a value from dictionary which is supposed to be str(bool).
        If key exists returns (value == 'True') else defval
    """
    if key in d:
        return d[key] == 'True'
    else:
        return defval


if __name__ == '__main__':
    a = [[1, 2, 4], [2, 4, 5], [2, 4, 2]]
    for [x, y, z] in a:
        print x, y
