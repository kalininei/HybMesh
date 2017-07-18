import os.path
import numbers
import inspect
from hybmeshpack.hmscript import flow, InvalidArgument


def icheck(iarg, cond):
    """ check input argument #iarg of the current frame function.
        cond - an object of any class defined in datachecks module.
    """
    args, _, _, vals = inspect.getargvalues(inspect.stack()[1][0])
    nm = args[iarg]
    value = vals[nm]
    try:
        cond.check_value(value)
    except InvalidArgument as e:
        f = inspect.stack()[1][3]
        s1 = "%s=" % nm + repr(value)
        s2 = "Wrong argument #%i for '%s' given as \n\t%s" % (iarg+1, f, s1)
        s3 = str(e)
        raise InvalidArgument('\n'.join([s2, s3]))


class Any(object):
    """any value"""
    def check_value(self, value):
        pass


class NoneOr(object):
    """ none or some value"""
    def __init__(self, ob):
        self.ob = ob

    def check_value(self, value):
        msg = ""
        if value is not None:
            try:
                self.ob.check_value(value)
                return
            except InvalidArgument as e:
                msg = str(e)
            raise InvalidArgument("should be None or " + msg)


class Tuple(object):
    """ tuple object with predefined length"""
    def __init__(self, *args):
        self.ob = args

    def check_value(self, value):
        if not isinstance(value, tuple):
            raise InvalidArgument("should be a tuple")

        if len(value) != len(self.ob):
            raise InvalidArgument("length should equal %i" % len(self.ob))

        for i, v in enumerate(value):
            try:
                self.ob[i].check_value(v)
            except InvalidArgument as e:
                raise InvalidArgument("tuple entry [%i] " % i + str(e))


class Bool(object):
    """bolean value"""
    def __init__(self, **kwargs):
        self.conditions = kwargs

    def check_value(self, value):
        if not isinstance(value, bool):
            raise InvalidArgument("should be boolean")


class UInt(object):
    """non-negative integer"""
    def __init__(self, **kwargs):
        """ minv: int
            maxv: int
        """
        self.conditions = kwargs

    def __minv(self, value):
        if value < self.conditions['minv']:
            raise InvalidArgument(
                "should be greater than %i" % (self.conditions['minv'] - 1))

    def __maxv(self, value):
        if value > self.conditions['maxv']:
            raise InvalidArgument(
                "should be lower than %i" % (self.conditions['maxv'] + 1))

    def check_value(self, value):
        if not isinstance(value, numbers.Integral) or value < 0:
            raise InvalidArgument("should be a non-negative integer")

        if 'maxv' in self.conditions:
            self.__maxv(value)
        if 'minv' in self.conditions:
            self.__minv(value)


class Float(object):
    """float number"""
    def __init__(self, **kwargs):
        """ grthan: float
            within: [float, float, "[]" or "()" or "(]" or "[)" ]
        """
        self.conditions = kwargs

    def __grthan(self, value):
        if value < self.conditions['grthan']:
            raise InvalidArgument(
                "should be greater than {:.16g}".format(
                    self.conditions['grthan']))

    def __within(self, value):
        c = self.conditions['within']
        leq = c[2][0] == "["
        req = c[2][1] == "]"
        bad = False
        if leq and value < c[0]:
            bad = True
        elif not leq and value <= c[0]:
            bad = True
        elif req and value > c[1]:
            bad = True
        elif not req and value >= c[1]:
            bad = True
        if bad:
            raise InvalidArgument(
                "should be within {}{:.16g}, {:.16g}{}".format(
                    c[2][0], c[0], c[1], c[2][1]))

    def check_value(self, value):
        if not isinstance(value, numbers.Real):
            raise InvalidArgument("should be a float number")

        if 'grthan' in self.conditions:
            self.__grthan(value)

        if 'within' in self.conditions:
            self.__within(value)


class List(object):
    """ List of objects"""
    def __init__(self, ob, **kwargs):
        """ llen: int,
            minlen: int,
            unique: bool
        """
        self.ob = ob
        self.conditions = kwargs

    def __unique(self, value):
        for i in range(len(value)):
            for j in range(i + 1, len(value)):
                if value[i] == value[j]:
                    raise InvalidArgument("list entries should be unique")

    def __minlen(self, value):
        ml = self.conditions['minlen']
        if len(value) < ml:
            raise InvalidArgument(
                "list length should be not less than %i" % ml)

    def __llen(self, value):
        ml = self.conditions['llen']
        if len(value) != ml:
            raise InvalidArgument(
                "list length should equal %i" % ml)

    def check_value(self, value):
        if not isinstance(value, list):
            raise InvalidArgument("should be a list")

        for i, v in enumerate(value):
            try:
                self.ob.check_value(v)
            except InvalidArgument as e:
                raise InvalidArgument("list entry [%i] " % i + str(e))

        if 'llen' in self.conditions:
            self.__llen(value)
        if 'minlen' in self.conditions:
            self.__minlen(value)
        if 'unique' in self.conditions:
            self.__unique(value)


class UList(List):
    """ List of unique objects"""
    def __init__(self, ob, **kwargs):
        super(UList, self).__init__(ob, **kwargs)
        self.conditions['unique'] = True


class CompoundList(object):
    """ Compound list of objects of different types"""
    def __init__(self, *args, **kwargs):
        """ minlen: int,
        """
        self.ob = args
        self.conditions = kwargs

    def __minlen(self, value):
        ml = self.conditions['minlen']
        if len(value) < ml:
            raise InvalidArgument(
                "list length should be not less than %i" % ml)

    def __llen(self, value):
        ml = self.conditions['llen']
        if len(value) != ml:
            raise InvalidArgument(
                "list length should equal %i" % ml)

    def check_value(self, value):
        if not isinstance(value, list):
            raise InvalidArgument("should be a list")
        if len(value) % len(self.ob) != 0:
            raise InvalidArgument(
                "list length should be divisible by %i" % len(self.ob))

        for i, v in enumerate(value):
            try:
                n = i % len(self.ob)
                self.ob[n].check_value(v)
            except InvalidArgument as e:
                raise InvalidArgument("list entry [%i] " % i + str(e))

        if 'minlen' in self.conditions:
            self.__minlen(value)
        if 'llen' in self.conditions:
            self.__llen(value)


class Or(object):
    """Choice between two conditions"""
    def __init__(self, c1, c2):
        self.c1 = c1
        self.c2 = c2

    def check_value(self, value):
        msg1, msg2 = None, None
        try:
            self.c1.check_value(value)
            return
        except InvalidArgument as e:
            msg1 = str(e)
        try:
            self.c2.check_value(value)
            return
        except InvalidArgument as e:
            msg2 = str(e)

        raise InvalidArgument("either " + msg1 + " or " + msg2)


class UListOr1(Or):
    """unique list of values or single value"""
    def __init__(self, c1, **kwargs):
        super(UListOr1, self).__init__(UList(c1, **kwargs), c1)


class ListOr1(Or):
    """list of values or single value"""
    def __init__(self, c1, **kwargs):
        super(ListOr1, self).__init__(List(c1, **kwargs), c1)


class Grid2D(object):
    """grid2d identifier"""
    def check_value(self, value):
        if value not in flow.receiver.get_grid2_names():
            raise InvalidArgument('should be a 2D grid identifier')


class Grid3D(object):
    """grid3d identifier"""
    def check_value(self, value):
        if value not in flow.receiver.get_grid3_names():
            raise InvalidArgument('should be a 3D grid identifier')


class Cont2D(object):
    """contour2D identifier"""
    def check_value(self, value):
        if value not in flow.receiver.get_contour2_names():
            raise InvalidArgument('should be a 2D contour identifier')


class Surf3D(object):
    """surface3d identifier"""
    def check_value(self, value):
        if value not in flow.receiver.get_surface3_names():
            raise InvalidArgument('should be a 3D surface identifier')


class ACont2D(object):
    """grid2d or contour2d identifier"""
    def check_value(self, value):
        if value not in flow.receiver.get_grid2_names() and\
                value not in flow.receiver.get_contour2_names():
            raise InvalidArgument('should be a 2D grid or contour identifier')


class ASurf3D(object):
    """grid3d or surface identifier"""
    def check_value(self, value):
        if value not in flow.receiver.get_grid3_names() and\
                value not in flow.receiver.get_surface3_names():
            raise InvalidArgument('should be a 3D grid or surface identifier')


class AObject(object):
    """geometrical object identifier"""
    def check_value(self, value):
        if value not in flow.receiver.get_names():
            raise InvalidArgument('should be a geometrical object identifier')


class OneOf(object):
    """choice between predefined values"""
    def __init__(self, *args):
        self.predef = args

    def check_value(self, value):
        if value not in self.predef:
            raise InvalidArgument("should be one of " + str(self.predef))


class Point2D(object):
    """2d point"""
    def __init__(self, **kwargs):
        """ grthan: [x, y]
            noteq: [[x, y], [x, y], ...]
        """
        self.conditions = kwargs

    def __grthan(self, value):
        p0 = self.conditions['grthan']
        if value[0] <= p0[0] or value[1] <= p0[1]:
            raise InvalidArgument("should be greater than " + str(p0))

    def __noteq(self, value):
        plist = self.conditions['noteq']
        for p in plist:
            if value == p:
                raise InvalidArgument("should not equal " + str(p))

    def check_value(self, value):
        if not isinstance(value, list) or len(value) != 2 or\
                not isinstance(value[0], numbers.Real) or\
                not isinstance(value[1], numbers.Real):
            raise InvalidArgument("should be a 2D point as [x, y]")

        if 'grthan' in self.conditions:
            self.__grthan(value)
        if 'noteq' in self.conditions:
            self.__noteq(value)


class Point3D(object):
    """3d point"""
    def __init__(self, **kwargs):
        """ noteq: [[x, y, z], [x, y, z], ...]
        """
        self.conditions = kwargs

    def __noteq(self, value):
        plist = self.conditions['noteq']
        for p in plist:
            if value == p:
                raise InvalidArgument("should not equal " + str(p))

    def check_value(self, value):
        if not isinstance(value, list) or len(value) != 3 or\
                not isinstance(value[0], numbers.Real) or\
                not isinstance(value[1], numbers.Real):
            raise InvalidArgument("should be a 3D point as [x, y, z]")

        if 'noteq' in self.conditions:
            self.__noteq(value)


class APoint(object):
    """2d or 3d point"""
    def __init__(self, **kwargs):
        """ noteq: [[x, y, (z)], [x, y, (z)], ...]
        """
        self.conditions = kwargs

    def __noteq(self, value):
        plist = self.conditions['noteq']
        if len(value) == 2:
            value = [value[0], value[1], 0.]
        for p in plist:
            p2 = p
            if len(p2) == 2:
                p2 = [p[0], p[1], 0.]
            if value == p2:
                raise InvalidArgument("should not equal " + str(p))

    def check_value(self, value):
        if not isinstance(value, list) or len(value) not in [2, 3] or\
                not isinstance(value[0], numbers.Real) or\
                not isinstance(value[1], numbers.Real) or\
                (len(value) == 3 and not isinstance(value[2], numbers.Real)):
            raise InvalidArgument("should be a point as [x, y] or [x, y, z]")

        if 'noteq' in self.conditions:
            self.__noteq(value)


class ZType(UInt):
    """ zone type identifier """
    def __init__(self):
        super(ZType, self).__init__()


class String(object):
    def check_value(self, value):
        if not isinstance(value, str):
            raise InvalidArgument("should be a string")


class ExistingFile(object):
    def check_value(self, value):
        try:
            r = os.path.isfile(value)
            if not r:
                raise
            else:
                return
        except:
            raise InvalidArgument("should be an existing file")


class IncList(List):
    def __init__(self, ob, **kwargs):
        super(IncList, self).__init__(ob, **kwargs)

    def __startfrom(self, value):
        v = self.conditions['startfrom']
        if value[0] != v:
            raise InvalidArgument("should start from " + repr(v))

    def check_value(self, value):
        super(IncList, self).check_value(value)

        for i in range(1, len(value)):
            if value[i] <= value[i - 1]:
                raise InvalidArgument("should be a list of increasing values")

        if len(value) >= 1 and "startfrom" in self.conditions:
            self.__startfrom(value)


class Func(object):
    def __init__(self, **kwargs):
        """ nargs: int
        """
        self.conditions = kwargs

    def __nargs(self, value):
        garg = -1
        try:
            garg = len(inspect.getargspec(value)[0])
        except:
            try:
                garg = len(inspect.getargspec(value.__call__)[0]) - 1
            except:
                pass
        if garg != self.conditions['nargs']:
            raise InvalidArgument("should be a function of %i arguments" %
                                  self.conditions['nargs'])

    def check_value(self, value):
        if not callable(value):
            raise InvalidArgument("should be a callable object")

        if "nargs" in self.conditions:
            self.__nargs(value)


class Dict(object):
    def __init__(self, cls1, cls2):
        self.cls1 = cls1
        self.cls2 = cls2

    def check_value(self, value):
        if not isinstance(value, dict):
            raise InvalidArgument("should be a dictionary")

        try:
            for k in value.keys():
                self.cls1.check_value(k)
        except InvalidArgument as e:
            raise InvalidArgument("each dict key " + str(e))

        try:
            for v in value.values():
                self.cls2.check_value(v)
        except InvalidArgument as e:
            raise InvalidArgument("each dict value " + str(e))
