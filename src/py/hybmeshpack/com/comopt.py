import hybmeshpack.basic as bp
import copy


# represent an option which has no default value
class _NoDefault:
    pass

no_default = _NoDefault()


class BasicOption(object):
    """ Represent option transformation before
        str(Command.option) will be called.
        For types like int, float, str, [int], {int: int}, ...
        Base class will be enough.
        For compound ones child class is needed.
        See example of child implementation for gdata.Point2 option below.
    """
    def __init__(self, tp=None, default=no_default):
        self.tp = tp
        self.default = default

    def has_default(self):
        return not isinstance(self.default, _NoDefault)

    def serial(self, v):
        'before writing'
        return v

    def unserial(self, s):
        'after reading'
        if self.tp is None:
            return s
        else:
            return self.tp(s)


class ListOfOptions(BasicOption):
    "list of BasicOption data"
    def __init__(self, tp, default=no_default):
        'tp - another BasicOption derived class'
        super(ListOfOptions, self).__init__(tp, default)

    def serial(self, v):
        return [self.tp.serial(x) for x in v]

    def unserial(self, v):
        return [self.tp.unserial(x) for x in v]


class NoneOr(BasicOption):
    "some data or None"
    def __init__(self, tp, default=no_default):
        'tp - another BasicOption object'
        super(NoneOr, self).__init__(tp, default)

    def serial(self, v):
        return None if v is None else self.tp.serial(v)

    def unserial(self, v):
        if v is None or v == 'None':
            return None
        else:
            return self.tp.unserial(v)


class SubDictOption(BasicOption):
    def __init__(self, **kwargs):
        """ Present dictionary of options:
            name1: value1,
            name2: value2,...

            kwargs are  name->BasicOption()
        """
        self.args = kwargs

    def serial(self, v):
        a = {}
        for k, val in v.items():
            if k in self.args:
                a[k] = self.args[k].serial(val)
        return a

    def unserial(self, v):
        a = {}
        for k, tp in self.args.items():
            if k in v:
                a[k] = tp.unserial(v[k])
            elif self.has_default():
                a[k] = tp.default
            else:
                raise KeyError("Mandatory option field is absent")
        return a


class Point2Option(BasicOption):
    def __init__(self, default):
        super(Point2Option, self).__init__(None, default)

    def serial(self, v):
        return (v[0], v[1])

    def unserial(self, v):
        return (float(v[0]), float(v[1]))


class Point3Option(BasicOption):
    def __init__(self, default):
        super(Point3Option, self).__init__(None, default)

    def serial(self, v):
        return (v[0], v[1], v[2])

    def unserial(self, v):
        return (float(v[0]), float(v[1]), float(v[2]))


class BoolOption(BasicOption):
    def __init__(self, default):
        super(BoolOption, self).__init__(bool, default)

    def serial(self, v):
        return v

    def unserial(self, v):
        if v in [0, False, '0', 'False', 'false']:
            return False
        else:
            return True


class ListCompressedOption(BasicOption):
    def __init__(self, default):
        super(ListCompressedOption, self).__init__(None, default)

    def serial(self, v):
        return bp.compress_int_list(v)

    def unserial(self, v):
        return bp.int_list_from_compress(v)


class BoundaryPickerOption(BasicOption):
    """
    {'name': contour name,
      bnd_type1: [list of edges indicies],
      bnd_type2: [list of edges indicies], ...}
    """
    def __init__(self):
        super(BoundaryPickerOption, self).__init__()

    def serial(self, v):
        'before writing'
        a = copy.deepcopy(v)
        for k, v in a.items():
            if k != 'name':
                a[k] = bp.compress_int_list(a[k])
        return a

    def unserial(self, s):
        'after reading'
        a = copy.deepcopy(s)
        for k, v in a.items():
            if k != 'name':
                a[k] = bp.int_list_from_compress(a[k])
        return a
