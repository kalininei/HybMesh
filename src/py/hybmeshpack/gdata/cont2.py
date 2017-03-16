import basic
from hybmeshpack.hmcore import c2 as c2core


class AbstractContour2(basic.GeomObject2):
    def __init__(self):
        super(AbstractContour2, self).__init__()

    def n_vertices(self):
        raise NotImplementedError

    def n_edges(self):
        raise NotImplementedError

    def length(self):
        raise NotImplementedError

    def contour2(self):
        "returns Contour2 object: self or deepcopied"
        raise NotImplementedError

    def end_points(self):
        """ -> p0, p1. Where p0, p1 are [x, y] points.
            returns None, None for compaund contours
            and p0 = p1 for closed contours
        """
        c = self.contour2()
        cont_type = c2core.contour_type(c.cdata)
        if cont_type in ["compound", "mdomain"]:
            return None, None
        else:
            if cont_type == "closed":
                p1 = c.point_at(0)
                p2 = p1
            else:
                p1, p2 = c2core.end_points(c.cdata)
            return p1, p2

    def raw_data(self, what):
        """ -> ctypes arrays
        what = 'vert' -> [x0, y0, x1, y1, ...]
        what = 'edge_vert' -> [v0, v1, v0, v1, ...]
        what = 'bt' -> [b0, b1, b2, .... ]
        what = 'edge_lengths' -> [len1, len2, len3, ....]
        """
        return self.contour2()._raw_data(what)

    # overriden from GeomObject
    def point_at(self, index):
        return self.contour2()._point_at(index)

    def closest_points(self, p0, proj):
        return self.contour2()._closest_points(p0, proj)


class Contour2(AbstractContour2):
    def __init__(self, cdata):
        import ctypes as ct
        super(Contour2, self).__init__()
        # pointer to data stored at c-side
        self.cdata = ct.cast(cdata, ct.c_void_p)

    def __del__(self):
        if self.cdata:
            c2core.free_cont2(self.cdata)

    @classmethod
    def empty(cls):
        p = c2core.build_from_points([], False, [])
        return cls(p)

    # overriden from GeomObject
    def deepcopy(self):
        c = c2core.deepcopy(self.cdata)
        return Contour2(c)

    def dims(self):
        return c2core.dims(self.cdata)

    def assign_boundary_type(self, bt):
        return c2core.assign_boundary_types(self.cdata, bt)

    # overriden from GeomObject2
    def move2(self, dx, dy):
        c2core.move(self.cdata, dx, dy)

    def scale2(self, xpc, ypc, x0, y0):
        c2core.scale(self.cdata, xpc, ypc, x0, y0)

    def reflect2(self, x0, y0, x1, y1):
        c2core.reflect(self.cdata, x0, y0, x1, y1)

    def rotate2(self, x0, y0, angle):
        c2core.rotate(self.cdata, x0, y0, angle)

    def area(self):
        return c2core.area(self.cdata)

    def _closest_points(self, pts, proj):
        "proj = 'vertex'/'edge'"
        return c2core.closest_points(self.cdata, pts, proj)

    def _raw_data(self, what):
        return c2core.raw_data(self.cdata, what)

    def _point_at(self, index):
        return c2core.point_by_index(self.cdata, index)

    # overriden from AbstractContour2
    def n_vertices(self):
        return c2core.dims(self.cdata)[0]

    def n_edges(self):
        return c2core.dims(self.cdata)[1]

    def contour2(self):
        return self

    def length(self):
        return c2core.length(self.cdata)


def _meas2cont(cont, p0):
    [p1] = cont.closest_points([p0], 'edge')
    p1[0] -= p0[0]
    p1[1] -= p0[1]
    return p1[0] * p1[0] + p1[1] * p1[1]


def closest_contour(contlist, pnt):
    if len(contlist) == 0:
        return None
    ret = contlist[0]
    mindist = _meas2cont(contlist[0], pnt)
    for i in range(1, len(contlist)):
        d = _meas2cont(contlist[i], pnt)
        if d < mindist:
            ret = contlist[i]
            mindist = d
    return ret
