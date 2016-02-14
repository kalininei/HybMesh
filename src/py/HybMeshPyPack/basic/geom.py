#!/usr/bin/env python
"basic geometry"
import copy
import math


class Point2(object):
    'point in (x,y) plane'
    def __init__(self, x, y):
        self.x = float(x)
        self.y = float(y)

    def __str__(self):
        return str(self.x) + " " + str(self.y)

    def dist(self, p1):
        "-> distance between self and p1 "
        x, y = p1.x - self.x, p1.y - self.y
        return math.sqrt(x * x + y * y)

    def dist0(self):
        " -> distance between self and (0, 0) "
        x, y = self.x, self.y
        return math.sqrt(x * x + y * y)

    @classmethod
    def fromstring(cls, s):
        s2 = s.split()
        return cls(float(s2[0]), float(s2[1]))


def angle_3pnt(p1, p2, p3):
    '-> [0, 2*pi]. Get the angle from 3 consecutive points'
    v1x, v1y = p1.x - p2.x, p1.y - p2.y
    v2x, v2y = p3.x - p2.x, p3.y - p2.y
    dot = v1x * v2x + v1y * v2y
    cross = v1x * v2y - v1y * v2x
    a = math.atan2(cross, dot)
    if a < 0:
        a += 2 * math.pi
    return a


class GeomStruct(object):
    """ Geometric structure class - parent for all complicated
        geometry objects
    """
    __num = 0

    @staticmethod
    def __new_struct_id():
        GeomStruct.__num += 1
        return GeomStruct.__num

    def __init__(self):
        ' sets a unique id for the structure '
        self.__id = self.__new_struct_id()
        #geometry change substribers: set of [function(void)]
        self.__subscribers_change_geom = set()

    def deepcopy(self):
        """ returns deepcopied GeomStruct object
            with unique id
            Geometry subscribers are not copied
        """
        __tmp = self.__subscribers_change_geom
        self.__subscribers_change_geom = set()
        ret = copy.deepcopy(self)
        ret.__id = GeomStruct.__new_struct_id()
        self.__subscribers_change_geom = __tmp
        return ret

    def get_id(self):
        "->int. Returns unique id of object"
        return self.__id

    def add_subscriber_change_geom(self, func):
        """ (void function(void) func) -> void
            Adds a function which will be called if
            geometry changes
        """
        self.__subscribers_change_geom.add(func)

    def remove_subsriber_change_geom(self, func):
        """ Removes subscription from grid geometry event """
        self.__subscribers_change_geom.remove(func)

    def clear_subscribers_change_geom(self):
        """ Clear geometry change subscribers """
        self.__subscribers_change_geom.clear()

    def _geom_changed(self):
        "send geometry change message to subscribers"
        [f() for f in self.__subscribers_change_geom]


def rotate_points(pnts, x0, y0, angle):
    """ rotates a set of points pnts around [x0, y0]
        on an angle (degrees).
        returns list of rotates points
    """
    angle = angle / 180.0 * math.pi
    sina, cosa = math.sin(angle), math.cos(angle)
    newpoints = map(lambda p: Point2(
                            (p.x - x0) * cosa - (p.y - y0) * sina + x0,
                            (p.x - x0) * sina + (p.y - y0) * cosa + y0),
                pnts)
    return newpoints


def scale_points(pnts, p0, xpc, ypc):
    """ scale points list using p0 as reference point
        and xpc% and ypc% as scaling procentages
    """
    sx, sy = xpc / 100.0, ypc / 100.0
    newpoints = map(lambda p: Point2((p.x - p0.x) * sx + p0.x,
        (p.y - p0.y) * sy + p0.y), pnts)
    return newpoints


def div_range(a, b, n, k=1):
    """ -> [a, ...., b] - list of floats

    Divides [a, b] range into n subsections.
    k is the refinement coefficient which works like:
        h[i] = k * h[i-1]
    """
    if a != 0:
        return [x + a for x in div_range(0, b - a, n, k)]
    #now a = 0, b = length of the section
    st = [k ** i for i in range(n)]
    a0 = float(b) / sum(st)
    ret = [0]
    for x in st:
        ret.append(ret[-1] + a0 * x)
    return ret


class Point2SetStruct(GeomStruct):
    "Geometry object defined as the set of 2D points"
    def __init__(self):
        super(Point2SetStruct, self).__init__()
        #points coords
        self.points = []

    def n_points(self):
        ' -> number of points '
        return len(self.points)

    def bounding_box(self):
        """ -> (Point2, Point2)

        returns bottom left and top right point of grid bounding box
        """
        x = [p.x for p in self.points]
        y = [p.y for p in self.points]
        return Point2(min(x), min(y)), Point2(max(x), max(y))

    # ---- Transformations
    def move(self, dx, dy):
        for p in self.points:
            p.x += dx
            p.y += dy
        self._geom_changed()

    def rotate(self, x0, y0, angle):
        pnew = rotate_points(self.points, x0, y0, angle)
        for p, p2 in zip(self.points, pnew):
            p.x, p.y = p2.x, p2.y
        self._geom_changed()

    def scale(self, p0, xpc, ypc):
        """ scale the grid using p0 as reference point
            and xpc% and ypc% as scaling procentages
        """
        pnew = scale_points(self.points, p0, xpc, ypc)
        for p, p2 in zip(self.points, pnew):
            p.x, p.y = p2.x, p2.y
        self._geom_changed()

    def unscale(self, p0, xpc, ypc):
        """ unscale the grid after scaling with p0 as
            reference point and xpc%, ypc% as scaling procentages
        """
        self.scale(p0, 10000.0 / xpc, 10000.0 / ypc)

    def points_to_str(self, separator=" "):
        "->str. Points to string"
        return separator.join(map(str, self.points))

    def fill_points_from_str(self, data):
        """fill points list from string "x0 y0 x1 y1 ..." """
        self.points = self.__str_to_pnts(data)

    def append_points_from_str(self, data):
        """append points to list from string "x0 y0 x1 y1 ..." """
        self.points.extend(self.__str_to_pnts(data))

    @staticmethod
    def __str_to_pnts(data):
        it = iter(map(float, data.split()))
        return map(Point2, it, it)

if __name__ == "__main__":
    print div_range(0, 5, 10, 1)
