import copy
import math


class GeomStruct(object):
    "Geometric structure class - parent for all structures"
    __num = 0

    @staticmethod
    def __new_struct_id():
        GeomStruct.__num += 1
        return GeomStruct.__num

    def deepcopy(self):
        """ returns deepcopied GeomStruct object
            with unique id
        """
        ret = copy.deepcopy(self)
        ret.__id = GeomStruct.__new_struct_id()
        return ret

    def __init__(self):
        ' sets a unique id for the structure '
        self.__id = self.__new_struct_id()

    def get_id(self):
        return self.__id


class Point2(object):
    ' point in (x,y) plane '
    def __init__(self, x, y):
        self.x = x
        self.y = y

    def __str__(self):
        return str(self.x) + " " + str(self.y)

    def dist(self, p1):
        " -> distance between self and p1 "
        x, y = p1.x - self.x, p1.y - self.y
        return math.sqrt(x * x + y * y)


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


def div_range(a, b, n, k=1):
    """ -> [a, ...., b] - list of floats

    Divides [a, b] range into n subsections.
    k is the refinement coefficient which works like:
        h[i] = k * h[i-1]
    """
    if a > 0:
        return [x + a for x in div_range(0, b - a, n, k)]
    #now a = 0, b = length of the section
    st = [k ** i for i in range(n)]
    a0 = b / sum(st)
    ret = [0]
    for x in st:
        ret.append(ret[-1] + a0 * x)
    return ret
