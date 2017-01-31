class BndType(object):
    def __init__(self, index):
        self.index = index

    def deepcopy(self):
        return BndType(self.index)


# abstract geometry object classes
class GeomObject(object):
    def __init__(self):
        pass

    def deepcopy(self):
        raise NotImplementedError

    def move(self, dx, dy, dz):
        raise NotImplementedError

    def scale(self, xpc, ypc, zpc, x0, y0, z0):
        raise NotImplementedError

    def point_at(self, index):
        raise NotImplementedError


class GeomObject2(GeomObject):
    def __init__(self):
        super(GeomObject2, self).__init__()

    def move2(self, dx, dy):
        raise NotImplementedError

    def scale2(self, xpc, ypc, x0, y0):
        raise NotImplementedError

    def reflect2(self, x0, y0, x1, y1):
        raise NotImplementedError

    def rotate2(self, x0, y0, angle):
        raise NotImplementedError

    def area(self):
        raise NotImplementedError

    def closest_points(self, pts, proj):
        " proj: 'vertex', 'edge' "
        raise NotImplementedError

    def move(self, dx, dy, dz):
        self.move2(dx, dy)

    def scale(self, xpc, ypc, zpc, x0, y0, z0):
        self.scale2(xpc, ypc, x0, y0)


class GeomObject3(GeomObject):
    def __init__(self):
        super(GeomObject3, self).__init__()

    def volume(self):
        raise NotImplementedError
