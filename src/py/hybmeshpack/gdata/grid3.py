import basic
import srf3
from hybmeshpack.hmcore import g3 as g3core


class Grid3(basic.GeomObject3):
    def __init__(self, cdata):
        super(Grid3, self).__init__()
        # pointer to data stored at c-side
        self.cdata = cdata

    def __del__(self):
        if self.cdata:
            g3core.free_grid3(self.cdata)

    def n_vertices(self):
        return g3core.dims(self.cdata)[0]

    def n_edges(self):
        return g3core.dims(self.cdata)[1]

    def n_faces(self):
        return g3core.dims(self.cdata)[2]

    def surface(self):
        return GridSurface(self.cdata)

    # overriden from GeomObject
    def deepcopy(self):
        c = g3core.deepcopy(self.cdata)
        return Grid3(c)

    def move(self, dx, dy, dz):
        g3core.move(self.cdata, dx, dy, dz)

    def scale(self, xpc, ypc, zpc, x0, y0, z0):
        g3core.scale(self.cdata, xpc, ypc, zpc, x0, y0, z0)

    def point_at(self, index):
        return g3core.point_by_index(index)


    def raw_data(self, what):
        """ -> ctypes arrays
        what = 'btypes' -> [b0, b1, b2, ...]
        """
        return g3core.raw_data(self.cdata, what)

    # overriden from GeomObject3
    def volume(self):
        return g3core.volume(self.cdata)


class GridSurface(srf3.AbstractSurface3):
    def __init__(self, g):
        super(GridSurface, self).__init__()
        self.cgrid = g.cdata

    def __del__(self):
        # cgrid is owned by Grid3 hence do nothing
        pass

    # overriden from GeomObject
    def deepcopy(self):
        return self.surface3()

    def move(self, dx, dy, dz):
        raise NotImplementedError

    def scale(self, xpc, ypc, zpc, x0, y0, z0):
        raise NotImplementedError

    # overriden from GeomObject3
    def volume(self):
        return g3core.volume(self.cdata)

    # overriden from AbstractSurface3
    def n_vertices(self):
        return g3core.bnd_dims(self.cdata)[0]

    def n_edges(self):
        return g3core.bnd_dims(self.cdata)[1]

    def n_faces(self):
        return g3core.bnd_dims(self.cdata)[2]

    def surface3(self):
        return g3core.extract_surface(self.cdata)

    def area(self):
        return g3core.bnd_area(self.cdata)
