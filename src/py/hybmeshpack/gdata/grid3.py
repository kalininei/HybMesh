import basic
import srf3
from hybmeshpack.hmcore import g3 as g3core


class Grid3(basic.GeomObject3):
    def __init__(self, cdata):
        import ctypes as ct
        super(Grid3, self).__init__()
        # pointer to data stored at c-side
        self.cdata = ct.cast(cdata, ct.c_void_p)

    def __del__(self):
        if self.cdata:
            g3core.free_grid3(self.cdata)

    def n_vertices(self):
        return g3core.dims(self.cdata)[0]

    def n_edges(self):
        return g3core.dims(self.cdata)[1]

    def n_faces(self):
        return g3core.dims(self.cdata)[2]

    def n_cells(self):
        return g3core.dims(self.cdata)[3]

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

    def dims(self):
        return g3core.dims(self.cdata)

    def raw_data(self, what):
        """ -> ctypes arrays
        what = 'vert' -> [x0, y0, z0, x1, ..]
        what = 'edge_vert' -> [e0v0, e0v1, ..]
        what = 'face_dim' -> [f0dim, f1dim, ...]
        what = 'face_edge' -> [f0e0, f0e1, ..., f1e0, ..]
        what = 'face_vert' -> [f0v0, f0v1, f0v2, .., f1v0]
        what = 'face_cell' -> [f0left, f0right, f1left, ...]
        what = 'cell_fdim' -> [c0nfaces, c1nfaces, ...]
        what = 'cell_vdim' -> [c0nvert, c1nvert, ...]
        what = 'cell_face' -> [c0f0, c0f1, c0f2, .., c1f0]
        what = 'cell_vert' -> [c0v0, c0v1, ..., c0v0,]
        what = 'bt'  -> [f0bt, f1bt, ....]
        what = 'bnd' -> [bf0, bf1, ...]
        what = 'bnd_bt' -> [bf0, bt0, bf1, bt1, ...]
        """
        return g3core.raw_data(self.cdata, what)

    def assign_boundary_type(self, bt):
        return g3core.assign_boundary_types(self.cdata, bt)

    # overriden from GeomObject3
    def volume(self):
        return g3core.volume(self.cdata)


class GridSurface(srf3.AbstractSurface3):
    def __init__(self, gcdata):
        super(GridSurface, self).__init__()
        self.cdata = gcdata

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
    def dims(self):
        return g3core.bnd_dims(self.cdata)

    def n_vertices(self):
        return g3core.bnd_dims(self.cdata)[0]

    def n_edges(self):
        return g3core.bnd_dims(self.cdata)[1]

    def n_faces(self):
        return g3core.bnd_dims(self.cdata)[2]

    def surface3(self):
        return srf3.Surface3(g3core.extract_surface(self.cdata))

    def area(self):
        return g3core.bnd_area(self.cdata)
