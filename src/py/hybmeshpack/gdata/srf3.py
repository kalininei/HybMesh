import basic
from hybmeshpack.hmcore import s3 as s3core


class AbstractSurface3(basic.GeomObject3):
    def __init__(self):
        super(AbstractSurface3, self).__init__()

    def n_vertices(self):
        raise NotImplementedError

    def n_edges(self):
        raise NotImplementedError

    def n_faces(self):
        raise NotImplementedError

    def raw_data(self, what):
        """ -> ctypes arrays
        what = 'vert' -> [x0, y0, z0, x1, y1, z1, ...]
        what = 'edge_vert' -> [e0v0, e0v1, e1v0, e1v1, ...]
        what = 'face_dim' -> [f0dim, f1dim]
        what = 'face_edge' -> [f0e0, f0e1, ..., f1e0,..]
        what = 'face_vert' -> [f0v0, f0v1, ..., f1v0,..]
        what = 'bt' -> [b0, b1, b2, ...]
        what = 'face_centers' -> [x0, y0, z0, x1, ...]
        what = 'face_areas' -> [area1, area2, ...]
        """
        return self.surface3()._raw_data(what)

    def surface3(self):
        "returns Surface3 object: self or deepcopied"
        raise NotImplementedError

    def area(self):
        raise NotImplementedError


class Surface3(AbstractSurface3):
    def __init__(self, cdata):
        import ctypes as ct
        super(Surface3, self).__init__()
        # pointer to data stored at c-side
        self.cdata = ct.cast(cdata, ct.c_void_p)

    def __del__(self):
        if self.cdata:
            s3core.free_surf3(self.cdata)

    def _raw_data(self, what):
        return s3core.raw_data(self.cdata, what)

    # overriden from GeomObject
    def deepcopy(self):
        c = s3core.deepcopy(self.cdata)
        return Surface3(c)

    def move(self, dx, dy, dz):
        s3core.move(self.cdata, dx, dy, dz)

    def scale(self, xpc, ypc, zpc, x0, y0, z0):
        s3core.scale(self.cdata, xpc, ypc, zpc, x0, y0, z0)

    def point_at(self, index):
        return s3core.point_by_index(index)

    def dims(self):
        return s3core.dims(self.cdata)

    def assign_boundary_type(self, bt):
        return s3core.assign_boundary_types(self.cdata, bt)

    # overriden from GeomObject3
    def volume(self):
        return s3core.volume(self.cdata)

    # overriden from AbstractSurface
    def n_vertices(self):
        return s3core.dims(self.cdata)[0]

    def n_edges(self):
        return s3core.dims(self.cdata)[1]

    def n_faces(self):
        return s3core.dims(self.cdata)[2]

    def surface3(self):
        return self

    def area(self):
        return s3core.area(self.cdata)
