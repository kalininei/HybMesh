import basic
import cont2
from hybmeshpack.hmcore import g2 as g2core


class Grid2(basic.GeomObject2):
    def __init__(self, cdata):
        import ctypes as ct
        super(Grid2, self).__init__()
        # pointer to data stored at c-side
        self.cdata = ct.cast(cdata, ct.c_void_p)

    def __del__(self):
        if self.cdata:
            g2core.free_grid2(self.cdata)

    def n_vertices(self):
        return g2core.dims(self.cdata)[0]

    def n_edges(self):
        return g2core.dims(self.cdata)[1]

    def n_cells(self):
        return g2core.dims(self.cdata)[2]

    def cell_types_info(self):
        sizes = self.raw_data("cell_dim")
        ret = {}
        for s in sizes:
            if s in ret:
                ret[s] += 1
            else:
                ret[s] = 1
        return ret

    def skewness(self, threshold):
        return g2core.skewness(self.cdata, threshold)

    def contour(self):
        return GridContour(self.cdata)

    def raw_data(self, what):
        """ returns ctypes arrays
        what = 'vert' -> [x0, y0, x1, y1, ...]
        what = 'edge_vert' -> [e0p0, e0p1, e1p0, e1p1, ...];
        what = 'edge_cell' -> [e0left, e0right, e1left, e1right, ...];
        what = 'cell_dim' -> [sz0, sz1, ... ]
        what = 'cell_vert' -> [v0, v1, .., vn, v0, ....]; need cellsizes
        what = 'cell_edge' -> [e0, e1, ..., ]; need cellsizes
        what = 'cell_center' -> [x0, y0, x1, y1, ....]
        what = 'bnd' -> [e0, e1, ...]
        what = 'bt' -> [b0, b1, ...]
        what = 'bnd_bt' -> [e0, b0, e1, b1, ...]
        """
        return g2core.raw_data(self.cdata, what)


    def assign_boundary_type(self, bt):
        return g2core.assign_boundary_types(self.cdata, bt)


    # overriden from GeomObject
    def dims(self):
        return g2core.dims(self.cdata)

    def deepcopy(self):
        return Grid2(g2core.deepcopy(self.cdata))

    def point_at(self, index):
        return g2core.point_by_index(index)

    # overriden from GeomObject2
    def move2(self, dx, dy):
        g2core.move(self.cdata, dx, dy)

    def scale2(self, xpc, ypc, x0, y0):
        g2core.scale(self.cdata, xpc, ypc, x0, y0)

    def reflect2(self, x0, y0, x1, y1):
        g2core.reflect(self.cdata, x0, y0, x1, y1)

    def rotate2(self, x0, y0, angle):
        g2core.rotate(self.cdata, x0, y0, angle)

    def area(self):
        return g2core.area(self.cdata)

    def closest_points(self, pts, proj):
        return g2core.closest_points(self.cdata, pts, proj)


class GridContour(cont2.AbstractContour2):
    def __init__(self, gcdata):
        super(GridContour, self).__init__()
        self.cdata = gcdata

    def __del__(self):
        # cgrid is owned by Grid3 hence do nothing
        pass

    # overriden from GeomObject
    def deepcopy(self):
        return self.contour2()

    # overriden from GeomObject2
    def move2(self, dx, dy):
        raise NotImplementedError

    def scale2(self, xpc, ypc, x0, y0):
        raise NotImplementedError

    def reflect2(self, x0, y0, x1, y1):
        raise NotImplementedError

    def rotate2(self, x0, y0, angle):
        raise NotImplementedError

    def area(self):
        return g2core.area(self.cdata)

    def dims(self):
        return g2core.bnd_dims(self.cdata)

    # overriden from AbstractContour2
    def n_points(self):
        return g2core.bnd_dims(self.cdata)[0]

    def n_edges(self):
        return g2core.bnd_dims(self.cdata)[1]

    def contour2(self):
        return cont2.Contour2(g2core.extract_contour(self.cdata))

    def length(self):
        return g2core.bnd_length(self.cdata)
