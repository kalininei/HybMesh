from hybmeshpack.hmcore import g3 as g3core


class Grid3(object):
    def __init__(self, cdata):
        # pointer to data stored at c-side
        self.cdata = cdata

    def __del__(self):
        if self.cdata:
            g3core.free_g3(self.cdata)

    def n_points(self):
        return g3core.dims(self.cdata)[0]

    def n_edges(self):
        return g3core.dims(self.cdata)[1]

    def n_faces(self):
        return g3core.dims(self.cdata)[2]

    def n_cells(self):
        return g3core.dims(self.cdata)[3]

    def surface(self):
        import srf3
        return srf3.GridSurface(self)

    def deepcopy(self):
        raise NotImplementedError
