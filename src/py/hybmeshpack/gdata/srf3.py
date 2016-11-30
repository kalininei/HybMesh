from hybmeshpack.hmcore import s3 as s3core


class _AbstractSurface3(object):
    def __init__(self):
        pass

    def n_points(self):
        raise NotImplementedError

    def n_edges(self):
        raise NotImplementedError

    def n_faces(self):
        raise NotImplementedError

    def btypes(self):
        raise NotImplementedError

    def deepcopy(self):
        raise NotImplementedError


class Surface3(_AbstractSurface3):
    def __init__(self, cdata):
        super(Surface3, self).__init__()
        # pointer to data stored at c-side
        self.cdata = cdata

    def n_points(self):
        return s3core.dims(self.cdata)[0]

    def n_edges(self):
        return s3core.dims(self.cdata)[1]

    def n_faces(self):
        return s3core.dims(self.cdata)[2]

    def btypes(self):
        return s3core.btypes(self.cdata)

    def __del__(self):
        if self.cdata:
            s3core.free_srf3(self.cdata)

    def shallow_separate(self):
        slist = s3core.extract_subsurfaces(self.cdata)
        ret = []
        for s in slist:
            ret.append(Surface3(s))
        return ret


class GridSurface(_AbstractSurface3):
    def __init__(self, g):
        super(GridSurface, self).__init__()
        self.cgrid = g.cdata

    def __del__(self):
        # cgrid belongs to Grid3 hence do nothing
        pass

    def n_points(self):
        return s3core.dims(grid=self.cgrid)[0]

    def n_edges(self):
        return s3core.dims(grid=self.cgrid)[1]

    def n_faces(self):
        return s3core.dims(grid=self.cgrid)[2]

    def btypes(self):
        return s3core.btypes(grid=self.cgrid)

    def deepcopy(self):
        cp = s3core.extract_grid3_surface(self.cgrid)
        return Surface3(cp)
