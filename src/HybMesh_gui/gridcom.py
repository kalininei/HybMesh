import ast
import copy
import bgeom
import grid2
import command
import objcom
from unite_grids import unite_grids, grid_excl_cont, setbc_from_conts


class NewGridCommand(objcom.AbstractAddRemove):
    "Command with a new grid addition"
    def __init__(self, argsdict):
        super(NewGridCommand, self).__init__(argsdict)
        try:
            self.name = argsdict["name"]
        except KeyError:
            self.name = ""
        if self.name == "":
            self.name = "NewGrid1"

    def _addrem_objects(self):
        #backup info
        g = self._build_grid()
        g.build_contour()
        return [(self.name, g)], [], [], []

    #function for overriding
    def _build_grid(self):
        '-> grid2.Grid2'
        raise NotImplementedError


class AddUnfRectGrid(NewGridCommand):
    "Add uniform rectangular grid "
    def __init__(self, argsdict):
        argsdict = copy.deepcopy(argsdict)
        super(AddUnfRectGrid, self).__init__(argsdict)
        self.p0, self.p1 = argsdict['p0'], argsdict['p1']
        self.Nx, self.Ny = argsdict['nx'], argsdict['ny']

    @classmethod
    def fromstring(cls, slist):
        arg = ast.literal_eval(slist)
        arg['p0'] = bgeom.Point2.fromstring(arg['p0'])
        arg['p1'] = bgeom.Point2.fromstring(arg['p1'])
        arg['nx'] = int(arg['nx'])
        arg['ny'] = int(arg['ny'])
        return cls(arg)

    def _build_grid(self):
        return grid2.UnfRectGrid(self.p0, self.p1, self.Nx, self.Ny)


class AddUnfCircGrid(NewGridCommand):
    "Add uniform circular grid "

    def __init__(self, argsdict):
        argsdict = copy.deepcopy(argsdict)
        super(AddUnfCircGrid, self).__init__(argsdict)
        self.p0 = argsdict['p0']
        self.rad = argsdict['rad']
        self.Na, self.Nr = argsdict['na'], argsdict['nr']
        try:
            self.coef = argsdict['coef']
        except KeyError:
            self.coef = 1.0
        try:
            self.is_trian = argsdict['is_trian']
        except KeyError:
            self.is_trian = True

    @classmethod
    def fromstring(cls, slist):
        a = ast.literal_eval(slist)
        a['p0'] = bgeom.Point2.fromstring(a['p0'])
        a['rad'] = float(a['rad'])
        a['na'], a['nr'] = int(a['na']), int(a['nr'])
        if 'coef' in a:
            a['coef'] = float(a['coef'])
        if 'is_trian' in a:
            a['is_trian'] = a['is_trian'] == 'True'
        return cls(a)

    def _build_grid(self):
        return grid2.UnfCircGrid(self.p0, self.rad, self.Na, self.Nr,
                self.coef, self.is_trian)


class AddUnfRingGrid(NewGridCommand):
    "Add uniform circular ring"

    def __init__(self, argsdict):
        argsdict = copy.deepcopy(argsdict)
        super(AddUnfRingGrid, self).__init__(argsdict)
        self.p0 = argsdict['p0']
        self.irad = argsdict['radinner']
        self.orad = argsdict['radouter']
        self.Na, self.Nr = argsdict['na'], argsdict['nr']
        try:
            self.coef = argsdict['coef']
        except KeyError:
            self.coef = 1.0
        if 'name' not in argsdict or argsdict['name'] == '':
            argsdict['name'] = 'CircularGrid1'

    @classmethod
    def fromstring(cls, slist):
        a = ast.literal_eval(slist)
        a['p0'] = bgeom.Point2.fromstring(a['p0'])
        a['radinner'] = float(a['radinner'])
        a['radouter'] = float(a['radouter'])
        a['na'], a['nr'] = int(a['na']), int(a['nr'])
        if 'coef' in a:
            a['coef'] = float(a['coef'])
        return cls(a)

    def _build_grid(self):
        return grid2.UnfRingGrid(self.p0, self.irad, self.orad,
                self.Na, self.Nr, self.coef)


class ExcludeContours(objcom.AbstractAddRemove):
    "Exclude contour area from grids"

    def __init__(self, name, grd, cnts, is_inner, keepg=None):
        """ (string name, string grd, [strings] cnts, bool is_inner,
                bool is_inner, bool keepg)

            grd, cnts - grid and a list of the contours names
            is_inner - exculde inner (if true) or outer (if false) area
            keepg - preserve source grid
        """
        if keepg is None:
            keepg = True
        a = {"name": name, "grdname": grd,
                "cntnames": ' '.join(cnts), "is_inner": is_inner,
                "keepg": keepg}
        super(ExcludeContours, self).__init__(a)

        #assign fields
        self.name = name
        self.grd, self.cnts = grd, cnts
        self.is_inner = is_inner
        self.keepg = keepg

    @classmethod
    def fromstring(cls, slist):
        a = ast.literal_eval(slist)
        c = a['cntnames'].split()
        isin = a['is_inner'] == 'True'
        if 'keepg' in a:
            keepg = a['keepg'] == 'True'
        else:
            keepg = None
        return cls(a['name'], a['grdname'], c, isin, keepg)

    #overriden from NewGridCommand keeping self.remove_com
    def __build_grid(self):
        _, _, base = self.receiver.get_grid(name=self.grd)
        for c in self.cnts:
            c = self.receiver.get_any_contour(c)
            base = grid_excl_cont(base, c, self.is_inner)
            if base is None:
                return None
        return base

    def _addrem_objects(self):
        g = self.__build_grid()
        if g is None:
            return [], [], [], []
        remg = [] if self.keepg else [self.grd]
        return [(self.name, g)], remg, [], []


class RenameGrid(command.Command):
    def __init__(self, oldname, newname):
        a = {"oldname": oldname, "newname": newname}
        super(RenameGrid, self).__init__(a)
        self.oldname, self.newname = oldname, newname
        #actual new grid name can differ from self.newname
        #due to unique names policy. Actual new name is stored
        #in self.backup_newname
        self.backup_newname = None

    #overriden from Command
    def doc(self):
        return "Rename grid: %s" % self.oldname

    @classmethod
    def _method_code(cls):
        return "RenameGrid"

    @classmethod
    def fromstring(cls, slist):
        a = ast.literal_eval(slist)
        return cls(a['oldname'], a['newname'])

    def _exec(self):
        i, _, _ = self.receiver.get_grid(name=self.oldname)
        self.receiver.grids2.change_key(self.oldname, self.newname)
        #backup new name as it can be different from self.newname
        _, self.backup_newname, _ = self.receiver.get_grid(ind=i)
        return True

    def _clear(self):
        pass

    def _undo(self):
        self.receiver.grids2.change_key(self.backup_newname, self.oldname)

    def _redo(self):
        self._exec()


class UniteOpts(object):
    ' Grids unification option: gridname + buffer size + density'

    def __init__(self, name, buf=0):
        self.name = name
        self.buf = buf

    def __str__(self):
        d = {'name': self.name, 'buf': self.buf}
        return str(d)

    @classmethod
    def fromstring(cls, slist):
        a = ast.literal_eval(slist)
        nm = a['name']
        b = a['buf'] if 'buf' in a else 0.0
        return cls(nm, b)


class UniteGrids(objcom.AbstractAddRemove):
    def __init__(self, **kwargs):
        """ args[name] - name of the new grid,
            args[fix_bnd] - whether to fix all boundary points (=False)
            args[keepsrc] - whether to remove source grids (=True)
            args[empty_holes] - keep all empty zone in 2nd grid (=False)
            args[s0] = UniteOpts,   <- list of grids. Keys start with 's'.
            args[s1] = UniteOpts,      No other args keys can start with 's'
            ....
        """
        super(UniteGrids, self).__init__(kwargs)
        self.grid_name = kwargs['name']
        self.source = [None for i in range(len(kwargs))]
        maxnum = 0
        for k, v in kwargs.items():
            if k[0] == 's':
                num = int(k[1:])
                if num > maxnum:
                    maxnum = num
                self.source[num] = v
        self.source = self.source[:maxnum + 1]

        self.fix_bnd = kwargs['fix_bnd'] if 'fix_bnd' in kwargs else False
        self.keepsrc = kwargs['keepsrc'] if 'keepsrc' in kwargs else True
        self.empty_holes = kwargs['empty_holes'] \
                if 'empty_holes' in kwargs else False

    def doc(self):
        ret = "Unite Grids: "
        for s in self.source:
            ret += s.name + " "
        return ret

    @classmethod
    def fromstring(cls, slist):
        a = ast.literal_eval(slist)
        for k, v in a.items():
            if k[0] == 's':
                a[k] = UniteOpts.fromstring(v)
            elif v in ['False', 'True']:
                a[k] = v == 'True'
        return cls(**a)

    def _get_grid(self, ind):
        """ returns (Grid, buffer size) from index """
        g = self.receiver.grids2[self.source[ind].name]
        b = self.source[ind].buf
        return g, b

    def _addrem_objects(self):
        #basic grid
        ret, _ = self._get_grid(0)

        #unification
        for i in range(1, len(self.source)):
            g, b = self._get_grid(i)
            ret = unite_grids(ret, g, b, self.fix_bnd, self.empty_holes)
            if (ret is None):
                return [], [], [], []

        #boundary conditions
        ret.build_contour()
        srccont = [self._get_grid(i)[0].cont for i in range(len(self.source))]
        setbc_from_conts(ret.cont, srccont)

        delgrd = [] if self.keepsrc else [n.name for n in self.source]
        return [(self.grid_name, ret)], delgrd, [], []


class BuildBLayer(objcom.AbstractAddRemove):
    def __init__(self, **kwargs):
        """ args[name] - name of the new grid,
            args[hfull]
            args[mult]
            args[h0]
            args[cont]
        """
        super(UniteGrids, self).__init__(kwargs)
        self.grid_name = kwargs['name']
        self.hfull = kwargs['hfull']
        self.mult = kwargs['mult']
        self.h0 = kwargs['h0']
        self.cont = kwargs['cont']

    @classmethod
    def fromstring(cls, slist):
        a = ast.literal_eval(slist)
        return cls(**a)


    def _addrem_objects(self):
        #basic grid
        ret, _ = self._get_grid(0)

        #unification
        for i in range(1, len(self.source)):
            g, b = self._get_grid(i)
            ret = unite_grids(ret, g, b, self.fix_bnd, self.empty_holes)
            if (ret is None):
                return [], [], [], []

        #boundary conditions
        ret.build_contour()
        srccont = [self._get_grid(i)[0].cont for i in range(len(self.source))]
        setbc_from_conts(ret.cont, srccont)

        delgrd = [] if self.keepsrc else [n.name for n in self.source]
        return [(self.grid_name, ret)], delgrd, [], []
