'commands for geometric objects'
import copy
import ast
import bgeom
import command


class RemoveGrid(command.Command):
    'remove grid/contour list'
    def __init__(self, names, contnames=[]):
        """ name - string name of the removing grid
            contnames - names of removing user contours
        """
        a = {}
        if len(names) > 0:
            a['names'] = ' '.join(names)
        if len(contnames) > 0:
            a['contnames'] = ' '.join(contnames)
        super(RemoveGrid, self).__init__(a)
        self.names = names
        self.contnames = contnames
        self._backup_grid = []
        self._backup_cont = []

    #overriden from Command
    @classmethod
    def _method_code(cls):
        return "RemoveGrid"

    @classmethod
    def fromstring(cls, slist):
        a = ast.literal_eval(slist)
        gnames, contnames = [], []
        if 'names' in a:
            gnames = a['names'].split()
        if 'contnames' in a:
            contnames = a['contnames'].split()
        return cls(gnames, contnames)

    def doc(self):
        n = self.names + self.contnames
        return "Remove objects: %s" % (', '.join(n))

    def _exec(self):
        self._backup_grid = []
        self._backup_cont = []
        for n in self.names:
            self._backup_grid.append(self.receiver.get_grid(name=n))
            self.receiver.remove_grid(n)
        for n in self.contnames:
            self._backup_cont.append(self.receiver.get_user_contour(name=n))
            self.receiver.remove_user_contour(n)
        return True

    def _clear(self):
        self._backup_grid = []
        self._backup_cont = []

    def _undo(self):
        for v in self._backup_grid:
            self.receiver.grids2.insert(*v)
        for v in self._backup_cont:
            self.receiver.contours2.insert(*v)

    def _redo(self):
        self._exec()


class MoveGrids(command.Command):
    "Move objects"

    def __init__(self, dx, dy, names, cnames):
        a = {'dx': dx, 'dy': dy, 'names': ' '.join(names),
                'contnames': ' '.join(cnames)}
        super(MoveGrids, self).__init__(a)
        self.dx, self.dy = dx, dy
        self.names = names
        self.cnames = cnames

    def doc(self):
        n = self.names + self.cnames
        return "Move objects: " + ", ".join(n)

    @classmethod
    def fromstring(cls, slist):
        a = ast.literal_eval(slist)
        names, contnames = [], []
        if 'names' in a:
            names = a['names'].split()
        if 'contnames' in a:
            contnames = a['contnames'].split()
        return cls(float(a['dx']), float(a['dy']), names, contnames)

    @classmethod
    def _method_code(cls):
        return "MoveGrids"

    def _exec(self):
        for g in self.names:
            self.receiver.grids2[g].move(self.dx, self.dy)
        for c in self.cnames:
            self.receiver.contours2[c].move(self.dx, self.dy)
        return True

    def _clear(self):
        pass

    def _undo(self):
        for g in self.names:
            self.receiver.grids2[g].move(-self.dx, -self.dy)
        for c in self.cnames:
            self.receiver.contours2[c].move(-self.dx, -self.dy)

    def _redo(self):
        self._exec()


class RotateGrids(command.Command):
    "Rotate objects"

    def __init__(self, p0, angle, names, cnames):
        p0 = copy.deepcopy(p0)
        a = {'p0': p0, 'angle': angle,
                'names': ' '.join(names), 'contnames': ' '.join(cnames)}
        super(RotateGrids, self).__init__(a)
        self.x0, self.y0, self.angle = p0.x, p0.y, angle
        self.names = names
        self.cnames = cnames

    def doc(self):
        n = self.names + self.cnames
        return "Rotate objects: " + ", ".join(n)

    @classmethod
    def fromstring(cls, slist):
        a = ast.literal_eval(slist)
        names, contnames = [], []
        if 'names' in a:
            names = a['names'].split()
        if 'contnames' in a:
            contnames = a['contnames'].split()
        return cls(bgeom.Point2.fromstring(a['p0']),
                float(a['angle']), names, contnames)

    @classmethod
    def _method_code(cls):
        return "RotateGrids"

    def _exec(self):
        for g in self.names:
            self.receiver.grids2[g].rotate(self.x0, self.y0, self.angle)
        for c in self.cnames:
            self.receiver.contours2[c].rotate(self.x0, self.y0, self.angle)
        return True

    def _clear(self):
        pass

    def _undo(self):
        for g in self.names:
            self.receiver.grids2[g].rotate(self.x0, self.y0, -self.angle)
        for c in self.cnames:
            self.receiver.contours2[c].rotate(self.x0, self.y0, -self.angle)

    def _redo(self):
        self._exec()


class ScaleGrids(command.Command):
    "Scale objects"

    def __init__(self, p0, xpc, ypc, names, cnames):
        p0 = copy.deepcopy(p0)
        a = {'p0': p0, 'xpc': xpc, 'ypc': ypc,
                'names': ' '.join(names), 'contnames': ' '.join(cnames)}
        super(ScaleGrids, self).__init__(a)
        self.p0 = p0
        self.xpc, self.ypc = xpc, ypc
        self.names = names
        self.cnames = cnames

    def doc(self):
        n = self.names + self.cnames
        return "Scale objects: " + ", ".join(n)

    @classmethod
    def fromstring(cls, slist):
        a = ast.literal_eval(slist)
        names, contnames = [], []
        if 'names' in a:
            names = a['names'].split()
        if 'contnames' in a:
            contnames = a['contnames'].split()
        return cls(bgeom.Point2.fromstring(a['p0']),
                float(a['xpc']), float(a['ypc']), names, contnames)

    @classmethod
    def _method_code(cls):
        return "ScaleGrids"

    def _exec(self):
        for g in self.names:
            self.receiver.grids2[g].scale(self.p0, self.xpc, self.ypc)
        for c in self.cnames:
            self.receiver.contours2[c].scale(self.p0, self.xpc, self.ypc)
        return True

    def _clear(self):
        pass

    def _undo(self):
        for g in self.names:
            self.receiver.grids2[g].unscale(self.p0, self.xpc, self.ypc)
        for c in self.cnames:
            self.receiver.contours2[c].unscale(self.p0, self.xpc, self.ypc)

    def _redo(self):
        self._exec()


class CopyGrid(command.Command):
    "Copy objects"

    def __init__(self, srcnames, copynames, cnames, copycnames):
        """ srcnames, copynames are the list of old and new grid names
            cnames, copycnames are the list of old and new contours names
        """
        a = {'srcnames': ' '.join(srcnames),
                'copynames': ' '.join(copynames),
                'contnames': ' '.join(cnames),
                'copycontnames': ' '.join(copycnames)}
        super(CopyGrid, self).__init__(a)
        self.srcnames = srcnames
        self.copynames = copynames
        self.cnames = cnames
        self.copycnames = copycnames
        #new grids as a tuple (ind, name, grid)
        self.created_grids = []
        self.created_conts = []

    def doc(self):
        n = self.srcnames + self.cnames
        return "Copy grids: " + ", ".join(n)

    @classmethod
    def fromstring(cls, slist):
        a = ast.literal_eval(slist)
        gnames, newgnames = [], []
        cnames, newcnames = [], []
        if 'srcnames' in a:
            gnames = a['srcnames'].split()
            newgnames = a['copynames'].split()
        if 'contnames' in a:
            cnames = a['contnames'].split()
            newcnames = a['copycontnames'].split()
        return cls(gnames, newgnames, cnames, newcnames)

    @classmethod
    def _method_code(cls):
        return "CopyGrid"

    def _exec(self):
        self.created_grids = []
        self.created_conts = []
        for name1, name2 in zip(self.srcnames, self.copynames):
            gold = self.receiver.get_grid(name=name1)[2]
            gnew = gold.deepcopy()
            self.receiver.add_grid(name2, gnew)
            self.created_grids.append(self.receiver.get_grid(grid=gnew))
        for name1, name2 in zip(self.cnames, self.copycnames):
            gold = self.receiver.get_user_contour(name=name1)[2]
            gnew = gold.deepcopy()
            self.receiver.add_user_contour(name2, gnew)
            self.created_conts.append(
                    self.receiver.get_user_contour(cont=gnew))
        return True

    def _clear(self):
        self.created_grids = []
        self.created_conts = []

    def _undo(self):
        for _, n, _ in self.created_grids:
            self.receiver.remove_grid(n)
        for _, n, _ in self.created_conts:
            self.receiver.remove_user_contour(n)

    def _redo(self):
        for v in self.created_grids:
            self.receiver.grids2.insert(*v)
        for v in self.created_conts:
            self.receiver.contours2.insert(*v)
