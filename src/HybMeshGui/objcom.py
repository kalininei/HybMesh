'commands for geometric objects'
import copy
import ast
import bgeom
import command


class _AddRemoveObjects(command.Command):
    """ subcommand to remove objects.
        Could be used only as the internal command for other user commands

        Order of execution:
            remove grids, remove contours, add grids, add contours
    """
    def __init__(self, addgrids, remgrids, addconts, remconts):
        """ addgrids -- [ (string name, Grid2 grid), ... ]
            remgrids -- [ string name, ... ]
            addconts -- [ (string name, Contour2 conts), ...]
            remconts -- [ string name, ...]
        """
        self.addgrids, self.addconts = addgrids, addconts
        self.remgrids, self.remconts = remgrids, remconts
        super(_AddRemoveObjects, self).__init__({})
        self._clear()

    @classmethod
    def fromstring(cls, slist):
        raise NotImplementedError

    def doc(self):
        raise NotImplementedError

    def do(self, receiver):
        self.receiver = receiver
        self._clear()
        if len(self.remgrids) + len(self.remconts) + len(self.addgrids) +\
                len(self.addconts) == 0:
            return False
        #remove grids
        for n in self.remgrids:
            self._bu_remgrids.append(receiver.get_grid(name=n))
            receiver.remove_grid(n)
        #remove contours
        for n in self.remconts:
            self._bu_remconts.append(receiver.get_user_contour(name=n))
            receiver.remove_user_contour(n)
        #add grids
        for v in self.addgrids:
            self._bu_addgrids.append(receiver.add_grid(*v))
        #add contours
        for v in self.addconts:
            self._bu_addconts.append(receiver.add_user_contour(*v))
        return True

    def _exec(self):
        return self.do(self.receiver)

    def _clear(self):
        self._bu_remgrids = []
        self._bu_remconts = []
        self._bu_addgrids = []
        self._bu_addconts = []

    def _undo(self):
        #everything in the reversed order
        #undo add contours
        for v in self._bu_addconts[::-1]:
            self.receiver.remove_user_contour(v[1])
        #undo add grids
        for v in self._bu_addgrids[::-1]:
            self.receiver.remove_grid(v[1])
        #undo remove contours
        for v in self._bu_remconts[::-1]:
            self.receiver.contours2.insert(*v)
        #undo remove grids
        for v in self._bu_remgrids[::-1]:
            self.receiver.grids2.insert(*v)

    def _redo(self):
        return self.do(self.receiver)


class AbstractAddRemove(command.Command):
    """ Abstract base for commands which end up
        with adding and/or removing objects """
    def __init__(self, argsdict):
        super(AbstractAddRemove, self).__init__(argsdict)
        self.__addrem = None

    def _exec(self):
        self.__addrem = _AddRemoveObjects(*self._addrem_objects())
        return self.__addrem.do(self.receiver)

    def _clear(self):
        self.__addrem = None

    def _undo(self):
        self.__addrem._undo()

    def _redo(self):
        self.__addrem._redo()

    #function for overriding
    def _addrem_objects(self):
        """ -> addgrids, remgrids, addconts, remconts.
            addgrids -- [ (string name, Grid2 grid), ... ]
            remgrids -- [ string name, ... ]
            addconts -- [ (string name, Contour2 conts), ...]
            remconts -- [ string name, ...]
        """
        raise NotImplementedError


class RemoveGrid(AbstractAddRemove):
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

    def _addrem_objects(self):
        return [], self.names, [], self.contnames


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


class CopyGrid(AbstractAddRemove):
    "Copy objects"

    def __init__(self, srcnames, copynames, cnames, copycnames,
                dx=0.0, dy=0.0):
        """ srcnames, copynames are the list of old and new grid names
            cnames, copycnames are the list of old and new contours names
            dx, dy - copied objects shift
        """
        a = {'srcnames': ' '.join(srcnames),
                'copynames': ' '.join(copynames),
                'contnames': ' '.join(cnames),
                'copycontnames': ' '.join(copycnames),
                'dx': dx, 'dy': dy}
        super(CopyGrid, self).__init__(a)
        self.srcnames = srcnames
        self.copynames = copynames
        self.cnames = cnames
        self.copycnames = copycnames
        if dx != 0 or dy != 0:
            self.shift = dx, dy
        else:
            self.shift = None

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
        if 'dx' in a:
            dx = float(a['dx'])
        else:
            dx = 0.0
        if 'dy' in a:
            dy = float(a['dy'])
        else:
            dy = 0.0
        return cls(gnames, newgnames, cnames, newcnames, dx, dy)

    def _addrem_objects(self):
        newg, newc = [], []
        #copy
        for name1, name2 in zip(self.srcnames, self.copynames):
            gold = self.receiver.get_grid(name=name1)[2]
            gnew = gold.deepcopy()
            newg.append((name2, gnew))
        for name1, name2 in zip(self.cnames, self.copycnames):
            gold = self.receiver.get_user_contour(name=name1)[2]
            gnew = gold.deepcopy()
            newc.append((name2, gnew))
        #shift
        if self.shift is not None:
            for g in newg:
                g[1].move(*self.shift)
            for c in newc:
                c[1].move(*self.shift)
        return newg, [], newc, []
