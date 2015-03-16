"contour commands"
import ast
import copy
import bp
import bgeom
import command
import contour2


class AddRectCont(command.Command):
    "Add rectangular contour"
    def __init__(self, name, p0, p1, bnds):
        self.name = name
        self.p0, self.p1 = copy.deepcopy(p0), copy.deepcopy(p1)
        self.bnds = copy.deepcopy(bnds)
        a = {"name": name, "p0": p0, "p1": p1,
                "bnds": ' '.join(map(str, bnds))}
        super(AddRectCont, self).__init__(a)
        self.cont = (None, None, None)

    @classmethod
    def fromstring(cls, slist):
        arg = ast.literal_eval(slist)
        p0 = bgeom.Point2.fromstring(arg['p0'])
        p1 = bgeom.Point2.fromstring(arg['p1'])
        bnds = map(int, arg['bnds'].split())
        return cls(arg['name'], p0, p1, bnds)

    @classmethod
    def _method_code(cls):
        return "AddRectCont"

    def _exec(self):
        c = contour2.ClosedContour2()
        p0, p1 = self.p0, bgeom.Point2(self.p1.x, self.p0.y)
        p2, p3 = self.p1, bgeom.Point2(self.p0.x, self.p1.y)
        c.append_points([p0, p1, p2, p3], self.bnds)
        self.cont = self.receiver.add_user_contour(self.name, c)
        return True

    def _clear(self):
        del self.cont

    def _undo(self):
        self.receiver.remove_user_contour(self.cont[1])

    def _redo(self):
        self.receiver.contours2.insert(*self.cont)


class EditBoundaryType(command.Command):
    "Add/edit/remove boundary type"
    def __init__(self, remindex, index, name, color):
        """ (int remindex, int index, str name, (char)*3 color)

            removes boundary with index remindex and
            sets new boundary (index, name, color).
            remindex=None if nothing to delete.
            index, name, color=None if nothing to set.
        """
        self.remindex = remindex
        self.index = index
        self.name = name
        self.color = copy.deepcopy(color)
        a = {}
        if self.name is not None:
            a['name'] = self.name
        if self.color is not None:
            a['r'] = self.color[0]
            a['g'] = self.color[1]
            a['b'] = self.color[2]
        if self.remindex is not None:
            a['remindex'] = self.remindex
        if self.index is not None:
            a['index'] = self.index
        super(EditBoundaryType, self).__init__(a)
        self._backup1 = None
        self._backup2 = None

    @classmethod
    def fromstring(cls, slist):
        a = ast.literal_eval(slist)
        remindex, index, name, color = None, None, None, None
        if 'name' in a:
            name = a['name']
        if 'r' in a:
            color = tuple(map(int, [a['r'], a['g'], a['b']]))
        if 'remindex' in a:
            remindex = int(a['remindex'])
        if 'index' in a:
            index = int(a['index'])
        return cls(remindex, index, name, color)

    @classmethod
    def _method_code(cls):
        return cls.__name__

    def doc(self):
        return "Edit boundary type"

    def _exec(self):
        #1 remove
        if self.remindex is not None:
            self._backup1 = copy.deepcopy(
                    self.receiver.boundary_types.get(index=self.remindex))
            self.receiver.boundary_types.rem_bnd(self.remindex)
        #2 set
        if self.index is not None:
            self._backup2 = copy.deepcopy(
                    self.receiver.boundary_types.get(index=self.index))
            self.receiver.boundary_types.set_bnd(self.index, self.name,
                    self.color)
        return True

    def _clear(self):
        self._backup1 = None
        self._backup2 = None

    def _undo(self):
        if self.index is not None:
            self.receiver.boundary_types.rem_bnd(self.index)
        if self._backup2 is not None:
            self.receiver.boundary_types.set_bnd(self._backup2.index,
                    self._backup2.name, self._backup2.color)
        if self._backup1 is not None:
            self.receiver.boundary_types.set_bnd(self._backup1.index,
                    self._backup1.name, self._backup1.color)
        self._clear()

    def _redo(self):
        self._exec()


class BTypePicker(object):
    def __init__(self, contname, bdir):
        """ bdir is the dictionary:
                {boundary index: [list of edges indicies]}
            cont - contour name
        """
        self.name = contname
        self.bdir = bdir

    def __str__(self):
        d = {}
        for bi, ei in self.bdir.items():
            ei.sort()
            d[bi] = bp.compress_int_list(ei)
        return str(d)

    @classmethod
    def fromstring(cls, cname, slist):
        arg = ast.literal_eval(slist)
        for k, v in arg.items():
            arg[k] = bp.int_list_from_compress(v)
        return cls(cname, arg)


class SetBTypeToContour(command.Command):
    'sets integer boundary types to contours edges'

    def __init__(self, conts_opts):
        """ (BTypePicker conts_opts)
        """
        self.conts_opts = conts_opts
        a = {}
        for c in conts_opts:
            a[c.name] = c
        super(SetBTypeToContour, self).__init__(a)
        self.backup_bnd = []
        self.new_bnd = []

    @classmethod
    def fromstring(cls, slist):
        a = ast.literal_eval(slist)
        conts_opts = []
        for k, v in a.items():
            conts_opts.append(BTypePicker.fromstring(k, v))
        return cls(conts_opts)

    def doc(self):
        return "Set boundary types to contours: " + ', '.join(
                [x.name for x in self.conts_opts])

    def _exec(self):
        for co in self.conts_opts:
            bu = {}
            new = {}
            cont = self.receiver.get_any_contour(co.name)
            #fill new and backup dicts
            for b, ei in co.bdir.items():
                for e in ei:
                    #backup
                    bu[e] = cont.edge_bnd(e)
                    #command
                    new[e] = b
            #evalution
            self.new_bnd.append(new)
            self.backup_bnd.append(bu)
        self._redo()
        return True

    def _clear(self):
        self.backup_bnd = []
        self.new_bnd = []

    def _undo(self):
        for co, bu in zip(self.conts_opts, self.backup_bnd):
            cont = self.receiver.get_any_contour(co.name)
            cont.add_edge_bnd(bu)

    def _redo(self):
        for co, new in zip(self.conts_opts, self.new_bnd):
            cont = self.receiver.get_any_contour(co.name)
            cont.add_edge_bnd(new)
