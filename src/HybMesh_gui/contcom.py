"contour commands"
import ast
import copy
import bp
import bgeom
import command
import objcom
import contour2


class AddRectCont(objcom.AbstractAddRemove):
    "Add rectangular contour"
    def __init__(self, name, p0, p1, bnds):
        self.name = name
        self.p0, self.p1 = copy.deepcopy(p0), copy.deepcopy(p1)
        self.bnds = copy.deepcopy(bnds)
        a = {"name": name, "p0": p0, "p1": p1,
                "bnds": ' '.join(map(str, bnds))}
        super(AddRectCont, self).__init__(a)

    def doc(self):
        return "Add rectangular contour"

    @classmethod
    def fromstring(cls, slist):
        arg = ast.literal_eval(slist)
        p0 = bgeom.Point2.fromstring(arg['p0'])
        p1 = bgeom.Point2.fromstring(arg['p1'])
        bnds = map(int, arg['bnds'].split())
        return cls(arg['name'], p0, p1, bnds)

    def _addrem_objects(self):
        c = contour2.ClosedContour2()
        p0, p1 = self.p0, bgeom.Point2(self.p1.x, self.p0.y)
        p2, p3 = self.p1, bgeom.Point2(self.p0.x, self.p1.y)
        c.append_points([p0, p1, p2, p3], self.bnds)
        return [], [], [(self.name, c)], []


class GridBndToContour(command.Command):
    'Copy grid boundary to user contour'
    def __init__(self, name, grd_name, simplify, separate):
        self.name = name
        self.grd_name = grd_name
        self.simplify = simplify
        self.separate = separate
        a = {"name": name, "grd_name": grd_name,
                "simplify": simplify, "separate": separate}
        super(GridBndToContour, self).__init__(a)
        self.cont = (None, None, None)
        self.simp_command = None

    def doc(self):
        return "Copy boundary from %s" % self.grd_name

    @classmethod
    def fromstring(cls, slist):
        arg = ast.literal_eval(slist)
        return cls(arg['name'], arg['grd_name'], arg['simplify'] == 'True',
                arg['separate'] == 'True')

    def _exec(self):
        #copy contour
        src = self.receiver.grids2[self.grd_name].cont
        cont = contour2.Contour2.create_from_abstract(src)
        self.cont = self.receiver.add_user_contour(self.name, cont)
        #simplify and separate
        if self.simplify or self.separate:
            self.simp_command = SimplifyContours(
                    self.separate, self.simplify, [self.cont[1]],
                    [self.cont[1]], False, 0.0)
            if not self.simp_command.do(self.receiver):
                self.simp_command = None
        return True

    def _clear(self):
        self.cont = None
        self.simp_command = None

    def _undo(self):
        if self.simp_command:
            self.simp_command._undo()
        self.receiver.remove_user_contour(self.cont[1])

    def _redo(self):
        self.receiver.contours2.insert(*self.cont)
        if self.simp_command:
            self.simp_command._redo()


class SimplifyContours(objcom.AbstractAddRemove):
    'Separate and simplify user contour'
    def __init__(self, sep, simp, conts, newnms, keepsrc, angle):
        '(separate, simplify, [cont_names], [new_cont_names], keepsrc, angle)'
        self.sep, self.simp = sep, simp
        self.conts = conts
        self.newnms, self.keepsrc = newnms, keepsrc
        self.angle = angle
        a = {"separate": sep, "simplify": simp, "names": ' '.join(conts)}
        if simp:
            a['angle'] = angle
        if sep:
            a['new_names'] = ' '.join(newnms)
            a['keepsrc'] = keepsrc
        super(SimplifyContours, self).__init__(a)
        if not self.sep:
            self.newnms = self.conts

    def doc(self):
        return "Simplify contours: %s" % ' '.join(self.conts)

    @classmethod
    def fromstring(cls, slist):
        arg = ast.literal_eval(slist)
        sep, simp = arg['separate'] == 'True', arg['simplify'] == 'True'
        conts = arg['names'].split()
        if sep:
            newnms = arg['new_names'].split()
            keepsrc = arg['keepsrc'] == 'True'
        else:
            newnms = None
            keepsrc = None
        if simp:
            angle = float(arg['angle'])
        else:
            angle = 0
        return cls(sep, simp, conts, newnms, keepsrc, angle)

    def _addrem_objects(self):
        added, removed = [], []
        #operate
        for i, nm in enumerate(self.conts):
            cont = self.receiver.get_user_contour(name=nm)
            res = []

            #simplify
            if self.simp:
                s = cont[2].simplify(self.angle)
                if s:
                    res = [s]

            #separate
            if self.sep:
                if res:
                    s = res[0].separate()
                    if s:
                        res = s
                else:
                    res = cont[2].separate()

            #continue if no result
            if not res:
                continue

            #delete
            # simplify + separate (keep src)
            if len(res) == 1 or not self.keepsrc:
                removed.append(cont[1])
            #add
            for r in res:
                added.append((self.newnms[i], r))

        return [], [], added, removed


class UniteContours(objcom.AbstractAddRemove):
    "unite contours to one contours with complicated connection"
    def __init__(self, name, srcnames, keepsrc):
        """ (string name, [string] srcnames, bool keepsrc)

            name - united contour name
            srcnames - names of contours to unite
            keepsrc - whether to keep or remove source contours
        """
        self.name = name
        self.srcnames = srcnames
        self.keepsrc = keepsrc
        a = {"name": name, "srcnames": ' '.join(srcnames), "keepsrc": keepsrc}
        super(UniteContours, self).__init__(a)

    def doc(self):
        return "Unite contours: %s" % ' '.join(self.srcnames)

    @classmethod
    def fromstring(cls, slist):
        arg = ast.literal_eval(slist)
        name = arg['name']
        conts = arg['srcnames'].split()
        keep = arg['keepsrc'] == 'True'
        return cls(name, conts, keep)

    def _addrem_objects(self):
        conts = [self.receiver.contours2[n] for n in self.srcnames]
        #operate
        newcont = contour2.Contour2.create_from_abstract(conts[0])
        for c in conts[1:]:
            newcont.add_from_abstract(c)
        ac = [(self.name, newcont)]
        rc = [] if self.keepsrc else self.srcnames
        return [], [], ac, rc


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


class RenameContour(command.Command):
    def __init__(self, oldname, newname):
        a = {"oldname": oldname, "newname": newname}
        super(RenameContour, self).__init__(a)
        self.oldname, self.newname = oldname, newname
        #actual new grid name can differ from self.newname
        #due to unique names policy. Actual new name is stored
        #in self.backup_newname
        self.backup_newname = None

    #overriden from Command
    def doc(self):
        return "Rename contour: %s" % self.oldname

    @classmethod
    def fromstring(cls, slist):
        a = ast.literal_eval(slist)
        return cls(a['oldname'], a['newname'])

    def _exec(self):
        i, _, _ = self.receiver.get_user_contour(name=self.oldname)
        self.receiver.contours2.change_key(self.oldname, self.newname)
        #backup new name as it can be different from self.newname
        _, self.backup_newname, _ = self.receiver.get_user_contour(ind=i)
        return True

    def _clear(self):
        pass

    def _undo(self):
        self.receiver.contours2.change_key(self.backup_newname, self.oldname)

    def _redo(self):
        self._exec()


class CreateContour(objcom.AbstractAddRemove):
    'Manually enter contour points'
    def __init__(self, name, pts, bnds, closed):
        'str name, [Points], [int bnd], bool closed'
        self.name = name
        self.pts, self.bnds, self.is_closed = pts, bnds, closed
        a = {'name': name, 'points': ', '.join(map(str, pts)),
                'bnd': bnds, 'closed': closed}
        super(CreateContour, self).__init__(a)

    def doc(self):
        return "Create contour"

    @classmethod
    def fromstring(cls, slist):
        arg = ast.literal_eval(slist)
        pts = map(bgeom.Point2.fromstring, arg['points'].split(','))
        b = arg['bnd']
        cl = arg['closed'] == 'True'
        return cls(arg['name'], pts, b, cl)

    def _addrem_objects(self):
        eds = [[i, i + 1] for i in range(len(self.pts))]
        if self.is_closed:
            eds[-1][1] = 0
        else:
            eds = eds[:-1]
        cont = contour2.Contour2.create_from_point_set(self.pts,
                eds, self.bnds)
        return [], [], [cont], []
