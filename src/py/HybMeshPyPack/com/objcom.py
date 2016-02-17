'commands for geometric objects'
import command


class RenameGeom(command.Command):
    "Rename"

    def __init__(self, argsdict):
        super(RenameGeom, self).__init__(argsdict)
        self.backup_newname = None
        self.renamed_object = None  # 'grid' or 'ucontour'

    #overriden from Command
    def doc(self):
        return "Rename object: %s" % self.oldname

    @classmethod
    def _arguments_types(cls):
        return {'newname': command.BasicOption(str),
                'oldname': command.BasicOption(str),
                }

    def _exec(self):
        on = self.options['oldname']
        nn = self.options['newname']
        if on in self.receiver.get_grid_names():
            self.renamed_object = 'grid'
            i, _, _ = self.receiver.get_grid(name=on)
            self.receiver.grids2.change_key(on, nn)
            #backup new name as it can be different from self.newname
            _, self.backup_newname, _ = self.receiver.get_grid(ind=i)
        elif on in self.receiver.get_ucontour_names():
            self.renamed_object = 'ucontour'
            i, _, _ = self.receiver.get_ucontour(on)
            self.receiver.contours2.change_key(on, nn)
            #backup new name as it can be different from self.newname
            _, self.backup_newname, _ = self.receiver.get_grid(ind=i)
        else:
            raise command.ExecutionError("Object %s doesn't exist" % on, self)
        return True

    def _clear(self):
        pass

    def _undo(self):
        if self.renamed_object == 'grid':
            self.receiver.grids2.change_key(
                self.backup_newname, self.options['oldname'])
        elif self.renamed_object == 'ucontour':
            self.receiver.contours2.change_key(
                self.backup_newname, self.options['oldname'])

    def _redo(self):
        self._exec()


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
            self._bu_remconts.append(receiver.get_ucontour(name=n))
            receiver.remove_ucontour(n)
        #add grids
        for v in self.addgrids:
            self._bu_addgrids.append(receiver.add_grid(*v))
        #add contours
        for v in self.addconts:
            self._bu_addconts.append(receiver.add_ucontour(*v))
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
            self.receiver.remove_ucontour(v[1])
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

    def _get_added_names(self):
        gnms, cnms = [], []
        if self.__addrem is not None:
            for g in self.__addrem._bu_addgrids:
                gnms.append(g[1])
            for c in self.__addrem._bu_addconts:
                cnms.append(c[1])
        return gnms, cnms

    #function for overriding
    def _addrem_objects(self):
        """ -> addgrids, remgrids, addconts, remconts.
            addgrids -- [ (string name, Grid2 grid), ... ]
            remgrids -- [ string name, ... ]
            addconts -- [ (string name, Contour2 conts), ...]
            remconts -- [ string name, ...]
        """
        raise NotImplementedError


class RemoveGeom(AbstractAddRemove):
    'remove grid/contour list'
    def __init__(self, argsdict):
        super(RemoveGeom, self).__init__(argsdict)

    def doc(self):
        return "Remove objects: %s" % (', '.join(self.options['names']))

    def _addrem_objects(self):
        gridnames, contnames = [], []
        for n in self.options['names']:
            if n in self.receiver.get_grid_names():
                gridnames.append(n)
            elif n in self.receiver.get_ucontour_names():
                contnames.append(n)
            else:
                raise command.ExecutionError('Object %s not found' % n, self)

        return [], gridnames, [], contnames

    @classmethod
    def _arguments_types(cls):
        return {'names': command.ListOfOptions(command.BasicOption(str))}


class MoveGeom(command.Command):
    "Move objects"

    def __init__(self, argsdict):
        super(MoveGeom, self).__init__(argsdict)

    @classmethod
    def _arguments_types(cls):
        return {'names': command.ListOfOptions(command.BasicOption(str)),
                'dx': command.BasicOption(float),
                'dy': command.BasicOption(float),
                }

    def doc(self):
        return "Move objects: " + ", ".join(self.options['names'])

    def _move(self, px, py):
        for n in self.options['names']:
            if n in self.receiver.get_grid_names():
                self.receiver.grids2[n].move(px, py)
            elif n in self.receiver.get_ucontour_names():
                self.receiver.contours2[n].move(px, py)
            else:
                raise command.ExecutionError('Object %s not found' % n, self)

    def _exec(self):
        p = (self.options['dx'], self.options['dy'])
        self._move(*p)
        return True

    def _undo(self):
        p = (-self.options['dx'], -self.options['dy'])
        self._move(*p)

    def _clear(self):
        pass

    def _redo(self):
        self._exec()


class RotateGeom(command.Command):
    "Rotate objects"

    def __init__(self, argsdict):
        super(RotateGeom, self).__init__(argsdict)

    def doc(self):
        return "Rotate objects: " + ", ".join(self.options['names'])

    @classmethod
    def _arguments_types(cls):
        'p0 - center of rotation, angle (deg) - rotation angle'
        return {'names': command.ListOfOptions(command.BasicOption(str)),
                'p0': command.Point2Option(),
                'angle': command.BasicOption(float)
                }

    def _rot(self, p0, a):
        for n in self.options['names']:
            if n in self.receiver.get_grid_names():
                self.receiver.grids2[n].rotate(p0.x, p0.y, a)
            elif n in self.receiver.get_ucontour_names():
                self.receiver.contour2[n].rotate(p0.x, p0.y, a)
            else:
                raise command.ExecutionError('Object %s not found' % n, self)

    def _exec(self):
        p0 = self.options['p0']
        a = self.options['angle']
        self._rot(p0, a)
        return True

    def _undo(self):
        p0 = self.options['p0']
        a = self.options['angle']
        self._rot(p0, a)

    def _clear(self):
        pass

    def _redo(self):
        self._exec()


class ScaleGeom(command.Command):
    "Scale objects"

    def __init__(self, argsdict):
        if 'xpc' not in argsdict:
            argsdict['xpc'] = 100
        if 'ypc' not in argsdict:
            argsdict['ypc'] = 100
        super(ScaleGeom, self).__init__(argsdict)

    def doc(self):
        return "Scale objects: " + ", ".join(self.options['names'])

    @classmethod
    def _arguments_types(cls):
        """ p0 - reference point,
            xpc, ypc - [0-100] percentages of scaling in x/y direction
        """
        return {'names': command.ListOfOptions(command.BasicOption(str)),
                'p0': command.Point2Option(),
                'xpc': command.BasicOption(float),
                'ypc': command.BasicOption(float),
                }

    def _exec(self):
        arg = (self.options['p0'], self.options['xpc'], self.options['ypc'])
        for n in self.options['names']:
            if n in self.receiver.get_grid_names():
                self.receiver.grids2[n].scale(*arg)
            elif n in self.receiver.get_ucontour_names():
                self.receiver.contours2[n].scale(*arg)
            else:
                raise command.ExecutionError('Object %s not found' % n, self)
        return True

    def _clear(self):
        pass

    def _undo(self):
        arg = (self.options['p0'], self.options['xpc'], self.options['ypc'])
        for n in self.options['names']:
            if n in self.receiver.get_grid_names():
                self.receiver.grids2[n].unscale(*arg)
            elif n in self.receiver.get_ucontour_names():
                self.receiver.contours2[n].unscale(*arg)
            else:
                raise command.ExecutionError('Object %s not found' % n, self)

    def _redo(self):
        self._exec()


class CopyGeom(AbstractAddRemove):
    "Copy objects"

    def __init__(self, argsdict):
        if 'dx' not in argsdict:
            argsdict['dx'] = 0
        if 'dy' not in argsdict:
            argsdict['dy'] = 0
        super(CopyGeom, self).__init__(argsdict)

    @classmethod
    def _arguments_types(cls):
        """ names - list of sources
            newnames - list of new names for copies objects
            dx, dy - shift for copied objects
        """
        return {'names': command.ListOfOptions(command.BasicOption(str)),
                'newnames': command.ListOfOptions(command.BasicOption(str)),
                'dx': command.BasicOption(float),
                'dy': command.BasicOption(float),
                }

    def doc(self):
        return "Copy grids: " + ", ".join(self.options['names'])

    def _addrem_objects(self):
        #get separate grid and contours names
        gnames, contnames = [], []
        for n, nn in zip(self.options['names'], self.options['newnames']):
            if n in self.receiver.get_grid_names():
                gnames.append((n, nn))
            elif n in self.receiver.get_ucontour_names():
                contnames.append((n, nn))
            else:
                raise command.ExecutionError('Object %s not found' % n, self)

        newg, newc = [], []
        #copy
        for name1, name2 in gnames:
            gold = self.receiver.get_grid(name=name1)[2]
            gnew = gold.deepcopy()
            newg.append((name2, gnew))
        for name1, name2 in contnames:
            gold = self.receiver.get_ucontour(name=name1)[2]
            gnew = gold.deepcopy()
            newc.append((name2, gnew))

        #shift
        if self.options['dx'] != 0 or self.options['dy'] != 0:
            for g in newg:
                g[1].move(self.options['dx'], self.options['dy'])
            for c in newc:
                c[1].move(self.options['dx'], self.options['dy'])
        return newg, [], newc, []
