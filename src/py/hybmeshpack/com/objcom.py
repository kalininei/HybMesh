'commands for geometric objects'
import command
import copy


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
            remove grids, remove contours, ..., add grids, add contours, ...
    """
    def __init__(self, addg2=[], remg2=[],
                 addc2=[], remc2=[],
                 addg3=[], remg3=[],
                 adds3=[], rems3=[]):
        """ addg2 -- [ (string name, Grid2 grid), ... ]
            remg2 -- [ string name, ... ]
            addc2 -- [ (string name, Contour2 conts), ...]
            remc2 -- [ string name, ...]
            addg3    -- [ (string name, Grid3 grid), ... ]
            remg3    -- [ string name, ...]
            adds3    -- [ (string name, Surface srf), ... ]
            rems3    -- [ string name, ...]
        """
        self.addg2, self.addc2 = addg2, addc2
        self.remg2, self.remc2 = remg2, remc2
        self.addg3, self.remg3 = addg3, remg3
        self.adds3, self.rems3 = adds3, rems3
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
        if len(self.addg2) + len(self.remg2) +\
                len(self.addc2) + len(self.remc2) +\
                len(self.addg3) + len(self.remg3) +\
                len(self.adds3) + len(self.rems3) == 0:
            return False

        #remove grids
        for n in self.remg2:
            self._bu_remg2.append(receiver.get_grid(name=n))
            receiver.remove_grid(n)
        #remove contours
        for n in self.remc2:
            self._bu_remc2.append(receiver.get_ucontour(name=n))
            receiver.remove_ucontour(n)
        #remove g3
        for n in self.remg3:
            self._bu_remg3.append(receiver.get_grid3(name=n))
            receiver.remove_grid3(n)
        #remove s3
        for n in self.rems3:
            self._bu_rems3.append(receiver.get_usurface(name=n))
            receiver.remove_usurface(n)
        #add grids
        for v in self.addg2:
            self._bu_addg2.append(receiver.add_grid(*v))
        #add contours
        for v in self.addc2:
            self._bu_addc2.append(receiver.add_ucontour(*v))
        #add grid3
        for v in self.addg3:
            self._bu_addg3.append(receiver.add_grid3(*v))
        #add surfacee3
        for v in self.adds3:
            self._bu_adds3.append(receiver.add_usurface(*v))
        return True

    def _exec(self):
        return self.do(self.receiver)

    def _clear(self):
        self._bu_addg2 = []
        self._bu_remg2 = []
        self._bu_addc2 = []
        self._bu_remc2 = []
        self._bu_addg3 = []
        self._bu_remg3 = []
        self._bu_adds3 = []
        self._bu_rems3 = []

    def _undo(self):
        #everything in the reversed order
        #undo add s3
        for v in self._bu_adds3[::-1]:
            self.receiver.remove_usurface(v[1])
        #undo add g3
        for v in self._bu_addg3[::-1]:
            self.receiver.remove_grid3(v[1])
        #undo add contours
        for v in self._bu_addc2[::-1]:
            self.receiver.remove_ucontour(v[1])
        #undo add grids
        for v in self._bu_add2[::-1]:
            self.receiver.remove_grid(v[1])
        #undo remove surface
        for v in self._bu_rems3[::-1]:
            self.receiver.surfaces3.insert(*v)
        #undo remove grid3
        for v in self._bu_remg3[::-1]:
            self.receiver.grids3.insert(*v)
        #undo remove contours
        for v in self._bu_remc2[::-1]:
            self.receiver.contours2.insert(*v)
        #undo remove grids
        for v in self._bu_remg2[::-1]:
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
        "[[g2 names], [c2 names], [g3 names], [s3 names]]"
        gnms, cnms, g3nms, s3nms = [], [], [], []
        if self.__addrem is not None:
            for g in self.__addrem._bu_addg2:
                gnms.append(g[1])
            for c in self.__addrem._bu_addc2:
                cnms.append(c[1])
            for g in self.__addrem._bu_addg3:
                g3nms.append(g[1])
            for s in self.__addrem._bu_adds3:
                s3nms.append(s[1])
        return gnms, cnms, g3nms, s3nms

    #function for overriding
    def _addrem_objects(self):
        """ -> addg2, remg2, addc2, remc2, addg3, remg3, adds3, rems3
            addg2 -- [ (string name, Grid2 grid), ... ]
            remg2 -- [ string name, ... ]
            addc2 -- [ (string name, Contour2 conts), ...]
            remc2 -- [ string name, ...]
            addg3 -- [ (string name, Grid3 grid), ... ]
            remg3 -- [ string name, ...]
            adds3 --
            rems3 --
        """
        raise NotImplementedError


class RemoveGeom(AbstractAddRemove):
    'remove grid/contour list'
    def __init__(self, argsdict):
        super(RemoveGeom, self).__init__(argsdict)

    def doc(self):
        return "Remove objects: %s" % (', '.join(self.options['names']))

    def _addrem_objects(self):
        gridnames, contnames, g3names, s3names = [], [], [], []
        for n in self.options['names']:
            if n in self.receiver.get_grid_names():
                gridnames.append(n)
            elif n in self.receiver.get_ucontour_names():
                contnames.append(n)
            elif n in self.receiver.get_grid3_names():
                g3names.append(n)
            elif n in self.receiver.get_usurface_names():
                s3names.append(n)
            else:
                raise command.ExecutionError('Object %s not found' % n, self)

        return [], gridnames, [], contnames, [], g3names, [], s3names

    @classmethod
    def _arguments_types(cls):
        return {'names': command.ListOfOptions(command.BasicOption(str))}


class RemoveAll(AbstractAddRemove):
    'remove grids, contours, grids3d, boundary types list'
    def __init__(self, argsdict):
        super(RemoveAll, self).__init__(argsdict)
        self._btypes_backup = None

    def doc(self):
        return "Remove all objects"

    @classmethod
    def _arguments_types(cls):
        return {}

    def _addrem_objects(self):
        gridnames = self.receiver.get_grid_names()
        contnames = self.receiver.get_ucontour_names()
        g3names = self.receiver.get_grid3_names()
        s3names = self.receiver.get_usurface_names()
        return [], gridnames, [], contnames, [], g3names, [], s3names

    def __remove_btypes(self):
        self._btypes_backup = []
        srbt = self.receiver.boundary_types
        for bt in srbt._data:
            self._btypes_backup.append(copy.deepcopy(bt))
        srbt.clear()

    def _exec(self):
        self.__remove_btypes()
        return super(RemoveAll, self)._exec()

    def _clear(self):
        self._btypes_backup = None
        super(RemoveAll, self)._clear()

    def _undo(self):
        if self._btypes_backup is not None:
            for bt in self._btypes_backup:
                self.receiver.boundary_types.set_bnd(
                    bt.index, bt.name, bt.color)
        self._btypes_backup = None
        super(RemoveAll, self)._undo()

    def _redo(self):
        self.__remove_btypes()
        super(RemoveAll, self)._redo()


class MoveGeom(command.Command):
    "Move 2d objects"

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
    "Rotate 2d objects"

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
                self.receiver.contours2[n].rotate(p0.x, p0.y, a)
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
    "Scale 2d objects"

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


class ReflectGeom(command.Command):
    "Reflect 2d objects"

    def __init__(self, argsdict):
        if (argsdict['p1'].x == argsdict['p2'].x and
                argsdict['p1'].y == argsdict['p2'].y):
            raise ValueError("Reflection over a zero vector is imposible")

        super(ReflectGeom, self).__init__(argsdict)

    @classmethod
    def _arguments_types(cls):
        return {'names': command.ListOfOptions(command.BasicOption(str)),
                'p1': command.Point2Option(),
                'p2': command.Point2Option(),
                }

    def doc(self):
        return "Reflect objects: " + ", ".join(self.options['names'])

    def _exec(self):
        for n in self.options['names']:
            if n in self.receiver.get_grid_names():
                ob = self.receiver.grids2[n]
            elif n in self.receiver.get_ucontour_names():
                ob = self.receiver.contours2[n]
            else:
                raise command.ExecutionError('Object %s not found' % n, self)
            ob.reflect(self.options['p1'], self.options['p2'])
        return True

    def _undo(self):
        self._exec()

    def _clear(self):
        pass

    def _redo(self):
        self._exec()


class CopyGeom(AbstractAddRemove):
    "Copy 2d objects"

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
        return newg, [], newc, [], [], [], [], []
