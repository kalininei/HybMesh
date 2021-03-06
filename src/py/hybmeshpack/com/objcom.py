'commands for geometric objects'
import command
import comopt
import addremove
from hybmeshpack import gdata


class RenameGeom(command.Command):
    "Rename object"

    def __init__(self, argsdict):
        super(RenameGeom, self).__init__(argsdict)
        self.renamed_object = None

    # overriden from Command
    def doc(self):
        return "Rename object: %s" % self.oldname

    @classmethod
    def _arguments_types(cls):
        return {'newname': comopt.BasicOption(str),
                'oldname': comopt.BasicOption(str),
                }

    def _change_name(self, nn):
        obj = self.renamed_object
        try:
            self.receiver.grids2.change(oldobj=obj, newname=nn)
            return True
        except KeyError:
            pass
        try:
            self.receiver.contours2.change(oldobj=obj, newname=nn)
            return True
        except KeyError:
            pass
        try:
            self.receiver.grids3.change(oldobj=obj, newname=nn)
            return True
        except KeyError:
            pass
        try:
            self.receiver.surfaces3.change(oldobj=obj, newname=nn)
            return True
        except KeyError:
            pass
        return False

    def _exec(self):
        self.renamed_object = self.any_object_by_name(self.option('oldname'))
        return self._change_name(self.get_option('newname'))

    def _clear(self):
        self.renamed_object = None

    def _undo(self):
        return self._change_name(self.get_option('oldname'))

    def _redo(self):
        self._exec()


class RemoveGeom(addremove.AbstractAddRemove):
    'remove grid/contour list'
    def __init__(self, argsdict):
        super(RemoveGeom, self).__init__(argsdict)

    @classmethod
    def _arguments_types(cls):
        return {'names': comopt.ListOfOptions(comopt.BasicOption(str))}

    def doc(self):
        return "Remove objects: %s" % (', '.join(self.options['names']))

    def _addrem_objects(self):
        g2names, c2names, g3names, s3names = [], [], [], []
        for n in self.get_option('names'):
            ob = self.any_object_by_name(n)
            if isinstance(ob, gdata.cont2.Contour2):
                c2names.append(n)
            elif isinstance(ob, gdata.srf3.Surface3):
                s3names.append(n)
            elif isinstance(ob, gdata.grid2.Grid2):
                g2names.append(n)
            elif isinstance(ob, gdata.grid3.Grid3):
                g3names.append(n)
            else:
                raise ValueError("Unknown object %" % str(ob), self)

        return [], g2names, [], c2names, [], g3names, [], s3names


class RemoveAll(command.Command):
    'removes all geometry objects and boundary types'
    def __init__(self, argsdict):
        super(RemoveAll, self).__init__(argsdict)
        self.bu = None

    def doc(self):
        return "Remove all objects"

    @classmethod
    def _arguments_types(cls):
        return {}

    def _exec(self):
        self.bu = self.receiver.backup_copy()
        self.receiver.to_zero_state()
        return True

    def _clear(self):
        self.bu = None

    def _undo(self):
        self.receiver.restore_from_backup(self.bu)
        self.bu = None

    def _redo(self):
        self._exec()


class MoveGeom(command.Command):
    "Move geometrical objects"

    def __init__(self, argsdict):
        super(MoveGeom, self).__init__(argsdict)

    @classmethod
    def _arguments_types(cls):
        return {'names': comopt.ListOfOptions(comopt.BasicOption(str)),
                'dx': comopt.BasicOption(float, 0.),
                'dy': comopt.BasicOption(float, 0.),
                'dz': comopt.BasicOption(float, 0.),
                }

    def doc(self):
        return "Move objects: " + ", ".join(self.options['names'])

    def _move(self, dx, dy, dz):
        for n in self.get_option('names'):
            self.any_object_by_name(n).move(dx, dy, dz)

    def _exec(self):
        dx = self.get_option('dx')
        dy = self.get_option('dy')
        dz = self.get_option('dz')
        self._move(dx, dy, dz)
        return True

    def _undo(self):
        dx = -self.get_option('dx')
        dy = -self.get_option('dy')
        dz = -self.get_option('dz')
        self._move(dx, dy, dz)

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
        return {'names': comopt.ListOfOptions(comopt.BasicOption(str)),
                'p0': comopt.Point2Option(),
                'angle': comopt.BasicOption(float)
                }

    def _rot(self, p, a):
        for n in self.get_option('names'):
            self.receiver.get_object2(n).rotate2(p[0], p[1], a)

    def _exec(self):
        p0 = self.get_option('p0')
        a = self.get_option('angle')
        self._rot(p0, a)
        return True

    def _undo(self):
        p0 = self.get_option('p0')
        a = self.get_option('angle')
        self._rot(p0, -a)

    def _clear(self):
        pass

    def _redo(self):
        self._exec()


class ScaleGeom(command.Command):
    "Scale geom objects"

    def __init__(self, argsdict):
        super(ScaleGeom, self).__init__(argsdict)

    def doc(self):
        return "Scale objects: " + ", ".join(self.options['names'])

    @classmethod
    def _arguments_types(cls):
        return {'names': comopt.ListOfOptions(comopt.BasicOption(str)),
                'p0': comopt.Point3Option((0., 0., 0.)),
                'xpc': comopt.BasicOption(float, 100.),
                'ypc': comopt.BasicOption(float, 100.),
                'zpc': comopt.BasicOption(float, 100.),
                }

    def _scale(self, px, py, pz):
        p0 = self.get_option('p0')
        for n in self.get_option('names'):
            self.any_object_by_name(n).scale(px, py, pz, p0[0], p0[1], p0[2])
        return True

    def _exec(self):
        xpc = self.get_option('xpc')
        ypc = self.get_option('ypc')
        zpc = self.get_option('zpc')
        return self._scale(xpc, ypc, zpc)

    def _clear(self):
        pass

    def _undo(self):
        xpc = 10000. / self.get_option('xpc')
        ypc = 10000. / self.get_option('ypc')
        zpc = 10000. / self.get_option('zpc')
        return self._scale(xpc, ypc, zpc)

    def _redo(self):
        self._exec()


class ReflectGeom(command.Command):
    "Reflect 2d objects"

    def __init__(self, argsdict):
        super(ReflectGeom, self).__init__(argsdict)

    @classmethod
    def _arguments_types(cls):
        return {'names': comopt.ListOfOptions(comopt.BasicOption(str)),
                'p1': comopt.Point2Option(),
                'p2': comopt.Point2Option(),
                }

    def doc(self):
        return "Reflect objects: " + ", ".join(self.options['names'])

    def _exec(self):
        (x0, y0) = self.get_option('p1')
        (x1, y1) = self.get_option('p2')
        for n in self.get_option('names'):
            ob = self.receiver.get_object2(n)
            ob.reflect2(x0, y0, x1, y1)
        return True

    def _undo(self):
        self._exec()

    def _clear(self):
        pass

    def _redo(self):
        self._exec()


class CopyGeom(addremove.AbstractAddRemove):
    "Copy geometrical objects"

    def __init__(self, argsdict):
        super(CopyGeom, self).__init__(argsdict)
        self.__odered = []

    @classmethod
    def _arguments_types(cls):
        """ names - list of sources
        """
        return {'names': comopt.ListOfOptions(comopt.BasicOption(str))}

    def doc(self):
        return "Copy objects: " + ", ".join(self.options['names'])

    def odered_output(self):
        g2, c2 = self.added_grids2(), self.added_contours2()
        g3, s3 = self.added_grids3(), self.added_surfaces3()
        ret = []
        for ob in self.__odered:
            if ob == 'c2':
                ret.append(c2[0])
                del c2[0]
            elif ob == 's3':
                ret.append(s3[0])
                del s3[0]
            elif ob == 'g2':
                ret.append(g2[0])
                del g2[0]
            elif ob == 'g3':
                ret.append(g3[0])
                del g3[0]
            else:
                raise ValueError("Unknown object", self)
        return ret

    def _addrem_objects(self):
        g2, c2, g3, s3 = [], [], [], []
        self.__odered = []
        for n in self.get_option('names'):
            ob = self.any_object_by_name(n)
            cp = (ob.deepcopy(), n + '_copy')
            if isinstance(ob, gdata.cont2.Contour2):
                self.__odered.append('c2')
                c2.append(cp)
            elif isinstance(ob, gdata.srf3.Surface3):
                self.__odered.append('s3')
                s3.append(cp)
            elif isinstance(ob, gdata.grid2.Grid2):
                self.__odered.append('g2')
                g2.append(cp)
            elif isinstance(ob, gdata.grid3.Grid3):
                self.__odered.append('g3')
                g3.append(cp)
            else:
                raise ValueError("Unknown object %" % str(ob), self)

        return g2, [], c2, [], g3, [], s3, []

    def _clear(self):
        self._odered = []
        return super(CopyGeom, self)._clear()


class SetBType(command.Command):
    'sets integer boundary types to object segments'

    def __init__(self, arg):
        super(SetBType, self).__init__(arg)
        self.backup_bnd = []
        self.new_bnd = []

    def doc(self):
        return "Set boundary types"

    @classmethod
    def _arguments_types(cls):
        """ name - string
            bdict - {bndtype: [list of edges indicies] }
            btypes - btype for the whole object
        """
        return {'name': comopt.BasicOption(str),
                'btypes': comopt.BoundaryPickerOption(),
                'whole': comopt.NoneOr(comopt.BasicOption(int), None),
                }

    def _exec(self):
        from hybmeshpack import hmcore
        self.inp = hmcore.proc.BndTypesDifference.from_pydata(
                self.get_option('whole'), self.get_option('btypes'))
        obj = self.any_object_by_name(self.get_option('name'))
        self.out = obj.assign_boundary_type(self.inp)
        return True

    def _clear(self):
        self.inp = None
        self.out = None

    def _undo(self):
        obj = self.any_object_by_name(self.get_option('name'))
        obj.assign_boundary_type(self.out)

    def _redo(self):
        obj = self.any_object_by_name(self.get_option('name'))
        obj.assign_boundary_type(self.inp)
