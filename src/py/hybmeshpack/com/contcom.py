"contour commands"
import command
import comopt
import addremove
from hybmeshpack.gdata.cont2 import Contour2
from hybmeshpack.hmcore import c2 as c2core
import math
import copy


class AddRectCont(addremove.AbstractAddRemove):
    "Add rectangular contour"
    def __init__(self, argsdict):
        super(AddRectCont, self).__init__(argsdict)

    @classmethod
    def _arguments_types(cls):
        return {'name': comopt.BasicOption(str, None),
                'p0': comopt.Point2Option(),
                'p1': comopt.Point2Option(),
                'bnds': comopt.ListOfOptions(comopt.BasicOption(int),
                                             [0, 0, 0, 0]),
                }

    def _addrem_contour2(self):
        p0, p1 = self.get_option('p0'), self.get_option('p1')
        pts = [[p0[0], p0[1]], [p1[0], p0[1]],
               [p1[0], p1[1]], [p0[0], p1[1]]]
        bnds = self.get_option('bnds')
        c = c2core.build_from_points(pts, True, bnds)
        return [(Contour2(c), self.get_option('name'))], []


class AddCircCont(addremove.AbstractAddRemove):
    "Add circular contour"
    def __init__(self, argsdict):
        super(AddCircCont, self).__init__(argsdict)

    @classmethod
    def _arguments_types(cls):
        return {'name': comopt.BasicOption(str, None),
                'p0': comopt.Point2Option(),
                'rad': comopt.BasicOption(float),
                'na': comopt.BasicOption(int),
                'bnd': comopt.BasicOption(int, 0),
                }

    def _addrem_contour2(self):
        op0 = self.get_option('p0')
        b = self.get_option('bnd')
        rad = self.get_option('rad')
        na = self.get_option('na')
        astep = 2 * math.pi / na
        pts = []
        for i in range(na):
            angle = i * astep
            pts.append([op0[0] + rad * math.cos(angle),
                        op0[1] + rad * math.sin(angle)])
        c = c2core.build_from_points(pts, True, [b])
        return [(Contour2(c), self.get_option('name'))], []


class GridBndToContour(addremove.AbstractAddRemove):
    'Copy grid boundary to user contour'

    def __init__(self, arg):
        super(GridBndToContour, self).__init__(arg)

    @classmethod
    def _arguments_types(cls):
        return {'grid_name': comopt.BasicOption(str),
                'cont_name': comopt.BasicOption(str, None),
                'simplify': comopt.BoolOption(True),
                'separate': comopt.BoolOption(False),
                }

    def _addrem_contour2(self):
        gn = self.get_option('grid_name')
        cn = self.get_option('cont_name')
        if cn is None:
            cn = gn + '_contour'
        simp = self.get_option('simplify')
        sep = self.get_option('separate')
        cc = self.grid2_by_name(gn).contour().contour2()
        if simp:
            c2core.simplify(cc.cdata, 0., True)
        if sep:
            cc = c2core.quick_separate(cc.cdata)
            for i, c in enumerate(cc):
                cc[i] = (Contour2(c), cn)
        else:
            cc = [(cc, cn)]
        return cc, []


class SimplifyContours(addremove.AbstractAddRemove):
    'Separate/simplify user contours'
    def __init__(self, arg):
        super(SimplifyContours, self).__init__(arg)

    @classmethod
    def _arguments_types(cls):
        """ [names] - list of sources contours
            [new_names] - list of names for created contours
            simplify - do simplification: make all nodes non-collinear
                       Edges will not be splitted if they have different
                       boundary types.
            angle - minimum acute angle between edges lines
                    after simplification (deg)
            separate - separate contour to ones with not-connected geometry
        """
        return {'names': comopt.ListOfOptions(comopt.BasicOption(str)),
                'new_names': comopt.ListOfOptions(comopt.BasicOption(str), []),
                'separate': comopt.BoolOption(False),
                'simplify': comopt.BoolOption(True),
                'angle': comopt.BasicOption(float, 1.0),
                }

    def _addrem_contour2(self):
        simp, sep, an, names, nnames = \
            [self.get_option(x) for x in ['simplify', 'separate',
                                          'angle', 'names', 'new_names']]
        if len(nnames) < names:
            for i, n in enumerate(names):
                if i >= len(nnames):
                    nnames.append(n)
        added = []
        # make copies
        for i, nm in enumerate(names):
            cont = self.any_acont_by_name(nm).deepcopy()
            added.append((cont, nnames[i]))

        # simplify
        if simp:
            for c in added:
                c2core.simplify(c[0].cdata, an, True)

        # separate
        if sep:
            newadded = []
            for c in added:
                for cc in c2core.quick_separate(c[0].cdata):
                    newadded.append((Contour2(cc), c[1]))
            added = newadded

        return added, []


class UniteContours(addremove.AbstractAddRemove):
    "Unite contours"
    def __init__(self, arg):
        super(UniteContours, self).__init__(arg)

    def doc(self):
        return "Unite contours: %s" % ' '.join(self.srcnames)

    @classmethod
    def _arguments_types(cls):
        return {'name': comopt.BasicOption(str, None),
                'sources': comopt.ListOfOptions(comopt.BasicOption(str)),
                }

    def _addrem_contour2(self):
        nm = self.get_option('name')
        src = self.get_option('sources')

        scont = []
        for s in src:
            scont.append(self.any_cont_by_name(s))
        pscont = [c.cdata for c in scont]
        newcont = c2core.unite_contours(pscont)

        return [(Contour2(newcont), nm)], []


class DecomposeContour(addremove.AbstractAddRemove):
    'Decompose contour'
    def __init__(self, arg):
        super(DecomposeContour, self).__init__(arg)

    @classmethod
    def _arguments_types(cls):
        return {'name': comopt.BasicOption(str, None),
                'source': comopt.BasicOption(str),
                }

    def _addrem_contour2(self):
        s = self.any_cont_by_name(self.get_option('source'))
        nm = self.get_option('name')
        nc = c2core.decompose_contour(s.cdata)
        ret = []
        for n in nc:
            ret.append((Contour2(n), nm))
        return ret, []


class EditBoundaryType(command.Command):
    "Add/edit/remove boundary type"
    def __init__(self, arg):
        super(EditBoundaryType, self).__init__(arg)
        self.bu = None

    @classmethod
    def _arguments_types(cls):
        """ (int remindex, int index, str name)

            removes boundary with index remindex and
            sets new boundary (index, name).
            remindex=None if nothing to delete.
            index=None, name=None if nothing to set.
            if index exists it will be overriden
            if name exists it will be changed
        """
        return {'remindex': comopt.BasicOption(int, None),
                'index': comopt.BasicOption(int, None),
                'name': comopt.BasicOption(str, None),
                }

    def doc(self):
        return "Edit boundary type"

    def _exec(self):
        self.bu = copy.deepcopy(self.receiver.zone_types)
        remi = self.get_option('remindex')
        addi = self.get_option('index')
        addn = self.get_option('name')
        if remi is not None:
            self.receiver.remove_zone_type(remi)
        if addi is not None or addn is not None:
            self.receiver.set_zone_type(addi, addn)
        return True

    def _clear(self):
        self.bu = None

    def _undo(self):
        self.receiver.zone_types = self.bu
        self.receiver._znames_pool = self.bu._names_pool
        self._clear()

    def _redo(self):
        self._exec()


class SetBTypeToContour(command.Command):
    'sets integer boundary types to contours edges'

    def __init__(self, arg):
        super(SetBTypeToContour, self).__init__(arg)
        self.backup_bnd = []
        self.new_bnd = []

    def doc(self):
        return "Set boundary types to contours"

    @classmethod
    def _arguments_types(cls):
        """ name - string
            btypes - {bndtype: [list of edges indicies] }
        """
        return {'name': comopt.BasicOption(str),
                'btypes': comopt.BoundaryPickerOption(),
                }

    def _exec(self):
        self._clear()

        cont = self.any_acont_by_name(self.get_option('name'))
        oldbnd = list(cont.raw_data('btypes'))
        newbnd = copy.deepcopy(oldbnd)
        di = self.get_option('btypes')
        for k, v in di.iteritems():
            for it in v:
                newbnd[it] = k

        self.backup_bnd = oldbnd
        self.new_bnd = newbnd

        self._redo()
        return True

    def _clear(self):
        self.backup_bnd = []
        self.new_bnd = []

    def _undo(self):
        cont = self.any_acont_by_name(self.get_option('name'))
        cont.set_bnd(self.backup_bnd)

    def _redo(self):
        cont = self.any_acont_by_name(self.get_option('name'))
        cont.set_bnd(self.new_bnd)


class CreateContour(addremove.AbstractAddRemove):
    'Manually enter contour points'
    def __init__(self, argsdict):
        super(CreateContour, self).__init__(argsdict)

    def doc(self):
        return "Create contour from sequence of points"

    @classmethod
    def _arguments_types(cls):
        return {'name': comopt.BasicOption(str, None),
                'points': comopt.ListOfOptions(comopt.Point2Option()),
                'bnds': comopt.ListOfOptions(comopt.BasicOption(int), []),
                }

    def _addrem_contour2(self):
        c = c2core.build_from_points(
            self.get_option('points'),
            False,
            self.get_option('bnds'))
        return [(Contour2(c), self.get_option('name'))], []


class CreateSpline(addremove.AbstractAddRemove):
    " Creates spline"

    def __init__(self, argsdict):
        super(CreateSpline, self).__init__(argsdict)

    def doc(self):
        return "Create spline from sequence of points"

    @classmethod
    def _arguments_types(cls):
        return {'name': comopt.BasicOption(str, None),
                'points': comopt.ListOfOptions(comopt.Point2Option()),
                'bnds': comopt.ListOfOptions(comopt.BasicOption(int), []),
                'nedges': comopt.BasicOption(int, 100),
                }

    def _addrem_contour2(self):
        c = c2core.spline(
            self.get_option('points'),
            self.get_option('bnds'),
            self.get_option('nedges'))
        return [(Contour2(c), self.get_option('name'))], []


class ClipDomain(addremove.AbstractAddRemove):
    'Domain clipping operation'

    def __init__(self, argsdict):
        super(ClipDomain, self).__init__(argsdict)

    def doc(self):
        return "Domain clipping operation"

    @classmethod
    def _arguments_types(cls):
        return {'name': comopt.BasicOption(str, None),
                'c1': comopt.BasicOption(str),
                'c2': comopt.BasicOption(str),
                'oper': comopt.BasicOption(str),
                'simplify': comopt.BoolOption(True),
                }

    def _addrem_contour2(self):
        c1 = self.any_cont_by_name(self.get_option('c1'))
        c2 = self.any_cont_by_name(self.get_option('c2'))
        op = self.get_option('oper')
        c3 = c2core.clip_domain(
            c1.cdata, c2.cdata, op,
            self.get_option('simplify'))
        return [(Contour2(c3), self.get_option('name'))], []


class PartitionContour(addremove.AbstractAddRemove):
    'Partition contour operation'

    def __init__(self, argsdict):
        super(PartitionContour, self).__init__(argsdict)

    def doc(self):
        return "Contour partition"

    @classmethod
    def _arguments_types(cls):
        return {'name': comopt.BasicOption(str, None),
                'algo': comopt.BasicOption(str),
                'step': comopt.ListOfOptions(comopt.BasicOption(float)),
                'angle0': comopt.BasicOption(float, 1.),
                'base': comopt.BasicOption(str),
                'keepbnd': comopt.BoolOption(True),
                'nedges': comopt.NoneOr(comopt.BasicOption(int), None),
                'crosses': comopt.ListOfOptions(comopt.BasicOption(str), []),
                'start': comopt.NoneOr(comopt.Point2Option(), None),
                'end': comopt.NoneOr(comopt.Point2Option(), None),
                'keep_pts': comopt.ListOfOptions(comopt.Point2Option(), [])
                }

    def _addrem_contour2(self):
        c = self.any_cont_by_name(self.get_option('base'))
        algo = self.get_option('algo')
        step = self.get_option('step')
        a0 = self.get_option('angle0')
        keepbnd = self.get_option('keepbnd')
        ned = self.get_option('nedges')
        crosses_cont = []
        for nm in self.get_option('crosses'):
            crosses_cont.append(self.any_cont_by_name(nm))
        crosses = []
        for cr in crosses_cont:
            crosses.append(cr.cdata)
        keeppts = self.get_option('keep_pts')
        start = self.get_option('start')
        end = self.get_option('end')
        ret = c2core.contour_partition(c.cdata, step, algo, a0, keepbnd,
                                       ned, crosses, keeppts, start, end)
        return [(Contour2(ret), self.get_option('name'))], []


class MatchedPartition(addremove.AbstractAddRemove):
    'Partition contour operation with condtitions'

    def __init__(self, argsdict):
        super(MatchedPartition, self).__init__(argsdict)

    def doc(self):
        return "Matched contour partition"

    @classmethod
    def _arguments_types(cls):
        return {'name': comopt.BasicOption(str, None),
                'base': comopt.BasicOption(str),
                'cconts': comopt.ListOfOptions(comopt.BasicOption(str), []),
                'step': comopt.BasicOption(float),
                'angle0': comopt.BasicOption(float, 1.),
                'infdist': comopt.BasicOption(str),
                'cpts': comopt.ListOfOptions(comopt.BasicOption(float), []),
                'power': comopt.BasicOption(float, 2.)
                }

    def _addrem_contour2(self):
        c = self.any_cont_by_name(self.get_option('base'))
        ccond, c_ccond = [], []
        for it in self.get_option('cconts'):
            ccond.append(self.any_cont_by_name(it))
            c_ccond.append(ccond[-1].cdata)
        cpts = self.get_option('cpts')
        step = self.get_option('step')
        infdist = self.get_option('infdist')
        power = self.get_option('power')
        a0 = self.get_option('angle0')
        c_ret = c2core.matched_partition(c.cdata, c_ccond, cpts, step,
                                         infdist, power, a0)
        return [(Contour2(c_ret), self.get_option('name'))], []


class ExtractSubcontours(addremove.AbstractAddRemove):
    "Extract subcontours"
    def __init__(self, arg):
        super(ExtractSubcontours, self).__init__(arg)

    @classmethod
    def _arguments_types(cls):
        return {'name': comopt.BasicOption(str, None),
                'src': comopt.BasicOption(str),
                'plist': comopt.ListOfOptions(comopt.Point2Option()),
                'project_to': comopt.BasicOption(str, 'vertex')
                }

    def _addrem_contour2(self):
        c = self.any_cont_by_name(self.get_option('src'))
        plist = self.get_option('plist')
        if self.get_option('project_to') == 'vertex':
            plist = c.closest_points(plist, 'vertex')
        elif self.get_option('project_to') == 'line':
            plist = c.closest_points(plist, 'edge')
        elif self.get_option('project_to') == 'corner':
            c2 = c.deepcopy()
            c2core.simplify(c2.cdata, 0., False)
            plist = c2.closest_points(plist, 'vertex')
        ret = c2core.extract_subcontours(c.cdata, plist)
        rr = [(Contour2(it), self.get_option('name')) for it in ret]
        return rr, []


class ConnectSubcontours(addremove.AbstractAddRemove):
    def __init__(self, arg):
        super(ConnectSubcontours, self).__init__(arg)

    @classmethod
    def _arguments_types(cls):
        return {'name': comopt.BasicOption(str, None),
                'src': comopt.ListOfOptions(comopt.BasicOption(str)),
                'fix': comopt.ListOfOptions(comopt.BasicOption(int), []),
                'close': comopt.BasicOption(str, 'no'),
                'shiftnext': comopt.BoolOption(False)
                }

    def _addrem_contour2(self):
        conts = map(self.any_cont_by_name, self.get_option('src'))
        cc = [c.cdata for c in conts]
        fx = self.get_option('fix')
        close = self.get_option('close')
        shift = self.get_option('shiftnext')
        ret = c2core.connect_subcontours(cc, fx, close, shift)
        return [(Contour2(ret), self.get_option('name'))], []
