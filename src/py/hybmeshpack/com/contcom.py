"contour commands"
import copy
import math
import hybmeshpack.basic.proc as bp
from hybmeshpack.basic.geom import Point2
import hybmeshpack.gdata.contour2 as contour2
from hybmeshpack.hmcore import c2 as c2core
from hybmeshpack import hmcore
from unite_grids import add_bc_from_cont
import command
import objcom


class AddRectCont(objcom.AbstractAddRemove):
    "Add rectangular contour"
    def __init__(self, argsdict):
        if 'name' not in argsdict:
            argsdict['name'] = "Contour1"
        if 'bnds' not in argsdict or len(argsdict['bnds']) == 0:
            argsdict['bnds'] = [0, 0, 0, 0]
        if len(argsdict['bnds']) < 4:
            b = argsdict['bnds'][-1]
            for i in range(4 - len(argsdict['bnds'])):
                argsdict['bnds'].append(b)
        super(AddRectCont, self).__init__(argsdict)

    @classmethod
    def _arguments_types(cls):
        return {'name': command.BasicOption(str),
                'p0': command.Point2Option(),
                'p1': command.Point2Option(),
                'bnds': command.ListOfOptions(command.BasicOption(int)),
                }

    def doc(self):
        return "Add rectangular contour"

    def _addrem_objects(self):
        c = contour2.ClosedContour2()
        op0, op1 = self.options['p0'], self.options['p1']
        b = self.options['bnds']
        p0, p1 = op0, Point2(op1.x, op0.y)
        p2, p3 = op1, Point2(op0.x, op1.y)
        c.append_points([p0, p1, p2, p3], b)
        return [], [], [(self.options['name'], c)]


class AddCircCont(objcom.AbstractAddRemove):
    "Add circular contour"
    def __init__(self, argsdict):
        if 'name' not in argsdict:
            argsdict['name'] = "Contour1"
        if 'bnd' not in argsdict:
            argsdict['bnd'] = 0
        super(AddCircCont, self).__init__(argsdict)

    @classmethod
    def _arguments_types(cls):
        return {'name': command.BasicOption(str),
                'p0': command.Point2Option(),
                'rad': command.BasicOption(float),
                'na': command.BasicOption(int),
                'bnd': command.BasicOption(int),
                }

    def doc(self):
        return "Add circlular contour"

    def _addrem_objects(self):
        c = contour2.ClosedContour2()
        op0 = self.options['p0']
        b = self.options['bnd']
        rad = self.options['rad']
        na = self.options['na']
        astep = 2 * math.pi / na
        pts = []
        for i in range(na):
            angle = i * astep
            pts.append(Point2(op0.x + rad * math.cos(angle),
                              op0.y + rad * math.sin(angle)))
        c.append_points(pts, [b] * len(pts))
        return [], [], [(self.options['name'], c)], []


class GridBndToContour(objcom.AbstractAddRemove):
    'Copy grid boundary to user contour'

    def __init__(self, arg):
        if "cont_name" not in arg:
            arg["cont_name"] = "Contour1"
        if "simplify" not in arg:
            arg['simplify'] = False
        if "separate" not in arg:
            arg['separate'] = False
        super(GridBndToContour, self).__init__(arg)
        self.simp_command = None

    def doc(self):
        return "Copy boundary from %s" % self.options['grid_name']

    @classmethod
    def _arguments_types(cls):
        return {'grid_name': command.BasicOption(str),
                'cont_name': command.BasicOption(str),
                'simplify': command.BoolOption(),
                'separate': command.BoolOption(),
                }

    def _addrem_objects(self):
        #copy contour
        try:
            g = self.receiver.get_grid(name=self.options['grid_name'])[2]
        except:
            raise command.ObjectNotFound(self.options['grid_name'])
        src = g.cont
        cont = contour2.Contour2.create_from_abstract(src)
        return [], [], [(self.options['cont_name'], cont)], []

    def _exec(self):
        r = super(GridBndToContour, self)._exec()
        nm = self._get_added_names()[1][0]
        #simplify and separate
        if self.options['simplify'] or self.options['separate']:
            sarg = {'separate': self.options['separate'],
                    'simplify': self.options['simplify'],
                    'names': [nm],
                    'new_names': [nm],
                    'keep_src': False,
                    'angle': 0.0}
            self.simp_command = SimplifyContours(sarg)
            if not self.simp_command.do(self.receiver):
                self.simp_command = None
        return r

    def _clear(self):
        self.simp_command = None
        super(GridBndToContour, self)._clear()

    def _undo(self):
        if self.simp_command:
            self.simp_command._undo()
        super(GridBndToContour, self)._undo()

    def _redo(self):
        super(GridBndToContour, self)._undo()
        if self.simp_command:
            self.simp_command._redo()


class SimplifyContours(objcom.AbstractAddRemove):
    'Separate and simplify user contour'
    def __init__(self, arg):
        if 'new_names' not in arg:
            arg['new_names'] = ['Contour1'] * len(arg['names'])
        if 'separate' not in arg:
            arg['separate'] = False
        if 'simplify' not in arg:
            arg['simplify'] = False
        if 'keep_src' not in arg:
            arg['keep_src'] = True
        if 'angle' not in arg:
            arg['angle'] = 1.0
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
            keep_src - keep source contours after procedure
        """
        return {'names': command.ListOfOptions(command.BasicOption(str)),
                'new_names': command.ListOfOptions(command.BasicOption(str)),
                'simplify': command.BoolOption(),
                'separate': command.BoolOption(),
                'angle': command.BasicOption(float),
                'keep_src': command.BoolOption(),
                }

    def doc(self):
        return "Simplify contours: %s" % ' '.join(self.conts)

    def _addrem_objects(self):
        simp, sep, an, keep, nnames = \
            [self.options[x] for x in
             ['simplify', 'separate', 'angle', 'keep_src', 'new_names']]
        added, removed = [], []

        #make copies
        for i, nm in enumerate(self.options['names']):
            try:
                cont = self.receiver.get_ucontour(name=nm)[2]
                if not keep:
                    removed.append(nm)
            except KeyError:
                cont = self.receiver.get_grid(name=nm)[2].cont

            added.append((nnames[i],
                          contour2.Contour2.create_from_abstract(cont)))

        #simplify
        if simp:
            for i in range(len(added)):
                cc = added[i][1].simplify(an)
                if cc is not None:
                    added[i] = (added[i][0], cc)

        #separate
        if sep:
            newadded = []
            for i in range(len(added)):
                sep = added[i][1].separate()
                if sep is None:
                    newadded.append(added[i])
                else:
                    for s in sep:
                        newadded.append((added[i][0], s))
            added = newadded

        return [], [], added, removed


class UniteContours(objcom.AbstractAddRemove):
    "unite contours to single contour with complicated connection"
    def __init__(self, arg):
        if 'keep_src' not in arg:
            arg['keep_src'] = True
        if "name" not in arg:
            arg["name"] = "Contour1"
        super(UniteContours, self).__init__(arg)

    def doc(self):
        return "Unite contours: %s" % ' '.join(self.srcnames)

    @classmethod
    def _arguments_types(cls):
        """
            name - united contour name
            sources - names of contours to unite
            keep_src - whether to keep or remove source contours
        """
        return {'name': command.BasicOption(str),
                'sources': command.ListOfOptions(command.BasicOption(str)),
                'keep_src': command.BoolOption(),
                }

    def _addrem_objects(self):
        conts = []
        for n in self.options['sources']:
            try:
                conts.append(self.receiver.contours2[n])
            except KeyError:
                conts.append(self.receiver.grids2[n].cont)

        #operate
        newcont = contour2.Contour2.create_from_abstract(conts[0])
        for c in conts[1:]:
            newcont.add_from_abstract(c)
        ac = [(self.options['name'], newcont)]
        rc = []
        if not self.options['keep_src']:
            for c in self.options['sources']:
                if c in self.receiver.contours2:
                    rc.append(c)
        return [], [], ac, rc


class EditBoundaryType(command.Command):
    "Add/edit/remove boundary type"
    def __init__(self, arg):
        a = copy.deepcopy(arg)
        if 'remindex' not in a:
            a['remindex'] = None
        if 'index' not in a:
            a['index'] = None
        if 'name' not in a:
            a['name'] = 'boundary1'
        if 'color' not in a:
            a['color'] = [0, 0, 0]
        super(EditBoundaryType, self).__init__(a)
        self._backup1 = None
        self._backup2 = None

    @classmethod
    def _arguments_types(cls):
        """ (int remindex, int index, str name, (char)*3 color)

            removes boundary with index remindex and
            sets new boundary (index, name, color).
            remindex=None if nothing to delete.
            index=None if nothing to set.
            name, color could be default
        """
        return {'remindex': command.BasicOption(int),
                'index': command.BasicOption(int),
                'name': command.BasicOption(str),
                'color': command.ListOfOptions(command.BasicOption(int))
                }

    def doc(self):
        return "Edit boundary type"

    def _exec(self):
        try:
            ri, i, n, c = [self.options[x] for x in ['remindex', 'index',
                                                     'name', 'color']]
            #1 remove
            if ri is not None:
                self._backup1 = copy.deepcopy(
                    self.receiver.boundary_types.get(index=ri))
                self.receiver.boundary_types.rem_bnd(ri)
            #2 set
            if i is not None:
                try:
                    self._backup2 = copy.deepcopy(
                        self.receiver.boundary_types.get(index=i))
                except KeyError:
                    self._backup2 = None
                self.receiver.boundary_types.set_bnd(i, n, c)
            return True
        except Exception as e:
            raise command.ExecutionError(str(e), self, e)

    def _clear(self):
        self._backup1 = None
        self._backup2 = None

    def _undo(self):
        if self.options['index'] is not None:
            self.receiver.boundary_types.rem_bnd(self.options['index'])
        if self._backup2 is not None:
            self.receiver.boundary_types.set_bnd(
                self._backup2.index, self._backup2.name, self._backup2.color)
        if self._backup1 is not None:
            self.receiver.boundary_types.set_bnd(
                self._backup1.index, self._backup1.name, self._backup1.color)
        self._clear()

    def _redo(self):
        self._exec()


class BoundaryPickerOption(command.BasicOption):
    def __init__(self):
        super(BoundaryPickerOption, self).__init__()

    def serial(self, v):
        'before writing'
        a = copy.deepcopy(v)
        for k, v in a.items():
            if k != 'name':
                a[k] = bp.compress_int_list(a[k])
        return a

    def unserial(self, s):
        'after reading'
        a = copy.deepcopy(s)
        for k, v in a.items():
            if k != 'name':
                a[k] = bp.int_list_from_compress(a[k])
        return a


class SetBTypeToContour(command.Command):
    'sets integer boundary types to contours edges'

    def __init__(self, arg):
        super(SetBTypeToContour, self).__init__(arg)
        self.backup_bnd = []
        self.new_bnd = []

    def doc(self):
        return "Set boundary types to contours: " + \
            ', '.join([x.name for x in self.conts_opts])

    @classmethod
    def _arguments_types(cls):
        """
            btypes: [
                {
                    'name': contour name,
                    bnd_type: [list of edges indicies],
                },
                {...}, ....
            ]
        """
        return {'btypes': command.ListOfOptions(BoundaryPickerOption())}

    def _exec(self):
        for co in self.options['btypes']:
            bu = {}
            new = {}
            try:
                cont = self.receiver.get_any_contour(co['name'])
            except:
                raise command.ObjectNotFound(co['name'])
            #fill new and backup dicts
            for b, ei in co.items():
                if b != 'name':
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
        for co, bu in zip(self.options['btypes'], self.backup_bnd):
            cont = self.receiver.get_any_contour(co['name'])
            cont.add_edge_bnd(bu)

    def _redo(self):
        for co, new in zip(self.options['btypes'], self.new_bnd):
            cont = self.receiver.get_any_contour(co['name'])
            cont.add_edge_bnd(new)


class CreateContour(objcom.AbstractAddRemove):
    'Manually enter contour points'
    def __init__(self, argsdict):
        if "name" not in argsdict:
            argsdict["name"] = "Contour1"
        if "bnds" not in argsdict:
            argsdict["bnds"] = 0
        if not isinstance(argsdict["bnds"], list):
            argsdict["bnds"] =\
                [argsdict["bnds"]] * (len(argsdict["points"]) - 1)
        super(CreateContour, self).__init__(argsdict)

    def doc(self):
        return "Create contour from sequence of points"

    @classmethod
    def _arguments_types(cls):
        return {'name': command.BasicOption(str),
                'points': command.ListOfOptions(command.Point2Option()),
                'bnds': command.ListOfOptions(command.BasicOption(int)),
                }

    def _addrem_objects(self):
        pt = self.options["points"]
        eds = [[i, i + 1] for i in range(len(pt) - 1)]
        p0 = self.options["points"][0]
        plast = self.options["points"][-1]
        if p0.x == plast.x and p0.y == plast.y:
            eds[-1][1] = 0
            pt = pt[:-1]
        cont = contour2.Contour2.create_from_point_set(
            pt, eds, self.options["bnds"])
        return [], [], [(self.options["name"], cont)], []


class CreateSpline(objcom.AbstractAddRemove):
    " Creates spline"

    def __init__(self, argsdict):
        if "name" not in argsdict:
            argsdict["name"] = "Contour1"
        if "bnds" not in argsdict:
            argsdict["bnds"] = 0
        super(CreateSpline, self).__init__(argsdict)

    def doc(self):
        return "Create spline from sequence of points"

    @classmethod
    def _arguments_types(cls):
        return {'name': command.BasicOption(str),
                'points': command.ListOfOptions(command.Point2Option()),
                'bnds': command.ListOfOptions(command.BasicOption(int)),
                'nedges': command.BasicOption(int),
                }

    def _addrem_objects(self):
        so = self.options
        # copy to c
        plist = []
        for p in so['points']:
            plist.extend([p.x, p.y])
        cplist = hmcore.list_to_c(plist, float)
        c_bt = hmcore.list_to_c(so['bnds'], int)
        # 4. call c procedure
        c_ret = 0
        try:
            c_ret, c_bnd = c2core.spline(cplist, c_bt, so['nedges'])
            ret = c2core.cont2_from_c(c_ret, c_bnd)
            return [], [], [(so["name"], ret)], []
        except Exception as e:
            raise command.ExecutionError("Spline builder error", self, e)
        finally:
            if c_ret != 0:
                c2core.free_cont2(c_ret)


class ClipDomain(objcom.AbstractAddRemove):
    'Domain clipping operation'

    def __init__(self, argsdict):
        if "name" not in argsdict:
            argsdict["name"] = "Contour1"
        if "simplify" not in argsdict:
            argsdict["simplify"] = True
        super(ClipDomain, self).__init__(argsdict)

    def doc(self):
        return "Domain clipping operation"

    @classmethod
    def _arguments_types(cls):
        return {'name': command.BasicOption(str),
                'c1': command.BasicOption(str),
                'c2': command.BasicOption(str),
                'oper': command.BasicOption(str),
                'simplify': command.BasicOption(bool),
                }

    def _addrem_objects(self):
        cont1 = self.receiver.get_any_contour(self.options['c1'])
        cont2 = self.receiver.get_any_contour(self.options['c2'])
        cc1 = c2core.cont2_to_c(cont1)
        cc2 = c2core.cont2_to_c(cont2)
        if self.options['oper'] == 'union':
            op = 1
        elif self.options['oper'] == 'difference':
            op = 2
        elif self.options['oper'] == 'intersection':
            op = 3
        elif self.options['oper'] == 'xor':
            op = 4
        else:
            raise Exception("Invalid operation %s" % str(self.options['oper']))
        cres = c2core.clip_domain(cc1, cc2, op, self.options['simplify'])
        if cres is not None:
            res = c2core.cont2_from_c(cres)
            add_bc_from_cont(res, cont1, cres, cc1)
            add_bc_from_cont(res, cont2, cres, cc2)
            c2core.free_cont2(cres)
        else:
            res = None
        c2core.free_cont2(cc1)
        c2core.free_cont2(cc2)

        if res is not None:
            return [], [], [(self.options["name"], res)], []
        else:
            return [], [], [], []


class PartitionContour(objcom.AbstractAddRemove):
    'Partition contour operation'

    def __init__(self, argsdict):
        if "name" not in argsdict:
            argsdict["name"] = "Contour1"
        super(PartitionContour, self).__init__(argsdict)

    def doc(self):
        return "Contour partition"

    @classmethod
    def _arguments_types(cls):
        return {'name': command.BasicOption(str),
                'algo': command.BasicOption(str),  # 'const'/'ref_points'
                'step': command.ListOfOptions(command.BasicOption(float)),
                'angle0': command.BasicOption(float),
                'base': command.BasicOption(str),
                'keepbnd': command.BoolOption(),
                'nedges': command.NoneOr(command.BasicOption(int)),
                'crosses': command.ListOfOptions(command.BasicOption(str)),
                }

    def _addrem_objects(self):
        so = self.options
        # 1. find contour
        cont = self.any_cont_by_name(so['base'])
        # 2. get boundary types array
        bt = []
        for i in range(cont.n_edges()):
            bt.append(cont.edge_bnd(i))
        # 3. copy to c
        c_cont = c2core.cont2_to_c(cont)
        c_bt = hmcore.list_to_c(bt, int)
        c_step = hmcore.list_to_c(so['step'], float)
        # 4. copy crosses conditions
        c_crosses = []
        for nm in so['crosses']:
            c = self.any_cont_by_name(nm)
            c_crosses.append(c2core.cont2_to_c(c))
        # 5. call c procedure
        try:
            c_ret, c_bnd = c2core.contour_partition(c_cont, c_bt, c_step,
                                                    so['algo'],
                                                    so['angle0'],
                                                    so['keepbnd'],
                                                    so['nedges'],
                                                    c_crosses)
            ret = c2core.cont2_from_c(c_ret, c_bnd)
            c2core.free_cont2(c_ret)
            return [], [], [(so["name"], ret)], []
        except Exception as e:
            raise command.ExecutionError("Contour partition error", self, e)
        finally:
            c2core.free_cont2(c_cont)
            for c in c_crosses:
                c2core.free_cont2(c)


class MatchedPartition(objcom.AbstractAddRemove):
    'Partition contour operation with condtitions'

    def __init__(self, argsdict):
        if "name" not in argsdict:
            argsdict["name"] = "Contour1"
        super(MatchedPartition, self).__init__(argsdict)

    def doc(self):
        return "Matched contour partition"

    @classmethod
    def _arguments_types(cls):
        return {'name': command.BasicOption(str),
                'base': command.BasicOption(str),
                'cconts': command.ListOfOptions(command.BasicOption(str)),
                'step': command.BasicOption(float),
                'angle0': command.BasicOption(float),
                'infdist': command.BasicOption(str),
                'cpts': command.ListOfOptions(command.BasicOption(float)),
                'power': command.BasicOption(float)
                }

    def _addrem_objects(self):
        so = self.options
        # 1. find contour
        cont = self.any_cont_by_name(so['base'])
        conds = []
        for s in so['cconts']:
            conds.append(self.any_cont_by_name(s))
        c_cont, c_ret, c_conds = 0, 0, []
        try:
            c_cont = c2core.cont2_to_c(cont)
            for c in conds:
                c_conds.append(c2core.cont2_to_c(c))
            c_rp = c2core.list_to_c(so['cpts'], float)
            c_ret = c2core.matched_partition(c_cont, c_conds, c_rp,
                                             so['step'],
                                             so['infdist'],
                                             so['power'],
                                             so['angle0'])
            ret = c2core.cont2_from_c(c_ret)
            add_bc_from_cont(ret, cont, c_ret, c_cont, force=3)

            return [], [], [(so["name"], ret)], []
        except Exception as e:
            raise command.ExecutionError("Contour partition error", self, e)
        finally:
            c2core.free_cont2(c_ret) if c_ret != 0 else None
            c2core.free_cont2(c_cont) if c_cont != 0 else None
            for c in c_conds:
                c2core.free_cont2(c)


class ExtractContour(objcom.AbstractAddRemove):
    def __init__(self, arg):
        if "name" not in arg:
            arg["name"] = "Contour1"
        super(ExtractContour, self).__init__(arg)

    def doc(self):
        return "Extract contour from source"

    @classmethod
    def _arguments_types(cls):
        return {'name': command.BasicOption(str),
                'src': command.BasicOption(str),
                'p0': command.Point2Option(),
                'p1': command.Point2Option(),
                'project_to': command.BasicOption(str)
                }

    def _addrem_objects(self):
        so = self.options
        cont = self.any_cont_by_name(so['src'])

        c_cont, c_ret = 0, 0
        try:
            c_cont = c2core.cont2_to_c(cont)
            c_ret = c2core.extract_contour(c_cont,
                                           so['p0'], so['p1'],
                                           so['project_to'])
            ret = c2core.cont2_from_c(c_ret)
            add_bc_from_cont(ret, cont, c_ret, c_cont, force=3)

            return [], [], [(so["name"], ret)], []
        except Exception as e:
            raise command.ExecutionError("ExtractContour error", self, e)
        finally:
            c2core.free_cont2(c_ret) if c_ret != 0 else None
            c2core.free_cont2(c_cont) if c_cont != 0 else None
