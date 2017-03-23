import copy
import math
from hybmeshpack.basic import geom as bgeom
from hybmeshpack.hmcore import g2 as g2core
from hybmeshpack.hmcore import c2 as c2core
from hybmeshpack.gdata.grid2 import Grid2
from hybmeshpack.gdata.cont2 import Contour2, closest_contour
import addremove
import comopt as co


class NewGridCommand(addremove.AbstractAddRemove):
    "Command with a new grid addition"
    def __init__(self, argsdict):
        super(NewGridCommand, self).__init__(argsdict)

    def _addrem_grid2(self):
        g = self._build_grid()
        if g is not None:
            return [(Grid2(g), self.get_option('name'))], []
        else:
            return [], []

    # function for overloading
    def _build_grid(self):
        '-> grid2.Grid2.cdata pointer'
        raise NotImplementedError


class AddUnfRectGrid(NewGridCommand):
    "Add uniform rectangular grid "
    def __init__(self, argsdict):
        super(AddUnfRectGrid, self).__init__(argsdict)

    @classmethod
    def _arguments_types(cls):
        return {'name': co.BasicOption(str, None),
                'p0': co.Point2Option((0., 0.)),
                'p1': co.Point2Option((1., 1.)),
                'nx': co.BasicOption(int, 3),
                'ny': co.BasicOption(int, 3),
                'custom_x': co.ListOfOptions(co.BasicOption(float), []),
                'custom_y': co.ListOfOptions(co.BasicOption(float), []),
                'bnds': co.ListOfOptions(co.BasicOption(int), [0, 0, 0, 0])
                }

    def _build_grid(self):
        xdata, ydata = self.get_option('custom_x'), self.get_option('custom_y')
        bnds = self.get_option('bnds')
        p0, p1 = self.get_option('p0'), self.get_option('p1')
        bgeom.build_seg(p0[0], p1[0], self.get_option('nx'), xdata)
        bgeom.build_seg(p0[1], p1[1], self.get_option('ny'), ydata)
        return g2core.build_rect_grid(xdata, ydata, bnds)


class AddCustomRectGrid(NewGridCommand):
    "Add custom rectangular grid "
    def __init__(self, argsdict):
        super(AddCustomRectGrid, self).__init__(argsdict)

    @classmethod
    def _arguments_types(cls):
        return {'name': co.BasicOption(str, None),
                'algo': co.BasicOption(str),
                'left': co.BasicOption(str),
                'bot': co.BasicOption(str),
                'right': co.NoneOr(co.BasicOption(str), None),
                'top': co.NoneOr(co.BasicOption(str), None),
                'her_w': co.ListOfOptions(co.BasicOption(float),
                                          [0., 0., 0., 0.]),
                'return_invalid': co.BoolOption(False),
                }

    def _build_grid(self):
        c = [self.get_option(n) for n in ['left', 'bot', 'right', 'top']]
        c = [self.any_cont_by_name(x) if x is not None else None
             for x in c]
        if c[2] is None:
            c[2] = c[0].deepcopy()
        if c[3] is None:
            c[3] = c[1].deepcopy()
        cc = [x.cdata if x is not None else None for x in c]
        cb = self.ask_for_callback()
        return g2core.custom_rectangular_grid(
            self.get_option('algo'), cc[0], cc[1], cc[2], cc[3],
            self.get_option('her_w'), self.get_option('return_invalid'), cb)


class AddCirc4Grid(NewGridCommand):
    "Add circular quadrangular grid "
    def __init__(self, argsdict):
        super(AddCirc4Grid, self).__init__(argsdict)

    @classmethod
    def _arguments_types(cls):
        return {'name': co.BasicOption(str, None),
                'algo': co.BasicOption(str, 'linear'),
                'p0': co.Point2Option((0., 0.)),
                'rad': co.BasicOption(float, 1.),
                'step': co.BasicOption(float, 0.1),
                'sqrside': co.BasicOption(float, 1.),
                'rcoef': co.BasicOption(float, 1.)
                }

    def _build_grid(self):
        return g2core.circ4grid(
            self.get_option('algo'), self.get_option('p0'),
            self.get_option('rad'), self.get_option('step'),
            self.get_option('sqrside'), self.get_option('rcoef'))


class AddUnfCircGrid(NewGridCommand):
    "Add uniform circular grid "

    def __init__(self, argsdict):
        super(AddUnfCircGrid, self).__init__(argsdict)

    @classmethod
    def _arguments_types(cls):
        return {'name': co.BasicOption(str, None),
                'p0': co.Point2Option(),
                'rad': co.BasicOption(float, 1.),
                'na': co.BasicOption(int, 8),
                'nr': co.BasicOption(int, 4),
                'coef': co.BasicOption(float, 1.),
                'is_trian': co.BoolOption(True),
                'custom_r': co.ListOfOptions(co.BasicOption(float), []),
                'custom_a': co.ListOfOptions(co.BasicOption(float), []),
                'bnd': co.BasicOption(int, 0),
                }

    def _build_grid(self):
        rdata, adata = self.get_option('custom_r'), self.get_option('custom_a')
        rad = self.get_option('rad')
        # radius segmentation calculation
        if len(rdata) == 0:
            rdata = bgeom.div_range(
                0, rad, self.get_option('nr'), self.get_option('coef'))
        else:
            if len(rdata) > 1 and rdata[0] != 0:
                rdata.insert(0, 0.0)
            bgeom.build_seg(0, rad, self.get_option('nr'), rdata)
        # arc segmentation calculation
        if len(adata) < 2:
            if len(adata) == 1:
                adata[0] = adata[0] / rdata[-1] / math.pi * 180.
            bgeom.build_seg(0., 360, self.get_option('na'), adata)
        else:
            amx = adata[-1] - adata[0]
            adata = [float(x) / amx * 360 for x in adata]

        return g2core.build_circ_grid(
            self.get_option('p0'), rdata, adata, self.get_option('is_trian'),
            self.get_option('bnd'))


class AddUnfHexGrid(NewGridCommand):
    "Add uniform circular ring"

    def __init__(self, argsdict):
        super(AddUnfHexGrid, self).__init__(argsdict)

    @classmethod
    def _arguments_types(cls):
        return {'name': co.BasicOption(str, None),
                'area': co.ListOfOptions(co.BasicOption(float)),
                'crad': co.BasicOption(float),
                'strict': co.BoolOption(False),
                }

    def _build_grid(self):
        return g2core.regular_hex_grid(self.get_option('area'),
                                       self.get_option('crad'),
                                       self.get_option('strict'))


class AddUnfRingGrid(NewGridCommand):
    "Add uniform circular ring"

    def __init__(self, argsdict):
        super(AddUnfRingGrid, self).__init__(argsdict)

    @classmethod
    def _arguments_types(cls):
        return {'name': co.BasicOption(str, None),
                'p0': co.Point2Option(),
                'radinner': co.BasicOption(float),
                'radouter': co.BasicOption(float),
                'na': co.BasicOption(int),
                'nr': co.BasicOption(int),
                'coef': co.BasicOption(float, 1.),
                'bnd': co.ListOfOptions(co.BasicOption(int), [0, 0])
                }

    def _build_grid(self):
        rdata = bgeom.div_range(
            self.get_option('radinner'), self.get_option('radouter'),
            self.get_option('nr'), self.get_option('coef'))
        adata = bgeom.div_range(0., 360., self.get_option('na'), 1.)
        return g2core.build_ring_grid(self.get_option('p0'), rdata, adata,
                                      self.get_option('bnd'))


class AddTriGrid(NewGridCommand):
    "Add grid in triangle area"
    def __init__(self, argsdict):
        super(AddTriGrid, self).__init__(argsdict)

    @classmethod
    def _arguments_types(cls):
        return {'name': co.BasicOption(str, None),
                'vertices': co.ListOfOptions(co.Point2Option()),
                'nedge': co.BasicOption(int),
                'bnd': co.ListOfOptions(co.BasicOption(int), [0, 0, 0])
                }

    def _build_grid(self):
        return g2core.build_tri_grid(
            self.get_option('vertices'), self.get_option('nedge'),
            self.get_option('bnd'))


class ExcludeContours(addremove.AbstractAddRemove):
    "Exclude contour area from grids"

    def __init__(self, arg):
        super(ExcludeContours, self).__init__(arg)

    @classmethod
    def _arguments_types(cls):
        """ name - new grid name,
            grid_name - source grid name,
            cont_names - contours names,
            is_inner - exclude inner (True) or outer (False) contour area
        """
        return {'name': co.BasicOption(str, None),
                'grid_name': co.BasicOption(str),
                'cont_names': co.ListOfOptions(co.BasicOption(str)),
                'is_inner': co.BoolOption(),
                }

    def _addrem_grid2(self):
        cbtot = self.ask_for_callback()
        base = self.grid2_by_name(self.get_option('grid_name'))

        clist = self.get_option('cont_names')
        for i, cname in enumerate(clist):
            cb = cbtot.subcallback(i, len(clist))
            cont = self.any_cont_by_name(cname)
            res = g2core.grid_excl_cont(
                base.cdata, cont.cdata, self.get_option('is_inner'), cb)
            base = Grid2(res)

        return [(base, self.get_option('name'))], []


class UniteGrids(addremove.AbstractAddRemove):
    def __init__(self, argsdict):
        super(UniteGrids, self).__init__(argsdict)

    def doc(self):
        ret = "Unite Grids: %s " % self.options['base']
        for s in self.options['plus']:
            ret += s['name'] + " "
        return ret[:-1]

    class Option(co.SubDictOption):
        'name, buf - dictionary'
        def __init__(self):
            super(UniteGrids.Option, self).__init__(
                name=co.BasicOption(str),
                buf=co.BasicOption(float),
            )

    @classmethod
    def _arguments_types(cls):
        """ name - name of the new grid,
            fix_bnd - whether to fix all boundary points (=False)
            keepsrc - whether to remove source grids (=True)
            empty_holes - keep all empty zone in 2nd grid (=False)
            base - name of base grid,
            plus - list of UniteOptions: entries define each next imposition
        """
        return {'name': co.BasicOption(str, None),
                'fix_bnd': co.BoolOption(False),
                'empty_holes': co.BoolOption(False),
                'base': co.BasicOption(str),
                'plus': co.ListOfOptions(UniteGrids.Option()),
                'angle0': co.BasicOption(float, 1.),
                'filler': co.BasicOption(str, '3')
                }

    def __get_grid(self, ind):
        """ returns (Grid, buffer size) by order of addition """
        if ind == 0:
            gname = self.get_option('base')
            b = 0
        else:
            gname = self.get_option('plus')[ind - 1]['name']
            b = self.get_option('plus')[ind - 1]['buf']
        return self.grid2_by_name(gname), b

    def _addrem_grid2(self):
        cbtotal = self.ask_for_callback()
        imax = len(self.get_option('plus'))
        wg, _ = self.__get_grid(0)

        for i in range(1, imax + 1):
            cb = cbtotal.subcallback(i - 1, imax)
            g, b = self.__get_grid(i)
            ret = g2core.unite_grids(
                wg.cdata, g.cdata, b, self.get_option('fix_bnd'),
                self.get_option('empty_holes'), self.get_option('angle0'),
                self.get_option('filler'), cb)
            wg = Grid2(ret)

        return [(wg, self.get_option('name'))], []


class BuildBoundaryGrid(NewGridCommand):
    def __init__(self, kwargs):
        super(BuildBoundaryGrid, self).__init__(kwargs)

    class Option(co.SubDictOption):
        """ {
            source: str - name of source contour
            partition: [0.0,...] -- layer partition starts from 0
            direction: int    - 1/-1 for positive/negative direction
            mesh_cont: 0/1/2/3/4 -- contour meshing strategy:
                0 - No meshing, use original contour nodes
                1 - keep all original contour nodes,
                2 - keep only non-collinear contour nodes,
                3 - ignore all existing contour nodes
                4 - incremental stepping
            mesh_cont_step: double -- contour mesh step if mesh_cont!=0
            step_start, step_end: contour mesh step for incremental
                                  stepping
            start: Point   -- Start Point of Section.
            end: Point     -- End Point of Section.
            force_conf: Bool -- force conformal mappings
            algo_acute/right/straight/reentr: float (deg) --
                 maximum angles for choosing corners treatment algo.
        }
        """
        def __init__(self):
            super(BuildBoundaryGrid.Option, self).__init__(
                source=co.BasicOption(str),
                partition=co.ListOfOptions(co.BasicOption(float)),
                direction=co.BasicOption(int),
                mesh_cont=co.BasicOption(int, 0),
                mesh_cont_step=co.BasicOption(float, 1.),
                step_start=co.BasicOption(float, 1.),
                step_end=co.BasicOption(float, 1.),
                algo_acute=co.BasicOption(float, 45.),
                algo_right=co.BasicOption(float, 120.),
                algo_straight=co.BasicOption(float, 240.),
                algo_reentr=co.BasicOption(float, 300.),
                start=co.NoneOr(co.Point2Option(), None),
                end=co.NoneOr(co.Point2Option(), None),
                force_conf=co.BoolOption(False)
            )

    @classmethod
    def _arguments_types(cls):
        return {
            'name': co.BasicOption(str, None),
            'opt': co.ListOfOptions(cls.Option())
        }

    def __adjust_arg1(self, op):
        """ changes 'source' from name to AbstractContour object
            if end points are None changes them to actual values
            and adds options for multiply connected domains
        """
        # separate contour
        sep = c2core.quick_separate(op['source'].cdata)
        sep = [Contour2(s) for s in sep]

        # if end points are None:
        #    transform them to real contour point
        #    create set of options for multiply connected domains
        newops = []
        if op['start'] is None or op['end'] is None:
            for s in sep:
                p0, p1 = s.end_points()
                if p0 is None:
                    continue
                addto = op
                if op['start'] is not None and op['end'] is not None:
                    _t, op['source'] = op['source'], None
                    newops.append(copy.deepcopy(op))
                    op['source'] = _t
                    addto = newops[-1]
                # addto['source'] = s
                addto['source'] = op['source']
                addto['start'] = copy.deepcopy(p0)
                addto['end'] = copy.deepcopy(p1)
            if op['start'] is None or op['end'] is None:
                raise Exception(
                    "cannot find correct source for boundary grid", self)

        # if points coincide find closest subcontour
        elif op['start'][0] == op['end'][0] and op['start'][1] == op['end'][1]:
            subcont = closest_contour(sep, op['start'])
            p0, p1 = subcont.end_points()
            if p0 != p1:
                op['start'] = p0
                op['end'] = p1

        return newops

    def __adjust_arg(self, arg):
        na = []
        for a in arg:
            na.extend(self.__adjust_arg1(a))
        arg.extend(na)

    def _build_grid(self):
        cb = self.ask_for_callback()
        arg = copy.deepcopy(self.get_option('opt'))
        usedconts = {}
        for op in arg:
            nm = op['source']
            if nm not in usedconts:
                usedconts[nm] = self.any_cont_by_name(nm)
            op['source'] = usedconts[nm]

        self.__adjust_arg(arg)
        cont = [a['source'] for a in arg]
        for a, b in zip(arg, cont):
            a['source'] = b.cdata
        return g2core.boundary_layer_grid(arg, cb)


class HealGrid(addremove.AbstractAddRemove):
    "Heal grid"

    def __init__(self, arg):
        super(HealGrid, self).__init__(arg)

    @classmethod
    def _arguments_types(cls):
        """ name - modified grid name,
            simp_bnd - angle (degree): is the maximum angles at which
                two adjacent boundary edges of the same cell will be merged
                (-1) ignores this procedure.
            convex - angle (degree). 0 - removes all concave cell segments,
                180 - only degenerate segments.
        """
        return {'name': co.BasicOption(str),
                'simp_bnd': co.BasicOption(float, -1.),
                'convex': co.BasicOption(float, -1.)
                }

    def __simplify_bnd(self, og, angle):
        ret = g2core.simplify_grid_boundary(og.cdata, angle)
        return Grid2(ret)

    def __convex_cells(self, og, angle):
        ret = g2core.convex_cells(og.cdata, angle)
        return Grid2(ret)

    def _addrem_grid2(self):
        g = self.grid2_by_name(self.get_option('name'))
        orig = g
        if self.get_option('simp_bnd') >= 0:
            g = self.__simplify_bnd(g, self.get_option('simp_bnd'))
        if self.get_option('convex') >= 0:
            g = self.__convex_cells(g, self.get_option('convex'))

        if g != orig:
            return [(g, self.get_option('name'))], [self.get_option('name')]
        else:
            return [], []


class MapGrid(NewGridCommand):
    "Maps grid area on given contour"

    def __init__(self, kwargs):
        super(MapGrid, self).__init__(kwargs)

    @classmethod
    def _arguments_types(cls):
        """ name - modified grid name,
            base - identifier of the basic grid
            target - identifier of contour
            base_points - points of the basic grid contour
            target_poitns - points of the target contour
            snap - snapping algo ("no", "add_vertices", "shift_vertices")
            btypes - source of boundary features ("from_grid", "from_contour")
            algo ('inverse_laplace', 'direct_laplace')
        """
        return {'name': co.BasicOption(str, None),
                'base': co.BasicOption(str),
                'target': co.BasicOption(str),
                'base_points': co.ListOfOptions(co.Point2Option()),
                'target_points': co.ListOfOptions(co.Point2Option()),
                'snap': co.BasicOption(str, 'no'),
                'btypes': co.BasicOption(str, 'from_contour'),
                'algo': co.BasicOption(str, 'inverse_laplace'),
                'is_reversed': co.BoolOption(False),
                'return_invalid': co.BoolOption(False),
                }

    def _build_grid(self):
        cb = self.ask_for_callback()
        g = self.grid2_by_name(self.get_option('base'))
        c = self.any_cont_by_name(self.get_option('target'))
        ret = g2core.map_grid(g.cdata, c.cdata,
                              self.get_option('base_points'),
                              self.get_option('target_points'),
                              self.get_option('snap'),
                              self.get_option('btypes') == 'from_contour',
                              self.get_option('algo'),
                              self.get_option('is_reversed'),
                              self.get_option('return_invalid'), cb)
        return ret


class UnstructuredFillArea(NewGridCommand):
    def __init__(self, argsdict):
        super(UnstructuredFillArea, self).__init__(argsdict)

    @classmethod
    def _arguments_types(cls):
        return {'name': co.BasicOption(str, None),
                'domain': co.ListOfOptions(co.BasicOption(str)),
                'constr': co.ListOfOptions(co.BasicOption(str), []),
                'pts': co.ListOfOptions(co.BasicOption(float), []),
                }

    def call_c_proc(self, c_dom, c_con, c_pts):
        raise NotImplementedError

    def _build_grid(self):
        # unite domain and constraints
        domc, domcc = [], []
        for n in self.get_option('domain'):
            domc.append(self.any_cont_by_name(n))
            domcc.append(domc[-1].cdata)
        fulldom = Contour2(c2core.unite_contours(domcc))

        conc, concc = [], []
        if len(self.get_option('constr')) > 0:
            for n in self.get_option('constr'):
                conc.append(self.any_cont_by_name(n))
                concc.append(conc[-1].cdata)
            fullconst = Contour2(c2core.unite_contours(concc))
        else:
            fullconst = Contour2.empty()

        # call builder
        pts = self.get_option('pts')
        return self.call_c_proc(fulldom.cdata, fullconst.cdata, pts)


class TriangulateArea(UnstructuredFillArea):
    def __init__(self, argsdict):
        super(TriangulateArea, self).__init__(argsdict)

    def call_c_proc(self, c_dom, c_con, pts):
        return g2core.unstructed_fill(c_dom, c_con, pts, '3')


class QuadrangulateArea(UnstructuredFillArea):
    def __init__(self, argsdict):
        super(QuadrangulateArea, self).__init__(argsdict)

    def call_c_proc(self, c_dom, c_con, pts):
        return g2core.unstructed_fill(c_dom, c_con, pts, '4')


class PebiFill(UnstructuredFillArea):
    def __init__(self, argsdict):
        super(PebiFill, self).__init__(argsdict)

    def call_c_proc(self, c_dom, c_con, pts):
        return g2core.unstructed_fill(c_dom, c_con, pts, 'pebi')


class StripeGrid(NewGridCommand):
    def __init__(self, kwargs):
        super(StripeGrid, self).__init__(kwargs)

    @classmethod
    def _arguments_types(cls):
        """ partition: increasing list of floats,
            tip: 'no', 'radial'
            bnd: boundary for left, bottom, right, top sections
        """
        return {'name': co.BasicOption(str, None),
                'source': co.BasicOption(str),
                'partition': co.ListOfOptions(co.BasicOption(float)),
                'tip': co.BasicOption(str, 'no'),
                'bnd': co.ListOfOptions(co.BasicOption(int), [0, 0, 0, 0])
                }

    def _build_grid(self):
        c = self.any_cont_by_name(self.get_option('source'))
        cb = self.ask_for_callback()
        return g2core.stripe_grid(c.cdata, self.get_option('partition'),
                                  self.get_option('tip'),
                                  self.get_option('bnd'), cb)


class SnapToContour(NewGridCommand):
    def __init__(self, kwargs):
        super(SnapToContour, self).__init__(kwargs)

    @classmethod
    def _arguments_types(cls):
        return {'name': co.BasicOption(str, None),
                'grid': co.BasicOption(str),
                'cont': co.BasicOption(str),
                'gp1': co.Point2Option(),
                'gp2': co.Point2Option(),
                'cp1': co.NoneOr(co.Point2Option(), None),
                'cp2': co.NoneOr(co.Point2Option(), None),
                'algo': co.BasicOption(str, "add"),
                }

    def _build_grid(self):
        g = self.grid2_by_name(self.get_option('grid'))
        c = self.any_cont_by_name(self.get_option('cont'))
        p1, p2 = self.get_option('gp1'), self.get_option('gp2')
        p3, p4 = self.get_option('cp1'), self.get_option('cp2')
        return g2core.snap_to_contour(
            g.cdata, c.cdata, p1, p2, p3, p4, self.get_option("algo"))
