import copy
import command
import objcom
from hybmeshpack import gdata, basic
import hybmeshpack.basic.geom as bgeom
from unite_grids import (setbc_from_conts, add_bc_from_cont,
                         boundary_layer_grid)
import hybmeshpack.gdata.contour2 as contour2
import math
from hybmeshpack.hmcore import g2 as g2core
from hybmeshpack.hmcore import c2 as c2core
from hybmeshpack.hmcore import libhmcport


class NewGridCommand(objcom.AbstractAddRemove):
    "Command with a new grid addition"
    def __init__(self, argsdict):
        if "name" not in argsdict:
            argsdict["name"] = "Grid1"
        super(NewGridCommand, self).__init__(argsdict)

    def _addrem_objects(self):
        #backup info
        g = self._build_grid()
        if g is not None:
            g.build_contour()
            return [(self.options['name'], g)], [], [], []
        else:
            return [], [], [], []

    #function for overloading
    def _build_grid(self):
        '-> grid2.Grid2'
        raise NotImplementedError


class AddUnfRectGrid(NewGridCommand):
    "Add uniform rectangular grid "
    def __init__(self, argsdict):
        super(AddUnfRectGrid, self).__init__(argsdict)

    @classmethod
    def _arguments_types(cls):
        return {'name': command.BasicOption(str),
                'p0': command.Point2Option(),
                'p1': command.Point2Option(),
                'nx': command.BasicOption(int),
                'ny': command.BasicOption(int),
                'custom_x': command.ListOfOptions(command.BasicOption(float)),
                'custom_y': command.ListOfOptions(command.BasicOption(float))
                }

    def _build_grid(self):
        so = self.options
        if len(so['custom_x']) == 0 and len(so['custom_y']) == 0:
            return gdata.grid2.UnfRectGrid(so['p0'], so['p1'],
                                           so['nx'], so['ny'])
        else:
            nx, ny, p0, p1 = so['nx'], so['ny'], so['p0'], so['p1']
            x, y = [], []
            if len(so['custom_x']) == 1:
                nx = int((p1.x - p0.x) / so['custom_x'][0]) + 1
            else:
                x = so['custom_x']
                nx = len(x) - 1
            if len(so['custom_y']) == 1:
                ny = int((p1.y - p0.y) / so['custom_x'][0]) + 1
            else:
                y = so['custom_y']
                ny = len(y) - 1
            ret = gdata.grid2.UnfRectGrid(p0, p1, nx, ny)

            if len(x) > 0:
                for j in range(ny + 1):
                    for i in range(nx + 1):
                        p = ret.points[j * (nx + 1) + i]
                        p.x = x[i]
            if len(y) > 0:
                for j in range(ny + 1):
                    for i in range(nx + 1):
                        p = ret.points[j * (nx + 1) + i]
                        p.y = y[i]
            return ret


class AddCustomRectGrid(NewGridCommand):
    "Add custom rectangular grid "
    def __init__(self, argsdict):
        super(AddCustomRectGrid, self).__init__(argsdict)

    @classmethod
    def _arguments_types(cls):
        return {'name': command.BasicOption(str),
                'algo': command.BasicOption(str),
                'left': command.BasicOption(str),
                'bot': command.BasicOption(str),
                'right': command.NoneOr(command.BasicOption(str)),
                'top': command.NoneOr(command.BasicOption(str)),
                'her_w': command.ListOfOptions(command.BasicOption(float))
                }

    def _build_grid(self):
        from unite_grids import add_bc_from_cont
        so = self.options
        # callback
        cb = self.ask_for_callback()
        # declare all c variables as 0 to invoke memory free at exit
        c_left, c_bot, c_right, c_top, c_ret, c_retcont =\
            0, 0, 0, 0, 0, 0
        try:
            # get contours
            left = self.any_cont_by_name(so['left'])
            bot = self.any_cont_by_name(so['bot'])
            # If no right/top use copy of left/bottom.
            # In custom_rectangular_grid procedure their points coordinates
            # will be properly changed at c-side.
            right = left if so['right'] is None\
                else self.any_cont_by_name(so['right'])
            top = bot if so['top'] is None\
                else self.any_cont_by_name(so['top'])
            # copy'em to c side
            c_left = c2core.cont2_to_c(left)
            c_bot = c2core.cont2_to_c(bot)
            c_right = c2core.cont2_to_c(right)
            c_top = c2core.cont2_to_c(top)
            # specific options
            c_herw = g2core.list_to_c(so['her_w'], float)

            # call c procedure
            c_ret = g2core.custom_rectangular_grid(
                so['algo'], c_left, c_bot, c_right, c_top, c_herw, cb)

            # build grid
            ret = g2core.grid_from_c(c_ret)
            ret.build_contour()
            retcont = ret.cont
            c_retcont = c2core.cont2_to_c(retcont)

            # boundary types
            # use c_* contours as sources since they were changed
            # in custom_rectangular_grid procedure to fit resulting grid
            # boundary
            add_bc_from_cont(retcont, left, c_retcont, c_left, force=2)
            add_bc_from_cont(retcont, bot, c_retcont, c_bot, force=2)
            add_bc_from_cont(retcont, right, c_retcont, c_right, force=2)
            add_bc_from_cont(retcont, top, c_retcont, c_top, force=2)

            # return
            return ret
        except Exception as e:
            raise command.ExecutionError(str(e), self, e)
        finally:
            for c in [c_left, c_bot, c_right, c_top, c_retcont]:
                if c != 0:
                    c2core.free_cont2(c)
            if c_ret != 0:
                g2core.free_c_grid(c_ret)


class AddCirc4Grid(NewGridCommand):
    "Add circular quadrangular grid "
    def __init__(self, argsdict):
        super(AddCirc4Grid, self).__init__(argsdict)

    @classmethod
    def _arguments_types(cls):
        return {'name': command.BasicOption(str),
                'algo': command.BasicOption(str),
                'p0': command.Point2Option(),
                'rad': command.BasicOption(float),
                'step': command.BasicOption(float),
                'sqrside': command.BasicOption(float),
                'rcoef': command.BasicOption(float)
                }

    def _build_grid(self):
        so = self.options
        c_ret = 0
        try:
            c_p0 = g2core.list_to_c([so['p0'].x, so['p0'].y], float)
            c_ret = g2core.circ4grid(
                so['algo'], c_p0, so['rad'], so['step'],
                so['sqrside'], so['rcoef'])
            return g2core.grid_from_c(c_ret)
        except Exception as e:
            raise command.ExecutionError(str(e), self, e)
        finally:
            if c_ret != 0:
                g2core.free_c_grid(c_ret)


class AddUnfCircGrid(NewGridCommand):
    "Add uniform circular grid "

    def __init__(self, argsdict):
        opt = copy.deepcopy(argsdict)
        if 'coef' not in opt:
            opt['coef'] = 1.0
        if 'is_trian' not in opt:
            opt['is_trian'] = True
        super(AddUnfCircGrid, self).__init__(argsdict)

    @classmethod
    def _arguments_types(cls):
        return {'name': command.BasicOption(str),
                'p0': command.Point2Option(),
                'rad': command.BasicOption(float),
                'na': command.BasicOption(int),
                'nr': command.BasicOption(int),
                'coef': command.BasicOption(float),
                'is_trian': command.BoolOption(),
                'custom_r': command.ListOfOptions(command.BasicOption(float)),
                'custom_a': command.ListOfOptions(command.BasicOption(float)),
                }

    def _build_grid(self):
        opt = self.options
        if len(opt['custom_r']) == 0 and len(opt['custom_a']) == 0:
            return gdata.grid2.UnfCircGrid(
                opt['p0'], opt['rad'], opt['na'], opt['nr'],
                opt['coef'], opt['is_trian'])
        else:
            # custom data
            rads, angles = [], []
            radius, na, nr = opt['rad'], opt['na'], opt['nr']

            # calculate rads
            if len(opt['custom_r']) == 1:
                nr = int(radius / opt['custom_r'][0]) + 1
            elif len(opt['custom_r']) > 1:
                rads = copy.deepcopy(opt['custom_r'])
                rads = rads[1:] if rads[0] == 0 else rads
                radius = rads[-1]
                nr = len(rads)

            # calculate angles
            if len(opt['custom_a']) == 1:
                arclen = 2 * math.pi * radius
                na = max(3, int(arclen / opt['custom_a'][0]))
            elif len(opt['custom_a']) > 1:
                angles = copy.deepcopy(opt['custom_a'])
                for i in range(len(angles)):
                    angles[i] = angles[i] * 2 * math.pi / angles[-1]
                angles.insert(-1, 0) if angles[0] != 0 else None
                angles = angles[:-1]
                na = len(angles)

            # build grid
            r = gdata.grid2.UnfCircGrid(
                bgeom.Point2(0, 0), opt['rad'], na, nr,
                opt['coef'], opt['is_trian'])

            # modify arch and radii
            for i in range(na):
                if len(angles) > 0:
                    sina, cosa = math.sin(angles[i]), math.cos(angles[i])
                else:
                    sina = r.points[i].y / opt['rad']
                    cosa = r.points[i].x / opt['rad']

                for j in range(nr):
                    p = r.points[na * j + i]
                    prad = rads[-1 - j] if len(rads) > 0 else p.dist0()
                    p.x, p.y = prad * cosa, prad * sina

            # move to pc
            r.move(opt['p0'].x, opt['p0'].y)
            return r


class AddUnfHexGrid(NewGridCommand):
    "Add uniform circular ring"

    def __init__(self, argsdict):
        super(AddUnfHexGrid, self).__init__(argsdict)

    @classmethod
    def _arguments_types(cls):
        return {'name': command.BasicOption(str),
                'area': command.ListOfOptions(command.BasicOption(float)),
                'crad': command.BasicOption(float),
                'strict': command.BoolOption(),
                }

    def _build_grid(self):
        c_ret = 0
        try:
            c_ret = g2core.regular_hex_grid(self.options['area'],
                                            self.options['crad'],
                                            self.options['strict'])
            ret = g2core.grid_from_c(c_ret)
            ret.build_contour()
            return ret
        except Exception as e:
            raise command.ExecutionError('AddUnfHexGrid: %s' % str(e), self, e)
        finally:
            g2core.free_c_grid(c_ret) if c_ret != 0 else None


class AddUnfRingGrid(NewGridCommand):
    "Add uniform circular ring"

    def __init__(self, argsdict):
        opt = copy.deepcopy(argsdict)
        if 'coef' not in opt:
            opt['coef'] = 1.0
        super(AddUnfRingGrid, self).__init__(argsdict)

    @classmethod
    def _arguments_types(cls):
        return {'name': command.BasicOption(str),
                'p0': command.Point2Option(),
                'radinner': command.BasicOption(float),
                'radouter': command.BasicOption(float),
                'na': command.BasicOption(int),
                'nr': command.BasicOption(int),
                'coef': command.BasicOption(float),
                }

    def _build_grid(self):
        opt = self.options
        return gdata.grid2.UnfRingGrid(
            opt['p0'], opt['radinner'], opt['radouter'],
            opt['na'], opt['nr'], opt['coef'])


class AddTriGrid(NewGridCommand):
    "Add grid in triangle area"
    def __init__(self, argsdict):
        super(AddTriGrid, self).__init__(argsdict)

    @classmethod
    def _arguments_types(cls):
        return {'name': command.BasicOption(str),
                'vertices': command.ListOfOptions(command.Point2Option()),
                'nedge': command.BasicOption(int)
                }

    def _build_grid(self):
        so = self.options
        ne = so['nedge']
        # check points order
        [p0, p1, p2] = map(copy.deepcopy, so['vertices'])
        x0, y0 = p0.x - p2.x, p0.y - p2.y
        x1, y1 = p1.x - p2.x, p1.y - p2.y
        if (x0 * y1 - y0 * x1) > 0:
            p0, p2 = p2, p0
        # build points
        plines = []
        for i in range(ne + 1):
            w1 = float(i) / ne
            x1, y1 = w1 * p2.x + (1 - w1) * p0.x, w1 * p2.y + (1 - w1) * p0.y
            x2, y2 = w1 * p2.x + (1 - w1) * p1.x, w1 * p2.y + (1 - w1) * p1.y
            plines.append([])
            for j in range(ne + 1 - i):
                w2 = float(j) / (ne - i) if ne != i else 1
                x3, y3 = (1 - w2) * x1 + w2 * x2, (1 - w2) * y1 + w2 * y2
                plines[-1].append(bgeom.Point2(x3, y3))

        allpoints = []
        for p in plines:
            allpoints.extend(p)
        # build cells
        pind = {}
        for i, p in enumerate(allpoints):
            pind[p] = i
        allcells = []
        for i in range(ne):
            p1 = pind[plines[i][0]]
            p2 = pind[plines[i + 1][0]]
            p3 = pind[plines[i][1]]
            allcells.append([p1, p2, p3])
            for j in range(ne - i - 1):
                p1 = pind[plines[i][j + 1]]
                p2 = pind[plines[i + 1][j]]
                p3 = pind[plines[i + 1][j + 1]]
                p4 = pind[plines[i][j + 2]]
                allcells.append([p1, p2, p3, p4])
        return gdata.grid2.Grid2.from_points_cells2(allpoints, allcells)


class ExcludeContours(objcom.AbstractAddRemove):
    "Exclude contour area from grids"

    def __init__(self, arg):
        if 'name' not in arg:
            arg['name'] = "Grid1"
        if 'keep_src' not in arg:
            arg['keep_src'] = True
        super(ExcludeContours, self).__init__(arg)

    @classmethod
    def _arguments_types(cls):
        """ name - new grid name,
            grid_name - source grid name,
            cont_names - contours names,
            is_inner - exclude inner (True) or outer (False) contour area
            keep_src - preserve source grid
        """
        return {'name': command.BasicOption(str),
                'grid_name': command.BasicOption(str),
                'cont_names': command.ListOfOptions(command.BasicOption(str)),
                'is_inner': command.BoolOption(),
                'keep_src': command.BoolOption()
                }

    #overriden from NewGridCommand keeping self.remove_com
    def __build_grid(self):
        base = self.grid_by_name(self.options['grid_name'])
        for cname in self.options['cont_names']:
            cb = self.ask_for_callback(basic.interf.Callback.CB_CANCEL2)
            cont = self.any_cont_by_name(cname)
            c_cont, c_grid, c_res = 0, 0, 0
            try:
                # copy to c-side
                c_grid = g2core.grid_to_c(base)
                c_cont = c2core.cont2_to_c(cont)

                # exclude
                c_res = g2core.grid_excl_cont(c_grid, c_cont,
                                              self.options['is_inner'], cb)

                # copy to py-side
                res = g2core.grid_from_c(c_res)
                res.build_contour()

                # boundary types
                bs = cont.bnd_types().union(base.cont.bnd_types())
                bs = bs.difference(set([0]))
                if len(bs) > 0:
                    res_cont = c2core.cont2_to_c(res.cont)
                    #1. from source grid
                    add_bc_from_cont(res.cont, base.cont, c_tar=res_cont)
                    #2. from contour
                    add_bc_from_cont(res.cont, cont, c_tar=res_cont,
                                     c_src=c_cont)
                    #3. free contour memory
                    c2core.free_cont2(res_cont)

                # assign to base
                base = res
            except Exception as e:
                raise command.ExecutionError('Contour exclusion error: %s' %
                                             str(e), self, e)
            finally:
                c2core.free_cont2(c_cont) if c_cont != 0 else None
                g2core.free_c_grid(c_grid) if c_grid != 0 else None
                g2core.free_c_grid(c_res) if c_res != 0 else None

        return base

    def _addrem_objects(self):
        g = self.__build_grid()
        if g is None:
            return [], [], [], []
        remg = [] if self.options['keep_src'] else [self.options['grid_name']]
        return [(self.options['name'], g)], remg, [], []


class UniteGrids(objcom.AbstractAddRemove):
    def __init__(self, argsdict):
        if 'name' not in argsdict:
            argsdict['name'] = 'Grid1'
        a = copy.deepcopy(argsdict)
        if 'fix_bnd' not in a:
            a['fix_bnd'] = False
        if 'keep_src' not in a:
            a['keep_src'] = True
        if 'empty_holes' not in a:
            a['empty_holes'] = False
        super(UniteGrids, self).__init__(a)

    def doc(self):
        ret = "Unite Grids: %s " % self.options['base']
        for s in self.options['plus']:
            ret += s['name'] + " "
        return ret[:-1]

    class Option(command.SubDictOption):
        'name, buf, den - dictionary'
        def __init__(self):
            super(UniteGrids.Option, self).__init__(
                name=command.BasicOption(str),
                buf=command.BasicOption(float),
                den=command.BasicOption(int),
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
        return {'name': command.BasicOption(str),
                'fix_bnd': command.BoolOption(),
                'empty_holes': command.BoolOption(),
                'keep_src': command.BoolOption(),
                'base': command.BasicOption(str),
                'plus': command.ListOfOptions(UniteGrids.Option()),
                'angle0': command.BasicOption(float),
                'filler': command.BasicOption(str)
                }

    def _get_grid(self, ind):
        """ returns (Grid, density, buffer size) by order of addition """
        if ind == 0:
            gname = self.options['base']
            b = 0
            d = 0
        else:
            gname = self.options['plus'][ind - 1]['name']
            b = self.options['plus'][ind - 1]['buf']
            d = self.options['plus'][ind - 1]['den']
        if gname not in self.receiver.get_grid_names():
            raise command.ExecutionError('Grid was not found: %s' % gname,
                                         self)
        _, _, g = self.receiver.get_grid(name=gname)
        return g, b, d

    def _addrem_objects(self):
        opt = self.options
        #basic grid
        ret, _, _ = self._get_grid(0)

        #unification
        for i in range(len(opt['plus'])):
            g, b, _ = self._get_grid(i + 1)
            #ask parent_flow for callback
            cb = self.ask_for_callback(basic.interf.Callback.CB_CANCEL2)
            c_g1, c_g2, c_ret = 0, 0, 0
            try:
                # copy grids on c-side
                c_g1 = g2core.grid_to_c(ret)
                c_g2 = g2core.grid_to_c(g)
                #execute
                c_ret = g2core.unite_grids(
                    c_g1, c_g2, b, opt['fix_bnd'], opt['empty_holes'],
                    opt['angle0'], opt['filler'], cb)
                ret = g2core.grid_from_c(c_ret)
            except Exception as e:
                raise command.ExecutionError('Unition Error: %s' % str(e),
                                             self, e)
            finally:
                g2core.free_c_grid(c_g1) if c_g1 != 0 else None
                g2core.free_c_grid(c_g2) if c_g2 != 0 else None
                g2core.free_c_grid(c_ret) if c_ret != 0 else None
            if (ret is None):
                return [], [], [], []

        #boundary conditions
        ret.build_contour()
        srccont = []
        for i in range(len(opt['plus']) + 1):
            srccont.append(self._get_grid(i)[0].cont)
        setbc_from_conts(ret.cont, srccont)

        #delete source
        if opt['keep_src']:
            delgrd = []
        else:
            delgrd = [opt['base']]
            for i in range(len(opt['plus'])):
                delgrd.append(opt['plus'][i]['name'])
        return [(opt['name'], ret)], delgrd, [], []


class BuildBoundaryGrid(NewGridCommand):
    def __init__(self, kwargs):
        if "name" not in kwargs:
            kwargs["name"] = "Grid1"
        for op in kwargs['opt']:
            if 'step_start' not in op:
                op['step_start'] = op['mesh_cont_step']
            if 'step_end' not in op:
                op['step_end'] = op['mesh_cont_step']
            if 'start' not in op or 'end' not in op:
                op['start'] = None
                op['end'] = None
        super(BuildBoundaryGrid, self).__init__(kwargs)

    class Option(command.SubDictOption):
        def __init__(self):
            super(BuildBoundaryGrid.Option, self).__init__(
                source=command.BasicOption(str),
                partition=command.ListOfOptions(command.BasicOption(float)),
                direction=command.BasicOption(int),
                mesh_cont=command.BasicOption(int),
                mesh_cont_step=command.BasicOption(float),
                step_start=command.BasicOption(float),
                step_end=command.BasicOption(float),
                algo_acute=command.BasicOption(float),
                algo_right=command.BasicOption(float),
                algo_straight=command.BasicOption(float),
                algo_reentr=command.BasicOption(float),
                start=command.NoneOr(command.Point2Option()),
                end=command.NoneOr(command.Point2Option()),
                force_conf=command.BoolOption()
            )

    @classmethod
    def _arguments_types(cls):
        """ arguments
                name: str  -- name of the new grid
                opt: [ {
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
                }, ...]
        """
        return {
            'name': command.BasicOption(str),
            'opt': command.ListOfOptions(cls.Option())
        }

    def __find_open_subcont(self, cont, pts):
        """ -> (p0, p1) or None.
        finds closest to pts subcontour in cont.
        if it is open -> returns its first and last points
        else returns None
        """
        # edges points
        ep = cont.edges_points()
        # array of subcontours
        edss = cont.sorted_edges()
        # subcontour closest to pts
        eds = None
        if len(edss) == 1:
            eds = edss[0]
        elif len(edss) == 0:
            return None
        else:
            pi = cont.closest_point_index(pts)
            for lst in edss:
                for e in lst:
                    if pi in ep[e][:2]:
                        eds = lst
                        break
                if eds is not None:
                    break

        # first/last point of subcontour
        if len(eds) > 1:
            e0, e1, em1, em2 = [ep[eds[i]][:2] for i in [0, 1, -1, -2]]
            p0 = e0[0] if e0[1] in e1 else e0[1]
            p1 = em1[1] if em1[0] in em2 else em1[1]
        elif len(eds) == 1:
            [p0, p1] = eds[0][:2]
        else:
            return None

        # return if subcontour is open
        return None if p0 == p1 else (cont.points[p0], cont.points[p1])

    def __adjust_arg(self, arg):
        """ changes 'source' from name to AbstractContour object
            if end points are None changes them to actual values
            and adds options for multiply connected domains
        """
        for op in arg:
            op['source'] = self.any_cont_by_name(op['source'])

        # if end points are None:
        #    transform them to real contour point
        #    create set of options for multiply connected domains
        newops = []
        for op in arg:
            if op['start'] is None or op['end'] is None:
                eds = op['source'].sorted_edges()
                epoints = op['source'].edges_points()
                p0 = op['source'].points[epoints[eds[0][0]][0]]
                op['start'] = op['end'] = copy.deepcopy(p0)
                for e in eds[1:]:
                    newop = copy.deepcopy(op)
                    pi = op['source'].points[epoints[e[0]][0]]
                    newop['start'] = newop['end'] = copy.deepcopy(pi)
                    newops.append(newop)
        arg.extend(newops)

        # if contour is open and points coincide
        # write first/last point of the contour
        for op in arg:
            if op['start'].x == op['end'].x and op['start'].y == op['end'].y:
                subcont = self.__find_open_subcont(op['source'], op['start'])
                if subcont is not None:
                    op['start'] = subcont[0]
                    op['end'] = subcont[-1]

    def _build_grid(self):
        '-> grid2.Grid2'
        arg = copy.deepcopy(self.options['opt'])
        self.__adjust_arg(arg)

        try:
            cb = self.ask_for_callback(basic.interf.Callback.CB_CANCEL2)
            ret = boundary_layer_grid(arg, cb)
        except Exception as e:
            raise command.ExecutionError(str(e), self, e)

        #boundary conditions
        if ret is not None:
            ret.build_contour()
            srccont = [op['source'] for op in arg]
            setbc_from_conts(ret.cont, srccont, force=2)

        return ret


class HealGrid(objcom.AbstractAddRemove):
    "simplification of grid boundaries"

    def __init__(self, arg):
        if 'simp_bnd' not in arg:
            arg['simp_bnd'] = -1
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
        return {'name': command.BasicOption(str),
                'simp_bnd': command.BasicOption(float),
                'convex': command.BasicOption(float)
                }

    def __simplify_bnd(self, og):
        import ctypes as ct
        gc = g2core.grid_to_c(og)
        ret = libhmcport.simplify_grid_boundary(
            gc, ct.c_double(self.options['simp_bnd']))
        if ret != 0:
            g2core.free_c_grid(gc)
            raise command.ExecutionError('Error in boundary simplification',
                                         self)
        ret = g2core.grid_from_c(gc)
        ret.build_contour()
        g2core.free_c_grid(gc)
        # write boundary types from og
        add_bc_from_cont(ret.cont, og.cont, force=3)
        return ret

    def __convex_cells(self, og):
        c_grid, c_ret = 0, 0
        try:
            c_grid = g2core.grid_to_c(og)
            c_ret = g2core.convex_cells(c_grid, self.options['convex'])
            ret = g2core.grid_from_c(c_ret)
            ret.build_contour()
            add_bc_from_cont(ret.cont, og.cont, force=3)
        except Exception as e:
            raise command.ExecutionError(str(e), self, e)
        finally:
            g2core.free_c_grid(c_ret) if c_ret != 0 else None
            g2core.free_c_grid(c_grid) if c_grid != 0 else None
        return ret

    def _addrem_objects(self):
        try:
            _, _, orig = self.receiver.get_grid(name=self.options['name'])
        except:
            raise command.ObjectNotFound(self.options['name'])
        g = orig
        if self.options['simp_bnd'] >= 0:
            g = self.__simplify_bnd(g)
        if self.options['convex'] >= 0:
            g = self.__convex_cells(g)

        if g != orig:
            return [(self.options['name'], g)], [self.options['name']], [], []


class MapGrid(NewGridCommand):
    "Maps grid area on given contour"

    def __init__(self, kwargs):
        if "name" not in kwargs:
            kwargs["name"] = "Grid1"
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
        return {'name': command.BasicOption(str),
                'base': command.BasicOption(str),
                'target': command.BasicOption(str),
                'base_points': command.ListOfOptions(command.Point2Option()),
                'target_points': command.ListOfOptions(command.Point2Option()),
                'snap': command.BasicOption(str),
                'btypes': command.BasicOption(str),
                'algo': command.BasicOption(str),
                'return_invalid': command.BoolOption(),
                }

    def _treat_boundaries(self, grid, method):
        grid.build_contour()
        # treat boundaries
        if method == "from_contour":
            c = self.any_cont_by_name(self.options['target'])
            add_bc_from_cont(grid.cont, c, force=3)
        elif method == "from_grid":
            g = self.grid_by_name(self.options['base'])
            # 1. assemble boundary edges: edge index -> bnd feature
            ret_full_bnd = grid.boundary_contours()
            ret_bnd_edges = []
            for ecol in ret_full_bnd:
                ret_bnd_edges.extend(ecol)
            grid_bnd_edges = []
            for ecol in g.boundary_contours():
                grid_bnd_edges.extend(ecol)
            # 2. whether boundary edge has its sibling in grid
            is_old_edge = []
            oldptsnum = g.n_points()
            for ei in ret_bnd_edges:
                e = grid.edges[ei]
                is_old_edge.append(e[0] < oldptsnum and e[1] < oldptsnum)

            # 3. build dictionary for boundary edges in grid
            def etostr(e):
                return str(min(e[0], e[1])) + '-' + str(max(e[0], e[1]))
            old_ebt = {}
            for ei in grid_bnd_edges:
                bt = 0
                if ei in g.bt:
                    bt = g.bt[ei]
                old_ebt[etostr(g.edges[ei])] = bt

            # 4. restore bfeatures for old edges
            for i, ei in enumerate(ret_bnd_edges):
                if is_old_edge[i]:
                    bt = old_ebt[etostr(grid.edges[ei])]
                    if bt > 0:
                        grid.bt[ei] = bt

            # 5. boundary types for new edges
            for c in ret_full_bnd:
                # find any old edge
                for istart, ei in enumerate(c):
                    if is_old_edge[ret_bnd_edges.index(ei)]:
                        break
                else:
                    # no old edges in current contour
                    # this should not happen
                    continue
                # starting from old edge fill all new ones
                bcur = grid.get_edge_bnd(c[istart])
                for i in range(len(c)):
                    icur = i + istart
                    if icur >= len(c):
                        icur -= len(c)
                    ecur = c[icur]
                    if is_old_edge[ret_bnd_edges.index(ecur)]:
                        bcur = grid.get_edge_bnd(c[icur])
                    else:
                        if bcur > 0:
                            grid.bt[c[icur]] = bcur
        else:
            raise ValueError("Unknown btypes = " + str(bt))

    def _build_grid(self):
        '-> grid2.Grid2'
        g = self.grid_by_name(self.options['base'])
        c = self.any_cont_by_name(self.options['target'])
        cb = self.ask_for_callback()

        c_grid, c_cont, c_ret = 0, 0, 0
        try:
            # build c grid, contour
            c_grid = g2core.grid_to_c(g)
            c_cont = c2core.cont2_to_c(c)

            # build points array
            bp, tp = self.options['base_points'], self.options['target_points']
            p1, p2 = [], []
            for i in range(min(len(bp), len(tp))):
                p1.append(bp[i].x)
                p1.append(bp[i].y)
                p2.append(tp[i].x)
                p2.append(tp[i].y)
            p1 = g2core.list_to_c(p1, float)
            p2 = g2core.list_to_c(p2, float)

            # mapping procedure
            c_ret = g2core.map_grid(c_grid, c_cont, p1, p2,
                                    self.options['snap'],
                                    self.options['algo'],
                                    self.options['return_invalid'], cb)

            # copy from c
            ret = g2core.grid_from_c(c_ret)
            self._treat_boundaries(ret, self.options['btypes'])
            return ret
        except Exception as e:
            raise command.ExecutionError('Mapping error: %s' % str(e),
                                         self, e)
        finally:
            g2core.free_c_grid(c_grid) if c_grid != 0 else None
            c2core.free_cont2(c_cont) if c_cont != 0 else None
            g2core.free_c_grid(c_ret) if c_ret != 0 else None


class UnstructuredFillArea(NewGridCommand):
    def __init__(self, argsdict):
        super(UnstructuredFillArea, self).__init__(argsdict)

    @classmethod
    def _arguments_types(cls):
        return {'name': command.BasicOption(str),
                'domain': command.ListOfOptions(command.BasicOption(str)),
                'constr': command.ListOfOptions(command.BasicOption(str)),
                'pts': command.ListOfOptions(command.BasicOption(float)),
                }

    def call_c_proc(self, c_dom, c_con, c_pts):
        raise NotImplementedError

    def _build_grid(self):
        dom = self.options['domain']
        con = self.options['constr']
        pts = self.options['pts']
        c_dom, c_con, c_ret = 0, 0, 0
        try:
            #unite domain
            c1 = self.any_cont_by_name(dom[0])
            fulldom = contour2.Contour2.create_from_abstract(c1)
            for c in dom[1:]:
                c1 = self.any_cont_by_name(c)
                fulldom.add_from_abstract(c1)

            #unite constraints
            if len(con) == 0:
                fullcon = None
            else:
                c1 = self.any_cont_by_name(con[0])
                fullcon = contour2.Contour2.create_from_abstract(c1)
                for c in con[1:]:
                    c1 = self.any_cont_by_name(c)
                    fullcon.add_from_abstract(c1)

            c_dom = c2core.cont2_to_c(fulldom)
            if fullcon is not None:
                c_con = c2core.cont2_to_c(fullcon)
            else:
                c_con = 0

            c_pts = g2core.list_to_c(pts, float) if pts is not None else 0
            c_ret = self.call_c_proc(c_dom, c_con, c_pts)
            ret = g2core.grid_from_c(c_ret)

            #boundary conditions
            ret.build_contour()
            add_bc_from_cont(ret.cont, fulldom, c_src=c_dom, force=3)
            return ret
        except Exception as e:
            raise command.ExecutionError(str(e), self, e)
        finally:
            g2core.free_c_grid(c_ret) if c_ret != 0 else None
            c2core.free_cont2(c_dom) if c_dom != 0 else None
            c2core.free_cont2(c_con) if c_con != 0 else None


class TriangulateArea(UnstructuredFillArea):
    def __init__(self, argsdict):
        super(TriangulateArea, self).__init__(argsdict)

    def call_c_proc(self, c_dom, c_con, c_pts):
        return g2core.unstructed_fill(c_dom, c_con, c_pts, '3')


class QuadrangulateArea(UnstructuredFillArea):
    def __init__(self, argsdict):
        super(QuadrangulateArea, self).__init__(argsdict)

    def call_c_proc(self, c_dom, c_con, c_pts):
        return g2core.unstructed_fill(c_dom, c_con, c_pts, '4')


class PebiFill(UnstructuredFillArea):
    def __init__(self, argsdict):
        super(PebiFill, self).__init__(argsdict)

    def call_c_proc(self, c_dom, c_con, c_pts):
        return g2core.unstructed_fill(c_dom, c_con, c_pts, 'pebi')


class StripeGrid(NewGridCommand):
    def __init__(self, kwargs):
        if "name" not in kwargs:
            kwargs["name"] = "Grid1"
        super(StripeGrid, self).__init__(kwargs)

    @classmethod
    def _arguments_types(cls):
        return {'name': command.BasicOption(str),
                'source': command.BasicOption(str),
                'target': command.BasicOption(str),
                'partition': command.ListOfOptions(command.BasicOption(float)),
                'tip': command.BasicOption(str),
                'bnd': command.ListOfOptions(command.BasicOption(int))
                }

    def _build_grid(self):
        '-> grid2.Grid2'
        so = self.options
        c = self.any_cont_by_name(so['source'])
        cb = self.ask_for_callback()

        c_ret, c_cont, c_cout = 0, 0, [0, 0, 0, 0]
        try:
            # build c grid, contour
            c_cont = c2core.cont2_to_c(c)

            c_p = g2core.list_to_c(self.options['partition'], float)

            # mapping procedure
            c_ret, c_cout = g2core.stripe_grid(c_cont, c_p, so['tip'], cb)

            # copy from c
            ret = g2core.grid_from_c(c_ret)

            #boundary conditions
            def bc_from_cont(c_cnt, b):
                c = c2core.cont2_from_c(c_cnt)
                c.bnds = {i: b for i in range(c.n_edges())}
                add_bc_from_cont(ret.cont, c, c_src=c_cnt, force=1)

            ret.build_contour()
            for i in range(4):
                bc_from_cont(c_cout[i], so['bnd'][i])

            return ret
        except Exception as e:
            raise command.ExecutionError('StripGrid: %s' % str(e), self, e)
        finally:
            c2core.free_cont2(c_cont) if c_cont != 0 else None
            g2core.free_c_grid(c_ret) if c_ret != 0 else None
            for i in range(4):
                c2core.free_cont2(c_cout[i]) if c_cout[i] != 0 else None
