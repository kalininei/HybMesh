import copy
import command
import objcom
from hybmeshpack import gdata, basic
from unite_grids import (unite_grids, grid_excl_cont, setbc_from_conts,
                         boundary_layer_grid)


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
                'ny': command.BasicOption(int)
                }

    def _build_grid(self):
        so = self.options
        return gdata.grid2.UnfRectGrid(so['p0'], so['p1'],
                                       so['nx'], so['ny'])


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
                }

    def _build_grid(self):
        opt = self.options
        return gdata.grid2.UnfCircGrid(
            opt['p0'], opt['rad'], opt['na'], opt['nr'],
            opt['coef'], opt['is_trian'])


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
        try:
            _, _, base = self.receiver.get_grid(name=self.options['grid_name'])
        except:
            raise command.ObjectNotFound(self.options['grid_name'])
        for c in self.options['cont_names']:
            c = self.receiver.get_any_contour(c)
            #ask parent_flow for callback
            cb = self.ask_for_callback(basic.interf.Callback.CB_CANCEL2)
            #invoke exclusion
            base = grid_excl_cont(base, c, self.options['is_inner'], cb)
            if base is None:
                return None
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
                den=command.BasicOption(int)
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
            try:
                #ask parent_flow for callback
                cb = self.ask_for_callback(basic.interf.Callback.CB_CANCEL2)
                #execute
                ret = unite_grids(ret, g, b,
                                  opt['fix_bnd'], opt['empty_holes'], cb)
            except Exception as e:
                raise command.ExecutionError('Unition Error: %s' % str(e),
                                             self, e)
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
        """
        return {'name': command.BasicOption(str),
                'simp_bnd': command.BasicOption(float),
                }

    def __simplify_bnd(self, og):
        import cobj
        import ctypes as ct
        import unite_grids
        libfa = cobj.cport_lib()
        gc = cobj.grid_to_c(og)
        ret = libfa.simplify_grid_boundary(
            gc, ct.c_double(self.options['simp_bnd']))
        if ret != 0:
            cobj.free_c_grid(gc)
            raise command.ExecutionError('Error in boundary simplification',
                                         self)
        ret = cobj.grid_from_c(gc)
        ret.build_contour()
        cobj.free_c_grid(gc)
        # write boundary types from og
        unite_grids.add_bc_from_cont(ret.cont, og.cont, force=3)
        return ret

    def _addrem_objects(self):
        try:
            _, _, orig = self.receiver.get_grid(name=self.options['name'])
        except:
            raise command.ObjectNotFound(self.options['name'])
        g = orig
        if self.options['simp_bnd'] >= 0:
            g = self.__simplify_bnd(g)

        if g != orig:
            return [(self.options['name'], g)], [self.options['name']], [], []
