import copy
import command
import objcom
import gdata.grid2
from basic.geom import Point2
from unite_grids import (unite_grids, grid_excl_cont, setbc_from_conts,
    boundary_layer_grid)
import basic.proc as bp


class NewGridCommand(objcom.AbstractAddRemove):
    "Command with a new grid addition"
    def __init__(self, argsdict):
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
        return gdata.grid2.UnfCircGrid(opt['p0'], opt['rad'],
                opt['na'], opt['nr'], opt['coef'], opt['is_trian'])


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
        return gdata.grid2.UnfRingGrid(opt['p0'], opt['radinner'],
                opt['radouter'], opt['na'], opt['nr'], opt['coef'])


class ExcludeContours(objcom.AbstractAddRemove):
    "Exclude contour area from grids"

    def __init__(self, arg):
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
        import basic.interf
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
        import basic.interf
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
    CONT_MESH_NO = 0
    CONT_MESH_KEEP_ORIGIN = 1
    CONT_MESH_KEEP_SHAPE = 2
    CONT_MESH_IGNORE_ORIGIN = 3

    def __init__(self, kwargs):
        for op in kwargs['opt']:
            if 'start' not in op or 'end' not in op:
                op['start'] = op['end'] = Point2(0, 0)
            bp.set_if_no(op, 'round_off', False)
            bp.set_if_no(op, 'maxsharp', 120)
            bp.set_if_no(op, 'mesh_cont', 0)
            bp.set_if_no(op, 'mesh_cont_step', 0.1)
        super(BuildBoundaryGrid, self).__init__(kwargs)

    class Option(command.SubDictOption):
        def __init__(self):
            super(BuildBoundaryGrid.Option, self).__init__(
                source=command.BasicOption(str),
                tp=command.BasicOption(str),
                start=command.Point2Option(),
                end=command.Point2Option(),
                round_off=command.BoolOption(),
                maxsharp=command.BasicOption(float),
                mesh_cont=command.BasicOption(int),
                mesh_cont_step=command.BasicOption(float),
                partition=command.ListOfOptions(command.BasicOption(float))
            )

    @classmethod
    def _arguments_types(cls):
        """ arguments
                name: str  -- name of the new grid
                opt: [ {
                    source: str -- name of source contour
                    tp: str    -- inside/outside/left/right/around boundary
                    start: Point   -- Start Point of Section.
                    end: Point     -- End Point of Section.
                    round_off: bool -- round off sharp corners
                    maxsharp: 20.0-160.0 -- maximum angle which is sharp (deg)
                    mesh_cont: 0/1/2/3 -- contour meshing strategy:
                        0 - No meshing, use original contour nodes
                        1 - keep all original contour nodes,
                        2 - keep only non-collinear contour nodes,
                        3 - ignore all existing contour nodes
                    mesh_cont_step: double -- contour mesh step if mesh_cont!=0
                    partition: [0.0,...] -- layer partition starts from 0
                }, ...]
        """
        return {
                'name': command.BasicOption(str),
                'opt': command.ListOfOptions(cls.Option())
        }

    #function for overloading
    def _build_grid(self):
        '-> grid2.Grid2'
        import basic.interf
        arg = copy.deepcopy(self.options['opt'])
        for op in arg:
            op['source'] = self.any_cont_by_name(op['source'])
        cb = self.ask_for_callback(basic.interf.Callback.CB_CANCEL2)
        return boundary_layer_grid(arg, cb)
