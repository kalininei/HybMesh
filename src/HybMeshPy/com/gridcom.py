import ast
import copy
import command
import objcom
import gdata.grid2
from unite_grids import (unite_grids, grid_excl_cont, setbc_from_conts)


class NewGridCommand(objcom.AbstractAddRemove):
    "Command with a new grid addition"
    def __init__(self, argsdict):
        super(NewGridCommand, self).__init__(argsdict)
        try:
            self.name = argsdict["name"]
        except KeyError:
            self.name = ""
        if self.name == "":
            self.name = "NewGrid1"

    def _addrem_objects(self):
        #backup info
        g = self._build_grid()
        g.build_contour()
        return [(self.name, g)], [], [], []

    #function for overriding
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
        try:
            _, _, base = self.receiver.get_grid(name=self.options['grid_name'])
        except:
            raise command.ObjectNotFound(self.options['grid_name'])
        for c in self.options['cont_names']:
            c = self.receiver.get_any_contour(c)
            #ask parent_flow for callback
            try:
                intf = self.parent_flow.get_interface()
                cb = intf.ok_cancel_2pg_callback()
            except:
                cb = None
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


class UniteOption(command.BasicOption):
    'name, buf, den - dictionary'
    def __init__(self):
        super(UniteOption, self).__init__()

    def serial(self, v):
        'before writing'
        return v

    def unserial(self, s):
        'after reading'
        ret = s
        ret['buf'] = float(ret['buf'])
        ret['den'] = float(ret['den'])
        return ret


class UniteGrids(objcom.AbstractAddRemove):
    def __init__(self, argsdict):
        """ args[name] - name of the new grid,
            args[fix_bnd] - whether to fix all boundary points (=False)
            args[keepsrc] - whether to remove source grids (=True)
            args[empty_holes] - keep all empty zone in 2nd grid (=False)
            args[s0] = UniteOpts,   <- list of grids. Keys start with 's'.
            args[s1] = UniteOpts,      No other args keys can start with 's'
            ....
        """
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
                'plus': command.ListOfOptions(UniteOption()),
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
                try:
                    intf = self.parent_flow.get_interface()
                    cb = intf.ok_cancel_2pg_callback()
                except:
                    cb = None
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


class BuildBoundaryGrid(objcom.AbstractAddRemove):
    CONT_MESH_NO = 0
    CONT_MESH_KEEP_ORIGIN = 1
    CONT_MESH_KEEP_SHAPE = 2
    CONT_MESH_IGNORE_ORIGIN = 3

    def __init__(self, **kwargs):
        """ arguments
                name: str  -- name of the new grid
                opt: [ {
                    source: st -- name of source contour
                    tp: str    -- inside/outside/left/right/around boundary
                    start: double   -- normalized to [0, 1] with respect to
                        source length coordinate of section start
                    end: double     -- normalized to [0, 1] with respect to
                        source length coordinate of section end
                    round_off: 0/1 -- round off sharp corners
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
        super(BuildBoundaryGrid, self).__init__(kwargs)
        #temporary. This should be done in Command class consructor
        self._Command__comLine = self._method_code() + str(kwargs)
        self.opt = copy.deepcopy(kwargs)

    @classmethod
    def fromstring(cls, slist):
        #temporary. This should be done in Command class
        args = ast.literal_eval(slist)
        return cls(**args)

    def _addrem_objects(self):
        print str(self.opt)
        return [], [], [], []
